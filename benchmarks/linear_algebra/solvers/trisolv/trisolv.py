# Copyright 2021 Universidade da Coruña
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Authors:
#    Miguel Ángel Abella González <miguel.abella@udc.es>
#    Gabriel Rodríguez <gabriel.rodriguez@udc.es>
#
# Contact:
#    Gabriel Rodríguez <gabriel.rodriguez@udc.es>


"""<replace_with_module_description>"""

from benchmarks.polybench import PolyBench
from benchmarks.polybench_classes import ArrayImplementation
from benchmarks.polybench_classes import PolyBenchOptions, PolyBenchSpec
from numpy.core.multiarray import ndarray
import numpy as np


class Trisolv(PolyBench):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        implementation = options.POLYBENCH_ARRAY_IMPLEMENTATION
        if implementation == ArrayImplementation.LIST:
            return _StrategyList.__new__(_StrategyList, options, parameters)
        elif implementation == ArrayImplementation.LIST_PLUTO:
            return _StrategyListPluto.__new__(_StrategyListPluto, options, parameters)
        elif implementation == ArrayImplementation.LIST_FLATTENED:
            return _StrategyListFlattened.__new__(_StrategyListFlattened, options, parameters)
        elif implementation == ArrayImplementation.NUMPY:
            return _StrategyNumPy.__new__(_StrategyNumPy, options, parameters)
        elif implementation == ArrayImplementation.LIST_FLATTENED_PLUTO:
            return _StrategyListFlattenedPluto.__new__(_StrategyListFlattenedPluto, options, parameters)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        # The parameters hold the necessary information obtained from "polybench.spec" file
        params = parameters.DataSets.get(self.DATASET_SIZE)
        if not isinstance(params, dict):
            raise NotImplementedError(f'Dataset size "{self.DATASET_SIZE.name}" not implemented '
                                      f'for {parameters.Category}/{parameters.Name}.')

        # Set up problem size from the given parameters (adapt this part with appropriate parameters)
        self.N = params.get('N')

    def print_array_custom(self, x: list, name: str):
        for i in range(0, self.N):
            self.print_value(x[i])
            if i % 20 == 0:
                self.print_message('\n')

    def run_benchmark(self):
        # Create data structures (arrays, auxiliary variables, etc.)
        L = self.create_array(2, [self.N, self.N], self.DATA_TYPE(0))
        x = self.create_array(1, [self.N], self.DATA_TYPE(0))
        b = self.create_array(1, [self.N], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(L, x, b)

        # Benchmark the kernel
        self.time_kernel(L, x, b)

        # Return printable data as a list of tuples ('name', value).
        # Each tuple element must have the following format:
        #   (A: str, B: matrix)
        #     - A: a representative name for the data (this string will be printed out)
        #     - B: the actual data structure holding the computed result
        #
        # The syntax for the return statement would then be:
        #   - For single data structure results:
        #     return [('data_name', data)]
        #   - For multiple data structure results:
        #     return [('matrix1', m1), ('matrix2', m2), ... ]
        return [('x', x)]


class _StrategyList(Trisolv):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, L: list, x: list, b: list):
        for i in range(0, self.N):
            x[i] = - 999
            b[i] = i
            for j in range(0, i + 1):
                L[i][j] = self.DATA_TYPE(i + self.N - j + 1) * 2 / self.N

    def kernel(self, L: list, x: list, b: list):
# scop begin
        for i in range(0, self.N):
            x[i] = b[i]
            for j in range(0, i):
                x[i] -= L[i][j] * x[j]
            x[i] = x[i] / L[i][i]
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, L: list, x: list, b: list):
# scop begin
        for i in range(0, self.N):
            x[i] = b[i]
            for j in range(0, i):
                x[i] -= L[i][j] * x[j]
            x[i] = x[i] / L[i][i]
# scop end

class _StrategyListFlattened(Trisolv):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        if options.LOAD_ELIMINATION: self.kernel = self.kernel_le
        else: self.kernel = self.kernel_regular

    def initialize_array(self, L: list, x: list, b: list):
        for i in range(0, self.N):
            x[i] = - 999
            b[i] = i
            for j in range(0, i + 1):
                L[self.N * i + j] = self.DATA_TYPE(i+self.N-j+1) * 2 / self.N

    def kernel_regular(self, L: list, x: list, b: list):
# scop begin
        if((self.N-1>= 0)):
            for c1 in range ((self.N-1)+1):
                x[c1] = b[c1]
            x[0] = x[0] / L[self.N*0+0]
            for c1 in range (1 , (self.N-1)+1):
                for c2 in range ((c1-1)+1):
                    x[c1] -= L[self.N * c1 + c2] * x[c2]
                x[c1] = x[c1] / L[self.N * c1 + c1]
# scop end

    def kernel_le(self, L: list, x: list, b: list):
# scop begin
        if((self.N-1>= 0)):
            for c1 in range ((self.N-1)+1):
                x[c1] = b[c1]
            x[0] = x[0] / L[self.N*0+0]
            for c1 in range (1 , (self.N-1)+1):
                tmp = x[c1] # load elimination
                for c2 in range ((c1-1)+1):
                    tmp -= L[self.N * c1 + c2] * x[c2] # load elimination
                x[c1] = tmp / L[self.N * c1 + c1] # load elimination
# scop end

class _StrategyNumPy(Trisolv):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, L: ndarray, x: ndarray, b: ndarray):
        for i in range(0, self.N):
            x[i] = - 999
            b[i] = i
            for j in range(0, i + 1):
                L[i, j] = self.DATA_TYPE(i+self.N-j+1) * 2 / self.N

    def kernel(self, L: ndarray, x: ndarray, b: ndarray):
# scop begin
        x[0:self.N] = b[0:self.N]
        for i in range(0, self.N):
            x[i] -= np.dot( L[i,0:i], x[0:i] )
            x[i] = x[i] / L[i, i]
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        self.kernel_vectorizer = self.kernel_pluto
        self.kernel = getattr( self, "kernel_%s" % (options.POCC) )

    def kernel_pluto(self, L: list, x: list, b: list):
# scop begin
        for i in range(0, self.N):
            x[i] = b[i]
            for j in range(0, i):
                x[i] -= L[self.N*(i) + j] * x[j]
            x[i] = x[i] / L[self.N*(i) + i]
# scop end

    def kernel_maxfuse(self, L: list, x: list, b: list):
# --pluto-fuse maxfuse
# scop begin
        if((self.N-1>= 0)):
            x[0] = b[0]
            x[0] = x[0] / L[(0)*self.N + 0]
            for c0 in range (1 , (self.N-1)+1):
                x[c0] = b[c0]
                for c1 in range (c0 , (c0 * 2-1)+1):
                    x[c0] -= L[(c0)*self.N + (-1 * c0) + c1] * x[(-1 * c0) + c1]
                x[c0] = x[c0] / L[(c0)*self.N + c0]
# scop end
