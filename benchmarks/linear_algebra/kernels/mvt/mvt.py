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


class Mvt(PolyBench):

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

    def print_array_custom(self, array: list, name: str):
        # Although two arrays will be printed, they share the same format, so there is no need to check which one comes.
        for i in range(0, self.N):
            if i % 20 == 0:
                self.print_message('\n')
            self.print_value(array[i])

    def run_benchmark(self):
        # Create data structures (arrays, auxiliary variables, etc.)
        A = self.create_array(2, [self.N, self.N], self.DATA_TYPE(0))
        x1 = self.create_array(1, [self.N], self.DATA_TYPE(0))
        x2 = self.create_array(1, [self.N], self.DATA_TYPE(0))
        y_1 = self.create_array(1, [self.N], self.DATA_TYPE(0))
        y_2 = self.create_array(1, [self.N], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(x1, x2, y_1, y_2, A)

        # Benchmark the kernel
        self.time_kernel(x1, x2, y_1, y_2, A)

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
        return [('x1', x1), ('x2', x2)]


class _StrategyList(Mvt):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, x1: list, x2: list, y_1: list, y_2: list, A: list):
        for i in range(0, self.N):
            x1[i] = self.DATA_TYPE(i % self.N) / self.N
            x2[i] = self.DATA_TYPE((i + 1) % self.N) / self.N
            y_1[i] = self.DATA_TYPE((i + 3) % self.N) / self.N
            y_2[i] = self.DATA_TYPE((i + 4) % self.N) / self.N
            for j in range(0, self.N):
                A[i][j] = self.DATA_TYPE(i * j % self.N) / self.N

    def kernel(self, x1: list, x2: list, y_1: list, y_2: list, A: list):
# scop begin
        for i in range(0, self.N):
            for j in range(0, self.N):
                x1[i] = x1[i] + A[i][j] * y_1[j]

        for i in range(0, self.N):
            for j in range(0, self.N):
                x2[i] = x2[i] + A[j][i] * y_2[j]
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, x1: list, x2: list, y_1: list, y_2: list, A: list):
# scop begin
        if((self.N-1>= 0)):
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.N-1)+1):
                    x1[c1] = x1[c1] + A[c1][c2] * y_1[c2]
                    x2[c1] = x2[c1] + A[c2][c1] * y_2[c2]
# scop end

class _StrategyListFlattened(Mvt):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        if options.LOAD_ELIMINATION: self.kernel = self.kernel_le
        else: self.kernel = self.kernel_regular

    def initialize_array(self, x1: list, x2: list, y_1: list, y_2: list, A: list):
        for i in range(0, self.N):
            x1[i] = self.DATA_TYPE(i % self.N) / self.N
            x2[i] = self.DATA_TYPE((i + 1) % self.N) / self.N
            y_1[i] = self.DATA_TYPE((i + 3) % self.N) / self.N
            y_2[i] = self.DATA_TYPE((i + 4) % self.N) / self.N
            for j in range(0, self.N):
                A[self.N * i + j] = self.DATA_TYPE(i * j % self.N) / self.N

    def kernel_regular(self, x1: list, x2: list, y_1: list, y_2: list, A: list):
# scop begin
        for i in range(0, self.N):
            for j in range(0, self.N):
                x1[i] = x1[i] + A[self.N * i + j] * y_1[j]

        for i in range(0, self.N):
            for j in range(0, self.N):
                x2[i] = x2[i] + A[self.N * j + i] * y_2[j]
# scop end

    def kernel_le(self, x1: list, x2: list, y_1: list, y_2: list, A: list):
# scop begin
        for i in range(0, self.N):
            tmp = x1[i] # load elimination
            for j in range(0, self.N):
                tmp = tmp + A[self.N * i + j] * y_1[j] # load elimination
            x1[i] = tmp # load elimination

        for i in range(0, self.N):
            tmp = x2[i] # load elimination
            for j in range(0, self.N):
                tmp = tmp + A[self.N * j + i] * y_2[j] # load elimination
            x2[i] = tmp # load elimination
# scop end

class _StrategyNumPy(Mvt):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, x1: ndarray, x2: ndarray, y_1: ndarray, y_2: ndarray, A: ndarray):
        for i in range(0, self.N):
            x1[i] = self.DATA_TYPE(i % self.N) / self.N
            x2[i] = self.DATA_TYPE((i + 1) % self.N) / self.N
            y_1[i] = self.DATA_TYPE((i + 3) % self.N) / self.N
            y_2[i] = self.DATA_TYPE((i + 4) % self.N) / self.N
            for j in range(0, self.N):
                A[i, j] = self.DATA_TYPE(i * j % self.N) / self.N

    def kernel(self, x1: ndarray, x2: ndarray, y_1: ndarray, y_2: ndarray, A: ndarray):
# scop begin
        x1[0:self.N] = x1[0:self.N] + np.dot( A[0:self.N,0:self.N], y_1[0:self.N] )
        x2[0:self.N] = x2[0:self.N] + np.dot( A[0:self.N,0:self.N].T, y_2[0:self.N] )
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        self.kernel = getattr( self, "kernel_%s" % (options.POCC) )

    def kernel_pluto(self, x1: list, x2: list, y_1: list, y_2: list, A: list):
# --pluto
# scop begin
        if((self.N-1>= 0)):
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.N-1)+1):
                    x1[c1] = x1[c1] + A[self.N*(c1) + c2] * y_1[c2]
                    x2[c1] = x2[c1] + A[self.N*(c2) + c1] * y_2[c2]
# scop end

    def kernel_vectorizer(self, x1: list, x2: list, y_1: list, y_2: list, A: list):
# --pluto --pluto-prevector --vectorizer --pragmatizer
# scop begin
        if((self.N-1>= 0)):
            for c2 in range ((self.N-1)+1):
                for c1 in range ((self.N-1)+1):
                    x1[c1] = x1[c1] + A[self.N*(c1) + c2] * y_1[c2]
                    x2[c1] = x2[c1] + A[self.N*(c2) + c1] * y_2[c2]
# scop end

    def kernel_maxfuse(self, x1: list, x2: list, y_1: list, y_2: list, A: list):
# --pluto --pluto-fuse maxfuse
# scop begin
        if((self.N-1>= 0)):
            for c0 in range ((self.N-1)+1):
                for c1 in range ((self.N-1)+1):
                    x1[c0] = x1[c0] + A[(c0)*self.N + c1] * y_1[c1]
                    x2[c0] = x2[c0] + A[(c1)*self.N + c0] * y_2[c1]
# scop end
