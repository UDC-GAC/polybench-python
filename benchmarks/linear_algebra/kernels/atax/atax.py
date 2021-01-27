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


class Atax(PolyBench):

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
        self.M = params.get('M')
        self.N = params.get('N')

    def print_array_custom(self, y: list, name: str):
        for i in range(0, self.N):
            if i % 20 == 0:
                self.print_message('\n')
            self.print_value(y[i])

    def run_benchmark(self):
        # Create data structures (arrays, auxiliary variables, etc.)
        A = self.create_array(2, [self.M, self.N], self.DATA_TYPE(0))
        x = self.create_array(1, [self.N])
        y = self.create_array(1, [self.N])
        tmp = self.create_array(1, [self.M])

        # Initialize data structures
        self.initialize_array(A, x, y, tmp)

        # Benchmark the kernel
        self.time_kernel(A, x, y, tmp)

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
        return [('y', y)]


class _StrategyList(Atax):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: list, x: list, y: list, tmp: list):
        fn = self.DATA_TYPE(self.N)

        for i in range(0, self.N):
            x[i] = 1 + (i / fn)

        for i in range(0, self.M):
            for j in range(0, self.N):
                A[i][j] = self.DATA_TYPE((i + j) % self.N) / (5 * self.M)

    def kernel(self, A: list, x: list, y: list, tmp: list):
# scop begin
        for i in range(0, self.N):
            y[i] = 0

        for i in range(0, self.M):
            tmp[i] = 0.0
            for j in range(0, self.N):
                tmp[i] = tmp[i] + A[i][j] * x[j]

            for j in range(0, self.N):
                y[j] = y[j] + A[i][j] * tmp[i]
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, A: list, x: list, y: list, tmp: list):
# scop begin
        for c1 in range ((self.M-1)+1):
            tmp[c1] = 0.0
        if((self.N-1>= 0)):
            for c1 in range ((self.M-1)+1):
                for c2 in range ((self.N-1)+1):
                    tmp[c1] = tmp[c1] + A[c1][c2] * x[c2]
        for c1 in range ((self.N-1)+1):
            y[c1] = 0
        if((self.M-1>= 0)):
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.M-1)+1):
                    y[c1] = y[c1] + A[c2][c1] * tmp[c2]
# scop end

class _StrategyListFlattened(Atax):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: list, x: list, y: list, tmp: list):
        fn = self.DATA_TYPE(self.N)

        for i in range(0, self.N):
            x[i] = 1 + (i / fn)

        for i in range(0, self.M):
            for j in range(0, self.N):
                A[self.N * i + j] = self.DATA_TYPE((i + j) % self.N) / (5 * self.M)

    def kernel(self, A: list, x: list, y: list, tmp: list):
# scop begin
        for i in range(0, self.N):
            y[i] = 0

        for i in range(0, self.M):
            tmp[i] = 0.0
            for j in range(0, self.N):
                tmp[i] = tmp[i] + A[self.N * i + j] * x[j]

            for j in range(0, self.N):
                y[j] = y[j] + A[self.N * i + j] * tmp[i]
# scop end


class _StrategyNumPy(Atax):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: list, x: list, y: list, tmp: list):
        fn = self.DATA_TYPE(self.N)

        for i in range(0, self.N):
            x[i] = 1 + (i / fn)

        for i in range(0, self.M):
            for j in range(0, self.N):
                A[i, j] = self.DATA_TYPE((i + j) % self.N) / (5 * self.M)

    def kernel(self, A: ndarray, x: ndarray, y: ndarray, tmp: ndarray):
# scop begin
        tmp[0:self.M] = tmp[0:self.M] + np.dot( A[0:self.M,0:self.N], x[0:self.N] )
        y[0:self.N] = np.dot( A[0:self.M].T, tmp[0:self.M] )
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        self.kernel = getattr( self, "kernel_%s" % (options.POCC) )

    def kernel_pluto(self, A: list, x: list, y: list, tmp: list):
# --pluto
# scop begin
        for c1 in range ((self.M-1)+1):
            tmp[c1] = 0.0
        if((self.N-1>= 0)):
            for c1 in range ((self.M-1)+1):
                for c2 in range ((self.N-1)+1):
                    tmp[c1] = tmp[c1] + A[self.N*(c1) + c2] * x[c2]
        for c1 in range ((self.N-1)+1):
            y[c1] = 0
        if((self.M-1>= 0)):
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.M-1)+1):
                    y[c1] = y[c1] + A[self.N*(c2) + c1] * tmp[c2]
# scop end

    def kernel_vectorizer(self, A: list, x: list, y: list, tmp: list):
# --pluto --pluto-prevector --vectorizer --pragmatizer
# scop begin
        for c1 in range ((self.M-1)+1):
            tmp[c1] = 0.0
        if((self.N-1>= 0)):
            for c1 in range ((self.M-1)+1):
                for c2 in range ((self.N-1)+1):
                    tmp[c1] = tmp[c1] + A[self.N*(c1) + c2] * x[c2]
        for c1 in range ((self.N-1)+1):
            y[c1] = 0
        if((self.M-1>= 0)):
            for c2 in range ((self.M-1)+1):
                for c1 in range ((self.N-1)+1):
                    y[c1] = y[c1] + A[self.N*(c2) + c1] * tmp[c2]
# scop end

    def kernel_maxfuse(self, A: list, x: list, y: list, tmp: list):
# --pluto --pluto-fuse maxfuse
# scop begin
        for c0 in range (min((self.M-1)+1 , (self.N-1)+1)):
            tmp[c0] = 0.0
            for c3 in range ((self.N-1)+1):
                tmp[c0] = tmp[c0] + A[(c0)*self.N + c3] * x[c3]
            y[c0] = 0
            for c3 in range ((c0)+1):
                y[c0 + (-1 * c3)] = y[c0 + (-1 * c3)] + A[(c3)*self.N + c0 + (-1 * c3)] * tmp[c3]
        if((self.M-1>= 0)):
            for c0 in range (self.M , (self.N-1)+1):
                y[c0] = 0
                for c3 in range ((self.M-1)+1):
                    y[c0 + (-1 * c3)] = y[c0 + (-1 * c3)] + A[(c3)*self.N + c0 + (-1 * c3)] * tmp[c3]
        if((self.M*-1>= 0)):
            for c0 in range ((self.N-1)+1):
                y[c0] = 0
        if((self.N-1>= 0)):
            for c0 in range (self.N , (self.M-1)+1):
                tmp[c0] = 0.0
                for c3 in range ((self.N-1)+1):
                    tmp[c0] = tmp[c0] + A[(c0)*self.N + c3] * x[c3]
                for c3 in range (self.N * -1 + c0 + 1 , (c0)+1):
                    y[c0 + (-1 * c3)] = y[c0 + (-1 * c3)] + A[(c3)*self.N + c0 + (-1 * c3)] * tmp[c3]
        if((self.N*-1>= 0)):
            for c0 in range ((self.M-1)+1):
                tmp[c0] = 0.0
        for c0 in range (max(self.M , self.N) , (self.M + self.N-2)+1):
            for c3 in range (self.N * -1 + c0 + 1 , (self.M-1)+1):
                y[c0 + (-1 * c3)] = y[c0 + (-1 * c3)] + A[(c3)*self.N + c0 + (-1 * c3)] * tmp[c3]
# scop end
