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


class Jacobi_2d(PolyBench):

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
        self.TSTEPS = params.get('TSTEPS')
        self.N = params.get('N')

    def run_benchmark(self):
        # Create data structures (arrays, auxiliary variables, etc.)
        A = self.create_array(2, [self.N, self.N], self.DATA_TYPE(0))
        B = self.create_array(2, [self.N, self.N], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(A, B)

        # Benchmark the kernel
        self.time_kernel(A, B)

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
        return [('A', A)]


class _StrategyList(Jacobi_2d):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: list, B: list):
        for i in range(0, self.N):
            for j in range(0, self.N):
                A[i][j] = (self.DATA_TYPE(i) * (j + 2) + 2) / self.N
                B[i][j] = (self.DATA_TYPE(i) * (j + 3) + 3) / self.N

    def print_array_custom(self, A: list, name: str):
        for i in range(0, self.N):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(A[i][j])

    def kernel(self, A: list, B: list):
# scop begin
        for t in range(0, self.TSTEPS):
            for i in range(1, self.N - 1):
                for j in range(1, self.N - 1):
                    B[i][j] = 0.2 * (A[i][j] + A[i][j-1] + A[i][1+j] + A[1+i][j] + A[i-1][j])

            for i in range(1, self.N - 1):
                for j in range(1, self.N - 1):
                    A[i][j] = 0.2 * (B[i][j] + B[i][j-1] + B[i][1+j] + B[1+i][j] + B[i-1][j])
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, A: list, B: list):
# scop begin
        if((self.N-3>= 0) and (self.TSTEPS-1>= 0)):
            for c0 in range ((self.TSTEPS-1)+1):
                for c2 in range (c0 * 2 + 1 , (self.N + c0 * 2-2)+1):
                    B[1][(-2 * c0) + c2] = 0.2 * (A[1][(-2 * c0) + c2] + A[1][(-2 * c0) + c2-1] + A[1][1+(-2 * c0) + c2] + A[1+1][(-2 * c0) + c2] + A[1 -1][(-2 * c0) + c2])
                for c1 in range (c0 * 2 + 2 , (self.N + c0 * 2-2)+1):
                    B[(-2 * c0) + c1][1] = 0.2 * (A[(-2 * c0) + c1][1] + A[(-2 * c0) + c1][1 -1] + A[(-2 * c0) + c1][1+1] + A[1+(-2 * c0) + c1][1] + A[(-2 * c0) + c1-1][1])
                    for c2 in range (c0 * 2 + 2 , (self.N + c0 * 2-2)+1):
                        A[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1] = 0.2 * (B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1] + B[((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1 -1] + B[((-2 * c0) + c1) + -1][1+((-2 * c0) + c2) + -1] + B[1+((-2 * c0) + c1) + -1][((-2 * c0) + c2) + -1] + B[((-2 * c0) + c1) + -1 -1][((-2 * c0) + c2) + -1])
                        B[(-2 * c0) + c1][(-2 * c0) + c2] = 0.2 * (A[(-2 * c0) + c1][(-2 * c0) + c2] + A[(-2 * c0) + c1][(-2 * c0) + c2-1] + A[(-2 * c0) + c1][1+(-2 * c0) + c2] + A[1+(-2 * c0) + c1][(-2 * c0) + c2] + A[(-2 * c0) + c1-1][(-2 * c0) + c2])
                    A[((-2 * c0) + c1) + -1][self.N + -2] = 0.2 * (B[((-2 * c0) + c1) + -1][self.N + -2] + B[((-2 * c0) + c1) + -1][self.N + -2 -1] + B[((-2 * c0) + c1) + -1][1+self.N + -2] + B[1+((-2 * c0) + c1) + -1][self.N + -2] + B[((-2 * c0) + c1) + -1 -1][self.N + -2])
                for c2 in range (c0 * 2 + 2 , (self.N + c0 * 2-1)+1):
                    A[self.N + -2][((-2 * c0) + c2) + -1] = 0.2 * (B[self.N + -2][((-2 * c0) + c2) + -1] + B[self.N + -2][((-2 * c0) + c2) + -1 -1] + B[self.N + -2][1+((-2 * c0) + c2) + -1] + B[1+self.N + -2][((-2 * c0) + c2) + -1] + B[self.N + -2 -1][((-2 * c0) + c2) + -1])
# scop end

class _StrategyListFlattened(Jacobi_2d):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: list, B: list):
        for i in range(0, self.N):
            for j in range(0, self.N):
                A[self.N * i + j] = (self.DATA_TYPE(i) * (j+2) + 2) / self.N
                B[self.N * i + j] = (self.DATA_TYPE(i) * (j+3) + 3) / self.N

    def print_array_custom(self, A: list, name: str):
        for i in range(0, self.N):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(A[self.N * i + j])

    def kernel(self, A: list, B: list):
# scop begin
        for t in range(0, self.TSTEPS):
            for i in range(1, self.N - 1):
                for j in range(1, self.N - 1):
                    B[self.N * i + j] = 0.2 * (A[self.N * i + j] + A[self.N * i + j - 1] + A[self.N * i + 1 + j]
                                               + A[self.N * (1 + i) + j] + A[self.N * (i - 1) + j])

            for i in range(1, self.N - 1):
                for j in range(1, self.N - 1):
                    A[self.N * i + j] = 0.2 * (B[self.N * i + j] + B[self.N * i + j - 1] + B[self.N * i + 1 + j]
                                               + B[self.N * (1 + i) + j] + B[self.N * (i - 1) + j])
# scop end

class _StrategyListFlattenedPluto(Jacobi_2d):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: list, B: list):
        for i in range(0, self.N):
            for j in range(0, self.N):
                A[self.N * i + j] = (self.DATA_TYPE(i) * (j+2) + 2) / self.N
                B[self.N * i + j] = (self.DATA_TYPE(i) * (j+3) + 3) / self.N

    def print_array_custom(self, A: list, name: str):
        for i in range(0, self.N):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(A[self.N * i + j])

    def kernel(self, A: list, B: list):
# scop begin
        for t in range(0, self.TSTEPS):
            for i in range(1, self.N - 1):
                for j in range(1, self.N - 1):
                    B[self.N * i + j] = 0.2 * (A[self.N * i + j] + A[self.N * i + j - 1] + A[self.N * i + 1 + j]
                                               + A[self.N * (1 + i) + j] + A[self.N * (i - 1) + j])

            for i in range(1, self.N - 1):
                for j in range(1, self.N - 1):
                    A[self.N * i + j] = 0.2 * (B[self.N * i + j] + B[self.N * i + j - 1] + B[self.N * i + 1 + j]
                                               + B[self.N * (1 + i) + j] + B[self.N * (i - 1) + j])

            for c0 in range( 0, TSTEPS ):
                if((N+TSTEPS*-1+c0*-1-1>= 0) and (optimization_iterator_1*-1+TSTEPS*-1-1>= 0) and (N+TSTEPS*-1-1>= 0)):
                    B[1]=(0.33333*(((A[(1 +(-1*1))])+(A[1]))+(A[(1 +1)])))
                if((TSTEPS+c0*-1-1>= 0) and (N*-1-1>= 0) and (TSTEPS-1>= 0)):
                    for c1 in range (c0 * 2 + 2 , (N * -1 + c0 * 2)+1):
                        B[(-2 * c0) + c1]=(0.33333*(((A[((-2 * c0) + c1+(-1*1))])+(A[(-2 * c0) + c1]))+(A[((-2 * c0) + c1+1)])))
                        A[((-2 * c0) + c1) + -1]=(0.33333*(((B[(((-2 * c0) + c1) + -1 +(-1*1))])+(B[((-2 * c0) + c1) + -1]))+(B[(((-2 * c0) + c1) + -1 +1)])))
                if((N+TSTEPS*-1+c0*-1-1>= 0) and (optimization_iterator_1*-1+TSTEPS*-1-1>= 0) and (N+TSTEPS*-1-1>= 0)):
                    A[1-N + -1]=(0.33333*(((B[(1-N + -1 +(-1*1))])+(B[1-N + -1]))+(B[(1-N + -1 +1)])))

# scop end

class _StrategyNumPy(Jacobi_2d):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: ndarray, B: ndarray):
        for i in range(0, self.N):
            for j in range(0, self.N):
                A[i, j] = (self.DATA_TYPE(i) * (j+2) + 2) / self.N
                B[i, j] = (self.DATA_TYPE(i) * (j+3) + 3) / self.N

    def print_array_custom(self, A: ndarray, name: str):
        for i in range(0, self.N):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(A[i, j])

    def kernel(self, A: ndarray, B: ndarray):
# scop begin
        for t in range(0, self.TSTEPS):
            B[1:self.N-1, 1:self.N-1] = 0.2 * (A[1:self.N-1, 1:self.N-1] + A[1:self.N-1, 0:self.N-2] + A[1:self.N-1, 2:self.N] + A[2:self.N, 1:self.N-1] + A[0:self.N-2, 1:self.N-1])
            A[1:self.N-1, 1:self.N-1] = 0.2 * (B[1:self.N-1, 1:self.N-1] + B[1:self.N-1, 0:self.N-2] + B[1:self.N-1, 2:self.N] + B[2:self.N, 1:self.N-1] + B[0:self.N-2, 1:self.N-1])
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    # No variations between --pluto / maxfuse / vectorizer
    def kernel(self, A: list, B: list):
# scop begin
        if((self.N-3>= 0) and (self.TSTEPS-1>= 0)):
            for c0 in range ((self.TSTEPS-1)+1):
                for c2 in range (c0 * 2 + 1 , (self.N + c0 * 2-2)+1):
                    B[self.N*(1) + (-2 * c0) + c2] = 0.2 * (A[self.N*(1) + (-2 * c0) + c2] + A[self.N*(1) + (-2 * c0) + c2-1] + A[self.N*(1) + 1+(-2 * c0) + c2] + A[self.N*(1+1) + (-2 * c0) + c2] + A[self.N*(1 -1) + (-2 * c0) + c2])
                for c1 in range (c0 * 2 + 2 , (self.N + c0 * 2-2)+1):
                    B[self.N*((-2 * c0) + c1) + 1] = 0.2 * (A[self.N*((-2 * c0) + c1) + 1] + A[self.N*((-2 * c0) + c1) + 1 -1] + A[self.N*((-2 * c0) + c1) + 1+1] + A[self.N*(1+(-2 * c0) + c1) + 1] + A[self.N*((-2 * c0) + c1-1) + 1])
                    for c2 in range (c0 * 2 + 2 , (self.N + c0 * 2-2)+1):
                        A[self.N*(((-2 * c0) + c1) + -1) + ((-2 * c0) + c2) + -1] = 0.2 * (B[self.N*(((-2 * c0) + c1) + -1) + ((-2 * c0) + c2) + -1] + B[self.N*(((-2 * c0) + c1) + -1) + ((-2 * c0) + c2) + -1 -1] + B[self.N*(((-2 * c0) + c1) + -1) + 1+((-2 * c0) + c2) + -1] + B[self.N*(1+((-2 * c0) + c1) + -1) + ((-2 * c0) + c2) + -1] + B[self.N*(((-2 * c0) + c1) + -1 -1) + ((-2 * c0) + c2) + -1])
                        B[self.N*((-2 * c0) + c1) + (-2 * c0) + c2] = 0.2 * (A[self.N*((-2 * c0) + c1) + (-2 * c0) + c2] + A[self.N*((-2 * c0) + c1) + (-2 * c0) + c2-1] + A[self.N*((-2 * c0) + c1) + 1+(-2 * c0) + c2] + A[self.N*(1+(-2 * c0) + c1) + (-2 * c0) + c2] + A[self.N*((-2 * c0) + c1-1) + (-2 * c0) + c2])
                    A[self.N*(((-2 * c0) + c1) + -1) + self.N + -2] = 0.2 * (B[self.N*(((-2 * c0) + c1) + -1) + self.N + -2] + B[self.N*(((-2 * c0) + c1) + -1) + self.N + -2 -1] + B[self.N*(((-2 * c0) + c1) + -1) + 1+self.N + -2] + B[self.N*(1+((-2 * c0) + c1) + -1) + self.N + -2] + B[self.N*(((-2 * c0) + c1) + -1 -1) + self.N + -2])
                for c2 in range (c0 * 2 + 2 , (self.N + c0 * 2-1)+1):
                    A[self.N*(self.N + -2) + ((-2 * c0) + c2) + -1] = 0.2 * (B[self.N*(self.N + -2) + ((-2 * c0) + c2) + -1] + B[self.N*(self.N + -2) + ((-2 * c0) + c2) + -1 -1] + B[self.N*(self.N + -2) + 1+((-2 * c0) + c2) + -1] + B[self.N*(1+self.N + -2) + ((-2 * c0) + c2) + -1] + B[self.N*(self.N + -2 -1) + ((-2 * c0) + c2) + -1])
# scop end
