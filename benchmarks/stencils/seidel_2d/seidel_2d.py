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


class Seidel_2d(PolyBench):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        implementation = options.POLYBENCH_ARRAY_IMPLEMENTATION
        if implementation == ArrayImplementation.LIST:
            return _StrategyList.__new__(_StrategyList, options, parameters)
        if implementation == ArrayImplementation.LIST_PLUTO:
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

        # Initialize data structures
        self.initialize_array(A)

        # Benchmark the kernel
        self.time_kernel(A)

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


class _StrategyList(Seidel_2d):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: list):
        for i in range(0, self.N):
            for j in range(0, self.N):
                A[i][j] = (self.DATA_TYPE(i) * (j + 2) + 2) / self.N

    def print_array_custom(self, A: list, name: str):
        for i in range(0, self.N):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(A[i][j])

    def kernel(self, A: list):
# scop begin
        for t in range(0, self.TSTEPS - 1 + 1):
            for i in range(1, self.N - 2 + 1):
                for j in range(1, self.N - 2 + 1):
                    A[i][j] = (A[i - 1][j - 1] + A[i - 1][j] + A[i - 1][j + 1]
                               + A[i][j - 1] + A[i][j] + A[i][j + 1]
                               + A[i + 1][j - 1] + A[i + 1][j] + A[i + 1][j + 1]) / 9.0
#scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, A: list):
# scop begin
        if((self.N-3>= 0) and (self.TSTEPS-1>= 0)):
            for c0 in range ((self.TSTEPS-1)+1):
                for c1 in range (c0 + 1 , (self.N + c0-2)+1):
                    for c2 in range (c0 + c1 + 1 , (self.N + c0 + c1-2)+1):
                        A[(-1 * c0) + c1][((-1 * c0) + (-1 * c1)) + c2] = (A[(-1 * c0) + c1-1][((-1 * c0) + (-1 * c1)) + c2-1] + A[(-1 * c0) + c1-1][((-1 * c0) + (-1 * c1)) + c2] + A[(-1 * c0) + c1-1][((-1 * c0) + (-1 * c1)) + c2+1] + A[(-1 * c0) + c1][((-1 * c0) + (-1 * c1)) + c2-1] + A[(-1 * c0) + c1][((-1 * c0) + (-1 * c1)) + c2] + A[(-1 * c0) + c1][((-1 * c0) + (-1 * c1)) + c2+1] + A[(-1 * c0) + c1+1][((-1 * c0) + (-1 * c1)) + c2-1] + A[(-1 * c0) + c1+1][((-1 * c0) + (-1 * c1)) + c2] + A[(-1 * c0) + c1+1][((-1 * c0) + (-1 * c1)) + c2+1])/9.0
#scop end

class _StrategyListFlattened(Seidel_2d):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: list):
        for i in range(0, self.N):
            for j in range(0, self.N):
                A[self.N * i + j] = (self.DATA_TYPE(i)*(j+2) + 2) / self.N

    def print_array_custom(self, A: list, name: str):
        for i in range(0, self.N):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(A[self.N * i + j])

    def kernel(self, A: list):
# scop begin
        for t in range(0, self.TSTEPS - 1 +1):
            for i in range(1, self.N - 2 + 1):
                for j in range(1, self.N - 2 + 1):
                    A[self.N * i + j] = (A[self.N * (i - 1) + j - 1] + A[self.N * (i - 1) + j] + A[self.N * (i - 1) + j + 1]
                                         + A[self.N * i + j - 1] + A[self.N * i + j] + A[self.N * i + j + 1]
                                         + A[self.N * (i + 1) + j - 1] + A[self.N * (i + 1) + j] + A[self.N * (i + 1) + j + 1]) / 9.0
# scop end


class _StrategyNumPy(Seidel_2d):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: ndarray):
        for i in range(0, self.N):
            for j in range(0, self.N):
                A[i, j] = (self.DATA_TYPE(i)*(j+2) + 2) / self.N

    def print_array_custom(self, A: ndarray, name: str):
        for i in range(0, self.N):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(A[i, j])

    def kernel(self, A: ndarray):
# scop begin
#        for t in range(0, self.TSTEPS - 1 +1):
#            for i in range(1, self.N - 2 + 1):
#                for j in range(1, self.N - 2 + 1):
#                    A[i,j] = (A[i - 1,j - 1] + A[i - 1,j] + A[i - 1,j + 1]              # Dependences prevent vectorization: A[i,j-1] -> A[i,j]
#                               + A[i,j - 1] + A[i,j] + A[i,j + 1]
#                               + A[i + 1,j - 1] + A[i + 1,j] + A[i + 1,j + 1]) / 9.0
        A_flattened = A.ravel()
        self_N_1 = self.N-1
        self_Nx2 = self.N*2
        for t in range( self.TSTEPS ):
            for i in range( 1, self_N_1 ):
                
                end = i * self.N
                slice_NW = slice(i-1, end-1, self_N_1 )
                slice_N = slice(i, end, self_N_1 )
                slice_NE = slice(i+1, end+1, self_N_1 )

                start = self.N+i
                end += self.N
                slice_W = slice(start-1, end-1, self_N_1)
                slice_center = slice(start, end, self_N_1)
                slice_E = slice( start+1, end+1, self_N_1)

                start += self.N
                end += self.N
                slice_SW = slice( start-1, end-1, self_N_1)
                slice_S = slice( start, end, self_N_1)
                slice_SE = slice( start+1, end, self_N_1)

                A_flattened[slice_center] = ( A_flattened[slice_NW] + A_flattened[slice_N] + A_flattened[slice_NE]+ 
                                              A_flattened[slice_W] + A_flattened[slice_center] + A_flattened[slice_E] + 
                                              A_flattened[slice_SW] + A_flattened[slice_S] + A_flattened[slice_SE] ) / 9.0
# diagonal version
#                diag_i = np.arange(1,i+1)
#                diag_j = np.arange(i,0,-1)
#                A[diag_i,diag_j] = ( A[diag_i-1,diag_j-1] + A[diag_i-1,diag_j] + A[diag_i-1,diag_j+1]+ A[diag_i,diag_j-1] + A[diag_i,diag_j] + A[diag_i,diag_j+1] + A[diag_i+1,diag_j-1] + A[diag_i+1,diag_j] + A[diag_i+1,diag_j+1] ) / 9.0
# scalar version
#                for j in range( 1, i+1 ):
#                    A[j,i-j+1] = (A[j-1,i-j] + A[j-1,i-j+1] + A[j-1,i-j+2]
#                        + A[j,i-j] + A[j,i-j+1] + A[j,i-j+2]
#                        + A[j+1,i-j] + A[j+1,i-j+1] + A[j+1,i-j+2] ) / 9.0
            for i in range( 2, self.N ):
                start = self.N*i-2
                end = self.N**2 -3*self.N +i +1
                slice_NW = slice( start-1, end-1, self.N-1 )
                slice_N = slice( start, end, self.N-1 )
                slice_NE = slice( start+1, end+1, self.N-1 )

                start += self.N
                end += self.N
                slice_W = slice( start-1, end-1, self.N-1 )
                slice_center = slice( start, end, self.N-1 )
                slice_E = slice( start+1, end+1, self.N-1 )

                start += self.N
                end += self.N
                slice_SW = slice( start-1, end-1, self.N-1 )
                slice_S = slice( start, end, self.N-1 )
                slice_SE = slice( start+1, end+1, self.N-1 )
                A_flattened[slice_center] = ( A_flattened[slice_NW] + A_flattened[slice_N] + A_flattened[slice_NE]+ 
                                              A_flattened[slice_W] + A_flattened[slice_center] + A_flattened[slice_E] + 
                                              A_flattened[slice_SW] + A_flattened[slice_S] + A_flattened[slice_SE] ) / 9.0
# diagonal version
#                diag_i = np.arange( i, self.N-1 )
#                diag_j = np.arange( self.N-2, i-1, -1 )
#                A[diag_i,diag_j] = ( A[diag_i-1,diag_j-1] + A[diag_i-1,diag_j] + A[diag_i-1,diag_j+1]+ A[diag_i,diag_j-1] + A[diag_i,diag_j] + A[diag_i,diag_j+1] + A[diag_i+1,diag_j-1] + A[diag_i+1,diag_j] + A[diag_i+1,diag_j+1] ) / 9.0
# scalar version
#                for j in range( self.N-2, i-1, -1 ):
#                    A[i+self.N-j-2,j] = ( A[i+self.N-j-3,j-1] + A[i+self.N-j-3,j] + A[i+self.N-j-3,j+1] 
#                            + A[i+self.N-j-2,j-1] + A[i+self.N-j-2,j] + A[i+self.N-j-2,j+1] 
#                            + A[i+self.N-j-1,j-1] + A[i+self.N-j-1,j] + A[i+self.N-j-1,j+1] ) / 9.0
#scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    # No variations between --pluto / maxfuse / vectorizer
    def kernel(self, A: list):
# scop begin
        if((self.N-3>= 0) and (self.TSTEPS-1>= 0)):
            for c0 in range ((self.TSTEPS-1)+1):
                for c1 in range (c0 + 1 , (self.N + c0-2)+1):
                    for c2 in range (c0 + c1 + 1 , (self.N + c0 + c1-2)+1):
                        A[self.N*((-1 * c0) + c1) + ((-1 * c0) + (-1 * c1)) + c2] = (A[self.N*((-1 * c0) + c1-1) + ((-1 * c0) + (-1 * c1)) + c2-1] + A[self.N*((-1 * c0) + c1-1) + ((-1 * c0) + (-1 * c1)) + c2] + A[self.N*((-1 * c0) + c1-1) + ((-1 * c0) + (-1 * c1)) + c2+1] + A[self.N*((-1 * c0) + c1) + ((-1 * c0) + (-1 * c1)) + c2-1] + A[self.N*((-1 * c0) + c1) + ((-1 * c0) + (-1 * c1)) + c2] + A[self.N*((-1 * c0) + c1) + ((-1 * c0) + (-1 * c1)) + c2+1] + A[self.N*((-1 * c0) + c1+1) + ((-1 * c0) + (-1 * c1)) + c2-1] + A[self.N*((-1 * c0) + c1+1) + ((-1 * c0) + (-1 * c1)) + c2] + A[self.N*((-1 * c0) + c1+1) + ((-1 * c0) + (-1 * c1)) + c2+1])/9.0
#scop end
