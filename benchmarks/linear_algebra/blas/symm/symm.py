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

class Symm(PolyBench):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        implementation = options.POLYBENCH_ARRAY_IMPLEMENTATION
        if implementation == ArrayImplementation.LIST:
            return _StrategyList.__new__(_StrategyList, options, parameters)
        elif implementation == ArrayImplementation.LIST_PLUTO:
            return _StrategyListPluto.__new__(_StrategyList, options, parameters)
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

    def run_benchmark(self):
        # Create data structures (arrays, auxiliary variables, etc.)
        alpha = 1.5
        beta = 1.2

        C = self.create_array(2, [self.M, self.N], self.DATA_TYPE(0))
        A = self.create_array(2, [self.M, self.M], self.DATA_TYPE(0))
        B = self.create_array(2, [self.M, self.N], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(alpha, beta, C, A, B)

        # Benchmark the kernel
        self.time_kernel(alpha, beta, C, A, B)

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
        return [('C', C)]


class _StrategyList(Symm):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, alpha, beta, C: list, A: list, B: list):
        for i in range(0, self.M):
            for j in range(0, self.N):
                C[i][j] = self.DATA_TYPE((i+j) % 100) / self.M
                B[i][j] = self.DATA_TYPE((self.N+i-j) % 100) / self.M

        for i in range(0, self.M):
            for j in range(0, i + 1):
                A[i][j] = self.DATA_TYPE((i+j) % 100) / self.M
            for j in range(i + 1, self.M):
                A[i][j] = -999  # regions of arrays that should not be used

    def print_array_custom(self, C: list, name: str):
        for i in range(0, self.M):
            for j in range(0, self.N):
                if (i * self.M + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(C[i][j])

    def kernel(self, alpha, beta, C: list, A: list, B: list):
        # BLAS PARAMS
        # SIDE = 'L'
        # UPLO = 'L'
        #  =>  Form  C := alpha*A*B + beta*C
        #  A is MxM
        #  B is MxN
        #  C is MxN
        # note that due to Fortran array layout, the code below more closely resembles upper triangular case in BLAS
# scop begin
        for i in range(0, self.M):
            for j in range(0, self.N):
                temp2 = 0
                for k in range(0, i):
                    C[k][j] += alpha * B[i][j] * A[i][k]
                    temp2 += B[k][j] * A[i][k]
                C[i][j] = beta * C[i][j] + alpha * B[i][j] * A[i][i] + alpha * temp2
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, alpha, beta, C: list, A: list, B: list):
# scop begin
        if((self.M-1>= 0) and (self.N-1>= 0)):
            for c2 in range ((self.N-1)+1):
                temp2 = 0
                C[0][c2] = beta * C[0][c2] + alpha*B[0][c2] * A[0][0] + alpha * temp2
            for c1 in range (1 , (self.M-1)+1):
                for c2 in range ((self.N-1)+1):
                    C[0][c2] += alpha*B[c1][c2] * A[c1][0]
                    temp2 = 0
                    temp2 += B[0][c2] * A[c1][0]
                    for c3 in range (1 , (c1-1)+1):
                        C[c3][c2] += alpha*B[c1][c2] * A[c1][c3]
                        temp2 += B[c3][c2] * A[c1][c3]
                    C[c1][c2] = beta * C[c1][c2] + alpha*B[c1][c2] * A[c1][c1] + alpha * temp2
# scop end

class _StrategyListFlattened(Symm):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, alpha, beta, C: list, A: list, B: list):
        for i in range(0, self.M):
            for j in range(0, self.N):
                C[self.N * i + j] = self.DATA_TYPE((i + j) % 100) / self.M
                B[self.N * i + j] = self.DATA_TYPE((self.N + i - j) % 100) / self.M

        for i in range(0, self.M):
            for j in range(0, i + 1):
                A[self.M * i + j] = self.DATA_TYPE((i + j) % 100) / self.M
            for j in range(i + 1, self.M):
                A[self.M * i + j] = -999  # regions of arrays that should not be used

    def print_array_custom(self, C: list, name: str):
        for i in range(0, self.M):
            for j in range(0, self.N):
                if (i * self.M + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(C[self.N * i + j])

    def kernel(self, alpha, beta, C: list, A: list, B: list):
# scop begin
        for i in range(0, self.M):
            for j in range(0, self.N):
                temp2 = 0
                for k in range(0, i):
                    C[self.N * k + j] += alpha * B[self.N * i + j] * A[self.M * i + k]
                    temp2 += B[self.N * k + j] * A[self.M * i + k]
                C[self.N * i + j] = beta * C[self.N * i + j] + alpha * B[self.N * i + j] * A[self.M * i + i] + alpha * temp2
# scop end


class _StrategyNumPy(Symm):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, alpha, beta, C: list, A: list, B: list):
        for i in range(0, self.M):
            for j in range(0, self.N):
                C[i, j] = self.DATA_TYPE((i+j) % 100) / self.M
                B[i, j] = self.DATA_TYPE((self.N+i-j) % 100) / self.M

        for i in range(0, self.M):
            for j in range(0, i + 1):
                A[i, j] = self.DATA_TYPE((i+j) % 100) / self.M
            for j in range(i + 1, self.M):
                A[i, j] = -999  # regions of arrays that should not be used

    def print_array_custom(self, C: ndarray, name: str):
        for i in range(0, self.M):
            for j in range(0, self.N):
                if (i * self.M + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(C[i, j])

    def kernel(self, alpha, beta, C: ndarray, A: ndarray, B: ndarray):
        # BLAS PARAMS
        # SIDE = 'L'
        # UPLO = 'L'
        #  =>  Form  C := alpha*A*B + beta*C
        #  A is MxM
        #  B is MxN
        #  C is MxN
        # note that due to Fortran array layout, the code below more closely resembles upper triangular case in BLAS
# scop begin
        for i in range(0, self.M):
            C[0:i,0:self.N] += alpha * np.dot( A[i,0:i,np.newaxis], B[np.newaxis,i,0:self.N] )
            temp2 = np.dot( A[i,0:i], B[0:i,0:self.N] )
            C[i,0:self.N] = beta * C[i,0:self.N] + alpha * B[i,0:self.N] * A[i,i] + alpha * temp2[0:self.N]
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        self.kernel_vectorizer = self.kernel_pluto
        self.kernel = getattr( self, "kernel_%s" % (options.POCC) )

    def kernel_pluto(self, alpha, beta, C: list, A: list, B: list):
# --pluto / --pluto-prevector --vectorizer --pragmatizer
# scop begin
        if((self.M-1>= 0) and (self.N-1>= 0)):
            for c2 in range ((self.N-1)+1):
                temp2 = 0
                C[self.N*(0) + c2] = beta * C[self.N*(0) + c2] + alpha*B[self.N*(0) + c2] * A[self.M*(0) + 0] + alpha * temp2
            for c1 in range (1 , (self.M-1)+1):
                for c2 in range ((self.N-1)+1):
                    C[self.N*(0) + c2] += alpha*B[self.N*(c1) + c2] * A[self.M*(c1) + 0]
                    temp2 = 0
                    temp2 += B[self.N*(0) + c2] * A[self.M*(c1) + 0]
                    for c3 in range (1 , (c1-1)+1):
                        C[self.N*(c3) + c2] += alpha*B[self.N*(c1) + c2] * A[self.M*(c1) + c3]
                        temp2 += B[self.N*(c3) + c2] * A[self.M*(c1) + c3]
                    C[self.N*(c1) + c2] = beta * C[self.N*(c1) + c2] + alpha*B[self.N*(c1) + c2] * A[self.M*(c1) + c1] + alpha * temp2
# scop end

    def kernel_maxfuse(self, alpha, beta, C: list, A: list, B: list):
# --pluto --pluto-fuse maxfuse
        if((self.M-1>= 0) and (self.N-1>= 0)):
            for c1 in range ((self.N-1)+1):
                temp2 = 0
                C[(0)*self.N + c1] = beta * C[(0)*self.N + c1] + alpha*B[(0)*self.N + c1] * A[(0)*self.M + 0] + alpha * temp2
            for c0 in range (1 , (self.M-1)+1):
                for c1 in range ((self.N-1)+1):
                    C[(0)*self.N + c1] += alpha*B[(c0)*self.N + c1] * A[(c0)*self.M + 0]
                    temp2 = 0
                    temp2 += B[(0)*self.N + c1] * A[(c0)*self.M + 0]
                    for c2 in range (1 , (c0-1)+1):
                        C[(c2)*self.N + c1] += alpha*B[(c0)*self.N + c1] * A[(c0)*self.M + c2]
                        temp2 += B[(c2)*self.N + c1] * A[(c0)*self.M + c2]
                    C[(c0)*self.N + c1] = beta * C[(c0)*self.N + c1] + alpha*B[(c0)*self.N + c1] * A[(c0)*self.M + c0] + alpha * temp2
# scop end
