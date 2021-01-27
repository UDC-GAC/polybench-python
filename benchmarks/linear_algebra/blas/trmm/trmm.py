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


class Trmm(PolyBench):

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

    def run_benchmark(self):
        # Create data structures (arrays, auxiliary variables, etc.)
        alpha = 1.5

        A = self.create_array(2, [self.M, self.M], self.DATA_TYPE(0))
        B = self.create_array(2, [self.M, self.N], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(alpha, A, B)

        # Benchmark the kernel
        self.time_kernel(alpha, A, B)

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
        return [('B', B)]


class _StrategyList(Trmm):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, alpha, A: list, B: list):
        for i in range(0, self.M):
            for j in range(0, i):
                A[i][j] = self.DATA_TYPE((i + j) % self.M) / self.M

            A[i][i] = 1.0

            for j in range(0, self.N):
                B[i][j] = self.DATA_TYPE((self.N + (i - j)) % self.N) / self.N

    def print_array_custom(self, B: list, name: str):
        for i in range(0, self.M):
            for j in range(0, self.N):
                if (i * self.M + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(B[i][j])

    def kernel(self, alpha, A: list, B: list):
        # BLAS parameters
        # SIDE   = 'L'
        # UPLO   = 'L'
        # TRANSA = 'T'
        # DIAG   = 'U'
        # = > Form  B := alpha * A ** T * B.
        # A is MxM
        # B is MxN
# scrop begin
        for i in range(0, self.M):
            for j in range(0, self.N):
                for k in range(i + 1, self.M):
                    B[i][j] += A[k][i] * B[k][j]
                B[i][j] = alpha * B[i][j]
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, alpha, A: list, B: list):
# scop begin
        if((self.M-1>= 0) and (self.N-1>= 0)):
            if((self.M-2>= 0)):
                for c1 in range ((self.N-1)+1):
                    for c2 in range ((self.M-2)+1):
                        for c3 in range (c2 + 1 , (self.M-1)+1):
                            B[c2][c1] += A[c3][c2] * B[c3][c1]
            for c1 in range ((self.M-1)+1):
                for c2 in range ((self.N-1)+1):
                    B[c1][c2] = alpha * B[c1][c2]
# scop end

class _StrategyListFlattened(Trmm):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        if options.LOAD_ELIMINATION: self.kernel = self.kernel_le
        else: self.kernel = self.kernel_regular

    def initialize_array(self, alpha, A: list, B: list):
        for i in range(0, self.M):
            for j in range(0, i):
                A[self.M * i + j] = self.DATA_TYPE((i+j) % self.M) / self.M

            A[self.M * i + i] = 1.0

            for j in range(0, self.N):
                B[self.N * i + j] = self.DATA_TYPE((self.N+(i-j)) % self.N) / self.N

    def print_array_custom(self, B: list, name: str):
        for i in range(0, self.M):
            for j in range(0, self.N):
                if (i * self.M + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(B[self.N * i + j])

    # Regular version
    def kernel_regular(self, alpha, A: list, B: list):
# scrop begin
        for i in range(0, self.M):
            for j in range(0, self.N):
                for k in range(i + 1, self.M):
                    B[self.N * i + j] += A[self.M * k + i] * B[self.N * k + j]
                B[self.N * i + j] = alpha * B[self.N * i + j]
# scop end

    # Load elimination
    def kernel_le(self, alpha, A: list, B: list): 
# scrop begin
        for i in range(0, self.M):
            for j in range(0, self.N):
                tmp = B[self.N*i+j]
                for k in range(i + 1, self.M):
                    tmp += A[self.M * k + i] * B[self.N * k + j]
                B[self.N * i + j] = alpha * tmp
# scop end

class _StrategyNumPy(Trmm):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, alpha, A: list, B: list):
        for i in range(0, self.M):
            for j in range(0, i):
                A[i, j] = self.DATA_TYPE((i + j) % self.M) / self.M

            A[i, i] = 1.0

            for j in range(0, self.N):
                B[i, j] = self.DATA_TYPE((self.N + (i - j)) % self.N) / self.N

    def print_array_custom(self, B: ndarray, name: str):
        for i in range(0, self.M):
            for j in range(0, self.N):
                if (i * self.M + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(B[i, j])

    def kernel(self, alpha, A: ndarray, B: ndarray):
        # BLAS parameters
        # SIDE   = 'L'
        # UPLO   = 'L'
        # TRANSA = 'T'
        # DIAG   = 'U'
        # = > Form  B := alpha * A ** T * B.
        # A is MxM
        # B is MxN
# scop begin
        for i in range(0, self.M):
            B[i,0:self.N] += (A[i+1:self.M,i,np.newaxis] * B[i+1:self.M,0:self.N]).sum(axis=0)
            B[i,0:self.N] = alpha * B[i,0:self.N]
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        self.kernel_vectorizer = self.kernel_pluto
        self.kernel = getattr( self, "kernel_%s" % (options.POCC) )

    def kernel_pluto(self, alpha, A: list, B: list):
# --pluto
# scop begin
        if((self.M-1>= 0) and (self.N-1>= 0)):
            if((self.M-2>= 0)):
                for c1 in range ((self.N-1)+1):
                    for c2 in range ((self.M-2)+1):
                        for c3 in range (c2 + 1 , (self.M-1)+1):
                            B[self.N*(c2) + c1] += A[self.M*(c3) + c2] * B[self.N*(c3) + c1]
            for c1 in range ((self.M-1)+1):
                for c2 in range ((self.N-1)+1):
                    B[self.N*(c1) + c2] = alpha * B[self.N*(c1) + c2]
# scop end

    def kernel_maxfuse(self, alpha, A: list, B: list):
# --pluto --pluto-fuse maxfuse
# scop begin
        if((self.M-1>= 0) and (self.N-1>= 0)):
            if((self.M-2>= 0)):
                for c0 in range ((self.N-1)+1):
                    for c1 in range ((self.M-2)+1):
                        for c4 in range (c1 + 1 , (self.M-1)+1):
                            B[(c1)*self.N + c0] += A[(c4)*self.M + c1] * B[(c4)*self.N + c0]
                        B[(c1)*self.N + c0] = alpha * B[(c1)*self.N + c0]
                    B[(self.M + -1)*self.N + c0] = alpha * B[(self.M + -1)*self.N + c0]
            if self.M == 1:
                for c0 in range ((self.N-1)+1):
                    B[(0)*self.N + c0] = alpha * B[(0)*self.N + c0]
# scop end
