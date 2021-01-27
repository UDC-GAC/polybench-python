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


class Syr2k(PolyBench):

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
        beta = 1.2

        C = self.create_array(2, [self.N, self.N], self.DATA_TYPE(0))
        A = self.create_array(2, [self.N, self.M], self.DATA_TYPE(0))
        B = self.create_array(2, [self.N, self.M], self.DATA_TYPE(0))

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


class _StrategyList(Syr2k):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, alpha, beta, C: list, A: list, B: list):
        for i in range(0, self.N):
            for j in range(0, self.M):
                A[i][j] = self.DATA_TYPE((i * j + 1) % self.N) / self.N
                B[i][j] = self.DATA_TYPE((i * j + 2) % self.M) / self.M

        for i in range(0, self.N):
            for j in range(0, self.N):
                C[i][j] = self.DATA_TYPE((i * j + 3) % self.N) / self.M

    def print_array_custom(self, C: list, name: str):
        for i in range(0, self.N):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(C[i][j])

    def kernel(self, alpha, beta, C: list, A: list, B: list):
        # BLAS PARAMS
        # UPLO  = 'L'
        # TRANS = 'N'
        # A is NxM
        # B is NxM
        # C is NxN
# scop begin
        for i in range(0, self.N):
            for j in range(0, i + 1):
                C[i][j] *= beta

            for k in range(0, self.M):
                for j in range(0, i + 1):
                    C[i][j] += A[j][k]*alpha*B[i][k] + B[j][k]*alpha*A[i][k]
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, alpha, beta, C: list, A: list, B: list):
# scop begin
        if((self.N-1>= 0)):
            for c1 in range ((self.N-1)+1):
                for c2 in range ((c1)+1):
                    C[c1][c2] *= beta
            if((self.M-1>= 0)):
                for c1 in range ((self.N-1)+1):
                    for c2 in range ((c1)+1):
                        for c3 in range ((self.M-1)+1):
                            C[c1][c2] += A[c2][c3]*alpha*B[c1][c3] + B[c2][c3]*alpha*A[c1][c3]
# scop end

class _StrategyListFlattened(Syr2k):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, alpha, beta, C: list, A: list, B: list):
        for i in range(0, self.N):
            for j in range(0, self.M):
                A[self.M * i + j] = self.DATA_TYPE((i * j + 1) % self.N) / self.N
                B[self.M * i + j] = self.DATA_TYPE((i * j + 2) % self.M) / self.M

        for i in range(0, self.N):
            for j in range(0, self.N):
                C[self.N * i + j] = self.DATA_TYPE((i * j + 3) % self.N) / self.M

    def print_array_custom(self, C: list, name: str):
        for i in range(0, self.N):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(C[self.N * i + j])

    def kernel(self, alpha, beta, C: list, A: list, B: list):
# scop begin
        for i in range(0, self.N):
            for j in range(0, i + 1):
                C[self.N * i + j] *= beta

            for k in range(0, self.M):
                for j in range(0, i + 1):
                    C[self.N * i + j] += A[self.M * j + k] * alpha * B[self.M * i + k] + B[self.M * j + k] * alpha * A[self.M * i + k]
# scop end


class _StrategyNumPy(Syr2k):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, alpha, beta, C: list, A: list, B: list):
        for i in range(0, self.N):
            for j in range(0, self.M):
                A[i, j] = self.DATA_TYPE((i * j + 1) % self.N) / self.N
                B[i, j] = self.DATA_TYPE((i * j + 2) % self.M) / self.M

        for i in range(0, self.N):
            for j in range(0, self.N):
                C[i, j] = self.DATA_TYPE((i * j + 3) % self.N) / self.M

    def print_array_custom(self, C: ndarray, name: str):
        for i in range(0, self.N):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(C[i, j])

    def kernel(self, alpha, beta, C: ndarray, A: ndarray, B: ndarray):
        # BLAS PARAMS
        # UPLO  = 'L'
        # TRANS = 'N'
        # A is NxM
        # B is NxM
        # C is NxN
# scop begin
        for i in range(0, self.N):
            C[i,0:i+1] *= beta
            C[i,0:i+1] += (A[0:i+1,0:self.M] * alpha * B[i,0:self.M]).sum(axis=1) + (B[0:i+1,0:self.M] * alpha * A[i,0:self.M]).sum(axis=1)
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        self.kernel_vectorizer = self.kernel_pluto
        self.kernel = getattr( self, "kernel_%s" % (options.POCC) )

    def kernel_pluto(self, alpha, beta, C: list, A: list, B: list):
# scop begin
        if((self.N-1>= 0)):
            for c1 in range ((self.N-1)+1):
                for c2 in range ((c1)+1):
                    C[self.N*(c1) + c2] *= beta
            if((self.M-1>= 0)):
                for c1 in range ((self.N-1)+1):
                    for c2 in range ((c1)+1):
                        for c3 in range ((self.M-1)+1):
                            C[self.N*(c1) + c2] += A[self.M*(c2) + c3]*alpha*B[self.M*(c1) + c3] + B[self.M*(c2) + c3]*alpha*A[self.M*(c1) + c3]
# scop end

    def kernel_maxfuse(self, alpha, beta, C: list, A: list, B: list):
# --pluto --pluto-fuse maxfuse
        if((self.N-1>= 0)):
            if((self.M-1>= 0)):
                for c0 in range ((self.N-1)+1):
                    for c1 in range ((c0)+1):
                        C[(c0)*self.N + c1] *= beta
                        for c2 in range (c0 , (self.M + c0-1)+1):
                            C[(c0)*self.N + c1] += A[(c1)*self.M + (-1 * c0) + c2]*alpha*B[(c0)*self.M + (-1 * c0) + c2] + B[(c1)*self.M + (-1 * c0) + c2]*alpha*A[(c0)*self.M + (-1 * c0) + c2]
            if((self.M*-1>= 0)):
                for c0 in range ((self.N-1)+1):
                    for c1 in range ((c0)+1):
                        C[(c0)*self.N + c1] *= beta
# scop end
