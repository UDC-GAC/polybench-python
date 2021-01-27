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
import math
import numpy as np


class Cholesky(PolyBench):

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


class _StrategyList(Cholesky):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: list):
        for i in range(0, self.N):
            for j in range(0, i + 1):
                A[i][j] = -self.DATA_TYPE(j % self.N) / self.N + 1

            for j in range(i + 1, self.N):
                A[i][j] = self.DATA_TYPE(0)
            A[i][i] = self.DATA_TYPE(1)

        # Make the matrix positive semi-definite.
        B = self.create_array(2, [self.N], self.DATA_TYPE(0))

        for t in range(0, self.N):
            for r in range(0, self.N):
                for s in range(0, self.N):
                    B[r][s] += A[r][t] * A[s][t]

        for r in range(0, self.N):
            for s in range(0, self.N):
                A[r][s] = B[r][s]

    def print_array_custom(self, A: list, name: str):
        for i in range(0, self.N):
            for j in range(0, i + 1):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(A[i][j])

    def kernel(self, A: list):
# scop begin
        for i in range(0, self.N):
            # j < i
            for j in range(0, i):
                for k in range(0, j):
                    A[i][j] -= A[i][k] * A[j][k]
                A[i][j] /= A[j][j]

            # i == j case
            for k in range(0, i):
                A[i][i] -= A[i][k] * A[i][k]

            A[i][i] = math.sqrt(A[i][i])
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, A: list):
# scop begin
        if((self.N-1>= 0)):
            A[0][0] = math.sqrt(A[0][0])
            if((self.N-2>= 0)):
                A[1][0] /= A[0][0]
            if((self.N-2>= 0)):
                A[1][1] -= A[1][0] * A[1][0]
            if((self.N-2>= 0)):
                A[1][1] = math.sqrt(A[1][1])
            for c0 in range (2 , (self.N-1)+1):
                A[c0][0] /= A[0][0]
                A[c0][c0] -= A[c0][0] * A[c0][0]
                for c1 in range (1 , (c0-1)+1):
                    for c2 in range ((c1-1)+1):
                        A[c0][c1] -= A[c0][c2] * A[c1][c2]
                    A[c0][c1] /= A[c1][c1]
                    A[c0][c0] -= A[c0][c1] * A[c0][c1]
                A[c0][c0] = math.sqrt(A[c0][c0])
# scop end

class _StrategyListFlattened(Cholesky):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        if options.LOAD_ELIMINATION: self.kernel = self.kernel_le
        else: self.kernel = self.kernel_regular

    def initialize_array(self, A: list):
        for i in range(0, self.N):
            for j in range(0, i + 1):
                A[self.N * i + j] = -self.DATA_TYPE(j % self.N) / self.N + 1

            for j in range(i + 1, self.N):
                A[self.N * i + j] = self.DATA_TYPE(0)
            A[self.N * i + i] = self.DATA_TYPE(1)

        # Make the matrix positive semi-definite.
        B = self.create_array(1, [self.N * self.N], self.DATA_TYPE(0))

        for t in range(0, self.N):
            for r in range(0, self.N):
                for s in range(0, self.N):
                    B[self.N * r + s] += A[self.N * r + t] * A[self.N * s + t]

        for r in range(0, self.N):
            for s in range(0, self.N):
                A[self.N * r + s] = B[self.N * r + s]

    def print_array_custom(self, A: list, name: str):
        for i in range(0, self.N):
            for j in range(0, i + 1):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(A[self.N * i + j])

    def kernel_regular(self, A: list):
# scop begin
        for i in range(0, self.N):
            # j < i
            for j in range(0, i):
                for k in range(0, j):
                    A[self.N * i + j] -= A[self.N * i + k] * A[self.N * j + k]
                A[self.N * i + j] /= A[self.N * j + j]

            # i == j case
            for k in range(0, i):
                A[self.N * i + i] -= A[self.N * i + k] * A[self.N * i + k]

            A[self.N * i + i] = math.sqrt(A[self.N * i + i])
# scop end

    def kernel_le(self, A: list):
# scop begin
        for i in range(0, self.N):
            # j < i
            for j in range(0, i):
                tmp = A[self.N*i+j] # load elimination
                for k in range(0, j):
                    tmp -= A[self.N * i + k] * A[self.N * j + k] # load elimination
                A[self.N * i + j] = tmp / A[self.N * j + j] # load elimination

            # i == j case
            tmp = A[self.N*i+i] # load elimination
            for k in range(0, i):
                tmp -= A[self.N * i + k] * A[self.N * i + k] # load elimination

            A[self.N * i + i] = math.sqrt(tmp) # load elimination
# scop end

class _StrategyNumPy(Cholesky):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: ndarray):
        for i in range(0, self.N):
            for j in range(0, i + 1):
                A[i, j] = -self.DATA_TYPE(j % self.N) / self.N + 1

            for j in range(i + 1, self.N):
                A[i, j] = self.DATA_TYPE(0)
            A[i, i] = self.DATA_TYPE(1)

        # Make the matrix positive semi-definite.
        B = self.create_array(2, [self.N], self.DATA_TYPE(0))

        for t in range(0, self.N):
            for r in range(0, self.N):
                B[r,0:self.N] += A[r,t] * A[0:self.N,t]
#                for s in range(0, self.N):
#                    B[r, s] += A[r, t] * A[s, t]

        for r in range(0, self.N):
            for s in range(0, self.N):
                A[r, s] = B[r, s]

    def print_array_custom(self, A: ndarray, name: str):
        for i in range(0, self.N):
            for j in range(0, i + 1):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(A[i, j])

    def kernel(self, A: ndarray):
# scop begin
        for i in range(0, self.N):
            # j < i
            for j in range(0, i):
                A[i,j] -= np.dot( A[i,0:j], A[j, 0:j] )
                A[i, j] /= A[j, j]

            # i == j case
            A[i,i] -= np.dot( A[i,0:i], A[i,0:i] )
        A[np.diag_indices(self.N)] = np.sqrt( np.diag(A[0:self.N,0:self.N]) )
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        self.kernel_vectorizer = self.kernel_pluto
        self.kernel = getattr( self, "kernel_%s" % (options.POCC) )

    def kernel_pluto(self, A: list):
# --pluto
# scop begin
        if((self.N-1>= 0)):
            A[self.N*(0) + 0] = math.sqrt(A[self.N*(0) + 0])
            if((self.N-2>= 0)):
                A[self.N*(1) + 0] /= A[self.N*(0) + 0]
            if((self.N-2>= 0)):
                A[self.N*(1) + 1] -= A[self.N*(1) + 0] * A[self.N*(1) + 0]
            if((self.N-2>= 0)):
                A[self.N*(1) + 1] = math.sqrt(A[self.N*(1) + 1])
            for c0 in range (2 , (self.N-1)+1):
                A[self.N*(c0) + 0] /= A[self.N*(0) + 0]
                A[self.N*(c0) + c0] -= A[self.N*(c0) + 0] * A[self.N*(c0) + 0]
                for c1 in range (1 , (c0-1)+1):
                    for c2 in range ((c1-1)+1):
                        A[self.N*(c0) + c1] -= A[self.N*(c0) + c2] * A[self.N*(c1) + c2]
                    A[self.N*(c0) + c1] /= A[self.N*(c1) + c1]
                    A[self.N*(c0) + c0] -= A[self.N*(c0) + c1] * A[self.N*(c0) + c1]
                A[self.N*(c0) + c0] = math.sqrt(A[self.N*(c0) + c0])

    def kernel_maxfuse(self, A: list):
# --pluto --pluto-fuse maxfuse
# scop begin
        if((self.N-1>= 0)):
            A[(0)*self.N + 0] = math.sqrt(A[(0)*self.N + 0])
            if((self.N-2>= 0)):
                A[(1)*self.N + 0] /= A[(0)*self.N + 0]
            if((self.N-2>= 0)):
                A[(1)*self.N + 1] -= A[(1)*self.N + 0] * A[(1)*self.N + 0]
            if((self.N-2>= 0)):
                A[(1)*self.N + 1] = math.sqrt(A[(1)*self.N + 1])
            for c0 in range (2 , (self.N-1)+1):
                A[(c0)*self.N + 0] /= A[(0)*self.N + 0]
                A[(c0)*self.N + c0] -= A[(c0)*self.N + 0] * A[(c0)*self.N + 0]
                for c1 in range (1 , (c0-1)+1):
                    for c2 in range ((c1-1)+1):
                        A[(c0)*self.N + c1] -= A[(c0)*self.N + c2] * A[(c1)*self.N + c2]
                    A[(c0)*self.N + c1] /= A[(c1)*self.N + c1]
                    A[(c0)*self.N + c0] -= A[(c0)*self.N + c1] * A[(c0)*self.N + c1]
                A[(c0)*self.N + c0] = math.sqrt(A[(c0)*self.N + c0])
# scop end
