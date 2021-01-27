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

class Gemver(PolyBench):

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

    def print_array_custom(self, w: list, name: str):
        for i in range(0, self.N):
            if i % 20 == 0:
                self.print_message('\n')
            self.print_value(w[i])

    def run_benchmark(self):
        # Create data structures (arrays, auxiliary variables, etc.)
        alpha = self.DATA_TYPE(1.5)
        beta = self.DATA_TYPE(1.2)

        A = self.create_array(2, [self.N, self.N], self.DATA_TYPE(0))
        u1 = self.create_array(1, [self.N], self.DATA_TYPE(0))
        v1 = self.create_array(1, [self.N], self.DATA_TYPE(0))
        u2 = self.create_array(1, [self.N], self.DATA_TYPE(0))
        v2 = self.create_array(1, [self.N], self.DATA_TYPE(0))
        w = self.create_array(1, [self.N], self.DATA_TYPE(0))
        x = self.create_array(1, [self.N], self.DATA_TYPE(0))
        y = self.create_array(1, [self.N], self.DATA_TYPE(0))
        z = self.create_array(1, [self.N], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(alpha, beta, A, u1, v1, u2, v2, w, x, y, z)

        # Benchmark the kernel
        self.time_kernel(alpha, beta, A, u1, v1, u2, v2, w, x, y, z)

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
        return [('w', w)]


class _StrategyList(Gemver):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, alpha, beta, A: list, u1: list, v1: list, u2: list, v2: list, w: list, x: list, y: list, z: list):
        fn = self.DATA_TYPE(self.N)

        for i in range(0, self.N):
            u1[i] = i
            u2[i] = ((i + 1) / fn) / 2.0
            v1[i] = ((i + 1) / fn) / 4.0
            v2[i] = ((i + 1) / fn) / 6.0
            y[i] = ((i + 1) / fn) / 8.0
            z[i] = ((i + 1) / fn) / 9.0
            x[i] = 0.0
            w[i] = 0.0
            for j in range(0, self.N):
                A[i][j] = self.DATA_TYPE(i * j % self.N) / self.N

    def kernel(self, alpha, beta, A: list, u1: list, v1: list, u2: list, v2: list, w: list, x: list, y: list, z: list):
# scop begin
        for i in range(0, self.N):
            for j in range(0, self.N):
                A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j]

        for i in range(0, self.N):
            for j in range(0, self.N):
                x[i] = x[i] + beta * A[j][i] * y[j]

        for i in range(0, self.N):
            x[i] = x[i] + z[i]

        for i in range(0, self.N):
            for j in range(0, self.N):
                w[i] = w[i] + alpha * A[i][j] * x[j]
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, alpha, beta, A: list, u1: list, v1: list, u2: list, v2: list, w: list, x: list, y: list, z: list):
# scop begin
        if((self.N-1>= 0)):
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.N-1)+1):
                    A[c2][c1] = A[c2][c1] + u1[c2] * v1[c1] + u2[c2] * v2[c1]
                    x[c1] = x[c1] + beta * A[c2][c1] * y[c2]
            for c1 in range ((self.N-1)+1):
                x[c1] = x[c1] + z[c1]
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.N-1)+1):
                    w[c1] = w[c1] + alpha * A[c1][c2] * x[c2]
# scop end

class _StrategyListFlattened(Gemver):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, alpha, beta, A: list, u1: list, v1: list, u2: list, v2: list, w: list, x: list, y: list, z: list):
        fn = self.DATA_TYPE(self.N)

        for i in range(0, self.N):
            u1[i] = i
            u2[i] = ((i + 1) / fn) / 2.0
            v1[i] = ((i + 1) / fn) / 4.0
            v2[i] = ((i + 1) / fn) / 6.0
            y[i] = ((i + 1) / fn) / 8.0
            z[i] = ((i + 1) / fn) / 9.0
            x[i] = 0.0
            w[i] = 0.0
            for j in range(0, self.N):
                A[self.N * i + j] = self.DATA_TYPE(i * j % self.N) / self.N

    def kernel(self, alpha, beta, A: list, u1: list, v1: list, u2: list, v2: list, w: list, x: list, y: list, z: list):
# scop begin
        for i in range(0, self.N):
            for j in range(0, self.N):
                A[self.N * i + j] = A[self.N * i + j] + u1[i] * v1[j] + u2[i] * v2[j]

        for i in range(0, self.N):
            for j in range(0, self.N):
                x[i] = x[i] + beta * A[self.N * j + i] * y[j]

        for i in range(0, self.N):
            x[i] = x[i] + z[i]

        for i in range(0, self.N):
            for j in range(0, self.N):
                w[i] = w[i] + alpha * A[self.N * i + j] * x[j]
# scop end


class _StrategyNumPy(Gemver):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, alpha, beta, A: list, u1: list, v1: list, u2: list, v2: list, w: list, x: list, y: list, z: list):
        fn = self.DATA_TYPE(self.N)

        for i in range(0, self.N):
            u1[i] = i
            u2[i] = ((i + 1) / fn) / 2.0
            v1[i] = ((i + 1) / fn) / 4.0
            v2[i] = ((i + 1) / fn) / 6.0
            y[i] = ((i + 1) / fn) / 8.0
            z[i] = ((i + 1) / fn) / 9.0
            x[i] = 0.0
            w[i] = 0.0
            for j in range(0, self.N):
                A[i, j] = self.DATA_TYPE(i * j % self.N) / self.N

    def kernel(self, alpha, beta, A: ndarray, u1: ndarray, v1: ndarray, u2: ndarray, v2: ndarray,
               w: ndarray, x: ndarray, y: ndarray, z: ndarray):
# scop begin
        A[0:self.N,0:self.N] = A[0:self.N,0:self.N] + np.dot( u1[:,np.newaxis], v1[np.newaxis,:] ) + np.dot( u2[:,np.newaxis], v2[np.newaxis,:] )

        x[0:self.N] = x[0:self.N] + beta * np.dot( A[0:self.N,0:self.N].T, y[0:self.N] )

        x[0:self.N] = x[0:self.N] + z[0:self.N]

        w[0:self.N] = w[0:self.N] + alpha * np.dot( A[0:self.N,0:self.N], x[0:self.N] )
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        self.kernel = getattr( self, "kernel_%s" % (options.POCC) )

    def kernel_pluto(self, alpha, beta, A: list, u1: list, v1: list, u2: list, v2: list, w: list, x: list, y: list, z: list):
# scop begin
        if((self.N-1>= 0)):
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.N-1)+1):
                    A[self.N*(c2) + c1] = A[self.N*(c2) + c1] + u1[c2] * v1[c1] + u2[c2] * v2[c1]
                    x[c1] = x[c1] + beta * A[self.N*(c2) + c1] * y[c2]
            for c1 in range ((self.N-1)+1):
                x[c1] = x[c1] + z[c1]
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.N-1)+1):
                    w[c1] = w[c1] + alpha * A[self.N*(c1) + c2] * x[c2]
# scop end

    def kernel_vectorizer(self, alpha, beta, A: list, u1: list, v1: list, u2: list, v2: list, w: list, x: list, y: list, z: list):
# --pluto --pluto-prevector --vectorizer --pragmatizer
# scop begin
        if((self.N-1>= 0)):
            for c2 in range ((self.N-1)+1):
                for c1 in range ((self.N-1)+1):
                    A[self.N*(c2) + c1] = A[self.N*(c2) + c1] + u1[c2] * v1[c1] + u2[c2] * v2[c1]
                    x[c1] = x[c1] + beta * A[self.N*(c2) + c1] * y[c2]
            for c1 in range ((self.N-1)+1):
                x[c1] = x[c1] + z[c1]
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.N-1)+1):
                    w[c1] = w[c1] + alpha * A[self.N*(c1) + c2] * x[c2]
# scop end

    def kernel_maxfuse(self, alpha, beta, A: list, u1: list, v1: list, u2: list, v2: list, w: list, x: list, y: list, z: list):
# --pluto --pluto-fuse maxfuse
# scop begin
        if((self.N-1>= 0)):
            for c0 in range ((self.N-1)+1):
                for c3 in range ((self.N-1)+1):
                    A[(c3)*self.N + c0] = A[(c3)*self.N + c0] + u1[c3] * v1[c0] + u2[c3] * v2[c0]
                    x[c0] = x[c0] + beta * A[(c3)*self.N + c0] * y[c3]
                x[c0] = x[c0] + z[c0]
                for c3 in range ((self.N-1)+1):
                    w[c3] = w[c3] + alpha * A[(c3)*self.N + c0] * x[c0]
# scop end
