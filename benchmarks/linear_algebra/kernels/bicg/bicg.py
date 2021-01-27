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


class Bicg(PolyBench):

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

    def print_array_custom(self, array: list, name: list):
        if name == 's':
            loop_bound = self.M
        else:
            loop_bound = self.N

        for i in range(0, loop_bound):
            if i % 20 == 0:
                self.print_message('\n')
            self.print_value(array[i])

    def run_benchmark(self):
        # Create data structures (arrays, auxiliary variables, etc.)
        A = self.create_array(2, [self.N, self.M], self.DATA_TYPE(0))
        s = self.create_array(1, [self.M], self.DATA_TYPE(0))
        q = self.create_array(1, [self.N], self.DATA_TYPE(0))
        p = self.create_array(1, [self.M], self.DATA_TYPE(0))
        r = self.create_array(1, [self.N], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(A, s, q, p, r)

        # Benchmark the kernel
        self.time_kernel(A, s, q, p, r)

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
        return [('s', s), ('q', q)]


class _StrategyList(Bicg):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: list, s: list, q: list, p: list, r: list):
        for i in range(0, self.M):
            p[i] = self.DATA_TYPE(i % self.M) / self.M

        for i in range(0, self.N):
            r[i] = self.DATA_TYPE(i % self.N) / self.N
            for j in range(0, self.M):
                A[i][j] = self.DATA_TYPE(i * (j+1) % self.N) / self.N

    def kernel(self, A: list, s: list, q: list, p: list, r: list):
# scop begin
        for i in range(0, self.M):
            s[i] = 0

        for i in range(0, self.N):
            q[i] = 0.0
            for j in range(0, self.M):
                s[j] = s[j] + r[i] * A[i][j]
                q[i] = q[i] + A[i][j] * p[j]
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, A: list, s: list, q: list, p: list, r: list):
# scop begin
        for c1 in range ((self.N-1)+1):
            q[c1] = 0.0
        if((self.M-1>= 0)):
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.M-1)+1):
                    q[c1] = q[c1] + A[c1][c2] * p[c2]
        for c1 in range ((self.M-1)+1):
            s[c1] = 0
        if((self.N-1>= 0)):
            for c1 in range ((self.M-1)+1):
                for c2 in range ((self.N-1)+1):
                    s[c1] = s[c1] + r[c2] * A[c2][c1]
# scop end

class _StrategyListFlattened(Bicg):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: list, s: list, q: list, p: list, r: list):
        for i in range(0, self.M):
            p[i] = self.DATA_TYPE(i % self.M) / self.M

        for i in range(0, self.N):
            r[i] = self.DATA_TYPE(i % self.N) / self.N
            for j in range(0, self.M):
                A[self.M * i + j] = self.DATA_TYPE(i * (j + 1) % self.N) / self.N

    def kernel(self, A: list, s: list, q: list, p: list, r: list):
# scop begin
        for i in range(0, self.M):
            s[i] = 0

        for i in range(0, self.N):
            q[i] = 0.0
            for j in range(0, self.M):
                s[j] = s[j] + r[i] * A[self.M * i + j]
                q[i] = q[i] + A[self.M * i + j] * p[j]
# scop end


class _StrategyNumPy(Bicg):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: list, s: list, q: list, p: list, r: list):
        for i in range(0, self.M):
            p[i] = self.DATA_TYPE(i % self.M) / self.M

        for i in range(0, self.N):
            r[i] = self.DATA_TYPE(i % self.N) / self.N
            for j in range(0, self.M):
                A[i, j] = self.DATA_TYPE(i * (j + 1) % self.N) / self.N

    def kernel(self, A: ndarray, s: ndarray, q: ndarray, p: ndarray, r: ndarray):
# scop begin
        s[0:self.M] = np.dot( r[0:self.N], A[0:self.N,0:self.M] )
        q[0:self.N] = np.dot( A[0:self.N,0:self.M], p[0:self.M] )
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        self.kernel = getattr( self, "kernel_%s" % (options.POCC) )

    def kernel_pluto(self, A: list, s: list, q: list, p: list, r: list):
# --pluto
# scop begin
        for c1 in range ((self.N-1)+1):
            q[c1] = 0.0
        if((self.M-1>= 0)):
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.M-1)+1):
                    q[c1] = q[c1] + A[self.M*(c1) + c2] * p[c2]
        for c1 in range ((self.M-1)+1):
            s[c1] = 0
        if((self.N-1>= 0)):
            for c1 in range ((self.M-1)+1):
                for c2 in range ((self.N-1)+1):
                    s[c1] = s[c1] + r[c2] * A[self.M*(c2) + c1]
# scop end

    def kernel_vectorizer(self, A: list, s: list, q: list, p: list, r: list):
# --pluto --pluto-prevector --pragmatizer --vectorizer
# scop begin
        for c1 in range ((self.N-1)+1):
            q[c1] = 0.0
        if((self.M-1>= 0)):
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.M-1)+1):
                    q[c1] = q[c1] + A[self.M*(c1) + c2] * p[c2]
        for c1 in range ((self.M-1)+1):
            s[c1] = 0
        if((self.N-1>= 0)):
            for c2 in range ((self.N-1)+1):
                for c1 in range ((self.M-1)+1):
                    s[c1] = s[c1] + r[c2] * A[self.M*(c2) + c1]
# scop end

    def kernel_maxfuse(self, A: list, s: list, q: list, p: list, r: list):
# --pluto --pluto-fuse maxfuse
# scop begin
        for c0 in range (min((self.M-1)+1 , (self.N-1)+1)):
            q[c0] = 0.0
            q[c0] = q[c0] + A[(c0)*self.M + 0] * p[0]
            s[c0] = 0
            s[c0] = s[c0] + r[0] * A[(0)*self.M + c0]
            for c1 in range (c0 + 1 , min((self.M + c0-1)+1 , (self.N + c0-1)+1)):
                q[c0] = q[c0] + A[(c0)*self.M + (-1 * c0) + c1] * p[(-1 * c0) + c1]
                s[c0] = s[c0] + r[(-1 * c0) + c1] * A[((-1 * c0) + c1)*self.M + c0]
            for c1 in range (self.M + c0 , (self.N + c0-1)+1):
                s[c0] = s[c0] + r[(-1 * c0) + c1] * A[((-1 * c0) + c1)*self.M + c0]
            for c1 in range (self.N + c0 , (self.M + c0-1)+1):
                q[c0] = q[c0] + A[(c0)*self.M + (-1 * c0) + c1] * p[(-1 * c0) + c1]
        if((self.N-1>= 0)):
            for c0 in range (self.N , (self.M-1)+1):
                s[c0] = 0
                for c1 in range (c0 , (self.N + c0-1)+1):
                    s[c0] = s[c0] + r[(-1 * c0) + c1] * A[((-1 * c0) + c1)*self.M + c0]
        if((self.N*-1>= 0)):
            for c0 in range ((self.M-1)+1):
                s[c0] = 0
        if((self.M-1>= 0)):
            for c0 in range (self.M , (self.N-1)+1):
                q[c0] = 0.0
                for c1 in range (c0 , (self.M + c0-1)+1):
                    q[c0] = q[c0] + A[(c0)*self.M + (-1 * c0) + c1] * p[(-1 * c0) + c1]
        if((self.M*-1>= 0)):
            for c0 in range ((self.N-1)+1):
                q[c0] = 0.0
# scop end
