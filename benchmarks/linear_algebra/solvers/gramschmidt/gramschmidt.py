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


class Gramschmidt(PolyBench):

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
        A = self.create_array(2, [self.M, self.N], self.DATA_TYPE(0))
        R = self.create_array(2, [self.N, self.N], self.DATA_TYPE(0))
        Q = self.create_array(2, [self.M, self.N], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(A, R, Q)

        # Benchmark the kernel
        self.time_kernel(A, R, Q)

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
        return [('R', R), ('Q', Q)]


class _StrategyList(Gramschmidt):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: list, R: list, Q: list):
        for i in range(0, self.M):
            for j in range(0, self.N):
                A[i][j] = ((self.DATA_TYPE((i * j) % self.M) / self.M) * 100) + 10
                Q[i][j] = 0.0

        for i in range(0, self.N):
            for j in range(0, self.N):
                R[i][j] = 0.0

    def print_array_custom(self, array: list, name: str):
        if name == 'R':
            loop_bound = self.N
        else:
            loop_bound = self.M

        for i in range(0, loop_bound):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(array[i][j])

    def kernel(self, A: list, R: list, Q: list):
# scop begin
        for k in range(0, self.N):
            nrm = 0.0
            for i in range(0, self.M):
                nrm += A[i][k] * A[i][k]
            R[k][k] = math.sqrt(nrm)

            for i in range(0, self.M):
                Q[i][k] = A[i][k] / R[k][k]

            for j in range(k + 1, self.N):
                R[k][j] = 0.0
                for i in range(0, self.M):
                    R[k][j] += Q[i][k] * A[i][j]

                for i in range(0, self.M):
                    A[i][j] = A[i][j] - Q[i][k] * R[k][j]
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, A: list, R: list, Q: list):
# scop begin
        if((self.N-1>= 0)):
            for c1 in range ((self.N-2)+1):
                for c3 in range (c1 + 1 , (self.N-1)+1):
                    R[c1][c3] = 0.0
            if((self.M-1>= 0)):
                for c1 in range ((self.N-2)+1):
                    nrm = 0.0
                    for c3 in range ((self.M-1)+1):
                        nrm += A[c3][c1] * A[c3][c1]
                    R[c1][c1] = math.sqrt(nrm)
                    for c3 in range ((self.M-1)+1):
                        Q[c3][c1] = A[c3][c1] / R[c1][c1]
                    for c3 in range (c1 + 1 , (self.N-1)+1):
                        for c6 in range ((self.M-1)+1):
                            R[c1][c3] += Q[c6][c1] * A[c6][c3]
                        for c6 in range ((self.M-1)+1):
                            A[c6][c3] = A[c6][c3] - Q[c6][c1] * R[c1][c3]
            if((self.M-1>= 0)):
                nrm = 0.0
            if((self.M-1>= 0)):
                for c3 in range ((self.M-1)+1):
                    nrm += A[c3][self.N + -1] * A[c3][self.N + -1]
            if((self.M-1>= 0)):
                R[self.N + -1][self.N + -1] = math.sqrt(nrm)
            if((self.M-1>= 0)):
                for c3 in range ((self.M-1)+1):
                    Q[c3][self.N + -1] = A[c3][self.N + -1] / R[self.N + -1][self.N + -1]
            if((self.M*-1>= 0)):
                for c1 in range ((self.N-1)+1):
                    nrm = 0.0
                    R[c1][c1] = math.sqrt(nrm)
# scop end

class _StrategyListFlattened(Gramschmidt):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: list, R: list, Q: list):
        for i in range(0, self.M):
            for j in range(0, self.N):
                A[self.N * i + j] = ((self.DATA_TYPE((i * j) % self.M) / self.M) * 100) + 10
                Q[self.N * i + j] = 0.0

        for i in range(0, self.N):
            for j in range(0, self.N):
                R[self.N * i + j] = 0.0

    def print_array_custom(self, array: list, name: str):
        if name == 'R':
            loop_bound = self.N
        else:
            loop_bound = self.M

        for i in range(0, loop_bound):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(array[self.N * i + j])

    def kernel(self, A: list, R: list, Q: list):
# scop begin
        for k in range(0, self.N):
            nrm = 0.0
            for i in range(0, self.M):
                nrm += A[self.N * i + k] * A[self.N * i + k]
            R[self.N * k + k] = math.sqrt(nrm)

            for i in range(0, self.M):
                Q[self.N * i + k] = A[self.N * i + k] / R[self.N * k + k]

            for j in range(k + 1, self.N):
                R[self.N * k + j] = 0.0
                for i in range(0, self.M):
                    R[self.N * k + j] += Q[self.N * i + k] * A[self.N * i + j]

                for i in range(0, self.M):
                    A[self.N * i + j] = A[self.N * i + j] - Q[self.N * i + k] * R[self.N * k + j]
# scop end


class _StrategyNumPy(Gramschmidt):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, A: ndarray, R: ndarray, Q: ndarray):
        for i in range(0, self.M):
            for j in range(0, self.N):
                A[i, j] = ((self.DATA_TYPE((i * j) % self.M) / self.M) * 100) + 10
                Q[i, j] = 0.0

        for i in range(0, self.N):
            for j in range(0, self.N):
                R[i, j] = 0.0

    def print_array_custom(self, array: ndarray, name: str):
        if name == 'R':
            loop_bound = self.N
        else:
            loop_bound = self.M

        for i in range(0, loop_bound):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(array[i, j])

    def kernel(self, A: ndarray, R: ndarray, Q: ndarray):
# scop begin
        for k in range(0, self.N):
            nrm = np.dot( A[0:self.M,k], A[0:self.M,k] )
            R[k, k] = math.sqrt(nrm)
            Q[0:self.M,k] = A[0:self.M,k] / R[k,k]
            R[k, k+1:self.N] = np.dot( Q[0:self.M, k].T, A[0:self.M,k+1:self.N] )
            A[0:self.M,k+1:self.N] = A[0:self.M,k+1:self.N] - np.dot(Q[0:self.M,k,np.newaxis], R[np.newaxis,k,k+1:self.N])
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        self.kernel = getattr( self, "kernel_%s" % (options.POCC) )

    def kernel_pluto(self, A: list, R: list, Q: list):
# --pluto
# scop begin
        if((self.N-1>= 0)):
            for c1 in range ((self.N-2)+1):
                for c3 in range (c1 + 1 , (self.N-1)+1):
                    R[self.N*(c1) + c3] = 0.0
            if((self.M-1>= 0)):
                for c1 in range ((self.N-2)+1):
                    nrm = 0.0
                    for c3 in range ((self.M-1)+1):
                        nrm += A[self.N*(c3) + c1] * A[self.N*(c3) + c1]
                    R[self.N*(c1) + c1] = math.sqrt(nrm)
                    for c3 in range ((self.M-1)+1):
                        Q[self.N*(c3) + c1] = A[self.N*(c3) + c1] / R[self.N*(c1) + c1]
                    for c3 in range (c1 + 1 , (self.N-1)+1):
                        for c6 in range ((self.M-1)+1):
                            R[self.N*(c1) + c3] += Q[self.N*(c6) + c1] * A[self.N*(c6) + c3]
                        for c6 in range ((self.M-1)+1):
                            A[self.N*(c6) + c3] = A[self.N*(c6) + c3] - Q[self.N*(c6) + c1] * R[self.N*(c1) + c3]
            if((self.M-1>= 0)):
                nrm = 0.0
            if((self.M-1>= 0)):
                for c3 in range ((self.M-1)+1):
                    nrm += A[self.N*(c3) + self.N + -1] * A[self.N*(c3) + self.N + -1]
            if((self.M-1>= 0)):
                R[self.N*(self.N + -1) + self.N + -1] = math.sqrt(nrm)
            if((self.M-1>= 0)):
                for c3 in range ((self.M-1)+1):
                    Q[self.N*(c3) + self.N + -1] = A[self.N*(c3) + self.N + -1] / R[self.N*(self.N + -1) + self.N + -1]
            if((self.M*-1>= 0)):
                for c1 in range ((self.N-1)+1):
                    nrm = 0.0
                    R[self.N*(c1) + c1] = math.sqrt(nrm)
# scop end

    def kernel_vectorizer(self, A: list, R: list, Q: list):
# --pluto --pluto-prevector --vectorizer --pragmatizer
# scop begin
        if (self.N >= 1):
          ub1 = (self.N + -2);
          for c1 in range( ub1+1 ):
              for c3 in range( c1+1, self.N ):
                  R[(c1)*self.N + c3] = 0.0;
          if (self.M >= 1):
            for c1 in range( self.N-1 ):
                nrm = 0.0;
                for c3 in range( self.M ):
                    nrm += A[(c3)*self.N + c1] * A[(c3)*self.N + c1];
                R[(c1)*self.N + c1] = math.sqrt(nrm);
                for c3 in range( self.M ):
                    Q[(c3)*self.N + c1] = A[(c3)*self.N + c1] / R[(c1)*self.N + c1];
                lb1 = (c1 + 1);
                ub1 = (self.N + -1);
                for c3 in range( lb1, ub1+1 ):
                    for c6 in range( self.M ):
                        R[(c1)*self.N + c3] += Q[(c6)*self.N + c1] * A[(c6)*self.N + c3];
                    for c6 in range( self.M ):
                        A[(c6)*self.N + c3] = A[(c6)*self.N + c3] - Q[(c6)*self.N + c1] * R[(c1)*self.N + c3];
          if (self.M >= 1):
            nrm = 0.0;
            for c3 in range( self.M ):
                nrm += A[(c3)*self.N + self.N + -1] * A[(c3)*self.N + self.N + -1];
            R[(self.N + -1)*self.N + self.N + -1] = math.sqrt(nrm);
            for c3 in range( self.M ):
                Q[(c3)*self.N + self.N + -1] = A[(c3)*self.N + self.N + -1] / R[(self.N + -1)*self.N + self.N + -1];
          if (self.M <= 0):
            for c1 in range( self.N ):
                nrm = 0.0;
                R[(c1)*self.N + c1] = math.sqrt(nrm);
# scop end

    def kernel_maxfuse(self, A: list, R: list, Q: list):
# --pluto --pluto-fuse maxfuse
# scop begin
        if((self.N-1>= 0)):
            if((self.M-1>= 0)):
                for c0 in range ((self.N-2)+1):
                    for c3 in range (c0 + 1 , (self.N-1)+1):
                        R[(c0)*self.N + c3] = 0.0
                    nrm = 0.0
                    for c3 in range ((self.M-1)+1):
                        nrm += A[(c3)*self.N + c0] * A[(c3)*self.N + c0]
                    R[(c0)*self.N + c0] = math.sqrt(nrm)
                    for c3 in range ((self.M-1)+1):
                        Q[(c3)*self.N + c0] = A[(c3)*self.N + c0] / R[(c0)*self.N + c0]
                    for c3 in range (c0 + 1 , (self.N-1)+1):
                        for c6 in range ((self.M-1)+1):
                            R[(c0)*self.N + c3] += Q[(c6)*self.N + c0] * A[(c6)*self.N + c3]
                        for c6 in range ((self.M-1)+1):
                            A[(c6)*self.N + c3] = A[(c6)*self.N + c3] - Q[(c6)*self.N + c0] * R[(c0)*self.N + c3]
            if((self.M-1>= 0)):
                nrm = 0.0
            if((self.M-1>= 0)):
                for c3 in range ((self.M-1)+1):
                    nrm += A[(c3)*self.N + self.N + -1] * A[(c3)*self.N + self.N + -1]
            if((self.M-1>= 0)):
                R[(self.N + -1)*self.N + self.N + -1] = math.sqrt(nrm)
            if((self.M-1>= 0)):
                for c3 in range ((self.M-1)+1):
                    Q[(c3)*self.N + self.N + -1] = A[(c3)*self.N + self.N + -1] / R[(self.N + -1)*self.N + self.N + -1]
            if((self.M*-1>= 0)):
                for c0 in range ((self.N-2)+1):
                    for c3 in range (c0 + 1 , (self.N-1)+1):
                        R[(c0)*self.N + c3] = 0.0
                    nrm = 0.0
                    R[(c0)*self.N + c0] = math.sqrt(nrm)
            if((self.M*-1>= 0)):
                nrm = 0.0
            if((self.M*-1>= 0)):
                R[(self.N + -1)*self.N + self.N + -1] = math.sqrt(nrm)
# scop end
