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


class _3mm(PolyBench):

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
        self.NI = params.get('NI')
        self.NJ = params.get('NJ')
        self.NK = params.get('NK')
        self.NL = params.get('NL')
        self.NM = params.get('NM')

    def run_benchmark(self):
        # Create data structures (arrays, auxiliary variables, etc.)
        E = self.create_array(2, [self.NI, self.NJ], self.DATA_TYPE(0))
        A = self.create_array(2, [self.NI, self.NK], self.DATA_TYPE(0))
        B = self.create_array(2, [self.NK, self.NJ], self.DATA_TYPE(0))
        F = self.create_array(2, [self.NJ, self.NL], self.DATA_TYPE(0))
        C = self.create_array(2, [self.NJ, self.NM], self.DATA_TYPE(0))
        D = self.create_array(2, [self.NM, self.NL], self.DATA_TYPE(0))
        G = self.create_array(2, [self.NI, self.NL], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(E, A, B, F, C, D, G)

        # Benchmark the kernel
        self.time_kernel(E, A, B, F, C, D, G)

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
        return [('G', G)]


class _StrategyList(_3mm):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, E: list, A: list, B: list, F: list, C: list, D: list, G: list):
        for i in range(0, self.NI):
            for j in range(0, self.NK):
                A[i][j] = self.DATA_TYPE((i * j + 1) % self.NI) / (5 * self.NI)

        for i in range(0, self.NK):
            for j in range(0, self.NJ):
                B[i][j] = self.DATA_TYPE((i * (j + 1) + 2) % self.NJ) / (5 * self.NJ)

        for i in range(0, self.NJ):
            for j in range(0, self.NM):
                C[i][j] = self.DATA_TYPE(i * (j + 3) % self.NL) / (5 * self.NL)

        for i in range(0, self.NM):
            for j in range(0, self.NL):
                D[i][j] = self.DATA_TYPE((i * (j + 2) + 2) % self.NK) / (5 * self.NK)

    def print_array_custom(self, G: list, name: str):
        for i in range(0, self.NI):
            for j in range(0, self.NL):
                if (i * self.NI + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(G[i][j])

    def kernel(self, E: list, A: list, B: list, F: list, C: list, D: list, G: list):
# scop begin
        # E := A * B
        for i in range(0, self.NI):
            for j in range(0, self.NJ):
                E[i][j] = 0.0
                for k in range(0, self.NK):
                    E[i][j] += A[i][k] * B[k][j]

        # F := C * D
        for i in range(0, self.NJ):
            for j in range(0, self.NL):
                F[i][j] = 0.0
                for k in range(0, self.NM):
                    F[i][j] += C[i][k] * D[k][j]

        # G := E * F
        for i in range(0, self.NI):
            for j in range(0, self.NL):
                G[i][j] = 0.0
                for k in range(0, self.NJ):
                    G[i][j] += E[i][k] * F[k][j]
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, E: list, A: list, B: list, F: list, C: list, D: list, G: list):
# scop begin
        if((self.NL-1>= 0)):
            for c1 in range (min((self.NI-1)+1 , (self.NJ-1)+1)):
                for c2 in range ((self.NL-1)+1):
                    G[c1][c2] = 0.0
                    F[c1][c2] = 0.0
        if((self.NL-1>= 0)):
            for c1 in range (max(0 , self.NI) , (self.NJ-1)+1):
                for c2 in range ((self.NL-1)+1):
                    F[c1][c2] = 0.0
        if((self.NL-1>= 0)):
            for c1 in range (max(0 , self.NJ) , (self.NI-1)+1):
                for c2 in range ((self.NL-1)+1):
                    G[c1][c2] = 0.0
        if((self.NL-1>= 0) and (self.NM-1>= 0)):
            for c1 in range ((self.NJ-1)+1):
                for c2 in range ((self.NL-1)+1):
                    for c5 in range ((self.NM-1)+1):
                        F[c1][c2] += C[c1][c5] * D[c5][c2]
        if((self.NJ-1>= 0)):
            for c1 in range ((self.NI-1)+1):
                for c2 in range ((self.NJ-1)+1):
                    E[c1][c2] = 0.0
        if((self.NJ-1>= 0) and (self.NK-1>= 0) and (self.NL-1>= 0)):
            for c1 in range ((self.NI-1)+1):
                for c2 in range ((self.NJ-1)+1):
                    for c5 in range ((self.NK-1)+1):
                        E[c1][c2] += A[c1][c5] * B[c5][c2]
                    for c5 in range ((self.NL-1)+1):
                        G[c1][c5] += E[c1][c2] * F[c2][c5]
        if((self.NJ-1>= 0) and (self.NK-1>= 0) and (self.NL*-1>= 0)):
            for c1 in range ((self.NI-1)+1):
                for c2 in range ((self.NJ-1)+1):
                    for c5 in range ((self.NK-1)+1):
                        E[c1][c2] += A[c1][c5] * B[c5][c2]
        if((self.NJ-1>= 0) and (self.NK*-1>= 0) and (self.NL-1>= 0)):
            for c1 in range ((self.NI-1)+1):
                for c2 in range ((self.NJ-1)+1):
                    for c5 in range ((self.NL-1)+1):
                        G[c1][c5] += E[c1][c2] * F[c2][c5]
# scop end

class _StrategyListFlattened(_3mm):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        if options.LOAD_ELIMINATION: self.kernel = self.kernel_le
        else: self.kernel = self.kernel_regular

    def initialize_array(self, E: list, A: list, B: list, F: list, C: list, D: list, G: list):
        for i in range(0, self.NI):
            for j in range(0, self.NK):
                A[self.NK * i + j] = self.DATA_TYPE((i * j + 1) % self.NI) / (5 * self.NI)

        for i in range(0, self.NK):
            for j in range(0, self.NJ):
                B[self.NJ * i + j] = self.DATA_TYPE((i * (j + 1) + 2) % self.NJ) / (5 * self.NJ)

        for i in range(0, self.NJ):
            for j in range(0, self.NM):
                C[self.NM * i + j] = self.DATA_TYPE(i * (j + 3) % self.NL) / (5 * self.NL)

        for i in range(0, self.NM):
            for j in range(0, self.NL):
                D[self.NL * i + j] = self.DATA_TYPE((i * (j + 2) + 2) % self.NK) / (5 * self.NK)

    def print_array_custom(self, G: list, name: str):
        for i in range(0, self.NI):
            for j in range(0, self.NL):
                if (i * self.NI + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(G[self.NL * i + j])

    def kernel_regular(self, E: list, A: list, B: list, F: list, C: list, D: list, G: list):
# scop begin
        # E := A * B
        for i in range(0, self.NI):
            for j in range(0, self.NJ):
                E[self.NJ * i + j] = 0.0
                for k in range(0, self.NK):
                    E[self.NJ * i + j] += A[self.NK * i + k] * B[self.NJ * k + j]

        # F := C * D
        for i in range(0, self.NJ):
            for j in range(0, self.NL):
                F[self.NL * i + j] = 0.0
                for k in range(0, self.NM):
                    F[self.NL * i + j] += C[self.NM * i + k] * D[self.NL * k + j]

        # G := E * F
        for i in range(0, self.NI):
            for j in range(0, self.NL):
                G[self.NL * i + j] = 0.0
                for k in range(0, self.NJ):
                    G[self.NL * i + j] += E[self.NJ * i + k] * F[self.NL * k + j]
# scop end

    def kernel_le(self, E: list, A: list, B: list, F: list, C: list, D: list, G: list):
# scop begin
        # E := A * B
        for i in range(0, self.NI):
            for j in range(0, self.NJ):
                tmp = 0.0 # load elimination
                for k in range(0, self.NK):
                    tmp += A[self.NK * i + k] * B[self.NJ * k + j] # load elimination
                E[self.NJ * i + j] = tmp # load elimination

        # F := C * D
        for i in range(0, self.NJ):
            for j in range(0, self.NL):
                tmp = 0.0 # load elimination
                for k in range(0, self.NM):
                    tmp += C[self.NM * i + k] * D[self.NL * k + j] # load elimination
                F[self.NL * i + j] = tmp # load elimination

        # G := E * F
        for i in range(0, self.NI):
            for j in range(0, self.NL):
                tmp = 0.0 # load elimination
                for k in range(0, self.NJ):
                    tmp += E[self.NJ * i + k] * F[self.NL * k + j] # load elimination
                G[self.NL * i + j] = tmp # load elimination
# scop end

class _StrategyNumPy(_3mm):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, E: list, A: list, B: list, F: list, C: list, D: list, G: list):
        for i in range(0, self.NI):
            for j in range(0, self.NK):
                A[i, j] = self.DATA_TYPE((i * j + 1) % self.NI) / (5 * self.NI)

        for i in range(0, self.NK):
            for j in range(0, self.NJ):
                B[i, j] = self.DATA_TYPE((i * (j + 1) + 2) % self.NJ) / (5 * self.NJ)

        for i in range(0, self.NJ):
            for j in range(0, self.NM):
                C[i, j] = self.DATA_TYPE(i * (j + 3) % self.NL) / (5 * self.NL)

        for i in range(0, self.NM):
            for j in range(0, self.NL):
                D[i, j] = self.DATA_TYPE((i * (j + 2) + 2) % self.NK) / (5 * self.NK)

    def print_array_custom(self, G: ndarray, name: str):
        for i in range(0, self.NI):
            for j in range(0, self.NL):
                if (i * self.NI + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(G[i, j])

    def kernel(self, E: ndarray, A: ndarray, B: ndarray, F: ndarray, C: ndarray, D: ndarray, G: ndarray):
# scop begin
        # E := A * B
        E[0:self.NI, 0:self.NJ] = np.dot( A[0:self.NI,0:self.NK], B[0:self.NK, 0:self.NJ] )

        # F := C * D
        F[0:self.NJ,0:self.NL] = np.dot( C[0:self.NJ, 0:self.NM], D[0:self.NM,0:self.NL] )

        # G := E * F
        G[0:self.NI,0:self.NL] = np.dot( E[0:self.NI,0:self.NJ], F[0:self.NJ,0:self.NL] )
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        self.kernel = getattr( self, "kernel_%s" % (options.POCC) )

    def kernel_pluto(self, E: list, A: list, B: list, F: list, C: list, D: list, G: list):
# scop begin
        if((self.NL-1>= 0)):
            for c1 in range (min((self.NI-1)+1 , (self.NJ-1)+1)):
                for c2 in range ((self.NL-1)+1):
                    G[self.NL*(c1) + c2] = 0.0
                    F[self.NL*(c1) + c2] = 0.0
        if((self.NL-1>= 0)):
            for c1 in range (max(0 , self.NI) , (self.NJ-1)+1):
                for c2 in range ((self.NL-1)+1):
                    F[self.NL*(c1) + c2] = 0.0
        if((self.NL-1>= 0)):
            for c1 in range (max(0 , self.NJ) , (self.NI-1)+1):
                for c2 in range ((self.NL-1)+1):
                    G[self.NL*(c1) + c2] = 0.0
        if((self.NL-1>= 0) and (self.NM-1>= 0)):
            for c1 in range ((self.NJ-1)+1):
                for c2 in range ((self.NL-1)+1):
                    for c5 in range ((self.NM-1)+1):
                        F[self.NL*(c1) + c2] += C[self.NM*(c1) + c5] * D[self.NL*(c5) + c2]
        if((self.NJ-1>= 0)):
            for c1 in range ((self.NI-1)+1):
                for c2 in range ((self.NJ-1)+1):
                    E[self.NJ*(c1) + c2] = 0.0
        if((self.NJ-1>= 0) and (self.NK-1>= 0) and (self.NL-1>= 0)):
            for c1 in range ((self.NI-1)+1):
                for c2 in range ((self.NJ-1)+1):
                    for c5 in range ((self.NK-1)+1):
                        E[self.NJ*(c1) + c2] += A[self.NK*(c1) + c5] * B[self.NJ*(c5) + c2]
                    for c5 in range ((self.NL-1)+1):
                        G[self.NL*(c1) + c5] += E[self.NJ*(c1) + c2] * F[self.NL*(c2) + c5]
        if((self.NJ-1>= 0) and (self.NK-1>= 0) and (self.NL*-1>= 0)):
            for c1 in range ((self.NI-1)+1):
                for c2 in range ((self.NJ-1)+1):
                    for c5 in range ((self.NK-1)+1):
                        E[self.NJ*(c1) + c2] += A[self.NK*(c1) + c5] * B[self.NJ*(c5) + c2]
        if((self.NJ-1>= 0) and (self.NK*-1>= 0) and (self.NL-1>= 0)):
            for c1 in range ((self.NI-1)+1):
                for c2 in range ((self.NJ-1)+1):
                    for c5 in range ((self.NL-1)+1):
                        G[self.NL*(c1) + c5] += E[self.NJ*(c1) + c2] * F[self.NL*(c2) + c5]
# scop end

    def kernel_vectorizer(self, E: list, A: list, B: list, F: list, C: list, D: list, G: list):
# --pluto --pluto-prevector --vectorizer --pragmatizer
# scop begin
        if((self.NL-1>= 0)):
            for c1 in range (min((self.NI-1)+1 , (self.NJ-1)+1)):
                for c2 in range ((self.NL-1)+1):
                    G[self.NL*(c1) + c2] = 0.0
                    F[self.NL*(c1) + c2] = 0.0
        if((self.NL-1>= 0)):
            for c1 in range (max(0 , self.NI) , (self.NJ-1)+1):
                for c2 in range ((self.NL-1)+1):
                    F[self.NL*(c1) + c2] = 0.0
        if((self.NL-1>= 0)):
            for c1 in range (max(0 , self.NJ) , (self.NI-1)+1):
                for c2 in range ((self.NL-1)+1):
                    G[self.NL*(c1) + c2] = 0.0
        if((self.NL-1>= 0) and (self.NM-1>= 0)):
            for c1 in range ((self.NJ-1)+1):
                for c5 in range ((self.NM-1)+1):
                    for c2 in range ((self.NL-1)+1):
                        F[self.NL*(c1) + c2] += C[self.NM*(c1) + c5] * D[self.NL*(c5) + c2]
        if((self.NJ-1>= 0)):
            for c1 in range ((self.NI-1)+1):
                for c2 in range ((self.NJ-1)+1):
                    E[self.NJ*(c1) + c2] = 0.0
        if((self.NJ-1>= 0) and (self.NK-1>= 0) and (self.NL-1>= 0)):
            for c1 in range ((self.NI-1)+1):
                for c2 in range ((self.NJ-1)+1):
                    for c5 in range ((self.NK-1)+1):
                        E[self.NJ*(c1) + c2] += A[self.NK*(c1) + c5] * B[self.NJ*(c5) + c2]
                    for c5 in range ((self.NL-1)+1):
                        G[self.NL*(c1) + c5] += E[self.NJ*(c1) + c2] * F[self.NL*(c2) + c5]
        if((self.NJ-1>= 0) and (self.NK-1>= 0) and (self.NL*-1>= 0)):
            for c1 in range ((self.NI-1)+1):
                for c5 in range ((self.NK-1)+1):
                    for c2 in range ((self.NJ-1)+1):
                        E[self.NJ*(c1) + c2] += A[self.NK*(c1) + c5] * B[self.NJ*(c5) + c2]
        if((self.NJ-1>= 0) and (self.NK*-1>= 0) and (self.NL-1>= 0)):
            for c1 in range ((self.NI-1)+1):
                for c2 in range ((self.NJ-1)+1):
                    for c5 in range ((self.NL-1)+1):
                        G[self.NL*(c1) + c5] += E[self.NJ*(c1) + c2] * F[self.NL*(c2) + c5]
# scop end

    def kernel_maxfuse(self, E: list, A: list, B: list, F: list, C: list, D: list, G: list):
# --pluto --pluto-fuse maxfuse
# scop begin
        if((self.NI-1>= 0) and (self.NK-1>= 0) and (self.NM-1>= 0)):
            for c0 in range (min((self.NJ-1)+1 , (self.NL-1)+1)):
                for c1 in range (min((self.NI-1)+1 , (self.NL-1)+1)):
                    G[(c1)*self.NL + c0] = 0.0
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range ((self.NM-1)+1):
                        F[(c0)*self.NL + c1] += C[(c0)*self.NM + c6] * D[(c6)*self.NL + c1]
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range ((self.NK-1)+1):
                        E[(c1)*self.NJ + c0] += A[(c1)*self.NK + c6] * B[(c6)*self.NJ + c0]
                    for c6 in range (max(0 , c0 * -1 + c1) , (c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NL , (self.NI-1)+1):
                    G[(c1)*self.NL + c0] = 0.0
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range ((self.NK-1)+1):
                        E[(c1)*self.NJ + c0] += A[(c1)*self.NK + c6] * B[(c6)*self.NJ + c0]
                    for c6 in range (c0 * -1 + c1 , (c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NI , min((self.NL-1)+1 , (self.NI + c0-1)+1)):
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range ((self.NM-1)+1):
                        F[(c0)*self.NL + c1] += C[(c0)*self.NM + c6] * D[(c6)*self.NL + c1]
                    for c6 in range (max(0 , c0 * -1 + c1) , (self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NI + c0 , (self.NL-1)+1):
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range ((self.NM-1)+1):
                        F[(c0)*self.NL + c1] += C[(c0)*self.NM + c6] * D[(c6)*self.NL + c1]
                for c1 in range (max(self.NI , self.NL) , (self.NI + c0-1)+1):
                    for c6 in range (c0 * -1 + c1 , (self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
        if((self.NI-1>= 0) and (self.NK-1>= 0) and (self.NL-1>= 0) and (self.NM-1>= 0)):
            for c0 in range (self.NL , (self.NJ-1)+1):
                for c1 in range (min((self.NI-1)+1 , (self.NL-1)+1)):
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range ((self.NM-1)+1):
                        F[(c0)*self.NL + c1] += C[(c0)*self.NM + c6] * D[(c6)*self.NL + c1]
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range ((self.NK-1)+1):
                        E[(c1)*self.NJ + c0] += A[(c1)*self.NK + c6] * B[(c6)*self.NJ + c0]
                    for c6 in range ((c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NL , (self.NI-1)+1):
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range ((self.NK-1)+1):
                        E[(c1)*self.NJ + c0] += A[(c1)*self.NK + c6] * B[(c6)*self.NJ + c0]
                    for c6 in range (self.NL * -1 + c1 + 1 , (c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NI , (self.NL-1)+1):
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range ((self.NM-1)+1):
                        F[(c0)*self.NL + c1] += C[(c0)*self.NM + c6] * D[(c6)*self.NL + c1]
                    for c6 in range ((self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (max(self.NI , self.NL) , (self.NI + self.NL-2)+1):
                    for c6 in range (self.NL * -1 + c1 + 1 , (self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
        if((self.NI-1>= 0) and (self.NK-1>= 0) and (self.NM*-1>= 0)):
            for c0 in range (min((self.NJ-1)+1 , (self.NL-1)+1)):
                for c1 in range (min((self.NI-1)+1 , (self.NL-1)+1)):
                    G[(c1)*self.NL + c0] = 0.0
                    F[(c0)*self.NL + c1] = 0.0
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range ((self.NK-1)+1):
                        E[(c1)*self.NJ + c0] += A[(c1)*self.NK + c6] * B[(c6)*self.NJ + c0]
                    for c6 in range (max(0 , c0 * -1 + c1) , (c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NL , (self.NI-1)+1):
                    G[(c1)*self.NL + c0] = 0.0
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range ((self.NK-1)+1):
                        E[(c1)*self.NJ + c0] += A[(c1)*self.NK + c6] * B[(c6)*self.NJ + c0]
                    for c6 in range (c0 * -1 + c1 , (c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NI , min((self.NL-1)+1 , (self.NI + c0-1)+1)):
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range (max(0 , c0 * -1 + c1) , (self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NI + c0 , (self.NL-1)+1):
                    F[(c0)*self.NL + c1] = 0.0
                for c1 in range (max(self.NI , self.NL) , (self.NI + c0-1)+1):
                    for c6 in range (c0 * -1 + c1 , (self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
        if((self.NI-1>= 0) and (self.NK-1>= 0) and (self.NL-1>= 0) and (self.NM*-1>= 0)):
            for c0 in range (self.NL , (self.NJ-1)+1):
                for c1 in range (min((self.NI-1)+1 , (self.NL-1)+1)):
                    F[(c0)*self.NL + c1] = 0.0
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range ((self.NK-1)+1):
                        E[(c1)*self.NJ + c0] += A[(c1)*self.NK + c6] * B[(c6)*self.NJ + c0]
                    for c6 in range ((c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NL , (self.NI-1)+1):
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range ((self.NK-1)+1):
                        E[(c1)*self.NJ + c0] += A[(c1)*self.NK + c6] * B[(c6)*self.NJ + c0]
                    for c6 in range (self.NL * -1 + c1 + 1 , (c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NI , (self.NL-1)+1):
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range ((self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (max(self.NI , self.NL) , (self.NI + self.NL-2)+1):
                    for c6 in range (self.NL * -1 + c1 + 1 , (self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
        if((self.NI-1>= 0) and (self.NK-1>= 0) and (self.NL*-1>= 0)):
            for c0 in range ((self.NJ-1)+1):
                for c1 in range ((self.NI-1)+1):
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range ((self.NK-1)+1):
                        E[(c1)*self.NJ + c0] += A[(c1)*self.NK + c6] * B[(c6)*self.NJ + c0]
        if((self.NI-1>= 0) and (self.NK*-1>= 0) and (self.NM-1>= 0)):
            for c0 in range (min((self.NJ-1)+1 , (self.NL-1)+1)):
                for c1 in range (min((self.NI-1)+1 , (self.NL-1)+1)):
                    G[(c1)*self.NL + c0] = 0.0
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range ((self.NM-1)+1):
                        F[(c0)*self.NL + c1] += C[(c0)*self.NM + c6] * D[(c6)*self.NL + c1]
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range (max(0 , c0 * -1 + c1) , (c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NL , (self.NI-1)+1):
                    G[(c1)*self.NL + c0] = 0.0
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range (c0 * -1 + c1 , (c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NI , min((self.NL-1)+1 , (self.NI + c0-1)+1)):
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range ((self.NM-1)+1):
                        F[(c0)*self.NL + c1] += C[(c0)*self.NM + c6] * D[(c6)*self.NL + c1]
                    for c6 in range (max(0 , c0 * -1 + c1) , (self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NI + c0 , (self.NL-1)+1):
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range ((self.NM-1)+1):
                        F[(c0)*self.NL + c1] += C[(c0)*self.NM + c6] * D[(c6)*self.NL + c1]
                for c1 in range (max(self.NI , self.NL) , (self.NI + c0-1)+1):
                    for c6 in range (c0 * -1 + c1 , (self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
        if((self.NI-1>= 0) and (self.NK*-1>= 0) and (self.NL-1>= 0) and (self.NM-1>= 0)):
            for c0 in range (self.NL , (self.NJ-1)+1):
                for c1 in range (min((self.NI-1)+1 , (self.NL-1)+1)):
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range ((self.NM-1)+1):
                        F[(c0)*self.NL + c1] += C[(c0)*self.NM + c6] * D[(c6)*self.NL + c1]
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range ((c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NL , (self.NI-1)+1):
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range (self.NL * -1 + c1 + 1 , (c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NI , (self.NL-1)+1):
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range ((self.NM-1)+1):
                        F[(c0)*self.NL + c1] += C[(c0)*self.NM + c6] * D[(c6)*self.NL + c1]
                    for c6 in range ((self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (max(self.NI , self.NL) , (self.NI + self.NL-2)+1):
                    for c6 in range (self.NL * -1 + c1 + 1 , (self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
        if((self.NI-1>= 0) and (self.NK*-1>= 0) and (self.NM*-1>= 0)):
            for c0 in range (min((self.NJ-1)+1 , (self.NL-1)+1)):
                for c1 in range (min((self.NI-1)+1 , (self.NL-1)+1)):
                    G[(c1)*self.NL + c0] = 0.0
                    F[(c0)*self.NL + c1] = 0.0
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range (max(0 , c0 * -1 + c1) , (c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NL , (self.NI-1)+1):
                    G[(c1)*self.NL + c0] = 0.0
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range (c0 * -1 + c1 , (c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NI , min((self.NL-1)+1 , (self.NI + c0-1)+1)):
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range (max(0 , c0 * -1 + c1) , (self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NI + c0 , (self.NL-1)+1):
                    F[(c0)*self.NL + c1] = 0.0
                for c1 in range (max(self.NI , self.NL) , (self.NI + c0-1)+1):
                    for c6 in range (c0 * -1 + c1 , (self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
        if((self.NI-1>= 0) and (self.NK*-1>= 0) and (self.NL-1>= 0) and (self.NM*-1>= 0)):
            for c0 in range (self.NL , (self.NJ-1)+1):
                for c1 in range (min((self.NI-1)+1 , (self.NL-1)+1)):
                    F[(c0)*self.NL + c1] = 0.0
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range ((c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NL , (self.NI-1)+1):
                    E[(c1)*self.NJ + c0] = 0.0
                    for c6 in range (self.NL * -1 + c1 + 1 , (c1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (self.NI , (self.NL-1)+1):
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range ((self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (max(self.NI , self.NL) , (self.NI + self.NL-2)+1):
                    for c6 in range (self.NL * -1 + c1 + 1 , (self.NI-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
        if((self.NI-1>= 0) and (self.NK*-1>= 0) and (self.NL*-1>= 0)):
            for c0 in range ((self.NJ-1)+1):
                for c1 in range ((self.NI-1)+1):
                    E[(c1)*self.NJ + c0] = 0.0
        if((self.NI*-1>= 0) and (self.NL-1>= 0) and (self.NM-1>= 0)):
            for c0 in range ((self.NJ-1)+1):
                for c1 in range ((self.NL-1)+1):
                    F[(c0)*self.NL + c1] = 0.0
                    for c6 in range ((self.NM-1)+1):
                        F[(c0)*self.NL + c1] += C[(c0)*self.NM + c6] * D[(c6)*self.NL + c1]
        if((self.NI*-1>= 0) and (self.NL-1>= 0) and (self.NM*-1>= 0)):
            for c0 in range ((self.NJ-1)+1):
                for c1 in range ((self.NL-1)+1):
                    F[(c0)*self.NL + c1] = 0.0
        if((self.NI-1>= 0) and (self.NJ-1>= 0)):
            for c0 in range (self.NJ , (self.NL-1)+1):
                for c1 in range (min((self.NI-1)+1 , (self.NJ * -1 + c0)+1)):
                    G[(c1)*self.NL + c0] = 0.0
                for c1 in range (self.NJ * -1 + c0 + 1 , (self.NI-1)+1):
                    G[(c1)*self.NL + c0] = 0.0
                    for c6 in range (max(0 , c0 * -1 + c1) , (self.NJ + c0 * -1 + c1-1)+1):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
                for c1 in range (max(self.NI , self.NJ * -1 + c0 + 1) , (self.NI + c0-1)+1):
                    for c6 in range (max(0 , c0 * -1 + c1) , min((self.NI-1)+1 , (self.NJ + c0 * -1 + c1-1)+1)):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
        if((self.NI-1>= 0) and (self.NJ*-1>= 0)):
            for c0 in range ((self.NL-1)+1):
                for c1 in range ((self.NI-1)+1):
                    G[(c1)*self.NL + c0] = 0.0
        if((self.NI-1>= 0)):
            for c0 in range (max(self.NJ , self.NL) , (self.NJ + self.NL-2)+1):
                for c1 in range (self.NJ * -1 + c0 + 1 , (self.NI + self.NL-2)+1):
                    for c6 in range (max(0 , self.NL * -1 + c1 + 1) , min((self.NI-1)+1 , (self.NJ + c0 * -1 + c1-1)+1)):
                        G[(c6)*self.NL + c1 + (-1 * c6)] += E[(c6)*self.NJ + (c0 + (-1 * c1)) + c6] * F[((c0 + (-1 * c1)) + c6)*self.NL + c1 + (-1 * c6)]
# scop end
