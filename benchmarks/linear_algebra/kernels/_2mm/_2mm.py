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

class _2mm(PolyBench):

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

    def run_benchmark(self):
        # Create data structures (arrays, auxiliary variables, etc.)
        alpha = 1.5
        beta = 1.2

        tmp = self.create_array(2, [self.NI, self.NJ], self.DATA_TYPE(0))
        A = self.create_array(2, [self.NI, self.NK], self.DATA_TYPE(0))
        B = self.create_array(2, [self.NK, self.NJ], self.DATA_TYPE(0))
        C = self.create_array(2, [self.NJ, self.NL], self.DATA_TYPE(0))
        D = self.create_array(2, [self.NI, self.NL], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(alpha, beta, tmp, A, B, C, D)

        # Benchmark the kernel
        self.time_kernel(alpha, beta, tmp, A, B, C, D)

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
        return [('D', D)]


class _StrategyList(_2mm):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, alpha, beta, tmp: list, A: list, B: list, C: list, D: list):
        for i in range(0, self.NI):
            for j in range(0, self.NK):
                A[i][j] = self.DATA_TYPE((i * j + 1) % self.NI) / self.NI

        for i in range(0, self.NK):
            for j in range(0, self.NJ):
                B[i][j] = self.DATA_TYPE(i * (j + 1) % self.NJ) / self.NJ

        for i in range(0, self.NJ):
            for j in range(0, self.NL):
                C[i][j] = self.DATA_TYPE((i * (j + 3) + 1) % self.NL) / self.NL

        for i in range(0, self.NI):
            for j in range(0, self.NL):
                D[i][j] = self.DATA_TYPE(i * (j + 2) % self.NK) / self.NK

    def print_array_custom(self, D: list, name: str):
        for i in range(0, self.NI):
            for j in range(0, self.NL):
                if (i * self.NI + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(D[i][j])

    def kernel(self, alpha, beta, tmp: list, A: list, B: list, C: list, D: list):
# scop begin
        # D := alpha * A * B * C + beta * D
        for i in range(self.NI):
            for j in range(self.NJ):
                tmp[i][j] = 0.0
                for k in range(0, self.NK):
                    tmp[i][j] += alpha * A[i][k] * B[k][j]

        for i in range(0, self.NI):
            for j in range(0, self.NL):
                D[i][j] *= beta
                for k in range(0, self.NJ):
                    D[i][j] += tmp[i][k] * C[k][j]
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, alpha, beta, tmp: list, A: list, B: list, C: list, D: list):
# scop begin
        if((self.NI-1>= 0)):
            if((self.NJ-1>= 0) and (self.NL-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range (min((self.NJ-1)+1 , (self.NL-1)+1)):
                        D[c1][c2] *= beta
                        tmp[c1][c2] = 0.0
                    for c2 in range (self.NL , (self.NJ-1)+1):
                        tmp[c1][c2] = 0.0
                    for c2 in range (self.NJ , (self.NL-1)+1):
                        D[c1][c2] *= beta
            if((self.NJ-1>= 0) and (self.NL*-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range ((self.NJ-1)+1):
                        tmp[c1][c2] = 0.0
            if((self.NJ*-1>= 0) and (self.NL-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range ((self.NL-1)+1):
                        D[c1][c2] *= beta
            if((self.NJ-1>= 0) and (self.NK-1>= 0) and (self.NL-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range ((self.NJ-1)+1):
                        for c5 in range ((self.NK-1)+1):
                            tmp[c1][c2] += alpha * A[c1][c5] * B[c5][c2]
                        for c5 in range ((self.NL-1)+1):
                            D[c1][c5] += tmp[c1][c2] * C[c2][c5]
            if((self.NJ-1>= 0) and (self.NK-1>= 0) and (self.NL*-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range ((self.NJ-1)+1):
                        for c5 in range ((self.NK-1)+1):
                            tmp[c1][c2] += alpha * A[c1][c5] * B[c5][c2]
            if((self.NJ-1>= 0) and (self.NK*-1>= 0) and (self.NL-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range ((self.NJ-1)+1):
                        for c5 in range ((self.NL-1)+1):
                            D[c1][c5] += tmp[c1][c2] * C[c2][c5]
# scop end


class _StrategyListFlattened(_2mm):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        if options.LOAD_ELIMINATION: self.kernel = self.kernel_le
        else: self.kernel = self.kernel_regular

    def initialize_array(self, alpha, beta, tmp: list, A: list, B: list, C: list, D: list):
        for i in range(0, self.NI):
            for j in range(0, self.NK):
                A[self.NK * i + j] = self.DATA_TYPE((i * j + 1) % self.NI) / self.NI

        for i in range(0, self.NK):
            for j in range(0, self.NJ):
                B[self.NJ * i + j] = self.DATA_TYPE(i * (j + 1) % self.NJ) / self.NJ

        for i in range(0, self.NJ):
            for j in range(0, self.NL):
                C[self.NL * i + j] = self.DATA_TYPE((i * (j + 3) + 1) % self.NL) / self.NL

        for i in range(0, self.NI):
            for j in range(0, self.NL):
                D[self.NL * i + j] = self.DATA_TYPE(i * (j + 2) % self.NK) / self.NK

    def print_array_custom(self, D: list, name: str):
        for i in range(0, self.NI):
            for j in range(0, self.NL):
                if (i * self.NI + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(D[self.NL * i + j])

    def kernel_regular(self, alpha, beta, tmp: list, A: list, B: list, C: list, D: list):
# scop begin
        # D := alpha * A * B * C + beta * D
        for i in range(self.NI):
            for j in range(self.NJ):
                tmp[self.NJ * i + j] = 0.0
                for k in range(0, self.NK):
                    tmp[self.NJ * i + j] += alpha * A[self.NK * i + k] * B[self.NJ * k + j]

        for i in range(0, self.NI):
            for j in range(0, self.NL):
                D[self.NL * i + j] *= beta
                for k in range(0, self.NJ):
                    D[self.NL * i + j] += tmp[self.NJ * i + k] * C[self.NL * k + j]
# scop end

    def kernel_le(self, alpha, beta, tmp: list, A: list, B: list, C: list, D: list):
# scop begin
        # D := alpha * A * B * C + beta * D
        for i in range(self.NI):
            for j in range(self.NJ):
                tmp2 = 0.0 # load elimination
                for k in range(0, self.NK):
                    tmp2 += alpha * A[self.NK * i + k] * B[self.NJ * k + j] # load elimination
                tmp[self.NJ * i + j] = tmp2 # load elimination

        for i in range(0, self.NI):
            for j in range(0, self.NL):
                tmp2 = D[self.NL * i + j] * beta # load elimination
                for k in range(0, self.NJ):
                    tmp2 += tmp[self.NJ * i + k] * C[self.NL * k + j] # load elimination
                D[self.NL * i + j] = tmp2 #load elimination
# scop end


class _StrategyNumPy(_2mm):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, alpha, beta, tmp: list, A: list, B: list, C: list, D: list):
        for i in range(0, self.NI):
            for j in range(0, self.NK):
                A[i, j] = self.DATA_TYPE((i * j + 1) % self.NI) / self.NI

        for i in range(0, self.NK):
            for j in range(0, self.NJ):
                B[i, j] = self.DATA_TYPE(i * (j + 1) % self.NJ) / self.NJ

        for i in range(0, self.NJ):
            for j in range(0, self.NL):
                C[i, j] = self.DATA_TYPE((i * (j + 3) + 1) % self.NL) / self.NL

        for i in range(0, self.NI):
            for j in range(0, self.NL):
                D[i, j] = self.DATA_TYPE(i * (j + 2) % self.NK) / self.NK

    def print_array_custom(self, D: ndarray, name: str):
        for i in range(0, self.NI):
            for j in range(0, self.NL):
                if (i * self.NI + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(D[i, j])

    def kernel(self, alpha, beta, tmp: ndarray, A: ndarray, B: ndarray, C: ndarray, D: ndarray):
# scop begin
        # D := alpha * A * B * C + beta * D
        tmp[0:self.NI, 0:self.NJ] = alpha * np.dot( A[0:self.NI, 0:self.NK], B[0:self.NK, 0:self.NK] )

        D[0:self.NI,0:self.NL] *= beta
        D[0:self.NI,0:self.NL] += np.dot( tmp[0:self.NI,0:self.NJ], C[0:self.NJ,0:self.NL] )
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        self.kernel = getattr( self, "kernel_%s" % (options.POCC) )

    def kernel_pluto(self, alpha, beta, tmp: list, A: list, B: list, C: list, D: list):
# --pluto
# scop begin
        if((self.NI-1>= 0)):
            if((self.NJ-1>= 0) and (self.NL-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range (min((self.NJ-1)+1 , (self.NL-1)+1)):
                        D[self.NL*(c1) + c2] *= beta
                        tmp[self.NJ*(c1) + c2] = 0.0
                    for c2 in range (self.NL , (self.NJ-1)+1):
                        tmp[self.NJ*(c1) + c2] = 0.0
                    for c2 in range (self.NJ , (self.NL-1)+1):
                        D[self.NL*(c1) + c2] *= beta
            if((self.NJ-1>= 0) and (self.NL*-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range ((self.NJ-1)+1):
                        tmp[self.NJ*(c1) + c2] = 0.0
            if((self.NJ*-1>= 0) and (self.NL-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range ((self.NL-1)+1):
                        D[self.NL*(c1) + c2] *= beta
            if((self.NJ-1>= 0) and (self.NK-1>= 0) and (self.NL-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range ((self.NJ-1)+1):
                        for c5 in range ((self.NK-1)+1):
                            tmp[self.NJ*(c1) + c2] += alpha * A[self.NK*(c1) + c5] * B[self.NJ*(c5) + c2]
                        for c5 in range ((self.NL-1)+1):
                            D[self.NL*(c1) + c5] += tmp[self.NJ*(c1) + c2] * C[self.NL*(c2) + c5]
            if((self.NJ-1>= 0) and (self.NK-1>= 0) and (self.NL*-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range ((self.NJ-1)+1):
                        for c5 in range ((self.NK-1)+1):
                            tmp[self.NJ*(c1) + c2] += alpha * A[self.NK*(c1) + c5] * B[self.NJ*(c5) + c2]
            if((self.NJ-1>= 0) and (self.NK*-1>= 0) and (self.NL-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range ((self.NJ-1)+1):
                        for c5 in range ((self.NL-1)+1):
                            D[self.NL*(c1) + c5] += tmp[self.NJ*(c1) + c2] * C[self.NL*(c2) + c5]
# scop end

    def kernel_vectorizer(self, alpha, beta, tmp: list, A: list, B: list, C: list, D: list):
# --pluto --pluto-prevector --vectorizer --pragmatizer
# scop begin
        if((self.NI-1>= 0)):
            if((self.NJ-1>= 0) and (self.NL-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range (min((self.NJ-1)+1 , (self.NL-1)+1)):
                        D[self.NL*(c1) + c2] *= beta
                        tmp[self.NJ*(c1) + c2] = 0.0
                    for c2 in range (self.NL , (self.NJ-1)+1):
                        tmp[self.NJ*(c1) + c2] = 0.0
                    for c2 in range (self.NJ , (self.NL-1)+1):
                        D[self.NL*(c1) + c2] *= beta
            if((self.NJ-1>= 0) and (self.NL*-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range ((self.NJ-1)+1):
                        tmp[self.NJ*(c1) + c2] = 0.0
            if((self.NJ*-1>= 0) and (self.NL-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range ((self.NL-1)+1):
                        D[self.NL*(c1) + c2] *= beta
            if((self.NJ-1>= 0) and (self.NK-1>= 0) and (self.NL-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range ((self.NJ-1)+1):
                        for c5 in range ((self.NK-1)+1):
                            tmp[self.NJ*(c1) + c2] += alpha * A[self.NK*(c1) + c5] * B[self.NJ*(c5) + c2]
                        for c5 in range ((self.NL-1)+1):
                            D[self.NL*(c1) + c5] += tmp[self.NJ*(c1) + c2] * C[self.NL*(c2) + c5]
            if((self.NJ-1>= 0) and (self.NK-1>= 0) and (self.NL*-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c5 in range ((self.NK-1)+1):
                        for c2 in range ((self.NJ-1)+1):
                            tmp[self.NJ*(c1) + c2] += alpha * A[self.NK*(c1) + c5] * B[self.NJ*(c5) + c2]
            if((self.NJ-1>= 0) and (self.NK*-1>= 0) and (self.NL-1>= 0)):
                for c1 in range ((self.NI-1)+1):
                    for c2 in range ((self.NJ-1)+1):
                        for c5 in range ((self.NL-1)+1):
                            D[self.NL*(c1) + c5] += tmp[self.NJ*(c1) + c2] * C[self.NL*(c2) + c5]
# scop end

    def kernel_maxfuse(self, alpha, beta, tmp: list, A: list, B: list, C: list, D: list):
# --pluto --pluto-fuse maxfuse
# scop begin
        if((self.NI-1>= 0)):
            if((self.NJ-1>= 0) and (self.NK-1>= 0) and (self.NL-1>= 0)):
                for c0 in range ((self.NI-1)+1):
                    for c1 in range (min((self.NJ-1)+1 , (self.NL-1)+1)):
                        D[(c0)*self.NL + c1] *= beta
                        tmp[(c0)*self.NJ + c1] = 0.0
                        for c6 in range ((self.NK-1)+1):
                            tmp[(c0)*self.NJ + c1] += alpha * A[(c0)*self.NK + c6] * B[(c6)*self.NJ + c1]
                        for c6 in range ((c1)+1):
                            D[(c0)*self.NL + c6] += tmp[(c0)*self.NJ + c1 + (-1 * c6)] * C[(c1 + (-1 * c6))*self.NL + c6]
                    for c1 in range (self.NL , (self.NJ-1)+1):
                        tmp[(c0)*self.NJ + c1] = 0.0
                        for c6 in range ((self.NK-1)+1):
                            tmp[(c0)*self.NJ + c1] += alpha * A[(c0)*self.NK + c6] * B[(c6)*self.NJ + c1]
                        for c6 in range ((self.NL-1)+1):
                            D[(c0)*self.NL + c6] += tmp[(c0)*self.NJ + c1 + (-1 * c6)] * C[(c1 + (-1 * c6))*self.NL + c6]
                    for c1 in range (self.NJ , (self.NL-1)+1):
                        D[(c0)*self.NL + c1] *= beta
                        for c6 in range (self.NJ * -1 + c1 + 1 , (c1)+1):
                            D[(c0)*self.NL + c6] += tmp[(c0)*self.NJ + c1 + (-1 * c6)] * C[(c1 + (-1 * c6))*self.NL + c6]
                    for c1 in range (max(self.NJ , self.NL) , (self.NJ + self.NL-2)+1):
                        for c6 in range (self.NJ * -1 + c1 + 1 , (self.NL-1)+1):
                            D[(c0)*self.NL + c6] += tmp[(c0)*self.NJ + c1 + (-1 * c6)] * C[(c1 + (-1 * c6))*self.NL + c6]
            if((self.NJ-1>= 0) and (self.NK-1>= 0) and (self.NL*-1>= 0)):
                for c0 in range ((self.NI-1)+1):
                    for c1 in range ((self.NJ-1)+1):
                        tmp[(c0)*self.NJ + c1] = 0.0
                        for c6 in range ((self.NK-1)+1):
                            tmp[(c0)*self.NJ + c1] += alpha * A[(c0)*self.NK + c6] * B[(c6)*self.NJ + c1]
            if((self.NJ-1>= 0) and (self.NK*-1>= 0) and (self.NL-1>= 0)):
                for c0 in range ((self.NI-1)+1):
                    for c1 in range (min((self.NJ-1)+1 , (self.NL-1)+1)):
                        D[(c0)*self.NL + c1] *= beta
                        tmp[(c0)*self.NJ + c1] = 0.0
                        for c6 in range ((c1)+1):
                            D[(c0)*self.NL + c6] += tmp[(c0)*self.NJ + c1 + (-1 * c6)] * C[(c1 + (-1 * c6))*self.NL + c6]
                    for c1 in range (self.NL , (self.NJ-1)+1):
                        tmp[(c0)*self.NJ + c1] = 0.0
                        for c6 in range ((self.NL-1)+1):
                            D[(c0)*self.NL + c6] += tmp[(c0)*self.NJ + c1 + (-1 * c6)] * C[(c1 + (-1 * c6))*self.NL + c6]
                    for c1 in range (self.NJ , (self.NL-1)+1):
                        D[(c0)*self.NL + c1] *= beta
                        for c6 in range (self.NJ * -1 + c1 + 1 , (c1)+1):
                            D[(c0)*self.NL + c6] += tmp[(c0)*self.NJ + c1 + (-1 * c6)] * C[(c1 + (-1 * c6))*self.NL + c6]
                    for c1 in range (max(self.NJ , self.NL) , (self.NJ + self.NL-2)+1):
                        for c6 in range (self.NJ * -1 + c1 + 1 , (self.NL-1)+1):
                            D[(c0)*self.NL + c6] += tmp[(c0)*self.NJ + c1 + (-1 * c6)] * C[(c1 + (-1 * c6))*self.NL + c6]
            if((self.NJ-1>= 0) and (self.NK*-1>= 0) and (self.NL*-1>= 0)):
                for c0 in range ((self.NI-1)+1):
                    for c1 in range ((self.NJ-1)+1):
                        tmp[(c0)*self.NJ + c1] = 0.0
            if((self.NJ*-1>= 0) and (self.NL-1>= 0)):
                for c0 in range ((self.NI-1)+1):
                    for c1 in range ((self.NL-1)+1):
                        D[(c0)*self.NL + c1] *= beta
# scop end
