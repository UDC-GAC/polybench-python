# Copyright 2019 Miguel Angel Abella Gonzalez <miguel.abella@udc.es>
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

"""<replace_with_module_description>"""

from benchmarks.polybench import PolyBench
from benchmarks.polybench_classes import ArrayImplementation
from benchmarks.polybench_classes import PolyBenchOptions, PolyBenchSpec
from numpy.core.multiarray import ndarray
import numpy as np


class Nussinov(PolyBench):

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
        seq = self.create_array(1, [self.N], int(0))  # base type = char; = int in Python
        table = self.create_array(2, [self.N, self.N], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(seq, table)

        # Benchmark the kernel
        self.time_kernel(seq, table)

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
        return [('table', table)]


class _StrategyList(Nussinov):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, seq: list, table: list):
        for i in range(0, self.N):
            seq[i] = (i + 1) % 4  # right side is AGCT/0..3

        for i in range(0, self.N):
            for j in range(0, self.N):
                table[i][j] = self.DATA_TYPE(0)

    def print_array_custom(self, table: list, name: str):
        t = 0
        for i in range(0, self.N):
            for j in range(i, self.N):
                if t % 20 == 0:
                    self.print_message('\n')
                self.print_value(table[i][j])
                t += 1

    def kernel(self, seq: list, table: list):
        # def max_score(s1, s2):
        #     if s1 >= s2:
        #         return s1
        #     else:
        #         return s2
        #
        def match(b1, b2): return b1+b2==3
# scop begin
        for i in range ((self.N-1)+1):
            for j in range (self.N + i * -1 , (self.N-1)+1):
                if((j-1>= 0)):
                    table[self.N-1-i][j] = max(table[self.N-1-i][j], table[self.N-1-i][j-1])
                if((i-1>= 0)):
                    table[self.N-1-i][j] = max(table[self.N-1-i][j], table[self.N-1-i+1][j])
                if((i-1>= 0) and (j-1>= 0) and (self.N*-1+i+j-1>= 0)):
                    table[self.N-1-i][j] = max(table[self.N-1-i][j], table[self.N-1-i+1][j-1]+match(seq[self.N-1-i], seq[j]))
                if((i-1>= 0) and (j-1>= 0) and (self.N+i*-1+j*-1>= 0)):
                    table[self.N-1-i][j] = max(table[self.N-1-i][j], table[self.N-1-i+1][j-1])
                for k in range (self.N + i * -1 , (j-1)+1):
                    table[self.N-1-i][j] = max(table[self.N-1-i][j], table[self.N-1-i][k] + table[k+1][j])
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, seq: list, table: list):
        def match(b1,b2):
            return b1+b2 == 3
# scop begin
        if((self.N-2>= 0)):
            table[self.N-1-1][self.N + -1] = max(table[self.N-1-1][self.N + -1], table[self.N-1-1][self.N + -1 -1])
            table[self.N-1-1][self.N + -1] = max(table[self.N-1-1][self.N + -1], table[self.N-1-1 +1][self.N + -1])
            table[self.N-1-1][self.N + -1] = max(table[self.N-1-1][self.N + -1], table[self.N-1-1 +1][self.N + -1 -1])
            for c0 in range (2 , (self.N-1)+1):
                table[self.N-1-c0][(-1 * c0) + self.N] = max(table[self.N-1-c0][(-1 * c0) + self.N], table[self.N-1-c0][(-1 * c0) + self.N-1])
                table[self.N-1-c0][(-1 * c0) + self.N] = max(table[self.N-1-c0][(-1 * c0) + self.N], table[self.N-1-c0+1][(-1 * c0) + self.N])
                table[self.N-1-c0][(-1 * c0) + self.N] = max(table[self.N-1-c0][(-1 * c0) + self.N], table[self.N-1-c0+1][(-1 * c0) + self.N-1])
                for c1 in range (self.N + c0 * -1 + 1 , (self.N-1)+1):
                    table[self.N-1-c0][c1] = max(table[self.N-1-c0][c1], table[self.N-1-c0][c1-1])
                    table[self.N-1-c0][c1] = max(table[self.N-1-c0][c1], table[self.N-1-c0+1][c1])
                    table[self.N-1-c0][c1] = max(table[self.N-1-c0][c1], table[self.N-1-c0+1][c1-1]+match(seq[self.N-1-c0], seq[c1]))
                    for c2 in range (self.N + c0 * -1 , (c1-1)+1):
                        table[self.N-1-c0][c1] = max(table[self.N-1-c0][c1], table[self.N-1-c0][c2] + table[c2+1][c1])
# scop end

class _StrategyListFlattened(Nussinov):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, seq: list, table: list):
        for i in range(0, self.N):
            seq[i] = (i + 1) % 4  # right side is AGCT/0..3

        for i in range(0, self.N):
            for j in range(0, self.N):
                table[self.N * i + j] = self.DATA_TYPE(0)

    def print_array_custom(self, table: list, name: str):
        t = 0
        for i in range(0, self.N):
            for j in range(i, self.N):
                if t % 20 == 0:
                    self.print_message('\n')
                self.print_value(table[self.N * i + j])
                t += 1

    def kernel(self, seq: list, table: list):
        def match(b1,b2): return b1+b2==3
# scop begin
        for i in range ((self.N-1)+1):
            for j in range (self.N + i * -1 , (self.N-1)+1):
                if((j-1>= 0)):
                    table[self.N*(self.N-1-i) + j] = max(table[self.N*(self.N-1-i)+j], table[self.N*(self.N-1-i)+j-1])
                if((i-1>= 0)):
                    table[self.N*(self.N-1-i) + j] = max(table[self.N*(self.N-1-i)+j], table[self.N*(self.N-1-i+1)+j])
                if((i-1>= 0) and (j-1>= 0) and (self.N*-1+i+j-1>= 0)):
                    table[self.N*(self.N-1-i)+j] = max(table[self.N*(self.N-1-i)+j], table[self.N*(self.N-1-i+1)+j-1]+match(seq[self.N-1-i], seq[j]))
                if((i-1>= 0) and (j-1>= 0) and (self.N+i*-1+j*-1>= 0)):
                    table[self.N*(self.N-1-i)+j] = max(table[self.N*(self.N-1-i)+j], table[self.N*(self.N-1-i+1)+j-1])
                for k in range (self.N + i * -1 , (j-1)+1):
                    table[self.N*(self.N-1-i)+j] = max(table[self.N*(self.N-1-i)+j], table[self.N*(self.N-1-i)+k] + table[self.N*(k+1)+j])
# scop end


class _StrategyNumPy(Nussinov):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, seq: ndarray, table: ndarray):
        for i in range(0, self.N):
            seq[i] = (i + 1) % 4  # right side is AGCT/0..3

        for i in range(0, self.N):
            for j in range(0, self.N):
                table[i, j] = self.DATA_TYPE(0)

    def print_array_custom(self, table: ndarray, name: str):
        t = 0
        for i in range(0, self.N):
            for j in range(i, self.N):
                if t % 20 == 0:
                    self.print_message('\n')
                self.print_value(table[i, j])
                t += 1

    def kernel(self, seq: ndarray, table: ndarray):
        # def max_score(s1, s2):
        #     if s1 >= s2:
        #         return s1
        #     else:
        #         return s2
        #
        def match(b1, b2): return b1+b2 == 3
        table_flattened = table.ravel()
# scop begin
        for i in range(self.N):
            diag_i = np.arange( 0, self.N-i-1 )
            diag_j = np.arange( i+1, self.N )
            slice_W      = slice( i , i + self.N*(self.N-i-1) , self.N+1 )
            slice_center = slice( i+1 , i+1 + self.N*(self.N-i-1) , self.N+1 )
            slice_SW = slice( i+self.N , i + self.N*(self.N-i) , self.N+1 )
            slice_S = slice( i+1+self.N , i+1 + self.N*(self.N-i) , self.N+1 )
            table_flattened[slice_center] = np.maximum( table_flattened[slice_center], table_flattened[slice_W] )
            table_flattened[slice_center] = np.maximum( table_flattened[slice_center], table_flattened[slice_S] )
            if i > 0: table_flattened[ slice_center ] = np.maximum( table_flattened[slice_center], table_flattened[slice_SW] + match( seq[0:self.N-1-i], seq[i+1:self.N] ) )
            else: table_flattened[ slice_center] = np.maximum( table_flattened[slice_center], table_flattened[slice_SW] )
            if i < self.N-1: table_flattened[slice_center] = np.maximum( table_flattened[slice_center], np.max( table[0:self.N-i-1, i+1:self.N], axis=0 ) )

# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def kernel(self, seq: list, table: list):
        def match(b1,b2):
            return b1+b2 == 3
# scop begin
        if((self.N-2>= 0)):
            table[self.N*(self.N-1-1) + self.N + -1] = max(table[self.N*(self.N-1-1) + self.N + -1], table[self.N*(self.N-1-1) + self.N + -1 -1])
            table[self.N*(self.N-1-1) + self.N + -1] = max(table[self.N*(self.N-1-1) + self.N + -1], table[self.N*(self.N-1-1 +1) + self.N + -1])
            table[self.N*(self.N-1-1) + self.N + -1] = max(table[self.N*(self.N-1-1) + self.N + -1], table[self.N*(self.N-1-1 +1) + self.N + -1 -1])
            for c0 in range (2 , (self.N-1)+1):
                table[self.N*(self.N-1-c0) + (-1 * c0) + self.N] = max(table[self.N*(self.N-1-c0) + (-1 * c0) + self.N], table[self.N*(self.N-1-c0) + (-1 * c0) + self.N-1])
                table[self.N*(self.N-1-c0) + (-1 * c0) + self.N] = max(table[self.N*(self.N-1-c0) + (-1 * c0) + self.N], table[self.N*(self.N-1-c0+1) + (-1 * c0) + self.N])
                table[self.N*(self.N-1-c0) + (-1 * c0) + self.N] = max(table[self.N*(self.N-1-c0) + (-1 * c0) + self.N], table[self.N*(self.N-1-c0+1) + (-1 * c0) + self.N-1])
                for c1 in range (self.N + c0 * -1 + 1 , (self.N-1)+1):
                    table[self.N*(self.N-1-c0) + c1] = max(table[self.N*(self.N-1-c0) + c1], table[self.N*(self.N-1-c0) + c1-1])
                    table[self.N*(self.N-1-c0) + c1] = max(table[self.N*(self.N-1-c0) + c1], table[self.N*(self.N-1-c0+1) + c1])
                    table[self.N*(self.N-1-c0) + c1] = max(table[self.N*(self.N-1-c0) + c1], table[self.N*(self.N-1-c0+1) + c1-1]+match(seq[self.N-1-c0], seq[c1]))
                    for c2 in range (self.N + c0 * -1 , (c1-1)+1):
                        table[self.N*(self.N-1-c0) + c1] = max(table[self.N*(self.N-1-c0) + c1], table[self.N*(self.N-1-c0) + c2] + table[self.N*(c2+1) + c1])
# scop end
