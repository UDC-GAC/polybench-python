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

class Floyd_warshall(PolyBench):

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
        path = self.create_array(2, [self.N, self.N], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(path)

        # Benchmark the kernel
        self.time_kernel(path)

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
        return [('path', path)]


class _StrategyList(Floyd_warshall):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, path: list):
        for i in range(0, self.N):
            for j in range(0, self.N):
                path[i][j] = i * j % 7 + 1
                if (i + j) % 13 == 0 or (i + j) % 7 == 0 or (i + j) % 11 == 0:
                    path[i][j] = 999

    def print_array_custom(self, path: list, name: str):
        for i in range(0, self.N):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(path[i][j])

    def kernel(self, path: list):
# scop begin
        for k in range(0, self.N):
            for i in range(0, self.N):
                for j in range(0, self.N):
                    path[i][j] = path[i][j] if path[i][j] < path[i][k] + path[k][j] else path[i][k]+path[k][j]
#                    if path[i][j] < path[i][k] + path[k][j]:
#                        path[i][j] = path[i][j]
#                    else:
#                        path[i][j] = path[i][k] + path[k][j]
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, path: list):
# scop begin
        if((self.N-1>= 0)):
            for c0 in range ((self.N-1)+1):
                for c1 in range ((self.N-1)+1):
                    for c2 in range ((self.N-1)+1):
                        path[c1][c2] = path[c1][c2] if path[c1][c2] < path[c1][c0] + path[c0][c2] else path[c1][c0]+path[c0][c2]
# scop end

class _StrategyListFlattened(Floyd_warshall):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, path: list):
        for i in range(0, self.N):
            for j in range(0, self.N):
                path[self.N * i + j] = i * j % 7 + 1
                if (i + j) % 13 == 0 or (i + j) % 7 == 0 or (i + j) % 11 == 0:
                    path[self.N * i + j] = 999

    def print_array_custom(self, path: list, name: str):
        for i in range(0, self.N):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(path[self.N * i + j])

    def kernel(self, path: list):
# scop begin
        for k in range(0, self.N):
            for i in range(0, self.N):
                for j in range(0, self.N):
                    if path[self.N * i + j] < path[self.N * i + k] + path[self.N * k + j]:
                        path[self.N * i + j] = path[self.N * i + j]
                    else:
                        path[self.N * i + j] = path[self.N * i + k] + path[self.N * k + j]
# scop end


class _StrategyNumPy(Floyd_warshall):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, path: ndarray):
        for i in range(0, self.N):
            for j in range(0, self.N):
                path[i, j] = i * j % 7 + 1
                if (i + j) % 13 == 0 or (i + j) % 7 == 0 or (i + j) % 11 == 0:
                    path[i, j] = 999

    def print_array_custom(self, path: ndarray, name: str):
        for i in range(0, self.N):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(path[i, j])

    def kernel(self, path: ndarray):
# scop begin
        for k in range(0, self.N):
            path_add = path[0:k+1,k,np.newaxis] + path[k, 0:self.N]
            if k > 0: path[ 0:k, 0:self.N ] = np.minimum( path[0:k,0:self.N], path_add[:k] )

            path[k,0:self.N] = np.minimum( path[k,0:self.N], path_add[k] )

            path_add = path[k+1:,k,np.newaxis] + path[k, 0:self.N]
            if k < self.N-1: path[ k+1:, 0:self.N ] = np.minimum( path[k+1:,0:self.N], path_add )

# Without index set splitting
#            for i in range(0, self.N):
#                path[i,0:self.N] = np.where( path[i,0:self.N] < path[i,k] + path[k,0:self.N], path[i,0:self.N], path[i,k] + path[k,0:self.N] )
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def kernel(self, path: list):
# scop begin
        if((self.N-1>= 0)):
            for c0 in range ((self.N-1)+1):
                for c1 in range ((self.N-1)+1):
                    for c2 in range ((self.N-1)+1):
                        path[self.N*(c1) + c2] = path[self.N*(c1) + c2] if path[self.N*(c1) + c2] < path[self.N*(c1) + c0] + path[self.N*(c0) + c2] else path[self.N*(c1) + c0]+path[self.N*(c0) + c2]
# scop end
