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
from benchmarks.polybench_classes import PolyBenchOptions, PolyBenchSpec, ArrayImplementation
import numpy as np


class Durbin(PolyBench):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        implementation = options.POLYBENCH_ARRAY_IMPLEMENTATION
        if implementation == ArrayImplementation.LIST:
            return _StrategyList.__new__(_StrategyList, options, parameters)
        elif implementation == ArrayImplementation.LIST_PLUTO:
            return _StrategyListPluto.__new__(_StrategyListPluto, options, parameters)
        elif implementation == ArrayImplementation.LIST_FLATTENED:
            return _StrategyList.__new__(_StrategyList, options, parameters)
        elif implementation == ArrayImplementation.NUMPY:
            return _StrategyNumPy.__new__(_StrategyNumPy, options, parameters)
        elif implementation == ArrayImplementation.LIST_FLATTENED_PLUTO:
            return _StrategyListPluto.__new__(_StrategyListPluto, options, parameters)

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
        r = self.create_array(1, [self.N], self.DATA_TYPE(0))
        y = self.create_array(1, [self.N], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(r, y)

        # Benchmark the kernel
        self.time_kernel(r, y)

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
        return [('y', y)]

class _StrategyList(Durbin):
    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, r: list, y: list):
        for i in range(0, self.N):
            r[i] = self.N + 1 - i

    def print_array_custom(self, y: list, name: str):
        for i in range(0, self.N):
            if i % 20 == 0:
                self.print_message('\n')
            self.print_value(y[i])

    def kernel(self, r: list, y: list):
        z = self.create_array(1, [self.N], self.DATA_TYPE(0))
# scop begin
        y[0] = -r[0]
        beta = 1.0
        alpha = -r[0]

        for k in range(1, self.N):
            beta = (1-alpha * alpha) * beta
            summ = 0.0

            for i in range(0, k):
                summ += r[k-i-1] * y[i]

            alpha = -(r[k] + summ) / beta

            for i in range(0, k):
                z[i] = y[i] + alpha * y[k-i-1]

            for i in range(0, k):
                y[i] = z[i]

            y[k] = alpha
# scop end

class _StrategyListPluto(_StrategyList):
    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, r: list, y: list):
        z = self.create_array(1, [self.N], self.DATA_TYPE(0))
# scop begin
# --pluto: some kind of bug, results incorrect
#        y[0] = -r[0]
#        beta = 1.0
#        alpha = -r[0]
#        for c1 in range (1 , (self.N-1)+1):
#            summ = 0.0
#            summ += r[c1-0 -1]*y[0]
#            beta = (1-alpha*alpha)*beta
#            for c2 in range (1 , (c1-1)+1):
#                summ += r[c1-c2-1]*y[c2]
#            alpha = - (r[c1] + summ)/beta
#            y[c1] = alpha
#            for c2 in range (c1 , (c1 * 2-1)+1):
#                z[(-1 * c1) + c2] = y[(-1 * c1) + c2] + alpha*y[c2-((-1*c1)+c2)-1]
#            for c2 in range (c1 * 2 , (c1 * 3-1)+1):
#                y[(-2 * c1) + c2] = z[(-2 * c1) + c2]
# --pluto --pluto-noskew
        y[0] = -r[0]
        beta = 1.0
        alpha = -r[0]
        for c1 in range (1 , (self.N-1)+1):
            sum = 0.0
            for c8 in range ((c1-1)+1):
                sum += r[c1-c8-1]*y[c8]
            beta = (1-alpha*alpha)*beta
            alpha = - (r[c1] + sum)/beta
            y[c1] = alpha
            for c8 in range ((c1-1)+1):
                z[c8] = y[c8] + alpha*y[c1-c8-1]
            for c8 in range ((c1-1)+1):
                y[c8] = z[c8]
# scop end

class _StrategyNumPy(Durbin):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, r: list, y: list):
        for i in range(0, self.N):
            r[i] = self.N + 1 - i

    def print_array_custom(self, y: list, name: str):
        for i in range(0, self.N):
            if i % 20 == 0:
                self.print_message('\n')
            self.print_value(y[i])

    def kernel(self, r: list, y: list):
        z = self.create_array(1, [self.N], self.DATA_TYPE(0))
# scop begin
        y[0] = -r[0]
        beta = 1.0
        alpha = -r[0]
        for k in range(1, self.N):
            beta = (1-alpha * alpha) * beta
            summ = np.dot( r[k-1::-1], y[0:k] )
            alpha = -(r[k] + summ) / beta
            z[0:k] = y[0:k] + alpha * y[k-1::-1]
            y[0:k] = z[0:k]
            y[k] = alpha

# --pluto --pluto-fuse maxfuse
#        y[0] = -r[0]
#        beta = 1.0
#        alpha = -r[0]
#        for c0 in range (1 , (self.N-1)+1):
#            sum = 0.0
#            sum += r[c0-0 -1]*y[0]
#            beta = (1-alpha*alpha)*beta
#            for c1 in range (1 , (c0-1)+1):
#                sum += r[c0-c1-1]*y[c1]
#            alpha = - (r[c0] + sum)/beta
#            y[c0] = alpha
#            for c1 in range (c0 , (c0 * 2-1)+1):
#                z[(-1 * c0) + c1] = y[(-1 * c0) + c1] + alpha*y[c0-(-1 * c0) + c1-1]
#            for c1 in range (c0 * 2 , (c0 * 3-1)+1):
#                y[(-2 * c0) + c1] = z[(-2 * c0) + c1]
#
# scop end
