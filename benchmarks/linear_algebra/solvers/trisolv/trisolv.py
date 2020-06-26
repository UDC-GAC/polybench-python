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

from benchmarks.polybench import PolyBench, PolyBenchParameters


class Trisolv(PolyBench):

    def __init__(self, options: dict, parameters: PolyBenchParameters):
        super().__init__(options)

        # Validate inputs
        if not isinstance(parameters, PolyBenchParameters):
            raise AssertionError(f'Invalid parameter "parameters": "{parameters}"')

        # The parameters hold the necessary information obtained from "polybench.spec" file
        params = parameters.DataSets.get(self.DATASET_SIZE)
        if not isinstance(params, dict):
            raise NotImplementedError(f'Dataset size "{self.DATASET_SIZE.name}" not implemented '
                                      f'for {parameters.Category}/{parameters.Name}.')

        # Adjust the data type and print modifier according to the data type
        self.DATA_TYPE = parameters.DataType
        self.set_print_modifier(parameters.DataType)

        # Set up problem size from the given parameters (adapt this part with appropriate parameters)
        self.N = params.get('N')

    def initialize_array(self, L: list, x: list, b: list):
        for i in range(0, self.N):
            x[i] = - 999
            b[i] = i
            if self.POLYBENCH_FLATTEN_LISTS:
                for j in range(0, i + 1):
                    L[self.N * i + j] = self.DATA_TYPE(i+self.N-j+1) * 2 / self.N
            else:
                for j in range(0, i + 1):
                    L[i][j] = self.DATA_TYPE(i + self.N - j + 1) * 2 / self.N

    def print_array_custom(self, x: list, name: str):
        for i in range(0, self.N):
            self.print_value(x[i])
            if i % 20 == 0:
                self.print_message('\n')

    def kernel(self, L: list, x: list, b: list):
# scop begin
        for i in range(0, self.N):
            x[i] = b[i]
            for j in range(0, i):
                x[i] -= L[i][j] * x[j]
            x[i] = x[i] / L[i][i]
# scop end

    def kernel_flat(self, L: list, x: list, b: list):
# scop begin
        for i in range(0, self.N):
            x[i] = b[i]
            for j in range(0, i):
                x[i] -= L[self.N * i + j] * x[j]
            x[i] = x[i] / L[self.N * i + i]
# scop end

    def run_benchmark(self):
        # Create data structures (arrays, auxiliary variables, etc.)
        if self.POLYBENCH_FLATTEN_LISTS:
            L = self.create_array(1, [self.N * self.N], self.DATA_TYPE(0))
        else:
            L = self.create_array(2, [self.N, self.N], self.DATA_TYPE(0))
        x = self.create_array(1, [self.N], self.DATA_TYPE(0))
        b = self.create_array(1, [self.N], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(L, x, b)

        if self.POLYBENCH_FLATTEN_LISTS:
            # Start instruments
            self.start_instruments()
            # Run kernel
            self.kernel_flat(L, x, b)
            # Stop and print instruments
            self.stop_instruments()
        else:
            # Start instruments
            self.start_instruments()
            # Run kernel
            self.kernel(L, x, b)
            # Stop and print instruments
            self.stop_instruments()

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
        return [('x', x)]
