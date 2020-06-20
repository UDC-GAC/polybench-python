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


class Syr2k(PolyBench):

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
        self.M = params.get('M')
        self.N = params.get('N')

    def initialize_array(self, C: list, A: list, B: list):
        for i in range(0, self.N):
            for j in range(0, self.M):
                A[i][j] = self.DATA_TYPE((i * j+1) % self.N) / self.N;
                B[i][j] = self.DATA_TYPE((i * j+2) % self.M) / self.M;

        for i in range(0, self.N):
            for j in range(0, self.N):
                C[i][j] = self.DATA_TYPE((i * j+3) % self.N) / self.M;

    def print_array_custom(self, C: list):
        for i in range(0, self.N):
            for j in range(0, self.N):
                if (i * self.N + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(C[i][j])

    def kernel(self, alpha, beta, C: list, A: list, B: list):
        # BLAS PARAMS
        # UPLO  = 'L'
        # TRANS = 'N'
        # A is NxM
        # B is NxM
        # C is NxN
# scop begin
        for i in range(0, self.N):
            for j in range(0, i + 1):
                C[i][j] *= beta

            for k in range(0, self.M):
                for j in range(0, i + 1):
                    C[i][j] += A[j][k]*alpha*B[i][k] + B[j][k]*alpha*A[i][k]
# scop end

    def run_benchmark(self):
        # Create data structures (arrays, auxiliary variables, etc.)
        alpha = 1.5
        beta = 1.2

        C = self.create_array(2, [self.N, self.N], self.DATA_TYPE(0))
        A = self.create_array(2, [self.N, self.M], self.DATA_TYPE(0))
        B = self.create_array(2, [self.N, self.M], self.DATA_TYPE(0))

        # Initialize data structures
        self.initialize_array(C, A, B)

        # Start instruments
        self.start_instruments()

        # Run kernel
        self.kernel(alpha, beta, C, A, B)

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
        return [('C', C)]