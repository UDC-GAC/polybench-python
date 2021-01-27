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


class Covariance(PolyBench):

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

        # Set up problem size
        self.M = params.get('M')
        self.N = params.get('N')

    def run_benchmark(self):
        # Array creation
        float_n = float(self.N)  # we will need a floating point version of N

        data = self.create_array(2, [self.N, self.M], self.DATA_TYPE(0))
        cov = self.create_array(2, [self.M, self.M], self.DATA_TYPE(0))
        mean = self.create_array(1, [self.M], self.DATA_TYPE(0))

        # Initialize array(s)
        self.initialize_array( float_n, data, cov, mean )

        # Benchmark the kernel
        self.time_kernel(float_n, data, cov, mean)

        # Return printable data as a list of tuples ('name', value)
        return [('cov', cov)]


class _StrategyList(Covariance):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, float_n:float, data: list, cov:list, mean:list):
        for i in range(0, self.N):
            for j in range(0, self.M):
                data[i][j] = self.DATA_TYPE(i * j) / self.M

    def print_array_custom(self, cov: list, name: str):
            for i in range(0, self.M):
                for j in range(0, self.M):
                    if (i * self.M + j) % 20 == 0:
                        self.print_message('\n')
                    self.print_value(cov[i][j])

    def kernel(self, float_n: float, data: list, cov: list, mean: list):
# scop begin
        for j in range(0, self.M):
            mean[j] = 0.0
            for i in range(0, self.N):
                mean[j] += data[i][j]
            mean[j] /= float_n

        for i in range(0, self.N):
            for j in range(0, self.M):
                data[i][j] -= mean[j]

        for i in range(0, self.M):
            for j in range(0, self.M):
                cov[i][j] = 0.0
                for k in range(0, self.N):
                    cov[i][j] += data[k][i] * data[k][j]
                cov[i][j] /= float_n - 1.0
                cov[j][i] = cov[i][j]
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, float_n: float, data: list, cov: list, mean: list):
# scop begin
        if((self.M-1>= 0)):
            for c1 in range ((self.M-1)+1):
                for c2 in range (c1 , (self.M-1)+1):
                    cov[c1][c2] = 0.0
            for c1 in range ((self.M-1)+1):
                mean[c1] = 0.0
            if((self.N-1>= 0)):
                for c1 in range ((self.M-1)+1):
                    for c2 in range ((self.N-1)+1):
                        mean[c1] += data[c2][c1]
            for c1 in range ((self.M-1)+1):
                mean[c1] /= float_n
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.M-1)+1):
                    data[c1][c2] -= mean[c2]
            if((self.N-1>= 0)):
                for c1 in range ((self.M-1)+1):
                    for c2 in range (c1 , (self.M-1)+1):
                        for c3 in range ((self.N-1)+1):
                            cov[c1][c2] += data[c3][c1] * data[c3][c2]
            for c1 in range ((self.M-1)+1):
                for c2 in range (c1 , (self.M-1)+1):
                    cov[c1][c2] /= (float_n - 1.0)
                    cov[c2][c1] = cov[c1][c2]
# scop end

class _StrategyListFlattened(Covariance):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, float_n:float, data: list, cov:list, mean:list):
        for i in range(0, self.N):
            for j in range(0, self.M):
                data[self.M * i + j] = self.DATA_TYPE(i * j) / self.M

    def print_array_custom(self, cov: list, name: str):
        for i in range(0, self.M):
            for j in range(0, self.M):
                if (i * self.M + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(cov[self.M * i + j])

    def kernel(self, float_n: float, data: list, cov: list, mean: list):
# scop begin
        # mean calculation: reduction + division
        for j in range(0, self.M):
            mean[j] = 0.0
            for i in range(0, self.N):
                mean[j] += data[self.M * i + j]
            mean[j] /= float_n

        # M - V
        for i in range(0, self.N):
            for j in range(0, self.M):
                data[self.M * i + j] -= mean[j]

        # element-by-element 
        # 2xcolumn traversal + reduction + transpotision
        for i in range(0, self.M):
            for j in range(0, self.M):
                cov[self.M * i + j] = 0.0
                for k in range(0, self.N):
                    cov[self.M * i + j] += data[self.M * k + i] * data[self.M * k + j]
                cov[self.M * i + j] /= float_n - 1.0
                cov[self.M * j + i] = cov[self.M * i + j]
# scop end


class _StrategyNumPy(Covariance):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, float_n:float, data: list, cov:list, mean:list):
        for i in range(0, self.N):
            for j in range(0, self.M):
                data[i, j] = self.DATA_TYPE(i * j) / self.M

    def print_array_custom(self, cov: ndarray, name: str):
        for i in range(0, self.M):
            for j in range(0, self.M):
                if (i * self.M + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(cov[i, j])

    def kernel(self, float_n: float, data: ndarray, cov: ndarray, mean: ndarray):
# scop begin
        mean[0:self.M] = data[0:self.N,0:self.M].sum(axis=0)
        mean[0:self.M] /= float_n

        data[0:self.N,0:self.M] -= mean[0:self.M]

        cov[0:self.M,0:self.M] = 0.0
        for i in range(0, self.M):
            cov[i,0:self.M] = (data[0:self.N,i] * data[0:self.N,0:self.N].T).sum(axis=1)
        cov[0:self.M,0:self.M] /= float_n - 1.0
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        self.kernel = getattr( self, "kernel_%s" % (options.POCC) )

    def kernel_pluto(self, float_n: float, data: list, cov: list, mean: list):
# --pluto
# scop begin
        if((self.M-1>= 0)):
            for c1 in range ((self.M-1)+1):
                for c2 in range (c1 , (self.M-1)+1):
                    cov[self.M*(c1) + c2] = 0.0
            for c1 in range ((self.M-1)+1):
                mean[c1] = 0.0
            if((self.N-1>= 0)):
                for c1 in range ((self.M-1)+1):
                    for c2 in range ((self.N-1)+1):
                        mean[c1] += data[self.M*(c2) + c1]
            for c1 in range ((self.M-1)+1):
                mean[c1] /= float_n
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.M-1)+1):
                    data[self.M*(c1) + c2] -= mean[c2]
            if((self.N-1>= 0)):
                for c1 in range ((self.M-1)+1):
                    for c2 in range (c1 , (self.M-1)+1):
                        for c3 in range ((self.N-1)+1):
                            cov[self.M*(c1) + c2] += data[self.M*(c3) + c1] * data[self.M*(c3) + c2]
            for c1 in range ((self.M-1)+1):
                for c2 in range (c1 , (self.M-1)+1):
                    cov[self.M*(c1) + c2] /= (float_n - 1.0)
                    cov[self.M*(c2) + c1] = cov[self.M*(c1) + c2]
# scop end

    def kernel_maxfuse(self, float_n: float, data: list, cov: list, mean: list):
# --pluto-fuse maxfuse
# scop begin
        if (self.M >= 1):
          if (self.N >= 1):
            for c0 in range( self.M ):
                for c3 in range( self.M ):
                    cov[(c0)*self.M + c3] = 0.0;
                mean[c0] = 0.0;
                for c3 in range( self.N ):
                    mean[c0] += data[(c3)*self.M + c0];
                mean[c0] /= float_n;
                for c3 in range( self.N ):
                    data[(c3)*self.M + c0] -= mean[c0];
                for c3 in range( c0+1 ):
                    for c4 in range( self.N ):
                        cov[(c3)*self.M + c0] += data[(c4)*self.M + c3] * data[(c4)*self.M + c0];
                for c3 in range( c0+1 ):
                    cov[(c3)*self.M + c0] /= (float_n - 1.0);
                    cov[(c0)*self.M + c3] = cov[(c3)*self.M + c0];

          if (self.N <= 0):
            for c0 in range( self.M ):
                for c3 in range( self.M ):
                    cov[(c0)*self.M + c3] = 0.0;
                mean[c0] = 0.0;
                mean[c0] /= float_n;
                for c3 in range( c0+1 ):
                    cov[(c3)*self.M + c0] /= (float_n - 1.0);
                    cov[(c0)*self.M + c3] = cov[(c3)*self.M + c0];
# scop end

    def kernel_vectorizer(self, float_n: float, data: list, cov: list, mean: list):
# --pluto --pluto-prevector --vectorizer --pragmatizer
# scop begin
        if (self.M >= 1):
          ub1 = (self.M + -1);
          for c1 in range(ub1+1):
              for c2 in range( c1, self.M ):
                  cov[(c1)*self.M + c2] = 0.0;
          for c1 in range( self.M ):
              mean[c1] = 0.0;
          if (self.N >= 1):
              for c2 in range( self.N ):
                for c1 in range( self.M ):
                    mean[c1] += data[(c2)*self.M + c1];
          for c1 in range( self.M ):
              mean[c1] /= float_n;
          ub1 = (self.N + -1);
          for c1 in range( ub1+1 ):
              for c2 in range( self.M ):
                  data[(c1)*self.M + c2] -= mean[c2];
          if (self.N >= 1):
            ub1 = (self.M + -1);
            for c1 in range( ub1+1 ):
                  for c3 in range( self.N ):
                    for c2 in range( c1, self.M ):
                        cov[(c1)*self.M + c2] += data[(c3)*self.M + c1] * data[(c3)*self.M + c2];
          ub1 = (self.M + -1);
          for c1 in range( ub1+1 ):
              for c2 in range( c1, self.M ):
                  cov[(c1)*self.M + c2] /= (float_n - 1.0);
                  cov[(c2)*self.M + c1] = cov[(c1)*self.M + c2];
# scop end
