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


"""Implements the correlation kernel in a PolyBench class."""

from benchmarks.polybench import PolyBench
from benchmarks.polybench_classes import ArrayImplementation
from benchmarks.polybench_classes import PolyBenchOptions, PolyBenchSpec
from numpy.core.multiarray import ndarray
from math import sqrt
import numpy as np


class Correlation(PolyBench):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        implementation = options.POLYBENCH_ARRAY_IMPLEMENTATION
        if implementation == ArrayImplementation.LIST:
            return _StrategyList.__new__(_StrategyList, options, parameters)
        elif implementation == ArrayImplementation.LIST_FLATTENED:
            return _StrategyListFlattened.__new__(_StrategyListFlattened, options, parameters)
        elif implementation == ArrayImplementation.NUMPY:
            return _StrategyNumPy.__new__(_StrategyNumPy, options, parameters)
        elif implementation == ArrayImplementation.LIST_PLUTO:
            return _StrategyListPluto.__new__(_StrategyListPluto, options, parameters)
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
        corr = self.create_array(2, [self.M, self.M], self.DATA_TYPE(0))
        mean = self.create_array(1, [self.M], self.DATA_TYPE(0))
        stddev = self.create_array(1, [self.M], self.DATA_TYPE(0))

        # Initialize array(s)
        self.initialize_array(float_n, data, corr, mean, stddev)

        # Benchmark the kernel
        self.time_kernel(float_n, data, corr, mean, stddev)

        # Return printable data as a list of tuples ('name', value)
        return [('corr', corr)]


class _StrategyList(Correlation):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyList)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, float_n: float, data: list, corr: list, mean: list, stddev: list):
        for i in range(0, self.N):
            for j in range(0, self.M):
                data[i][j] = (self.DATA_TYPE(i * j) / self.M) + i

    def print_array_custom(self, corr: list, name: str):
        for i in range(0, self.M):
            for j in range(0, self.M):
                if (i * self.M + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(corr[i][j])

    def kernel(self, float_n: float, data: list, corr: list, mean: list, stddev: list):
        eps = 0.1

# scop begin
        for j in range(0, self.M):
            mean[j] = 0.0
            for i in range(0, self.N):
                mean[j] += data[i][j]
            mean[j] /= float_n

        for j in range(0, self.M):
            stddev[j] = 0.0
            for i in range(0, self.N):
                stddev[j] += (data[i][j] - mean[j]) * (data[i][j] - mean[j])
            stddev[j] /= float_n
            stddev[j] = sqrt(stddev[j])
            # The following in an elegant but usual way to handle near-zero std. dev. values, which below would cause a
            # zero divide.
            stddev[j] = 1.0 if stddev[j] <= eps else stddev[j]

        # Center and reduce the column vectors.
        for i in range(0, self.N):
            for j in range(0, self.M):
                data[i][j] -= mean[j]
                data[i][j] /= sqrt(float_n) * stddev[j]

        # Calculate the m*n correlation matrix.
        for i in range(0, self.M-1):
            corr[i][i] = 1.0
            for j in range(i+1, self.M):
                corr[i][j] = 0.0
                for k in range(0, self.N):
                    corr[i][j] += (data[k][i] * data[k][j])
                corr[j][i] = corr[i][j]
        corr[self.M-1][self.M-1] = 1.0
# scop end


class _StrategyListFlattened(Correlation):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattened)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, float_n: float, data: list, corr: list, mean: list, stddev: list):
        for i in range(0, self.N):
            for j in range(0, self.M):
                data[self.M * i + j] = (self.DATA_TYPE(i * j) / self.M) + i

    def print_array_custom(self, corr: list, name: str):
        for i in range(0, self.M):
            for j in range(0, self.M):
                if (i * self.M + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(corr[self.M * i + j])

    def kernel(self, float_n: float, data: list, corr: list, mean: list, stddev: list):
        eps = 0.1

# scop begin
        # mean calculation: reduction + division
        for j in range(0, self.M):
            mean[j] = 0.0
            for i in range(0, self.N):
                mean[j] += data[self.M * i + j]
            mean[j] /= float_n

        # (M - V) * (M-V) reduction
        for j in range(0, self.M):
            stddev[j] = 0.0
            for i in range(0, self.N):
                stddev[j] += (data[self.M * i + j] - mean[j]) * (data[self.M * i + j] - mean[j])
            stddev[j] /= float_n
            stddev[j] = sqrt(stddev[j])
            # The following in an elegant but usual way to handle near-zero std. dev. values, which below would cause a
            # zero divide.
            stddev[j] = 1.0 if stddev[j] <= eps else stddev[j]

        # Center and reduce the column vectors.
        # (M-V) / V
        for i in range(0, self.N):
            for j in range(0, self.M):
                data[self.M * i + j] -= mean[j]
                data[self.M * i + j] /= sqrt(float_n) * stddev[j]

        # Calculate the m*n correlation matrix.
        # 2x column traversal + reduction + transposition
        for i in range(0, self.M - 1):
            corr[self.M * i + i] = 1.0
            for j in range(i + 1, self.M):
                corr[self.M * i + j] = 0.0
                for k in range(0, self.N):
                    corr[self.M * i + j] += (data[self.M * k + i] * data[self.M * k + j])
                corr[self.M * j + i] = corr[self.M * i + j]
        corr[self.M * (self.M - 1) + (self.M - 1)] = 1.0
# scop end


class _StrategyNumPy(Correlation):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyNumPy)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

    def initialize_array(self, float_n: float, data: list, corr: list, mean: list, stddev: list):
        for i in range(0, self.N):
            for j in range(0, self.M):
                data[i, j] = (self.DATA_TYPE(i * j) / self.M) + i

    def print_array_custom(self, corr: ndarray, name: str):
        for i in range(0, self.M):
            for j in range(0, self.M):
                if (i * self.M + j) % 20 == 0:
                    self.print_message('\n')
                self.print_value(corr[i, j])

    def kernel(self, float_n: float, data: ndarray, corr: ndarray, mean: ndarray, stddev: ndarray):
        eps = 0.1

# scop begin
        mean[0:self.M] = data[0:self.N,0:self.M].sum(axis=0)
        mean[0:self.M] /= float_n

        stddev[0:self.M] = ((data[0:self.N,0:self.M]-mean[0:self.M]) * (data[0:self.N,0:self.M]-mean[0:self.M])).sum(axis=0)
        stddev[0:self.M] /= float_n
        stddev[0:self.M] = np.sqrt(stddev[0:self.M])
            # The following in an elegant but usual way to handle near-zero std. dev. values, which below would cause a
            # zero divide.
        stddev[0:self.M] = np.where( stddev[0:self.M] <= eps, 1.0, stddev[0:self.M] )

        # Center and reduce the column vectors.
        data[0:self.N, 0:self.M] -= mean[0:self.M]
        data[0:self.N, 0:self.M] /= sqrt(float_n) * stddev[0:self.M]

        # Calculate the m*n correlation matrix.
        corr[np.diag_indices(corr.shape[0])] = 1.0
        triu_indices = np.triu_indices( n=self.M, m=self.M, k=1 )
        corr[triu_indices] = 0.0
        for i in range(0, self.M - 1):
            corr[i,i+1:self.M] = (data[0:self.N,i] * data[0:self.N,i+1:self.M].T).sum(axis=1)
        tril_indices = np.tril_indices( n=self.M, m=self.M, k=-1 )
        corr[tril_indices] = corr[triu_indices]
        corr[self.M - 1, self.M - 1] = 1.0
# scop end

class _StrategyListPluto(_StrategyList):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListPluto)

    def kernel(self, float_n: float, data: list, corr: list, mean: list, stddev: list):
        eps = 0.1

# scop begin
        corr[self.M-1][self.M-1] = 1.0
        for c1 in range ((self.M-2)+1):
            for c2 in range (c1 + 1 , (self.M-1)+1):
                corr[c1][c2] = 0.0
        for c1 in range ((self.M-2)+1):
            corr[c1][c1] = 1.0
            stddev[c1] = 0.0
            mean[c1] = 0.0
        if((self.M-1>= 0)):
            stddev[self.M + -1] = 0.0
            mean[self.M + -1] = 0.0
        if((self.N-1>= 0)):
            for c1 in range ((self.M-1)+1):
                for c2 in range ((self.N-1)+1):
                    mean[c1] += data[c2][c1]
        for c1 in range ((self.M-1)+1):
            mean[c1] /= float_n
        if((self.N-1>= 0)):
            for c1 in range ((self.M-1)+1):
                for c2 in range ((self.N-1)+1):
                    stddev[c1] += (data[c2][c1] - mean[c1]) * (data[c2][c1] - mean[c1])
                    data[c2][c1] -= mean[c1]
        for c1 in range ((self.M-1)+1):
            stddev[c1] /= float_n
            stddev[c1] = sqrt(stddev[c1])
            stddev[c1] = 1.0 if stddev[c1] <= eps else stddev[c1]
        if((self.M-1>= 0)):
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.M-1)+1):
                    data[c1][c2] /= sqrt(float_n) * stddev[c2]
        if((self.N-1>= 0)):
            for c1 in range ((self.M-2)+1):
                for c2 in range (c1 + 1 , (self.M-1)+1):
                    for c3 in range ((self.N-1)+1):
                        corr[c1][c2] += (data[c3][c1] * data[c3][c2])
        for c1 in range ((self.M-2)+1):
            for c2 in range (c1 + 1 , (self.M-1)+1):
                corr[c2][c1] = corr[c1][c2]
# scop end

class _StrategyListFlattenedPluto(_StrategyListFlattened):

    def __new__(cls, options: PolyBenchOptions, parameters: PolyBenchSpec):
        return object.__new__(_StrategyListFlattenedPluto)

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        super().__init__(options, parameters)

        self.kernel = getattr( self, "kernel_%s" % (options.POCC) )

    def kernel_pluto(self, float_n: float, data: list, corr: list, mean: list, stddev: list):
        eps = 0.1
# scop begin
        corr[self.M*(self.M-1)+self.M-1] = 1.0
        for c1 in range ((self.M-2)+1):
            for c2 in range (c1 + 1 , (self.M-1)+1):
                corr[self.M*c1+c2] = 0.0
        for c1 in range ((self.M-2)+1):
            corr[self.M*c1+c1] = 1.0
            stddev[c1] = 0.0
            mean[c1] = 0.0
        if((self.M-1>= 0)):
            stddev[self.M + -1] = 0.0
            mean[self.M + -1] = 0.0
        if((self.N-1>= 0)):
            for c1 in range ((self.M-1)+1):
                for c2 in range ((self.N-1)+1):
                    mean[c1] += data[self.M*c2+c1]
        for c1 in range ((self.M-1)+1):
            mean[c1] /= float_n
        if((self.N-1>= 0)):
            for c1 in range ((self.M-1)+1):
                for c2 in range ((self.N-1)+1):
                    stddev[c1] += (data[self.M*c2+c1] - mean[c1]) * (data[self.M*c2+c1] - mean[c1])
                    data[self.M*c2+c1] -= mean[c1]
        for c1 in range ((self.M-1)+1):
            stddev[c1] /= float_n
            stddev[c1] = sqrt(stddev[c1])
            stddev[c1] = 1.0 if stddev[c1] <= eps else stddev[c1]
        if((self.M-1>= 0)):
            for c1 in range ((self.N-1)+1):
                for c2 in range ((self.M-1)+1):
                    data[self.M*c1+c2] /= sqrt(float_n) * stddev[c2]
        if((self.N-1>= 0)):
            for c1 in range ((self.M-2)+1):
                for c2 in range (c1 + 1 , (self.M-1)+1):
                    for c3 in range ((self.N-1)+1):
                        corr[self.M*c1+c2] += (data[self.M*c3+c1] * data[self.M*c3+c2])
        for c1 in range ((self.M-2)+1):
            for c2 in range (c1 + 1 , (self.M-1)+1):
                corr[self.M*c2+c1] = corr[self.M*c1+c2]
# scop end

    def kernel_maxfuse(self, float_n: float, data: list, corr: list, mean: list, stddev: list):
        eps = 0.1
# --pluto-fuse maxfuse
# scop begin
        corr[(self.M-1)*self.M + self.M-1] = 1.0;
        if ((self.M >= 2) and (self.N >= 1)):
          for c3 in range( 1, self.M ):
              corr[(0)*self.M + c3] = 0.0;
          corr[(0)*self.M + 0] = 1.0;
          stddev[0] = 0.0;
          mean[0] = 0.0;
          for c3 in range( self.N ):
              mean[0] += data[(c3)*self.M + 0];
          mean[0] /= float_n;
          for c3 in range( self.N ):
              stddev[0] += (data[(c3)*self.M + 0] - mean[0]) * (data[(c3)*self.M + 0] - mean[0]);
              data[(c3)*self.M + 0] -= mean[0];
          stddev[0] /= float_n;
          stddev[0] = sqrt(stddev[0]);
          stddev[0] = 1.0 if stddev[0] <= eps else stddev[0]
          for c3 in range( self.N ):
              data[(c3)*self.M + 0] /= sqrt(float_n) * stddev[0];

        if ((self.M == 1) and (self.N >= 1)):
          stddev[0] = 0.0;
          mean[0] = 0.0;
          for c3 in range( self.N ):
              mean[0] += data[(c3)*self.M + 0];
          mean[0] /= float_n;
          for c3 in range( self.N ):
              stddev[0] += (data[(c3)*self.M + 0] - mean[0]) * (data[(c3)*self.M + 0] - mean[0]);
              data[(c3)*self.M + 0] -= mean[0];
          stddev[0] /= float_n;
          stddev[0] = sqrt(stddev[0]);
          stddev[0] = 1.0 if stddev[0] <= eps else stddev[0]
          for c3 in range( self.N ):
              data[(c3)*self.M + 0] /= sqrt(float_n) * stddev[0];

        if ((self.M >= 2) and (self.N <= 0)):
          for c3 in range( 1, self.M ):
            corr[(0)*self.M + c3] = 0.0;
          corr[(0)*self.M + 0] = 1.0;
          stddev[0] = 0.0;
          mean[0] = 0.0;
          mean[0] /= float_n;
          stddev[0] /= float_n;
          stddev[0] = sqrt(stddev[0]);
          stddev[0] = 1.0 if sttdev[0] <= eps else stddev[0]

        if ((self.M == 1) and (self.N <= 0)):
          stddev[0] = 0.0;
          mean[0] = 0.0;
          mean[0] /= float_n;
          stddev[0] /= float_n;
          stddev[0] = sqrt(stddev[0]);
          stddev[0] = 1.0 if stddev[0] <= eps else stddev[0]
        if (self.N >= 1):
            for c0 in range( 1, self.M-1 ):
              for c3 in range( c0+1, self.M ):
                  corr[(c0)*self.M + c3] = 0.0;
              corr[(c0)*self.M + c0] = 1.0;
              stddev[c0] = 0.0;
              mean[c0] = 0.0;
              for c3 in range( self.N ):
                  mean[c0] += data[(c3)*self.M + c0];
              mean[c0] /= float_n;
              for c3 in range( self.N ):
                  stddev[c0] += (data[(c3)*self.M + c0] - mean[c0]) * (data[(c3)*self.M + c0] - mean[c0]);
                  data[(c3)*self.M + c0] -= mean[c0];
              stddev[c0] /= float_n;
              stddev[c0] = sqrt(stddev[c0]);
              stddev[c0] = 1.0 if stddev[c0] <= eps else stddev[c0]
              for c3 in range( self.N ):
                  data[(c3)*self.M + c0] /= sqrt(float_n) * stddev[c0];
              for c3 in range( c0 ):
                  for c4 in range( self.N ):
                      corr[(c3)*self.M + c0] += (data[(c4)*self.M + c3] * data[(c4)*self.M + c0]);
              for c3 in range( c0 ):
                  corr[(c0)*self.M + c3] = corr[(c3)*self.M + c0];

        if ((self.M >= 2) and (self.N >= 1)):
            stddev[self.M + -1] = 0.0;
            mean[self.M + -1] = 0.0;
            for c3 in range( self.N ):
                mean[self.M + -1] += data[(c3)*self.M + self.M + -1];
            mean[self.M + -1] /= float_n;
            for c3 in range( self.N ):
                stddev[self.M + -1] += (data[(c3)*self.M + self.M + -1] - mean[self.M + -1]) * (data[(c3)*self.M + self.M + -1] - mean[self.M + -1]);
                data[(c3)*self.M + self.M + -1] -= mean[self.M + -1];
            stddev[self.M + -1] /= float_n;
            stddev[self.M + -1] = sqrt(stddev[self.M + -1]);
            stddev[self.M + -1] = 1.0 if stddev[self.M + -1] <= eps else stddev[self.M + -1];
            for c3 in range( self.N ):
                data[(c3)*self.M + self.M + -1] /= sqrt(float_n) * stddev[self.M + -1];
            for c3 in range( self.M-1 ):
                for c4 in range( self.N ):
                    corr[(c3)*self.M + self.M + -1] += (data[(c4)*self.M + c3] * data[(c4)*self.M + self.M + -1]);
            for c3 in range( self.M-1 ):
                corr[(self.M + -1)*self.M + c3] = corr[(c3)*self.M + self.M + -1];

        if (self.N <= 0):
            for c0 in range( self.M-1 ):
              for c3 in range( c0+1, self.M ):
                  corr[(c0)*self.M + c3] = 0.0;
              corr[(c0)*self.M + c0] = 1.0;
              stddev[c0] = 0.0;
              mean[c0] = 0.0;
              mean[c0] /= float_n;
              stddev[c0] /= float_n;
              stddev[c0] = sqrt(stddev[c0]);
              stddev[c0] = 1.0 if stddev[c0] <= eps else stddev[c0];
              for c3 in range( c0 ):
                  corr[(c0)*self.M + c3] = corr[(c3)*self.M + c0];

        if ((self.M >= 2) and (self.N <= 0)):
          stddev[self.M + -1] = 0.0;
          mean[self.M + -1] = 0.0;
          mean[self.M + -1] /= float_n;
          stddev[self.M + -1] /= float_n;
          stddev[self.M + -1] = sqrt(stddev[self.M + -1]);
          stddev[self.M + -1] = 1.0 if stddev[self.M + -1] <= eps else stddev[self.M + -1];
          for c3 in range( self.M-1 ):
              corr[(self.M + -1)*self.M + c3] = corr[(c3)*self.M + self.M + -1];
# scop end

    def kernel_vectorizer(self, float_n: float, data: list, corr: list, mean: list, stddev: list):
        eps = 0.1
# --pluto --pluto-prevector --vectorizer --pragmatizer
# scop begin
        corr[(self.M-1)*self.M + self.M-1] = 1.0;
        ub1 = (self.M + -2);
        for c1 in range( ub1+1 ):
            for c2 in range( c1+1, self.M ):
                corr[(c1)*self.M + c2] = 0.0;
        for c1 in range( self.M -1 ):
            corr[(c1)*self.M + c1] = 1.0;
            stddev[c1] = 0.0;
            mean[c1] = 0.0;
        if (self.M >= 1):
          stddev[self.M + -1] = 0.0;
          mean[self.M + -1] = 0.0;
        if (self.N >= 1):
            for c2 in range( self.N ):
              for c1 in range( self.M ):
                  mean[c1] += data[(c2)*self.M + c1];
        for c1 in range( self.M ):
            mean[c1] /= float_n;
        if (self.N >= 1):
            for c2 in range( self.N ):
              for c1 in range( self.M ):
                  stddev[c1] += (data[(c2)*self.M + c1] - mean[c1]) * (data[(c2)*self.M + c1] - mean[c1]);
                  data[(c2)*self.M + c1] -= mean[c1];
        for c1 in range( self.M ):
            stddev[c1] /= float_n;
            stddev[c1] = sqrt(stddev[c1]);
            stddev[c1] = 1.0 if stddev[c1] <= eps else stddev[c1];
        if (self.M >= 1):
          ub1 = (self.N + -1);
          for c1 in range( ub1+1 ):
              for c2 in range( self.M ):
                  data[(c1)*self.M + c2] /= sqrt(float_n) * stddev[c2];
        if (self.N >= 1):
          ub1 = (self.M + -2);
          for c1 in range( ub1+1 ):
              for c3 in range( self.N ):
                  for c2 in range( c1+1, self.M ):
                      corr[(c1)*self.M + c2] += (data[(c3)*self.M + c1] * data[(c3)*self.M + c2]);
        ub1 = (self.M + -2);
        for c1 in range( ub1+1 ):
            for c2 in range( c1+1, self.M ):
                corr[(c2)*self.M + c1] = corr[(c1)*self.M + c2];
# scop end
