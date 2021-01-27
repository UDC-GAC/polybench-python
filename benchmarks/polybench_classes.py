#
# Copyright 2021 Universidade da Coruña
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


from sys import stderr
from enum import Enum, auto


class _CustomDict(dict):
    """This class implements a Python dict in order to provide dict-like attribute access to inheriting subclasses."""

    def __init__(self):
        super(_CustomDict, self).__init__()
        self.__dict__ = self

    #
    # The following group of methods allow to expose the elements of the __dict__ as class fields.
    #

    # def __getattr__(self, item):
    #     if item in self:
    #         return self[item]
    #     else:
    #         raise AttributeError(f"No such attribute: {item}")
    #
    # def __setattr__(self, key, value):
    #     self[key] = value
    #
    # def __delattr__(self, item):
    #     if item in self:
    #         del self[item]
    #     else:
    #         raise AttributeError(f"No such attribute: {item}")


class DataSetSize(Enum):
    """Define the possible values for selecting DataSetSize sizes.

    Instead of manually managing the values of this enumeration we let the Python interpreter initialize them.
    """
    MINI = auto()
    SMALL = auto()
    MEDIUM = auto()
    LARGE = auto()
    EXTRA_LARGE = auto()


class ArrayImplementation(Enum):
    """Defines the possible values for selecting array implementations."""
    LIST = auto()
    LIST_FLATTENED = auto()
    NUMPY = auto()
    LIST_PLUTO = auto()
    LIST_FLATTENED_PLUTO = auto()


class PolyBenchOptions(_CustomDict):
    """Defines all of the available PolyBench options for PolyBench/Python and initializes them to proper defaults.

    This class inherits from _CustomDict in order to allow dict-like attribute access."""

    def __init__(self):
        super(PolyBenchOptions, self).__init__()

        # Typical options
        self.POLYBENCH_TIME = False             # Print out execution time
        self.POLYBENCH_DUMP_ARRAYS = False      # Dump live-out arrays

        # Options that may lead to better performance
        self.POLYBENCH_PADDING_FACTOR = 0       # Pad all dimensions of arrays by this value

        # Timing/profiling options
        self.POLYBENCH_PAPI = False                     # Turn on PAPI timing
        self.POLYBENCH_CACHE_SIZE_KB = 32770            # Cache size to flush, in KiB (32+ MiB)
        self.POLYBENCH_NO_FLUSH_CACHE = False           # Don't flush the cache before calling the timer
        self.POLYBENCH_CYCLE_ACCURATE_TIMER = False     # Use Time Stamp Counter
        self.POLYBENCH_LINUX_FIFO_SCHEDULER = False     # Use FIFO scheduler (must run as root)

        # Other options (not present in the README file)
        self.POLYBENCH_DUMP_TARGET = stderr     # Dump user messages into stderr, as in Polybench/C
        self.POLYBENCH_GFLOPS = False           # Unused/not implemented
        self.POLYBENCH_PAPI_VERBOSE = False     # When printing PAPI values include a descriptive name

        # Custom definitions
        # Custom option defining the problem size. The value comes from the commandline option --dataset-size and its
        # possible values are the same as in PolyBench/C:
        #   MINI_DATASET, SMALL_DATASET, MEDIUM_DATASET, LARGE_DATASET and EXTRALARGE_DATASET.
        self.POLYBENCH_DATASET_SIZE = DataSetSize.LARGE

        # PolyBench/Python options
        self.POLYBENCH_ARRAY_IMPLEMENTATION = ArrayImplementation.LIST  # Dictates the underlying array implementation


class PolyBenchSpec(_CustomDict):
    """This class stores the parameters from the polybench.spec file for a given benchmark.

    This class inherits from _CustomDict in order to allow dict-like attribute access."""

    def __init__(self, parameters: dict):
        super(PolyBenchSpec, self).__init__()

        """Process the parameters dictionary and store its values on public class fields."""
        self.Name = parameters['kernel']
        self.Category = parameters['category']

        if parameters['datatype'] == 'float' or parameters['datatype'] == 'double':
            self.DataType = float
        else:
            self.DataType = int

        mini_dict = {}
        small_dict = {}
        medium_dict = {}
        large_dict = {}
        extra_large_dict = {}
        for i in range(0, len(parameters['params'])):
            mini_dict[parameters['params'][i]] = parameters['MINI'][i]
            small_dict[parameters['params'][i]] = parameters['SMALL'][i]
            medium_dict[parameters['params'][i]] = parameters['MEDIUM'][i]
            large_dict[parameters['params'][i]] = parameters['LARGE'][i]
            extra_large_dict[parameters['params'][i]] = parameters['EXTRALARGE'][i]

        self.DataSets = {
            DataSetSize.MINI: mini_dict,
            DataSetSize.SMALL: small_dict,
            DataSetSize.MEDIUM: medium_dict,
            DataSetSize.LARGE: large_dict,
            DataSetSize.EXTRA_LARGE: extra_large_dict
        }


class PolyBenchSpecFile:
    """A .spec file contains a table, each row representing a benchmark and the columns represent different aspects of
    the benchmark (name, category, data type, problem sizes, etc.).

    This class allows to parse the contents of a PolyBench .spec file and store it in memory as a list of
    PolyBenchParameters object."""

    def __init__(self, spec_file_name: str = 'polybench.spec'):
        self.specs = []

        # Parse the passed file as if it is the file "polybench.spec".
        with open(spec_file_name) as spec_file:
            # The spec file is expected to be relatively small in size (maybe some KiB).
            # Process it line by line.
            spec_file.readline()  # skip header line
            for line in spec_file:
                dictionary = {}
                elements = line.split('\t')
                dictionary['kernel'] = elements[0]
                dictionary['category'] = elements[1]
                dictionary['datatype'] = elements[2]
                dictionary['params'] = elements[3].split(' ')
                not_numbers = elements[4].split(' ')
                numbers = []
                for nn in not_numbers:
                    numbers.append(int(nn))
                dictionary['MINI'] = numbers

                not_numbers = elements[5].split(' ')
                numbers = []
                for nn in not_numbers:
                    numbers.append(int(nn))
                dictionary['SMALL'] = numbers

                not_numbers = elements[6].split(' ')
                numbers = []
                for nn in not_numbers:
                    numbers.append(int(nn))
                dictionary['MEDIUM'] = numbers

                not_numbers = elements[7].split(' ')
                numbers = []
                for nn in not_numbers:
                    numbers.append(int(nn))
                dictionary['LARGE'] = numbers

                not_numbers = elements[8].split(' ')
                numbers = []
                for nn in not_numbers:
                    numbers.append(int(nn))
                dictionary['EXTRALARGE'] = numbers

                spec = PolyBenchSpec(dictionary)
                self.specs.append(spec)
