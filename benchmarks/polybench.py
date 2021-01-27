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


"""This module offers the base Polybench class for implementing benchmarks."""

# Check if the current platform is Linux, otherwise running PolyBench makes no sense at all.
# Currently PolyBench/Python can only run on Linux due to some runtime dependencies.
from platform import system
if system() != 'Linux':
    raise NotImplementedError('PolyBench/Python is only available for Linux.')

#
# PolyBench requirements
#
from abc import abstractmethod
from benchmarks.polybench_classes import ArrayImplementation, DataSetSize
from benchmarks.polybench_classes import PolyBenchOptions, PolyBenchSpec

#
# Standard types and methods
#
from time import time
import os  # For controlling Linux scheduler

#
# External libraries
#
# NumPy
import numpy
# python_papi
from pypapi import events as papi_events, papi_high
# inlineasm
from inlineasm import assemble
from ctypes import c_ulonglong

#
# Workarounds
#
from platform import python_implementation  # Used to determine the interpreter


class PolyBench:
    """This class offers common methods for building new benchmarks.

    This class is not meant to be instantiated as it is an abstract class. Each benchmark must inherit from this class
    and override the following methods:
        - __init__()
        - kernel()
        - print_array_custom()
    """

    # Timing counters
    __polybench_program_total_flops = 0
    __polybench_timer_start = 0
    __polybench_timer_stop = 0

    # PAPI counters
    __papi_available_counters = []
    __papi_counters = []
    __papi_counters_result = []

    DATASET_SIZE = DataSetSize.LARGE  # The default dataset size for selecting bounds
    DATA_TYPE = int  # The data type used for the current benchmark (used for conversions and formatting)
    DATA_PRINT_MODIFIER = '{:d} '  # A default print modifier. Should be set up in run()

    def __init__(self, options: PolyBenchOptions, parameters: PolyBenchSpec):
        """Class constructor.

        Since this is an abstract class, this method prevents its instantiation by throwing a RuntimeError.
        This method **MUST** be overridden by subclasses.
        This method may be called from subclasses during their initialization, so we have to check whether the __init__
        call came from a subclass or not.
        """

        # Check whether __init__ is being called from a subclass.
        # The first check is for preventing the issubclass() call from returning True when directly instantiating
        # PolyBench. As the documentation states, issubclass(X, X) -> True
        if self.__class__ != PolyBench and issubclass(self.__class__, PolyBench):
            #
            # Validate inputs
            #
            if not isinstance(parameters, PolyBenchSpec):
                raise AssertionError(f'Invalid parameter "parameters": "{parameters}"')

            #
            # Set up benchmark parameters
            #
            # ... Adjust the data type and print modifier according to the data type
            self.DATA_TYPE = parameters.DataType
            # ... Adjust the print modifier to the data type
            if self.DATA_TYPE == int:
                self.DATA_PRINT_MODIFIER = '{:d} '
            elif self.DATA_TYPE == float:
                self.DATA_PRINT_MODIFIER = '{:0.2f} '
            else:
                raise NotImplementedError(f'Unknown print modifier for type {self.DATA_TYPE}')

            #
            # Set up PolyBench options
            #
            # The options dictionary is expected to have all possible options. Blindly assign values.
            # Typical options
            self.POLYBENCH_TIME = options.POLYBENCH_TIME
            self.POLYBENCH_DUMP_ARRAYS = options.POLYBENCH_DUMP_ARRAYS

            # Options that may lead to better performance
            self.POLYBENCH_PADDING_FACTOR = options.POLYBENCH_PADDING_FACTOR

            # Timing/profiling options
            self.POLYBENCH_PAPI = options.POLYBENCH_PAPI
            self.POLYBENCH_CACHE_SIZE_KB = options.POLYBENCH_CACHE_SIZE_KB
            self.POLYBENCH_NO_FLUSH_CACHE = options.POLYBENCH_NO_FLUSH_CACHE
            self.POLYBENCH_CYCLE_ACCURATE_TIMER = options.POLYBENCH_CYCLE_ACCURATE_TIMER
            self.POLYBENCH_LINUX_FIFO_SCHEDULER = options.POLYBENCH_LINUX_FIFO_SCHEDULER

            # Other options (not present in the README file)
            self.POLYBENCH_DUMP_TARGET = options.POLYBENCH_DUMP_TARGET
            self.POLYBENCH_GFLOPS = options.POLYBENCH_GFLOPS
            self.POLYBENCH_PAPI_VERBOSE = options.POLYBENCH_PAPI_VERBOSE

            # Redefine POLYBENCH_DATASET_SIZE for use in benchmarks as DATASET_SIZE
            self.DATASET_SIZE = options.POLYBENCH_DATASET_SIZE

            # PolyBench/Python options
            self.POLYBENCH_ARRAY_IMPLEMENTATION = options.POLYBENCH_ARRAY_IMPLEMENTATION

            #
            # Define in-line C functions for interpreters different than CPython
            #
            if python_implementation() != 'CPython':
                from inline import c
                # Linux scheduler code snippets taken from PolyBench/C
                linux_shedulers = c('''
                    #include <sched.h>
                    void polybench_linux_fifo_scheduler() {
                        struct sched_param schedParam;
                        schedParam.sched_priority = sched_get_priority_max (SCHED_FIFO);
                        sched_setscheduler (0, SCHED_FIFO, &schedParam);
                    }
                    void polybench_linux_standard_scheduler() {
                        struct sched_param schedParam;
                        schedParam.sched_priority = sched_get_priority_max (SCHED_OTHER);
                        sched_setscheduler (0, SCHED_OTHER, &schedParam);
                    }
                ''')
                self.__native_linux_fifo_scheduler = linux_shedulers.polybench_linux_fifo_scheduler
                self.__native_linux_standard_scheduler = linux_shedulers.polybench_linux_standard_scheduler

            #
            # Define the inline-assembly function _read_tsc()
            #
            asm_code = """
                 bits 64
                 RDTSC
                 sal     rdx, 32
                 mov     eax, eax
                 or      rax, rdx
                 ret
            """
            self._read_tsc = assemble(asm_code, c_ulonglong)
        else:
            raise RuntimeError('Abstract classes cannot be instantiated.')

    @abstractmethod
    def initialize_array(self, *args, **kwargs):
        """Implements the array initialization procedure.

        Implement this method when requiring a special array initialization.
        """
        raise NotImplementedError('Initialize array not implemented')

    @abstractmethod
    def print_array_custom(self, array: list, dump_message: str = ''):
        """Prints the benchmark array using the same format as in Polybench/C.

        If a different format is to be used, this method can be overridden.
        """
        raise NotImplementedError('Custom array print not implemented')

    @abstractmethod
    def kernel(self, *args, **kwargs):
        """Implements the kernel to be benchmarked."""
        raise NotImplementedError('Kernel not implemented')

    @abstractmethod
    def run_benchmark(self):
        """Implements the code required for running a specific benchmark.

        This method **MUST** be overridden by subclasses.

        This method should declare the data structures required for running the kernel in a similar manner as done in
        the main() function of Polybench/C.

        :returns: a list of tuples. For each tuple (X, Y), X represents the name of the output and Y is the actual
         output. The output(s) will be used by the print_array() method when required.
        :rtype: list[tuple]
        """
        raise NotImplementedError('Kernel not implemented')

    def __create_array_rec(self, dimensions: int, sizes: list, initialization_value: int = 0) -> list:
        """Auxiliary recursive method for creating a new array based upon Python lists.

        This method assumes that the parameters were previously validated (in the create_array method).

        :param int dimensions: the number of dimensions to create. One dimension creates a list, two a matrix and so on.
        :param list[int] sizes: a list of integers, each one representing the size of a dimension. The first element of
            the list represents the size of the first dimension, the second element the size of the second dimension and
            so on. If this list is smaller than the actual number of dimensions then the last size read is used for the
            remaining dimensions.
        :param int initialization_value: (optional; default = 0) the value to use for initializing the arrays during
            their creation.
        :return: a list representing an array of N dimensions.
        :rtype list:
        """
        if dimensions == 1:
            # Just create a list with as many zeros as specified in sizes[0]
            return [initialization_value for _ in range(sizes[0])]

        if len(sizes) == 1:
            # Generate lists of the same size per dimension
            return [self.__create_array_rec(dimensions - 1, sizes, initialization_value) for _ in range(sizes[0])]
        else:
            # Generate lists with unique sizes per dimension
            return [self.__create_array_rec(dimensions - 1, sizes[1:], initialization_value) for _ in range(sizes[0])]

    def create_array(self, dimensions: int, sizes: list, initialization_value: int = 0):
        """
        Create a new array with a specified size.

        :param int dimensions: specifies the number of dimensions of the array.
        :param list[int] sizes: allows to specify the number of elements for each dimension. If this parameter is a list
            with one element, this element represents the size of all dimensions, otherwise the list must have as many
            elements as specified by the dimensions parameter and each element represents the size of a dimension.
            The size of the first dimension is specified by the first element on the list, the size of the second
            dimension is represented by the second element of the list and so on.
        :param int initialization_value: (optional; default = 0) the value at which all array elements are set.
        :return: either a list representing an array of M dimensions or a NumPy array.
        """
        # Sanity check: "dimensions" must be of type integer.
        if not isinstance(dimensions, int):
            raise AssertionError('Invalid type for parameter "dimensions". '
                                 f'Expected "{type(0)}"; received "{type(dimensions)}"')

        # Sanity check: "sizes" must be a list (of integers).
        if not isinstance(sizes, list):
            raise AssertionError('Invalid type for parameter "sizes". '
                                 f'Expected "{type([])}"; received "{type(sizes)}"')

        # NOTE: type checking for "initialization_value" may not be interesting at all...
        # Sanity check: "initialization_value" must be of type integer.
        if not isinstance(initialization_value, (int, float, complex)):
            raise AssertionError('Invalid type for parameter "initialization_value". '
                                 f'Expected one of "[{type(0)}, {type(0.0)}, {type(0j)}]"; '
                                 f'received "{type(initialization_value)}"')

        # Sanity check: "dimensions" must be a non-zero positive integer.
        if dimensions < 1:
            raise AssertionError('Invalid value for parameter "dimensions". '
                                 f'Expected "dimensions > 0"; received "{dimensions}"')

        # Sanity check: "sizes" must not be an empty list.
        if not sizes:
            raise AssertionError('Invalid value for parameter "sizes". '
                                 f'Expected "non-empty list"; received "{sizes}"')

        # The following sanity checks will use Python's list comprehension for checking conditions over the "sizes" list
        # and returning non-empty lists on error with the offending elements.

        # Sanity check: "sizes" must only contain "int" type values.
        not_ints = [val for val in sizes if not isinstance(val, int)]
        if not_ints:
            raise AssertionError('Invalid value for parameter "sizes". '
                                 f'Expected "list of integer"; received "{not_ints}"')

        # Sanity check: "sizes" must only contain positive integer numbers, excluding zero.
        not_positives = [val for val in sizes if val < 1]
        if not_positives:
            raise AssertionError('Invalid value for parameter "sizes". '
                                 f'Expected "list of positive integer"; received "{not_positives}"')

        # Add post-padding to every array dimension
        new_sizes = [size + self.POLYBENCH_PADDING_FACTOR for size in sizes]

        # Expand the new_sizes list to match the number of dimensions
        while len(new_sizes) < dimensions:
            new_sizes.append(sizes[-1])

        # At this point it is safe to say that both dimensions and sizes are valid.
        # Use the appropriate "array" implementation.
        if self.POLYBENCH_ARRAY_IMPLEMENTATION == ArrayImplementation.LIST or self.POLYBENCH_ARRAY_IMPLEMENTATION == ArrayImplementation.LIST_PLUTO:
            return self.__create_array_rec(dimensions, new_sizes, initialization_value)
        elif self.POLYBENCH_ARRAY_IMPLEMENTATION == ArrayImplementation.LIST_FLATTENED or self.POLYBENCH_ARRAY_IMPLEMENTATION == ArrayImplementation.LIST_FLATTENED_PLUTO:
            # A flattened list only has one dimension, whose value is the product of all dimensions.
            dimension_size = 1
            for dim_size in new_sizes:
                dimension_size *= dim_size
            return self.__create_array_rec(1, [dimension_size], initialization_value)
        elif self.POLYBENCH_ARRAY_IMPLEMENTATION == ArrayImplementation.NUMPY:
            # Create an auxiliary list for creating an initialized NumPy array.
            list_array = self.__create_array_rec(dimensions, new_sizes, initialization_value)
            return numpy.array(list_array, self.DATA_TYPE)
        else:
            raise NotImplementedError(f'Unknown internal array implementation: "{self.POLYBENCH_ARRAY_IMPLEMENTATION}"')

    def print_array(self, array: list, native_style: bool = True, dump_message: str = ''):
        """
        Prints the benchmarked array.

        :param list array: the array to be printed.
        :param bool native_style: (optional; default = True) allows to switch between native Python list printing
        :param string dump_message: (optional; default = '') allows to set a custom message on begin and end dump
        (default) or a custom format defined in the print_array_custom() method.
        """
        self.print_message(f'begin dump: {dump_message}')
        if native_style:
            print(array)
        else:
            self.print_array_custom(array, dump_message)
        self.print_message(f'\nend   dump: {dump_message}\n')

    def print_message(self, *args, **kwargs):
        """
        Prints a user message into the configured output.
        This method also removes the newline from vanilla print(), so the user must use them manually.

        This method is inspired by: https://stackoverflow.com/a/14981125
        """
        print(*args, file=self.POLYBENCH_DUMP_TARGET, end='', **kwargs)

    def print_value(self, value):
        """
        Prints a data value using the configured data formatter.

        :param value: the value to be printed.
        """
        self.print_message(self.DATA_PRINT_MODIFIER.format(value))

    def run(self) -> dict:
        """Runs a benchmark and returns its results.

        **DO NOT OVERRIDE THIS METHOD UNLESS YOU KNOWN WHAT YOU ARE DOING!**

        This method is used for performing aditional tasks after the benchmark is run.
        One of the tasks is to control the print-out of benchmark data and the benchmarking results.
        This method is also responsible for returning the benchmarking results when either POLYBENCH_TIME or
        POLYBENCH_PAPI is enabled.
        """
        #
        # Run the benchmark
        #
        outputs = self.run_benchmark()

        #
        # Perform post-benchmark actions
        #
        self.__print_instruments()
        if self.POLYBENCH_DUMP_ARRAYS:
            self.print_message('==BEGIN DUMP_ARRAYS==\n')
            for out in outputs:
                self.print_array(out[1], False, out[0])
            self.print_message('==END   DUMP_ARRAYS==\n')
            self.POLYBENCH_DUMP_TARGET.flush()

        if self.POLYBENCH_TIME:
            # Return execution time
            return {"POLYBENCH_TIME": self.polybench_result}
        if self.POLYBENCH_PAPI:
            # Return PAPI counters
            return {"POLYBENCH_PAPI": self.polybench_result}
        return {}

    def time_kernel(self, *args, **kwargs):
        if self.POLYBENCH_TIME or self.POLYBENCH_GFLOPS:
            # Simple time measurement
            self.__timer_start()
            self.kernel(*args, **kwargs)
            self.__timer_stop()
        elif self.POLYBENCH_PAPI:
            # Measuring performance counters is a bit tricky. The API allows to monitor multiple counters at once, but
            # that is not accurate so we need to measure each counter independently within a loop to ensure proper
            # operation.
            i = 0
            self.__papi_init()  # Initializes self.__papi_counters and self.__papi_available_counters
            self.__prepare_instruments()
            # Information for the following loop:
            # * self.__papi_counters holds a list of available counter ids
            # * self.__papi_counters_result holds the actual counter return values
            for counter in self.__papi_counters:
                if i > 0:
                    self.initialize_array(*args, **kwargs)  # force initialization
                i += 1
                papi_high.start_counters([counter])  # requires a list of counters
                self.kernel(*args, **kwargs)
                self.__papi_counters_result.extend(papi_high.stop_counters())  # returns a list of counter results
        else:
            # Default kernel run
            self.__prepare_instruments()
            self.kernel(*args, **kwargs)

        # Something like stop_instruments()
        if self.POLYBENCH_LINUX_FIFO_SCHEDULER:
            self.__linux_standard_scheduler()

    def __print_instruments(self):
        """Print the state of the instruments."""
        if self.POLYBENCH_TIME or self.POLYBENCH_GFLOPS:
            self.__timer_print()
        elif self.POLYBENCH_PAPI:
            self.__papi_print()

    def __prepare_instruments(self):
        if not self.POLYBENCH_NO_FLUSH_CACHE:
            self.__flush_cache()
        if self.POLYBENCH_LINUX_FIFO_SCHEDULER:
            self.__linux_fifo_scheduler()

    def __timer_start(self):
        if not self.POLYBENCH_CYCLE_ACCURATE_TIMER:
            self.__prepare_instruments()
            self.__timer_start_t = time()
        else:
            self.__prepare_instruments()
            self.__timer_start_t = self._read_tsc()

    def __timer_stop(self):
        if not self.POLYBENCH_CYCLE_ACCURATE_TIMER:
            self.__timer_stop_t = time()
        else:
            self.__timer_stop_t = self._read_tsc()

    def __timer_print(self):
        self.polybench_result = self.__timer_stop_t - self.__timer_start_t
        if not self.POLYBENCH_CYCLE_ACCURATE_TIMER:
            print(f'{self.polybench_result:0.6f}')
        else:
            print(f'{self.polybench_result:d}')

    def __papi_init(self):
        """
        Performs initialization of the PAPI library.

        Since the library python-papi only exports standard events, it is necessary to know which ones are available and
        which ones are requested by the user to perform a validation phase.
        :return: None. Modifies self.__papi_counters.
        """
        def get_available_counters() -> list:
            """
            Gets a list of available counters. Each element of the list is a tuple containing the name of the event and
            the numerical value corresponding to the event. This value is masked.
            :return: A list with the available benchmarks
            :rtype: list[tuple]
            """
            def is_number(x):
                return isinstance(x, int) or isinstance(x, float) or isinstance(x, complex)
            # See: https://stackoverflow.com/a/9794849
            from inspect import getmembers
            return getmembers(papi_events, is_number)

        def parse_counters_file() -> list:
            result = []
            with open('papi_counters.list') as f:
                contents = f.read()
                # Remove both empty lines and whitespaces
                lines = contents.splitlines()
                lines = [line.strip() for line in lines]

                is_in_comment = False
                for line in lines:
                    if not is_in_comment:
                        if line.startswith('/*'):
                            is_in_comment = True
                            continue
                        elif line.startswith('//'):
                            continue
                        else:
                            result.append(line.strip('",'))  # store plain counter names
                    else:
                        if line.endswith('*/'):
                            is_in_comment = False
            return result

        self.__papi_counters.clear()
        self.__papi_counters_result.clear()
        self.__papi_available_counters = get_available_counters()
        user_counters = parse_counters_file()

        # Check if there are any user supplied counters after reading the file
        if len(user_counters) < 1:
            raise NotImplementedError('Must supply at least one counter name in file "papi_counters.list"')

        # Check if the user counters exist within the available standard set.
        for usr_counter in user_counters:
            is_available = False
            for avlbl_counter in self.__papi_available_counters:
                if usr_counter == avlbl_counter[0]:
                    self.__papi_counters.append(avlbl_counter[1])
                    is_available = True
                    break
            if not is_available:
                if usr_counter.startswith( "0x" ): self.__papi_counters.append( int(usr_counter,16) )
                else: self.__papi_counters.append( int(usr_counter) )
#                else: print(f'WARNING: counter "{usr_counter}" not available.')

    def __papi_print(self):
        def papi_counter_names():
            """
            Translate back the PAPI identifiers into user-readable ones.
            :return: a list of PAPI identifiers
            """
            result = []
            for usr_counter in self.__papi_counters:
                for avlbl_counter in self.__papi_available_counters:
                    if usr_counter == avlbl_counter[1]:
                        result.append(avlbl_counter[0])
                        break
                else:
                    result.append( "native: " + hex(usr_counter) )
            return result

        self.polybench_result = {}
        counter_names = papi_counter_names()
        for i in range(0, len(self.__papi_counters)):
            if self.POLYBENCH_PAPI_VERBOSE:
                print(f'{counter_names[i]}=', end='')
            print(f'{self.__papi_counters_result[i]} ', end='')
            if self.POLYBENCH_PAPI_VERBOSE:
                print()  # new line
            # Append key-value to result (name-value)
            self.polybench_result[counter_names[i]] = self.__papi_counters_result[i]
        print()  # new line

    def __flush_cache(self):
        """Thrashes the cache by generating a very large data structure."""
        cs = int(self.POLYBENCH_CACHE_SIZE_KB * 1024 / 8)  # divided by sizeof(double)
        flush = [0.0 for _ in range(cs)]  # 0.0 forces data type to be a float
        tmp = 0.0
        for i in range(cs):
            tmp += flush[i]
        assert tmp <= 10.0

    def __linux_fifo_scheduler(self):
        if python_implementation() == 'CPython':
            param = os.sched_param(os.SCHED_FIFO)
            os.sched_setscheduler(0, os.SCHED_FIFO, param)
        else:
            self.__native_linux_fifo_scheduler()

    def __linux_standard_scheduler(self):
        if python_implementation() == 'CPython':
            param = os.sched_param(os.SCHED_OTHER)
            os.sched_setscheduler(0, os.SCHED_OTHER, param)
        else:
            self.__native_linux_standard_scheduler()
