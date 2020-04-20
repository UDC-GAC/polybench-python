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

"""Polybench/Python is the reimplementation of the polyhedral benchmark Polybench/C in the Python programming language.

This module implements a main program which allows the user to run benchmarks easily without the burden of creating
makefiles or using configuration scripts for generating and using those.
This program allows to use only the Python runtime for everything the user should need when evaluating or implementing
the different kernels."""


# Import the basic elements for searching kernel implementations
from kernels import kernel_classes
from kernels.polybench import Polybench

# Using argparse for parsing commandline options. See: https://docs.python.org/3.7/library/argparse.html
import argparse
import logging
import os


if __name__ == '__main__':

    main_logger = logging.getLogger(__name__)
    if not (main_logger is None):  # I just want to initialize the logger in a different scope
        log_handler_console = logging.StreamHandler()
        log_formatter_console = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

        log_handler_console.setLevel(logging.DEBUG)
        log_handler_console.setFormatter(log_formatter_console)

        main_logger.addHandler(log_handler_console)

    def check_benchmark_availability() -> None:
        """
        Checks whether there are any benchmarks available or not.
        :return: None.
        :raise: NotImplementedError when there are no benchmarks available.
        """
        if len(kernel_classes) < 1:
            raise NotImplementedError("There are no available benchmarks to run.")


    def print_available_benchmarks() -> None:
        """
        Prints on screen the available benchmarks (if any)
        :return: information on screen (commandline)
        """
        check_benchmark_availability()
        print('List of available benchmarks:')
        for impl in kernel_classes:
            print(f'  {impl.__module__.replace(".", "/")}.py')


    def parse_command_line() -> {
        'benchmark': str
    }:
        """
        Parse command line arguments and generate normalized results.

        :return: A dictionary with the decoded parameters.
        :rtype: dict[str, str]

        The possible return values for the dictionary are listed below:
            - **benchmark**: the module name implementing the user selected benchmark.
        """
        parser = argparse.ArgumentParser(description='Runs a given benchmark without setting up a shell environment at '
                                                     'all.')
        parser.add_argument('benchmark', metavar='benchmark.py', nargs='?', default=None,
                            help='The path, relative to this script, to any file having a class implementing Polybench.'
                                 ' All implementations must reside somewhere inside the "kernels" folder.')
        parser.add_argument('-v', dest='verbose', action='count', default=0,
                            help='Prints more information on screen. The amount of information is determined by the '
                                 'number of "v" character appearance and is directly related to logging levels. The '
                                 'default behavior (without this flag set) is to just log error messages. Setting this '
                                 'flag once rises the logging to WARNING level, "-vv" sets INFO level and "-vvv" sets '
                                 'DEBUG level. Please note that this option does not affect the behavior of other '
                                 'messages which may be printed when setting other options (for instance, printing the '
                                 'algorithm results).')
        # Parse the commandline arguments. This process will fail on error
        args = parser.parse_args()

        # Process verbosity first, as this affects logging
        if args.verbose == 0:    # show errors (or worse) only
            main_logger.setLevel(logging.ERROR)
        elif args.verbose == 1:    # show warnings
            main_logger.setLevel(logging.WARNING)
        elif args.verbose == 2:  # show info
            main_logger.setLevel(logging.INFO)
        else:  # show all logs
            main_logger.setLevel(logging.DEBUG)

        # Process the "benchmark" argument
        if args.benchmark is None:
            print_available_benchmarks()
            exit(-1)

        result = {}

        # Blindly replace the directory separator character with a commonly supported forward slash.
        # Reason: user may input the benchmark file by using tab-completion from a command line. On Windows,
        # tab-completion only works with backslashes while on Linux, MacOS and BSDs it works with forward slashes.
        # We are not taking into consideration other systems where the forward slash separator may not work.
        benchmark_py = args.benchmark.replace(os.sep, '/')
        # Remove the "py" extension and change slashes into dots.
        # What we are left with should be compatible with Class.__module__
        result['benchmark'] = benchmark_py.split('.')[0].replace('/', '.')

        return result


    def run(module_name: str) -> None:
        # Search the module within available implementations
        instance = None
        for implementation in kernel_classes:
            if implementation.__module__ == module_name:
                # Module found! Instantiate a new class with it
                if main_logger.isEnabledFor(logging.INFO):
                    main_logger.info('Module "%s" found. Implementor class is "%s"',
                                     implementation.__module__, implementation.__name__)
                instance = implementation()

        # Check if the module was found
        if instance is None:
            module = module_name.replace(".", "/") + '.py'
            raise NotImplementedError(f'Module {module} not implemented.')

        # When the module was found, run it.
        if isinstance(instance, Polybench):
            if main_logger.isEnabledFor(logging.INFO):
                main_logger.info(f'Running "{instance.__class__.__name__}"')
            instance.run()


    # Parse the command line arguments first. We need at least one mandatory parameter.
    parameters = parse_command_line()

    benchmark = parameters['benchmark']
    run(benchmark)
