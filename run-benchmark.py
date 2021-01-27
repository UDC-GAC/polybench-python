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


"""Polybench/Python is the reimplementation of the polyhedral benchmark Polybench/C in the Python programming language.

This module implements a main program which allows the user to run benchmarks easily without the burden of creating
makefiles or using configuration scripts for generating and using those.
This program allows to use only the Python runtime for everything the user should need when evaluating or implementing
the different kernels."""


# Import the basic elements for searching benchmark implementations
from platform import python_implementation

from benchmarks import benchmark_classes
from benchmarks.polybench_classes import ArrayImplementation, DataSetSize
from benchmarks.polybench_classes import PolyBenchOptions, PolyBenchSpecFile

# Using argparse for parsing commandline options. See: https://docs.python.org/3.7/library/argparse.html
import argparse
import os
from pathlib import Path

from filecmp import cmp  # used when verifying if two files have the same content
from sys import stderr, stdout


if __name__ == '__main__':

    def check_benchmark_availability() -> None:
        """
        Checks whether there are any benchmarks available or not.
        :return: None.
        :raise: NotImplementedError when there are no benchmarks available.
        """
        if len(benchmark_classes) < 1:
            raise NotImplementedError("There are no available benchmarks to run.")


    def print_available_benchmarks() -> None:
        """
        Prints on screen the available benchmarks (if any)
        :return: information on screen (commandline)
        """
        check_benchmark_availability()
        print('List of available benchmarks:')
        for impl in benchmark_classes:
            print(f'  {impl.__module__.replace(".", "/")}.py')


    def parse_command_line() -> {
        'benchmark': str,
        'polybench_options': dict,
        'save_results': bool,
        'verify': dict,
    }:
        """
        Parse command line arguments and generate normalized results.

        :return: A dictionary with the decoded parameters.
        :rtype: dict[str, Any]
        """
        parser = argparse.ArgumentParser(description='Runs a given benchmark without setting up a shell environment at '
                                                     'all.')
        # Default parameter, mandatory
        parser.add_argument('benchmark', metavar='benchmark.py', nargs='?', default=None,
                            help='The path, relative to this script, to any file having a class implementing Polybench.'
                                 ' All implementations must reside somewhere inside the "benchmarks" folder. The '
                                 'keyword "all" allows to run all of the benchmarks sequentially.')
        # Optional parameters
        parser.add_argument('--polybench-options', dest='options', default=None,
                            help='A comma separated list of options passed to PolyBench. Available options can be found'
                                 ' in the README file. Usage: run.py --options '
                                 'POLYBENCH_PADDING_FACTOR=3,POLYBENCH_LINUX_FIFO_SCHEDULER')
        parser.add_argument('--dataset-size', dest='dataset_size', default=None,
                            help='Specify a working dataset size to use from "polybench.spec" file. Valid values are:'
                                 '"MINI", "SMALL", "MEDIUM", "LARGE", "EXTRALARGE".')
        parser.add_argument('--save-results', dest='save_results', action='store_true',
                            help='Saves execution results into an automatically named file next to the benchmark '
                                 'implementation.')
        parser.add_argument('--output-array', dest='output_array', default=None,
                            help='Alias for POLYBENCH_DUMP_TARGET. Also enables POLYBENCH_DUMP_ARRAYS. Prints the '
                                 'benchmark''s result into a file. In order to print into the console use either '
                                 '"stdout" or "stderr".')
        parser.add_argument('--verify-file', dest='verify_file_name', default=None,
                            help='Verify the results of the benchmark against the results stored in a file. This '
                                 'option enables --output-array and makes the output file to target a file name whose '
                                 'name matches the one passed by this argument, appending the .verify suffix. When the '
                                 'benchmark terminates, its result is compared and a message indicating the comparison '
                                 'result will be printed on "stdout".')
        parser.add_argument('--verify-polybench-path', dest='verify_polybench_path', default=None,
                            help='Combined with --verify-file, this parameter allows to specify the path where '
                                 'PolyBench/C is present to allow for automatic discovery of the appropriate output '
                                 'file. Please note that the files containing the results of PolyBench/C must already '
                                 'exist and must be next to where the actual implementation resides. This is the '
                                 'default behavior when PolyBench/C is run from Perl scripts.')
        parser.add_argument('--iterations', dest='iterations', default=1,
                            help='Performs N runs of the benchmark.')
        parser.add_argument('--array-implementation', dest='array_implementation', default=0,
                            help='Allows to select the internal array implementation in use. 0: Python List; 1: Python '
                                 'List with flattened indexes; 2: NumPy array. Default: 0.')
        parser.add_argument('--load-elimination', action="store_true", default=False,
                            help="Enables load elimination (where applicable). Default: False." )
        parser.add_argument('--pocc', dest="pocc_option", default="pluto",
                            help="Selects the PoCC optimization applied. Valid values are 'pluto', 'vectorizer', and 'maxfuse'" )
        # Parse the commandline arguments. This process will fail on error
        args = parser.parse_args()

        # Process the "benchmark" argument
        if args.benchmark is None:
            print_available_benchmarks()
            exit(-1)

        # Initialize the result dictionary
        result = {
            'benchmark': None,  # should hold a string
            'polybench_options': PolyBenchOptions(),  # Stores PolyBench specific options
            'save_results': False,  # Allows to save execution results into a file
            'verify': {
                'enabled': False,  # Controls whether to verify results or not
                'file': '',        # The file name to verify against
                'path': '',        # The path to PolyBench/C
                'full_path': ''
            },
        }

        # Blindly replace the directory separator character with a commonly supported forward slash.
        # Reason: user may input the benchmark file by using tab-completion from a command line. On Windows,
        # tab-completion only works with backslashes while on Linux, MacOS and BSDs it works with forward slashes.
        # We are not taking into consideration other systems where the forward slash separator may not work.
        benchmark_py = args.benchmark.replace(os.sep, '/')
        # Remove the "py" extension and change slashes into dots.
        # What we are left with should be compatible with Class.__module__
        result['benchmark'] = benchmark_py.split('.')[0].replace('/', '.')

        print_result = False  # By default, do not print anything

        # Define some auxiliary functions.
        # Some parameters assign values with the same meanings. This prevents typing errors from happening.
        def set_output(file_name: str):
            if file_name == 'stderr':
                handle = stderr
            elif file_name == 'stdout':
                handle = stdout
            else:
                handle = open(file_name, 'w', newline='\n')
            result['output_array'] = handle

        # Process the "output_array" argument
        if not (args.output_array is None):
            set_output(args.output_array)
            print_result = True
        else:
            set_output('stderr')  # Just setting a default value

        # Process the "verify" arguments. May alter 'output_array'
        result['verify'] = {
            'enabled': False
        }
        if not (args.verify_file_name is None) and result['benchmark'] != 'all':
            result['verify']['enabled'] = True
            result['verify']['file'] = str(args.verify_file_name)
            set_output(args.verify_file_name + '.verify')
            print_result = True

        if not (args.verify_polybench_path is None):
            result['verify']['path'] = str(args.verify_polybench_path).rstrip('/')

        # Check if the arguments passed to "verify*" are valid
        if result['verify']['enabled']:
            # Calculate the full file name (path + file)
            if result['verify']['path'] == '':
                result['verify']['full_path'] = result['verify']['file']
            else:
                # PolyBench/C path given. Search the appropriate category/benchmark
                # split_cat will contain, as a list, the category without "benchmarks/" nor "benchmark.py"
                split_cat = result['benchmark'].replace('benchmarks.', '').split('.')[:-1]
                category_name = ''
                for token in split_cat:
                    category_name += token + '/'
                result['verify']['full_path'] = result['verify']['path'] + '/' +\
                    category_name + result['verify']['file']

            file = Path(result['verify']['full_path'])
            if not file.is_file():
                # The first check failed. It may mean that the benchmark path was converted to meet Python naming
                # conventions. There are two things to revert back:
                # - Underscores at the beginning of a path token must be removed
                # - Underscores in the middle of a path token may need to be converted into hyphens
                tokenized_path = result['verify']['full_path'].split('/')
                validated_tokenized_path = []
                for token in tokenized_path:
                    if len(token) > 0:
                        fixed_token = token.lstrip('_')
                        fixed_token = fixed_token.replace('_', '-')
                        validated_tokenized_path.append(fixed_token)
                    else:
                        # Probably the path starts with "/". Example: "/home/user"
                        validated_tokenized_path.append(token)
                result['verify']['full_path'] = '/'.join(validated_tokenized_path)

                file = Path(result['verify']['full_path'])
                if not file.is_file():
                    raise RuntimeError(f'Validation file does not exist: "{result["verify"]["full_path"]}"')
                # Update the validated full_path
                result['verify']['full_path'] = result['verify']['full_path']

        # Process PolyBench options
        # First import default options into polybench_options
        if not (args.options is None):
            # Comma separated text -> split
            options = args.options.split(',')
            polybench_opts = result['polybench_options']
            for option in options:
                if option in polybench_opts:  # simple "exists" validation
                    # Only boolean options can pass this validation
                    polybench_opts[option] = True
                else:  # may not exist if the text does not match
                    # Check if it is of the form OPT=val
                    opval = list(option.split('='))
                    if len(opval) == 2:
                        # Ok. Key-value found. Only a few of these exist... check manually
                        # ... for numerical conversions (currently all integers)
                        if opval[1].isnumeric():
                            polybench_opts[opval[0]] = int(opval[1])

        # Custom command line options can override output printing (the verify option). Update polybench_options
        if print_result:
            result['polybench_options'].POLYBENCH_DUMP_ARRAYS = True
            result['polybench_options'].POLYBENCH_DUMP_TARGET = result['output_array']

        # Append the dataset size if required
        if not (args.dataset_size is None):
            # Try to set the enumeration value from user input. On error, an exception is raised.
            if args.dataset_size not in [DataSetSize.MINI.name, DataSetSize.SMALL.name, DataSetSize.MEDIUM.name,
                                         DataSetSize.LARGE.name, DataSetSize.EXTRA_LARGE.name]:
                raise RuntimeError(f'Invalid value for parameter --dataset-size: "{args.dataset_size}"')
            result['polybench_options'].POLYBENCH_DATASET_SIZE = DataSetSize[args.dataset_size]

        # Process the number of iterations.
        if not str(args.iterations).isnumeric() or int(args.iterations) < 1:
            raise RuntimeError(f'Invalid value for parameter --iterations: "{args.iterations}"')
        else:
            result['iterations'] = int(args.iterations)

        # Process array implementation
        if str(args.array_implementation).isnumeric():
            n = int(args.array_implementation)

            if n == 0:
                result['polybench_options'].POLYBENCH_ARRAY_IMPLEMENTATION = ArrayImplementation.LIST
            elif n == 1:
                result['polybench_options'].POLYBENCH_ARRAY_IMPLEMENTATION = ArrayImplementation.LIST_FLATTENED
            elif n == 2:
                result['polybench_options'].POLYBENCH_ARRAY_IMPLEMENTATION = ArrayImplementation.NUMPY
            elif n == 3:
                result['polybench_options'].POLYBENCH_ARRAY_IMPLEMENTATION = ArrayImplementation.LIST_PLUTO
            elif n == 4:
                result['polybench_options'].POLYBENCH_ARRAY_IMPLEMENTATION = ArrayImplementation.LIST_FLATTENED_PLUTO
                # Process PoCC options
                if args.pocc_option in ["pluto","vectorizer","maxfuse"]: result['polybench_options'].POCC = args.pocc_option
                else: raise AssertionError( 'Argument "pocc=%s" not valid.' %(args.pocc_option) )
            else:
                raise AssertionError( 'Argument "array-implementation=%d" not valid.' % (n) )
        else:
            raise AssertionError('Argument "array-implementation" must be a number.')

        # Process save results. Only save when results are available (POLYBENCH_TIME or POLYBENCH_PAPI)
        if args.save_results:
            opts = result['polybench_options']
            if opts.POLYBENCH_TIME or opts.POLYBENCH_PAPI:
                result['save_results'] = True

        # Process load elimination
        opts = result['polybench_options']
        opts.LOAD_ELIMINATION = bool(args.load_elimination)

        return result


    def validate_benchmark_results(options: dict):
        """Compare two files and report if they match or not.

        This validation is intentionally very simple. In case of error, the user must manually compare the results"""
        output_file_name = options['output_array'].name
        verify_file_name = options['verify']['full_path']

        print(f'Verifying if files "{output_file_name}" and "{verify_file_name}" match... ', end='')
        if cmp(output_file_name, verify_file_name):
            print('OK')
        else:
            print('FAIL')
            print('Please, check contents manually.')


    def get_output_file(module_name: str, options: dict):
        output_str = module_name

        # Append interpreter name
        output_str += '_' + python_implementation()

        # Append measurement type information
        if options['polybench_options'].POLYBENCH_TIME:
            if options['polybench_options'].POLYBENCH_CYCLE_ACCURATE_TIMER:
                output_str += '_timer-ca'
            else:
                output_str += '_timer'
        elif options['polybench_options'].POLYBENCH_PAPI:
            output_str += '_papi'

        # Append array type implementation
        if options['polybench_options'].POLYBENCH_ARRAY_IMPLEMENTATION == ArrayImplementation.LIST:
            output_str += '_array=list'
        elif options['polybench_options'].POLYBENCH_ARRAY_IMPLEMENTATION == ArrayImplementation.LIST_FLATTENED:
            output_str += '_array=flattenedlist'
        else:
            output_str += '_array=numpy'
        output_str += '.output'
        return open(output_str, 'w')


    def run(options: dict, spec_file: PolyBenchSpecFile) -> None:
        # Set up parameters which may modify execution behavior
        module_name = options['benchmark']
        # Parameters used in case of verification
        verify_result = options['verify']
        iterations = options['iterations']

        instance = None
        # Search the module within available implementations
        for implementation in benchmark_classes:
            if module_name == 'all' or implementation.__module__ == module_name:
                # Module found!

                # TODO: remove this debug messages
                if True:
                    print(f'Running {implementation.__module__}')
                    from datetime import datetime
                    print(f'  Start time: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
                    print(f'  Interpreter: {python_implementation()}')
                    print(f'  Options: ')
                    print(f'    (iterations, {iterations})')
                    ooo = options['polybench_options']
                    print(f'    {"POLYBENCH_TIME", ooo.POLYBENCH_TIME}')
                    print(f'    {"POLYBENCH_CYCLE_ACCURATE_TIMER", ooo.POLYBENCH_CYCLE_ACCURATE_TIMER}')
                    print(f'    {"POLYBENCH_PAPI", ooo.POLYBENCH_PAPI}')
                    print(f'    {"POLYBENCH_ARRAY_IMPLEMENTATION", ooo.POLYBENCH_ARRAY_IMPLEMENTATION}')

                # Retrieve the appropriate parameters for initializing the current benchmark class
                non_pythonic_benchmark = implementation.__module__.replace('_', '-')
                for spec in spec_file.specs:
                    if spec.Name in non_pythonic_benchmark:
                        bench_specs = spec
                        break

                if options['save_results']:
                    output_f = get_output_file(implementation.__module__.replace('.', '/'), options)

                first_run = True  # For printing available columns on PAPI result

                # Run the benchmark N times. N will be either 1 or a greater number passed by argument.
                for i in range(iterations):
                    # Instantiate a new class with it
                    instance = implementation(options['polybench_options'], bench_specs)  # creates a new instance
                    # Run the benchmark. The returned value is a dictionary.
                    polybench_result = instance.run()

                    if options['save_results']:
                        # Perform operations against the output data when the appropriate option is enabled.
                        if options['polybench_options'].POLYBENCH_TIME:
                            output_f.write(f'{polybench_result["POLYBENCH_TIME"]}\n')
                            output_f.flush()

                        if options['polybench_options'].POLYBENCH_PAPI:
                            if first_run:
                                # Print headers
                                for counter in polybench_result["POLYBENCH_PAPI"]:
                                    output_f.write(f'{counter}\t')
                                output_f.write('\n')
                            for counter in polybench_result["POLYBENCH_PAPI"]:
                                output_f.write(f'{polybench_result["POLYBENCH_PAPI"][counter]}\t')
                            output_f.write('\n')
                            output_f.flush()

                    first_run = False
                if options['save_results']:
                    output_f.close()

                # Verify benchmark's results against other implementation's results
                if verify_result['enabled']:
                    validate_benchmark_results(options)

                # Terminate the loop for single-benchmark run
                if module_name != 'all':
                    break

        # Check if the module was not found and report an error accordingly
        if instance is None:
            module = module_name.replace(".", "/") + '.py'
            raise NotImplementedError(f'Module {module} not implemented.')


    # Parse the command line arguments first. We need at least one mandatory parameter.
    opts = parse_command_line()
    # Parse the spec file for obtaining all of the benchmark's parameters
    spec_f = PolyBenchSpecFile('polybench.spec')
    # Run the benchmark (and other user options)
    run(opts, spec_f)
