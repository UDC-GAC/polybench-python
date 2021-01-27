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
# Author:
#    Miguel Ángel Abella González <miguel.abella@udc.es>
#
# Contact:
#    Gabriel Rodríguez <gabriel.rodriguez@udc.es>


"""This utility simplifies the task of creating a new benchmark.

Manually generating a benchmark implies creating a directory structure and
the required initialization files for the Python's package system to work.

This tool does all required steps for creating a package and a template file
for the new benchmark. Once created, this benchmark will be available when
running the launcher script (run-benchmark.py)"""

# Using argparse for parsing commandline options. See: https://docs.python.org/3.7/library/argparse.html
import argparse
# Using Path for creating folders representing categories
from pathlib import Path
import os.path
import shutil

from string import ascii_lowercase, digits


if __name__ == '__main__':

    # Define some "constant" references
    benchmark_path = 'benchmarks'
    template_benchmark_file = 'util/template-benchmark.py'
    template_init_file = 'util/template-package-init.py'

    def parse_command_line() -> {
        'benchmark_name': str,
        'module_name': str,
        'category': str
    }:
        """
        Parse command line arguments and normalize the results.

        :return: A dictionary with the decoded parameters.
        :rtype: dict[str, str]

        The possible return values for the dictionary are listed below:
            - **benchmark_name** [param: *--name*]: the name for the benchmark.
            - **category** [param: *--category*]: the name for the category. This value will be normalized to use
            slashes as separators.
        """
        parser = argparse.ArgumentParser(description='Generate a template for implementing a new benchmark.',
                                         epilog='All benchmarks are created under the "benchmarks" folder.')
        parser.add_argument('--name', '-N', type=str, required=True,
                            metavar='<BenchmarkName>',
                            help='The name of the benchmark using PascalCase.')
        parser.add_argument('--category', '-C', type=str,
                            metavar='<category[.subcategory] | category[/subcategory]>',
                            help='The category name. This will place the benchmark template into a package, generating '
                                 'all required intermediate files. Subcategories can be separated either with a dot or '
                                 'a slash. Example: datamining/correlation and datamining.correlation are valid. '
                                 'Category names are expected to be in lower case and will be internally converted.')
        # Parse the commandline arguments. This process will fail on error
        args = parser.parse_args()

        result = {}

        # Process the benchmark name
        bench_name = str(args.name).replace('-', '_')

        # Process the category name. Since the category is optional, we need to check if it exists first
        almost_category = None
        if not (args.category is None):
            # Convert "." into "/" if necessary and remove tailing "/" when present
            almost_category = str(args.category).replace('.', '/').lower().rstrip('/')

        # Process both benchmark_name and category for making them respect Python module and package naming conventions
        # Modules should be all lowercase and underscores.
        # We can relax this convention by also including digits.
        valid_mod_chars = ascii_lowercase + '_' + digits
        valid_pkg_chars = valid_mod_chars

        def fix_naming(input_str: str) -> str:
            # Initialize fn_result
            if len(input_str) > 0 and input_str[0].isdigit():
                fn_result = '_'  # Prepend an underscore
            else:
                fn_result = ''

            # Iterate the string and filter out invalid elements
            for char in input_str.lower():
                if char in valid_mod_chars:
                    fn_result += char
                elif char == '-':
                    fn_result += '_'
            return fn_result

        # Class, module and package names cannot start with a number, Must prepend an underscore.
        if bench_name[0].isdigit():
            bench_name = '_' + bench_name

        mod_name = fix_naming(bench_name)

        cat_name = ''
        if not (almost_category is None):
            # Category processing rules are similar to module requirements.
            # Split the path into smaller chunks so we can reuse code.
            tokens = almost_category.split('/')

            cat_tokens = []
            for token in tokens:
                cat_tokens.append(fix_naming(token))

            # Join the results in a new and valid path
            cat_name = '/'.join(cat_tokens)

        result['benchmark_name'] = bench_name
        result['module_name'] = mod_name
        if almost_category is None:
            result['category'] = None
        else:
            result['category'] = cat_name

        return result

    def create_package_structure(category: str) -> str:
        """
        Create a valid python package structure (of files and folders) given a category name.

        :param str category: The category name for which to create the package structure. If this parameter is set to
            "None" then no actions will be performed by this function.
        :return: The full path for the created package relative to the running script's location. If "None" was passed
            as argument, the returned path will be the directory specified by "benchmark_path".
        :rtype str:

        Notes
        -----
        The package structure for Polybench/Python is always relative to the path specified by "benchmark_path". This
        implies that package structures will always be inside that path.
        """
        target_path = benchmark_path
        if not (category is None):
            target_path += '/' + category
            # Create, if necessary, all of the directories in the package path
            Path(target_path).mkdir(parents=True, exist_ok=True)

            # Check the existence of "__init__.py" files and, if necessary, create them from the template
            init_file_locations = target_path.split('/')[1:-1]  # exclude "benchmarks" and the final dir from this list
            init_file_path = benchmark_path + '/'
            for check in init_file_locations:
                init_file_path += check + '/'
                init_file = init_file_path + '__init__.py'
                if not os.path.isfile(init_file):
                    shutil.copy2(template_init_file, init_file)

        return target_path

    def create_benchmark(directory: str, benchmark: str) -> None:
        """
        Create a new benchmark template and its associated files.

        While creating the benchmark, Python naming conventions will be respected. Class names must be capitalized and
        package names should be lower case.

        :param str directory: The directory where the benchmark will be created.
        :param str benchmark: The name of the benchmark. This name will be modified for fitting Python's naming
            conventions.
        """
        # Blindly capitalize the benchmark name...
        benchmark_capitalized = benchmark[0].upper() + benchmark[1:]    # ... for using as class name
        benchmark_lowered = benchmark.lower()   # ... for using in package and file names

        # Check benchmark location for the existence of an "__init__.py" file
        init_file = directory + '/__init__.py'
        if not os.path.isfile(init_file):
            # The init file does not exist. Read the template's contents into memory, update them for this benchmark and
            # save the results into the new init file.
            with open(template_init_file) as f:
                contents = f.read()
                # Append: from path.to.package.benchmark import Benchmark
                contents += f'\nfrom {directory.replace("/", ".")}.{benchmark_lowered} import {benchmark_capitalized}\n'
                # Save the new init file
                with open(init_file, 'x') as g:
                    g.write(contents)

        target_file = directory + '/' + benchmark_lowered + '.py'
        # Check benchmark location for the existence of the benchmark. If a previous benchmark exists it will be kept,
        # otherwise a new benchmark template will be created with the correct benchmark's name.
        if not os.path.isfile(target_file):
            with open(template_benchmark_file) as f:
                contents = f.read()
                with open(target_file, 'x') as g:
                    contents = contents.replace('TemplateClass', benchmark_capitalized)
                    g.write(contents)


    # Parse the commandline and obtain the required parameters
    parameters = parse_command_line()
    benchmark_name = parameters['benchmark_name']
    category_name = parameters['category']

    # Build the directory structure represented by the category
    package_path = create_package_structure(category_name)
    create_benchmark(package_path, benchmark_name)
