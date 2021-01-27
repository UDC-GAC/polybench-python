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

"""Perform package discovery and module importing for retrieving a list of all available benchmarks."""


import pkgutil
from benchmarks.polybench import PolyBench


def __build_module_list__() -> (set, set):
    """Builds the list of available modules and PolyBench subclasses."""
    from pathlib import Path
    # Prepare the function result.
    candidates = set()

    # For some reason pkgutil.walk_packages() does not recurse from the "benchmarks" package. Use the parent directory
    # to make it work as expected and filter out the results to include only the packages from the "benchmarks" package.
    from_path = [Path(__path__[0]).parent]

    # See: https://docs.python.org/3/library/pkgutil.html#pkgutil.ModuleInfo
    # See: https://docs.python.org/3/library/pkgutil.html#pkgutil.walk_packages
    for moduleInfo_finder, moduleInfo_name, moduleInfo_ispackage in pkgutil.walk_packages(from_path):
        # Exclude all names not starting with "benchmarks."
        if str(moduleInfo_name).startswith("benchmarks."):
            # Include only modules, not packages.
            if not moduleInfo_ispackage:
                candidates.add(moduleInfo_name)

    # The algorithm added the module "polybench" which defines the abstract class Polybench. Remove it from the modules
    # set as it is not interesting to have here.
    candidates.remove('benchmarks.polybench')
    # Exploit the side effect of "walk_packages()" (package and module loading) and invoke "__subclasses__()"
    # See: https://docs.python.org/3/library/stdtypes.html#class.__subclasses__
    return sorted(candidates), PolyBench.__subclasses__()


# WARNING: unused "benchmark_modules"
benchmark_modules, benchmark_classes = __build_module_list__()
