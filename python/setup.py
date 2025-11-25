# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 23:48:22 2015

@author: thorsten
"""

from setuptools import setup
from setuptools import Extension
from setuptools import __version__ as setuptools_version

import sys
import math
import os
import platform
import subprocess
import pathlib
import glob

# If users want to import local modules (for setup.py itself), they are
# recommended to explicitly add the current directory to sys.path
# https://github.com/pypa/setuptools/issues/3939
sys.path.append(str(pathlib.Path(__file__).parent))
from bootstrap.find_package import add_csxcad, add_h5py


LICENSE = "GPL-3.0-or-later"


def get_license():
    # mutually exclusive "license" and "license_expression"
    if int(setuptools_version.split(".")[0]) < 77:
        return {"license": LICENSE}
    else:
        return {"license_expression": LICENSE}


def get_fallback_version(pyproject_toml, fallback_file):
    # Always generate a fallback_version in case importlib.metadata
    # introspection is unsupported or fails.
    with open(pyproject_toml, 'r') as toml:
        matching_list = [
            line for line in toml if line.startswith("fallback_version")
        ]
        fallback_version_quoted = matching_list[0].split("=")[1]
        fallback_version = fallback_version_quoted.replace('"', "").strip()

    with open(fallback_file, "w+") as fallback_file:
       fallback_file.write("__fallback_version__ = '%s'" % fallback_version)

    try:
        # pyproject.toml is respected by new pip, which calls setuptools_scm
        # to inject a version automatically.
        import setuptools_scm
        from importlib.metadata import version
        if int(version("setuptools_scm").split(".")[0]) >= 8:
            return None
        else:
            raise ValueError("setuptools_scm version too low.")
    except (ImportError, ValueError):
        return fallback_version


def normalize_path_subdir(path_str, subdir):
    # We want to find the correct $CSXCAD_INSTALL_PATH even if the
    # user input is slightly "incorrect".
    #
    # /home/opt/physics/venv -> /home/opt/physics (subdir == "venv")
    # /home/opt/physics/bin -> /home/opt/physics (subdir == "bin")
    path = pathlib.Path(path_str)
    if path.stem == subdir and (path / "../include/openEMS").exists():
        return (path / "../").resolve()
    else:
        return path


def determine_build_options():
    cpu_bits = str(int(math.log2(sys.maxsize) + 1))

    local_prefix_list = []
    use_prefix = True
    user_entered_prefix = False

    if "OPENEMS_INSTALL_PATH_IGNORE" in os.environ:
        use_prefix = False

    if use_prefix and 'OPENEMS_INSTALL_PATH' in os.environ:
        # Because CSXCAD and openEMS are usually installed to a non-standard
        # location such as the user's home directory rather than a standard
        # directory, their paths must be told to the compilers explicitly.
        # Since a general solution is impossible (the actual installation path
        # is unpredictable), users should set OPENEMS_INSTALL_PATH to ensure a
        # successful installation.
        local_prefix_list.append(
            normalize_path_subdir(os.environ["OPENEMS_INSTALL_PATH"], "bin")
        )
        user_entered_prefix = True
    if 'VIRTUAL_ENV' in os.environ:
        # if Python venv used, our documentation recommends using the same
        # path as OPENEMS_INSTALL_PATH, so that the custom C++ and Python prefix
        # overlaps.
        local_prefix_list.append(
            normalize_path_subdir(os.environ["VIRTUAL_ENV"], "venv")
        )
        user_entered_prefix = True

    if use_prefix and (not user_entered_prefix):
        raise RuntimeError(
            "No environment variable OPENEMS_INSTALL_PATH or VIRTUAL_ENV found, "
            "installation may fail due to missing headers and libraries! "
            "Please set the environment variable OPENEMS_INSTALL_PATH to the path "
            "of CSXCAD/openEMS installation, check documentation for details. "
            "If you know what you're doing, set OPENEMS_INSTALL_PATH_IGNORE=1 to "
            "suppress this error."
        )

    if platform.system() == "Darwin":
        # In additional to libraries in $OPENEMS_INSTALL_PATH and $VIRTUAL_ENV, we
        # also need to list custom headers and libraries installed to the local
        # system but are not used by compilers by default (such as a custom Boost).
        # On macOS, the required prefix are -L $(brew --prefix)/include and
        # -R $(brew --prefix)/lib respectively. Hardcode it as a special treatment.
        path = None
        try:
            status = subprocess.run(["brew", "--prefix"], capture_output=True)
            path = status.stdout.decode("UTF-8").replace("\n", "")
            local_prefix_list.append(
                pathlib.Path(path)
            )
        except FileNotFoundError:
            pass

        # If pip build isolation is disbled, users must install
        # dependencies manually. Unfortunately Homebrew's Cython
        # is a keg-only internal package outside the search path.
        # Add Cython to Python path manually.
        if path and not any(["pip-build-env" in i for i in sys.path]):
            homebrew_prefix = pathlib.Path(path) / "Cellar" / "cython"
            python_ver = (sys.version_info.major, sys.version_info.minor)
            python_dir = "python%d.%d/site-packages" % python_ver

            find_candidate_cython_list = [
                i for i in homebrew_prefix.rglob(python_dir) if i.is_dir()
            ]
            if find_candidate_cython_list:
                sys.path.append(str(find_candidate_cython_list[0].resolve()))

    if os.name == "posix":
        # The path /usr/local is also too common on Unix systems, so we hardcode it
        # as a special treatment on Unix-like systems. For example, on CentOS, the
        # paths -L /usr/local/include and -R /usr/local/lib must be listed if a
        # custom Boost is installed here.
        local_prefix_list.append(pathlib.Path("/usr/local"))

    build_options = {
        "extra_compile_args": [],
        "include_dirs": [],
        "library_dirs": [],
        "runtime_library_dirs": [],
    }

    for prefix in local_prefix_list:
        prefix_path = pathlib.Path(prefix)
        if os.name == 'nt':
            build_options["library_dirs"] += [str(prefix_path)]

        build_options["include_dirs"] += [
            str(prefix_path / "include")
        ]
        build_options["library_dirs"] += [
            str(prefix_path / "lib"),
            str(prefix_path / "lib" / cpu_bits)
        ]
        build_options["runtime_library_dirs"] += [
            str(prefix_path / "lib"),
            str(prefix_path / "lib" / cpu_bits)
        ]

    # Strictly speaking we should detect compiler, not platform,
    # unfortunately there's no easy way to do so without implementing
    # full compiler detection logic. This is good enough for 90% of
    # use cases.
    if os.name == "posix":
        build_options["extra_compile_args"].append("-std=c++11")

    # Setting this will cause an exception during build on Windows platforms.
    if os.name != "posix":
        del build_options["runtime_library_dirs"]

    return build_options


def get_modules_list(module_prefix, path_glob_pattern, build_options):
    output_list = []

    if int(setuptools_version.split(".")[0]) < 18:
        from Cython.Build import cythonize
        output_list = cythonize(
            Extension(name="*", sources=[path_glob_pattern], **build_options)
        )
    else:
        # Above setuptools 18, setuptools contains a special case for `.pyx`
        # files if setup_requires=["cython"] is set, pyx sources are auto-
        # Cythonized. This is recommended since setup.py / pip still runs
        # even if Cython is not installed. But wildcards are not supported,
        # we collect modules by hand.

        # e.g. [openEMS/_nf2ff.pyx, openEMS/openEMS.pyx, ...]
        filepath_list = glob.glob(path_glob_pattern)
        for filepath in filepath_list:
            # e.g. openEMS/openEMS.pyx
            filename = os.path.basename(filepath)
            # e.g. openEMS.openEMS
            module_name = module_prefix + filename.replace(".pyx", "")

            output_list.append(
                Extension(name=module_name, sources=[filepath], **build_options)
            )
    return output_list


build_opt = determine_build_options()
build_opt["language"] = "c++"
build_opt["libraries"] = ["CSXCAD", "openEMS", "nf2ff"]

extensions = get_modules_list(
    module_prefix="openEMS.",
    path_glob_pattern="openEMS/*.pyx",
    build_options=build_opt
)

# DO add new run-time dependencies here.
install_requires = [
    # BSD 3-Clause (https://github.com/numpy/numpy/blob/main/LICENSE.txt)
    'numpy',
]

# BSD 3-Clause (https://github.com/h5py/h5py/blob/master/LICENSE)
install_requires += add_h5py()

try:
    # CSXCAD already installed, use a plain package name as dependency,
    # pip won't rebuild it. This ensures manual package management still
    # works, and pip won't redownload anything from git.
    import CSXCAD
    install_requires += ["CSXCAD"]
except ImportError:
    # CSXCAD not installed, use CSXCAD with file:// or git:// as dependency.
    # As a side-effect, pip always rebuilds CSXCAD from scratch regardless
    # of whether it's installed. We don't want this side-effect.
    install_requires += add_csxcad()

setup(
  name="openEMS",
  version=get_fallback_version(
    "pyproject.toml", "openEMS/__fallback_version__.py"
  ),
  packages=["openEMS", ],
  package_data={'openEMS': ['*.pxd']},
  # DO NOT add any new build-time dependency in setup_requires.
  # We should use pyproject.toml exclusively. The only item
  # "cython" is meant to activate auto-Cython feature in
  # setuptools 18.
  setup_requires=[
    'cython'
  ],
  install_requires=install_requires,
  ext_modules=extensions,
  **get_license()
)
