# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 23:48:22 2015

@author: thorsten
"""

from setuptools import setup
from setuptools import Extension
from Cython.Build import cythonize

import sys
import math
import os
import platform
import subprocess
import pathlib


ROOT_DIR = os.path.dirname(__file__)

sys.path.append(os.path.join(ROOT_DIR,'..','..','CSXCAD','python'))


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
        try:
            status = subprocess.run(["brew", "--prefix"], capture_output=True)
            path = status.stdout.decode("UTF-8").replace("\n", "")
            local_prefix_list.append(
                pathlib.Path(path)
            )
        except FileNotFoundError:
            pass
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


build_options = determine_build_options()
extensions = [
    Extension(
        name="*",
        sources=[os.path.join("openEMS", "*.pyx")],
        language="c++",
        libraries=['CSXCAD', 'openEMS', 'nf2ff'],
        **build_options
    ),
]

setup(
  name="openEMS",
  version = '0.0.36',
  description = "Python interface for the openEMS FDTD library",
  classifiers = [
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Developers',
      'Intended Audience :: Information Technology',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
      'Programming Language :: Python',
      'Programming Language :: Python :: 3',
      'Programming Language :: Python :: 3.9',
      'Programming Language :: Python :: 3.10',
      'Programming Language :: Python :: 3.11',
      'Programming Language :: Python :: 3.12',
      'Programming Language :: Python :: Implementation :: CPython',
      'Topic :: Scientific/Engineering',
      'Topic :: Software Development :: Libraries :: Python Modules',
      'Operating System :: POSIX :: Linux',
      'Operating System :: Microsoft :: Windows',
  ],
  author = 'Thorsten Liebig',
  author_email = 'Thorsten.Liebig@gmx.de',
  maintainer = 'Thorsten Liebig',
  maintainer_email = 'Thorsten.Liebig@gmx.de',
  url = 'https://openEMS.de',
  packages=["openEMS", ],
  package_data={'openEMS': ['*.pxd']},
  python_requires='>=3.9',
  install_requires=[
    'cython', # Apache License 2.0 (https://github.com/cython/cython/blob/master/LICENSE.txt)
    'h5py',   # BSD 3-Clause (https://github.com/h5py/h5py/blob/master/LICENSE)
    'numpy',  # BSD 3-Clause (https://github.com/numpy/numpy/blob/main/LICENSE.txt)
  ],
  ext_modules = cythonize(extensions, language_level = 3)
 )
