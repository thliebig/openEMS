from setuptools import Extension, setup
from Cython.Build import cythonize
from pathlib import Path
import os
import sys

# TODO remove hardcoded path:
PATH_TO_CSXCAD_PYTHON_PACKAGE = Path(__file__).parent.parent.parent/'CSXCAD/python'

# Define some functions to check that the packages are installed where we expect them to be installed, and raise meaningful errors otherwise:
def looks_like_openEMS_is_installed_there(path_to_some_folder:Path):
    path_to_some_folder = Path(path_to_some_folder)
    FILES_THAT_SHOULD_EXIST = [
        path_to_some_folder/'bin/openEMS',
        path_to_some_folder/'include/openEMS/openems.h',
    ]
    if any([not p.is_file() for p in FILES_THAT_SHOULD_EXIST]):
       return False
    return True

def looks_like_CSXCAD_python_package(path_to_some_folder:Path):
    path_to_some_folder = Path(path_to_some_folder)
    FILES_THAT_SHOULD_EXIST = [
        path_to_some_folder/'pyproject.toml',
        path_to_some_folder/'CSXCAD/__init__.py',
        path_to_some_folder/'CSXCAD/CSXCAD.pyx',
    ]
    if any([not p.is_file() for p in FILES_THAT_SHOULD_EXIST]):
        return False
    return True

if 'OPENEMS_INSTALL_PATH' not in os.environ:
    raise SystemExit('Please set the environment variable OPENEMS_INSTALL_PATH to point to where openEMS was installed. ')

path_to_openEMS_installation = Path(os.environ['OPENEMS_INSTALL_PATH']) # Environment variables are always strings, so this should never raise any error.

if not looks_like_openEMS_is_installed_there(path_to_openEMS_installation):
    raise SystemExit(f'I was expecting to find openEMS installed in {path_to_openEMS_installation}, but it does not look like it is installed there. ')

if not looks_like_CSXCAD_python_package(PATH_TO_CSXCAD_PYTHON_PACKAGE):
    raise SystemExit(f'Cannot find python module of CSXCAD in {PATH_TO_CSXCAD_PYTHON_PACKAGE}. Please modify the `PATH_TO_CSXCAD_PYTHON_PACKAGE` variable in this setup.py file, which is located in {Path(__file__).resolve()}, to point to wherever you have the CSXCAD python package. ')

sys.path.append(str(PATH_TO_CSXCAD_PYTHON_PACKAGE))

# Strictly speaking we should detect compiler, not platform,
# unfortunately there's no easy way to do so without implementing
# full compiler detection logic. This is good enough for 90% of
# use cases.
cxxflags = []
if os.name == "posix":
    cxxflags.append("-std=c++11")

# Strictly speaking we should detect compiler, not platform,
# unfortunately there's no easy way to do so without implementing
# full compiler detection logic. This is good enough for 90% of
# use cases.
cxxflags = []
if os.name == "posix":
    cxxflags.append("-std=c++11")

extensions = [
<<<<<<< ours
    Extension(
        '*',
        ['openEMS/*.pyx',],
        include_dirs = [str(path_to_openEMS_installation/'include')],
        library_dirs = [str(path_to_openEMS_installation/'lib')],
        runtime_library_dirs = [str(path_to_openEMS_installation/'lib')],
        language = 'c++',
        libraries = ['openEMS','nf2ff'],
        extra_compile_args = cxxflags,
    )
]

setup(
    packages = ['openEMS'],
    ext_modules = cythonize(extensions),
)
=======
    Extension("*", [os.path.join(os.path.dirname(__file__), "openEMS","*.pyx")],
        language="c++",             # generate C++ code
        libraries    = ['CSXCAD','openEMS', 'nf2ff'],
        extra_compile_args=cxxflags),
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
>>>>>>> theirs
