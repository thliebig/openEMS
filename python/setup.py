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

extensions = [
    Extension(
        '*', 
        ['openEMS/*.pyx',],
        include_dirs = [str(path_to_openEMS_installation/'include')],
        library_dirs = [str(path_to_openEMS_installation/'lib')],
        runtime_library_dirs = [str(path_to_openEMS_installation/'lib')],
        language = 'c++',
        libraries = ['openEMS','nf2ff'],
    )
]

setup(
    packages = ['openEMS'],
    ext_modules = cythonize(extensions),
)

