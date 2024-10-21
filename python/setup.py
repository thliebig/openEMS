from setuptools import Extension, setup
from Cython.Build import cythonize
from pathlib import Path
import sys

# Hardcoded paths that would be better to remove but I don't know how to:
PATH_TO_OPEN_EMS_INSTALLATION = Path(__file__).parent.parent.parent.parent
PATH_TO_CSXCAD_PYTHON_PACKAGE = Path(__file__).parent.parent.parent/'CSXCAD/python'

# Define some functions to check that the packages are installed where we expect them to be installed, and raise meaningful errors otherwise:
def check_that_folder_looks_like_openEMS_installation(path_to_some_folder:Path):
    path_to_some_folder = Path(path_to_some_folder)
    FILES_THAT_SHOULD_EXIST = [
        path_to_some_folder/'bin/openEMS',
        path_to_some_folder/'include/openEMS/openems.h',
    ]
    if any([not p.is_file() for p in FILES_THAT_SHOULD_EXIST]):
        raise RuntimeError(f'I was expecting to find the installation of openEMS in {path_to_some_folder}, but it does not look like it is installed there. Please provide the path to where you have it installed and write it in `PATH_TO_OPEN_EMS_INSTALLATION` in the file {Path(__file__).resolve()}. Should you have it installed there, please create an issue in https://github.com/thliebig/openEMS reporting this. ')

def check_that_folder_looks_like_CSXCAD_python_package(path_to_some_folder:Path):
    path_to_some_folder = Path(path_to_some_folder)
    FILES_THAT_SHOULD_EXIST = [
        path_to_some_folder/'pyproject.toml',
        path_to_some_folder/'CSXCAD/__init__.py',
        path_to_some_folder/'CSXCAD/CSXCAD.pyx',
    ]
    if any([not p.is_file() for p in FILES_THAT_SHOULD_EXIST]):
        raise RuntimeError(f'I was expecting to find the python package CSXCAD in {path_to_some_folder}, but it does not look like it is there. Please provide the path to where you have it installed and write it in `PATH_TO_CSXCAD_PYTHON_PACKAGE` in the file {Path(__file__).resolve()}. Should you have it installed there, please create an issue in https://github.com/thliebig/CSXCAD reporting this. ')

# Check that the packages are installed there:
check_that_folder_looks_like_openEMS_installation(PATH_TO_OPEN_EMS_INSTALLATION)
check_that_folder_looks_like_CSXCAD_python_package(PATH_TO_CSXCAD_PYTHON_PACKAGE)

# Perform the actual installation:

sys.path.append(str(PATH_TO_CSXCAD_PYTHON_PACKAGE))

extensions = [
    Extension(
        '*', 
        ['openEMS/*.pyx',],
        include_dirs = [str(Path(PATH_TO_OPEN_EMS_INSTALLATION)/'include')],
        library_dirs = [str(Path(PATH_TO_OPEN_EMS_INSTALLATION)/'lib')],
        runtime_library_dirs = [str(Path(PATH_TO_OPEN_EMS_INSTALLATION)/'lib')],
        language = 'c++',
        libraries = ['openEMS','nf2ff'],
    )
]

setup(
    packages = ['openEMS'],
    ext_modules = cythonize(extensions),
)

