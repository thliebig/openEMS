# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 23:48:22 2015

@author: thorsten
"""

from setuptools import setup
from setuptools import Extension
from Cython.Build import cythonize

import os, sys
ROOT_DIR = os.path.dirname(__file__)

sys.path.append(os.path.join(ROOT_DIR,'..','..','CSXCAD','python'))

extensions = [
    Extension("*", [os.path.join(os.path.dirname(__file__), "openEMS","*.pyx")],
        language="c++",             # generate C++ code
        libraries    = ['CSXCAD','openEMS', 'nf2ff']),
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
    'cython>=3.0.6',  # Apache License 2.0 (https://github.com/cython/cython/blob/master/LICENSE.txt)
    'h5py>=3.10.0',   # BSD 3-Clause (https://github.com/h5py/h5py/blob/master/LICENSE)
    'numpy>=1.26.2',  # BSD 3-Clause (https://github.com/numpy/numpy/blob/main/LICENSE.txt)
    'typing_extensions', # Used for `typing_extensions.deprecated` warning, will not be needed as of Python 3.13 because it will become part of the standard `warnings` library.
  ],
  ext_modules = cythonize(extensions, language_level = 3)
 )
