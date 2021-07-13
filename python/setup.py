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
  version = '0.0.33',
  description = "Python interface for the openEMS FDTD library",
  classifiers = [
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Developers',
      'Intended Audience :: Information Technology',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering',
      'Topic :: Software Development :: Libraries :: Python Modules',
      'Operating System :: POSIX :: Linux',
      'Operating System :: Microsoft :: Windows',
  ],
  author = 'Thorsten Liebig',
  author_email = 'Thorsten.Liebig@gmx.de',
  maintainer = 'Thorsten Liebig',
  maintainer_email = 'Thorsten.Liebig@gmx.de',
  url = 'http://openEMS.de',
  packages=["openEMS", ],
  package_data={'openEMS': ['*.pxd']},
  ext_modules = cythonize(extensions)
 )
