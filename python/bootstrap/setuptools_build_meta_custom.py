from find_package import add_csxcad
from setuptools import build_meta as _orig

try:
    from setuptools import build_meta as _orig
except:
    raise RuntimeError(
        "setuptools is not installed. If pip --no-build-isolation "
        "is used, install setuptools and cython manually before using "
        "pip!"
    )

# This "import all" is intentional, since we're selectively
# overriding setuptools while keeping most things unchanged,
# don't edit this.
from setuptools.build_meta import *


def get_requires_for_build_wheel(self, config_settings=None):
    orig_list = _orig.get_requires_for_build_wheel(config_settings)
    orig_list += add_csxcad()
    return orig_list


def get_requires_for_build_sdist(self, config_settings=None):
    orig_list = _orig.get_requires_for_build_sdist(config_settings)
    orig_list += add_csxcad()
    return orig_list
