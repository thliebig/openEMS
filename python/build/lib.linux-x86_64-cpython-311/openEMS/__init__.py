# -*- coding: utf-8 -*-
#
# Shortcut openEMS import
from __future__ import absolute_import

# try to import CSXCAD, e.g. make sure the (windows) dll path are set
import CSXCAD

try:
    from importlib.metadata import version
    __version__ = version("openEMS")
except ImportError:
    try:
        from openEMS.__fallback_version__ import __fallback_version__
        __version__ = __fallback_version__
    except ImportError:
        __version__ = "0.0.0"

from openEMS.openEMS import openEMS
