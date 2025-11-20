# This file is used to find and insert dependencies dynamically.
# It's used by both "setuptools_build_meta_custom.py" to find
# build-time dependencies, and also by "setup.py" to find run-
# time dependencies.
import os
from pathlib import Path


def is_project_root(path_prefix):
    # Can we find the following files?
    file_list = ["update_openEMS.sh", "CSXCAD/python"]

    for file in file_list:
        if not (path_prefix / Path(file)).exists():
            return False

    # If so, it means we're build openEMS/python from the openEMS-Project
    # directory, and so should point CSXCAD's path to:
    # openEMS-Project/CSXCAD/python
    return True


def pip_in_tree_build():
    # pip has two modes: in-tree build and out-of-tree build,
    # depending on its version.
    #
    # Out-of-tree was the default for pip < 21.3, the Python
    # source directory is copied to /tmp for a clean build, but
    # this caused problems if the project uses anything from the
    # parent directory (what we do here). For pip >= 21.3,
    # out-of-tree builds are used instead.
    #
    # A way to distinguish them is to see whether an expected file
    # from the parent directory exists.
    if Path("../openems.h").exists():
        return True
    else:
        return False


def add_csxcad(version_str=""):
    if "CSXCAD_PYSRC_PATH" in os.environ:
        # User provides an explicit source path, just use it without
        # trying anything else.
        if os.environ["CSXCAD_PYSRC_PATH"].startswith("git+"):
            return ["CSXCAD %s @ %s" % (version_str, os.environ["CSXCAD_PYSRC_PATH"])]
        else:
            path = Path(os.environ["CSXCAD_PYSRC_PATH"]).resolve()
            # "file://localhost/" is usually written out as just
            # "file:///", they're equivalent per RFC 8089. We use the
            # full form here because pypa/packaging doesn't correctly
            # handle "file:///" until v19.0.
            #
            # See:
            # https://github.com/pypa/packaging/commit/eb0243865079c7c2179730da93130447556ea020
            # https://github.com/pypa/pip/blob/7e49dca9277bf4e325b85cfb9ebe70401f194fb6/src/pip/_internal/utils/urls.py#L30
            return ["CSXCAD %s @ file://localhost/%s" % (version_str, str(path))]

    # User provides no path, let's try some heuristics now.
    path = None

    if is_project_root(Path("../../")):
        # Find path relative to setup.py. Do we see a project root?
        path = Path("../../CSXCAD/python/")
    elif (
        "PWD" in os.environ and
        is_project_root(Path(os.environ["PWD"]) / "../../")
    ):
        # Find path relative to user shell. It's not reliable (user's
        # shell must be in "openEMS/python" but pip can be invoked anywhere),
        # but it's our last chance (e.g. on pip < 21.3 with out-of-tree
        # builds).
        path = Path(os.environ["PWD"]) / "../../CSXCAD/python/"

    # found path
    if path:
        return ["CSXCAD %s @ file://localhost/%s" % (version_str, path.resolve())]

    # path not found, try fallbacks
    if pip_in_tree_build():
        # We're confident that the user has nothing but an isolated openEMS.git
        # repo, we instruct pip to download CSXCAD from scratch.
        path = "git+https://github.com/thliebig/CSXCAD.git#subdirectory=python"
        return ["CSXCAD %s @ %s" % (version_str, path)]
    else:
        # The pip version is new enough to support pyproject.toml and build
        # backends, but old enough so it uses out-of-tree builds, preventing
        # us from locating CSXCAD's source code. It's too ambiguous to make
        # any decision: Are we trying to build openEMS itself, or build openEMS
        # submodule within the openEMS-Project directory? Instead of returning
        # the latest git version of CSXCAD as dependency, which is potentially
        # incompatible, we raise an error.
        raise RuntimeError(
            "Unable to detect CSXCAD's Python source code path. "
            "You're likely using an old pip without 'in-tree build' "
            "support. You can pick one solution below: "
            "(1) Rerun pip with 'pip install . --no-build-isolation' "
            "if CSXCAD Python extension is already installed (recommended). "
            "(2) Provide the path via CSXCAD_PYSRC_PATH and rerun pip "
            "(e.g. 'export "
            "CSXCAD_PYSRC_PATH=/home/user/openEMS-Project/CSXCAD/python/ && "
            "pip install . '). "
            "(3) Upgrade to pip 21.3 or newer."
        )


def add_setuptool_scm():
    if "OPENEMS_NOSCM" in os.environ:
        return []
    else:
        return ['setuptools_scm >= 8; python_version >= "3.9"']
