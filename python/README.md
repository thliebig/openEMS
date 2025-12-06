# openEMS Python Interface

openEMS Python extension, providing a binding to the C++ openEMS
electromagnetic field solver.

## Latest documentation

This file may be outdated. For the latest documentation, check
the online web page.

* [CSXCAD/openEMS Python Interface Installation](https://openems.readthedocs.io/en/latest/python/install.html)

### Dependencies

#### Python version

Ensure the system's Python interpreter is officially supported by
Python, such as Python 3.9. Use unsupported Python versions at your
own risk.

Installation is allowed using practically all Python versions (Python
3.4+), but for testing purposes only. Use at your own risk. As time of
writing, Python 3.5 is known to partially work (i.e., run a trivial script),
while Python 3.6 is likely the lowest fully-functional version.
In Python 3.5 and lower versions, `SyntaxError` may be encountered.

#### Other Dependencies

This file assumes readers have already installed required C++
and Python dependencies (including compiling and installing
the C++ CSXCAD library and openEMS field solver into the system,
and installing the CSXCAD Python extension). If not, follow the
online manual.

* [Requirements of Building CSXCAD/openEMS](https://openems.readthedocs.io/en/latest/install/index.html)

The manual includes instructions for multiple systems (Alpine, AlmaLinux,
CentOS, Debian, Fedora, Ubuntu, FreeBSD, macOS), they're not repeated here
for brevity.

## Quick Start

If the CSXCAD library and openEMS field solver were both installed into
`~/opt/physics`, install this package with:

```bash
# create an isolated venv in ~/opt/physics/venv and activate it
python3 -m venv $HOME/opt/physics/venv
source $HOME/opt/physics/venv/bin/activate

# build openEMS Python extension
# (CSXCAD_INSTALL_PATH and OPENEMS_INSTALL_PATH must be set!)
export CSXCAD_INSTALL_PATH=$HOME/opt/physics
export OPENEMS_INSTALL_PATH=$HOME/opt/physics
pip3 install .
```

Replace `$HOME/opt/physics` with the path prefix to CSXCAD/openEMS.
Both projects should be installed to the same prefix (installation
to different prefixes are unsupported and untested).

Once installed, test openEMS from a neutral directory. Don't test
openEMS in the source code directory `openEMS/python` to avoid importing
local files.

```bash
$ cd /  # Important: always leave "openEMS/python" first.

$ cd / && python3 -c "import openEMS; print(openEMS.openEMS())"
<openEMS.openEMS.openEMS object at 0x7f47f8dffb20>

$ cd / && python3 -c "import openEMS; print(openEMS.__version__)"
'0.0.36.post1.dev115+gfbb03a107.d20251112'
```

## Environment Variables

The following environment variables control the behaviors of the
Python extension installation.

1. **(Required)** `CSXCAD_INSTALL_PATH`, `OPENEMS_INSTALL_PATH`:
path prefix of the CSXCAD and openEMS installation. Without these
variables, installation is terminated with an error.

2. **(Optional)** `CSXCAD_PYSRC_PATH`: path to the CSXCAD Python
source code. By default, it's auto-detected by checking a few files
in the directory structure `openEMS-Projects`. If they don't exist,
CSXCAD source is auto-downloaded from GitHub. It can be overridden
with a filesystem path or a `git+https://` URL if auto-detection fails.

3. **(Optional)** `VIRTUAL_ENV`: path prefix of the Python `venv`,
set automatically if a Python `venv` is activated. If `venv` exists
in the C++ path prefix's `/venv` subdirectory
(`VIRTUAL_ENV=$OPENEMS_INSTALL_PATH/venv`), or overlaps
(`VIRTUAL_ENV=$OPENEMS_INSTALL_PATH`),
both `CSXCAD_INSTALL_PATH` and `OPENEMS_INSTALL_PATH` can be omitted,
activating the venv is sufficient.

4. **(Optional)** `CSXCAD_INSTALL_PATH_IGNORE`,
`OPENEMS_INSTALL_PATH_IGNORE`: disable `CSXCAD_INSTALL_PATH`
and `OPENEMS_INSTALL_PATH` usages and error checking. Useful only
if their installation paths are specified manually through other
methods, such as `CXXFLAGS` or `LDFLAGS`.

5. **(Optional)** `CSXCAD_NOSCM` and `OPENEMS_NOSCM`: pip no
longer downloads `setuptools_scm`, git-based version numbers are
no longer generated.

If build isolation is disabled (see below), all `CSXCAD_` variables
are made optional, because CSXCAD won't be rebuilt by pip if it has
already been installed.

## Basic Install

### Step 1: Set `CSXCAD_INSTALL_PATH` and `OPENEMS_INSTALL_PATH`

The environment variable `CSXCAD_INSTALL_PATH` and `OPENEMS_INSTALL_PATH`
must be set to ensure a successful installation. If either of them is not
set, a `RuntimeError` is generated.

If CSXCAD/openEMS were installed into `~/opt/physics`, set
`CSXCAD_INSTALL_PATH` and `OPENEMS_INSTALL_PATH` with:

```bash
export CSXCAD_INSTALL_PATH="$HOME/opt/physics"
export OPENEMS_INSTALL_PATH="$HOME/opt/physics"
```

Replace `$HOME/opt/physics` with the path prefix to CSXCAD/openEMS.
Both projects should be installed to the same prefix (installation
to different prefixes are unsupported and untested).

You
should be able to find `/lib` and `/include` in this prefix:

```bash
$ ls $OPENEMS_INSTALL_PATH
bin  include  lib  lib64  share

$ ls $OPENEMS_INSTALL_PATH/include
CSXCAD  fparser.hh  openEMS

$ ls $OPENEMS_INSTALL_PATH/lib
libCSXCAD.so    libCSXCAD.so.0.6.3  libfparser.so.4      libnf2ff.so    libnf2ff.so.0.1.0  libopenEMS.so.0
libCSXCAD.so.0  libfparser.so       libfparser.so.4.5.1  libnf2ff.so.0  libopenEMS.so      libopenEMS.so.0.0.36
```

As a hardcoded special case, the path of the current Python venv
(`VIRTUAL_ENV`) is also considered as a search path prefix by default.
If your C++ and Python `venv` paths exactly overlap, one doesn't need
to set any environment variables if a Python venv is activated prior
to installation. We still use `$CSXCAD_INSTALL_PATH` and
`$OPENEMS_INSTALL_PATH` throughout the documentation for consistency.

### Step 2: venv

Installing Python packages into Python's default search paths
(such as `/usr/` in the base system, or `~/.local` in the home
directory) is discouraged by most operating systems, because
there's a risk of dependency conflicts between a system-supplied
and a user-installed Python package.

To ensure Python packages are installed in a conflict-free manner,
it's suggested by most systems to create an isolated environment
for Python packages, known as a *virtual environment* (`venv`).

You should already have created a `venv` while installing CSXCAD,
skip the following commands. But if not, create a `venv`
specifically for CSXCAD/openEMS now:

```bash
python3 -m venv $HOME/opt/physics/venv/
```

This creates the Python venv in a subdirectory of `$OPENEMS_INSTALL_PATH`.
But if you prefer separation, you can use a different path, such
as `~/venvs/openems`, or activate an existing `venv` you already
have.

Remember, if the Python extension has been installed to an
isolated `venv`, all Python scripts that use CSXCAD
can only be executed inside this venv while it's activated.
Likewise, only Python packages installed into the venv can be
seen.

The venv can be entered via:

```bash
source $HOME/opt/physics/venv/bin/activate

# leave the venv with "deactivate"
```

Once the venv is activated, follow the next steps.

### Step 3: pip

Assuming that the correct `CSXCAD_INSTALL_PATH` and `OPENEMS_INSTALL_PATH`
have already been set, and a `venv` has been activated, run:

```bash
pip3 install .
```

When installing packages inside a `venv`, avoid using `--user` because
it doesn't respect the activated `venv`, effectively undoing it.

## Advanced Install

This section is written for experienced users, sysadmins and
developers who need to troubleshoot or customize their installation.
Ordinary users can skip this section.

### Suppress `RuntimeError`

CSXCAD/openEMS is usually installed to a non-standard location such as
`~/opt/physics` in the user home directory. By default, system compilers
are unable to find necessary C++ libraries, because only global path
prefixes such as `/usr/` or `/usr/local` are considered.

A `RuntimeError` is generated if `CSXCAD_INSTALL_PATH` or
`OPENEMS_INSTALL_PATH` are not set. If you know what you're doing
(e.g., both libraries are already added to the compiler search paths
manually), you can bypass these errors with:

```bash
export CSXCAD_INSTALL_PATH_IGNORE=1
export OPENEMS_INSTALL_PATH_IGNORE=1
```

### Search Path Management

By default, all necessary search paths are configured automatically
by `CSXCAD_INSTALL_PATH` and `OPENEMS_INSTALL_PATH`, or `VIRTUAL_ENV`.
On all Unix-like system, `/usr/local` is also always added into the search
paths as another hardcoded special case, regardless of the system's
default search paths.  Likewise, on macOS, the prefix reported by
`brew --prefix` is automatically added to the search paths.

If manual control is needed, set `CXXFLAGS` and `LDFLAGS` instead of
`CSXCAD_INSTALL_PATH` or `OPENEMS_INSTALL_PATH`. These flags include
the following arguments:

* `-I`: header include path, including the `/include` suffix.
* `-L`: library linking path, including the `/lib` suffix.
* `-Wl,-rpath,`: library runtime path, including the `/lib` suffix.

The following example assumes the installation prefix is
`$HOME/opt/physics`, and some dependent libraries have been installed
to `/usr/local`.

```bash
export CSXCAD_INSTALL_PATH_IGNORE=1
export OPENEMS_INSTALL_PATH_IGNORE=1

export CXXFLAGS="-I$HOME/opt/physics/include -I/usr/local/include $CXXFLAGS"
export LDFLAGS="-L$HOME/opt/physics/lib -L/usr/local/lib -Wl,-rpath,$HOME/opt/physics/lib $LDFLAGS"
```

To use these options properly, one needs to understand the motivation
behind specifying them. Basically, building a Python module requires
headers and libraries from three distinct sources:

1. Standard global headers and libraries provided by the system,
   and used by compilers by default. Typical paths are `/usr/include`
   and `/usr/lib`.
   They paths *do not* need any special listing, since they're used
   by default. All dependencies installed by the system's package
   manager typically also belong to this category, without special
   treatment (but exceptions exist, such as macOS Homebrew).

2. Non-standard global headers and libraries installed by the user
   (usually dependencies such as a custom Boost or VTK newer than the
   system's own version). They're outside the system's control, and
   not used by compilers by default.
   For example, on CentOS, the paths
   `-L/usr/local/include` and `-Wl,-rpath,/usr/local/lib` *must be*
   listed if any custom packages are installed to `/usr/local`.
   Unlike other system-wide package managers, macOS's Homebrew
   also belong to this category, because it's a 3rd-party package
   manager, thus it requires
   `-L$(brew --prefix)/include` and `-Wl,-rpath,$(brew --prefix)/lib`.

3. Non-standard local CSXCAD/openEMS headers and libraries.
   These files are usually installed to an arbitrary prefix in the user's
   home directory, not used by any compilers by default, as such
   `-L$HOME/opt/physics/include` and `-Wl,-rpath,$HOME/opt/physics/lib`.
   These paths *must be* listed.

If multiple paths are needed, repeat the option for each path, and
separate each option by spaces.

### Expose System Packages in venv

In an isolated `venv`, only Python packages installed into the `venv` can
be seen. Optionally, one can expose external system-wide packages to a
`venv` via `--system-site-packages` during `venv` creation:

```bash
# create venv, expose system packages
# (run once during installation)
python3 -m venv --system-site-packages $HOME/opt/physics/venv/
```

In this `venv`, the packages within `venv` stays within
the `venv`, but system-wide packages are also available. Activation
is still needed prior to using CSXCAD/openEMS in Python.

```bash
source $HOME/opt/physics/venv/bin/activate
```

### Install Python Extension to Home Directory Instead of venv

It's possible to install a package via `pip` into the default
Python search path under a user's home directory via `--user`.

After 2021, this practice is deprecated on most systems by
[PEP 668](https://peps.python.org/pep-0668/), since it bypasses
a system's own package manager, risking dependency conflicts.
The new option `--break-system-packages` is required.

```bash
pip3 install . --user --break-system-packages
```

This is still the default behavior of CSXCAD and openEMS when using
the `./update_openEMS.sh` script from `openEMS-Project.git`. If this
is undesirable, activate a `venv` manually prior to calling it.
We plan to improve the script in the future.

As suggested by the option `--break-system-packages`, it has the
risk of creating dependency conflicts between the same package
from the system and from `pip`. Using `--break-system-packages` is
only considered safe if all Python dependencies are installed via
your system's package manager (e.g. `apt`, `dnf`), as recommended in
the [documentation](https://openems.readthedocs.io/en/latest/install/requirements.html),
prior to running `pip3 install .`. Otherwise, `pip` may attempt to
install dependent packages on its own, risking dependency conflicts
with system packages.

### Installation without Re-downloading

By default, `pip` prefers to ignore existing packages in the
system, aggressively redownloading them, either for constructing
a fresh user `venv`, or for constructing the isolated build
environment. This includes building CSXCAD twice. It can be
problematic if you want to manage packages external to `pip`,
or if Internet access is not always online.

It's possible to suppress most redownloading behaviors, making
it a useful solution for external package management.

#### Install Dependencies

To manage packages manually, ensure that all dependencies have
been installed via your system's package manager. In theory,
one can use a DVD as the software repository. A full list
of package manager dependencies on various systems can be
found in the [documentation](https://openems.readthedocs.io/en/latest/install/requirements.html):

#### Expose System Packages to `venv`

Create a `venv`, and expose existing system packages, so we don't
need to install anything into the venv if the system already has
them (only needed to run once during installation).

```bash
python3 -m venv --system-site-packages $HOME/opt/physics
```

#### Disable Build Isolation

During installation via `pip`, `pip` will redownload build-time
dependencies such as `setuptools`, or `cython`, even when those
dependencies are already available to the system or a `venv`.
In particular, CSXCAD is built twice for this reason. This is
why one has to set both `CSXCAD_` and `OPENEMS_` environment
variables, even if `CSXCAD_` has already been installed.

This is the result of the *build isolation* feature in `pip`.
When building packages, `pip` creates an internal `venv` for itself,
isolated from both the base system and a user's own `venv`. This way,
it allows users to install packages with conflicting build-time
dependencies.

By default, CSXCAD/openEMS also uses `setuptools_scm` to automatically
create a version number based on the current `git` history. Since a
fairly new version is required, pre-installation via a system's package
manager may be impractical.

Both behaviors can disabled via:

```bash
export OPENEMS_NOSCM=1
pip3 install . --no-build-isolation
```

The variable `OPENEMS_NOSCM` is specific to openEMS.

The `CSXCAD_` variables are *not* needed. Without build isolation,
CSXCAD won't be rebuilt again while installing openEMS, the
existing installation is used.

### Override `CSXCAD_PYSRC_PATH`

Note: If *Build Isolation* has been disabled, this section is not
needed. Without build isolation, the CSXCAD Python extension which
you've already installed will be used. Just ensure to install
`CSXCAD/python` before `openEMS/python`.

Internally, pip requires source code to CSXCAD's Python extension
while building openEMS. Due to build isolation, CSXCAD is built
twice by pip: once for pip's internal build use, once for actual
end-user use.

Its path is automatically detected by assuming the following
directory structure with respect to `CSXCAD/python`.

    ├── CSXCAD
    │   ├── python
    ├── openEMS
    │   ├── openems.h
    │   ├── python
    └── update_openEMS.sh

If detection fails, we assume `openEMS` is used in isolation,
so a fallback URL `git+https://github.com/thliebig/CSXCAD.git#subdirectory=python`
is used instead.

The auto-detection and the fallback GitHub URL can be overridden
simultaneously by the environment variable `CSXCAD_PYSRC_PATH`,
such as:

    # use a local copy of CSXCAD
    export CSXCAD_PYSRC_PATH="$HOME/CSXCAD/python"

    # use a different Git repo
    export CSXCAD_PYSRC_PATH="git+https://example.com/CSXCAD.git#subdirectory=python"

### Legacy Installation via `setup.py`

`setup.py` method was the traditional way of building and installing
Python extensions. It has been deprecated by Python developers in favor
of `pip`. Follow this section only for the purpose of debugging a build.

Assuming that the correct `OPENEMS_INSTALL_PATH` have already been set
(or have been bypassed via `OPENEMS_INSTALL_PATH_IGNORE`), this extension
can be built manually via:

```bash
export OPENEMS_INSTALL_PATH=$HOME/opt/physics

python setup.py build_ext

# install to user's home directory, equivalent to
# pip3 install . --user --break-system-packages
python setup.py install --user

# if using a venv, remove --user so the venv path is respected
# python setup.py install
```

The `CSXCAD_` variables are *not* needed. Without build
isolation, CSXCAD won't be rebuilt again while installing
openEMS, the existing installation is used.

#### Advanced: setup.py search path management

On Unix-like systems, one can use the standard `CXXFLAGS` and
`LDFLAGS` to control compiler headers and libraries paths for
both `pip` and `setup.py` in the beginning of this tutorial.

In `setup.py`, it also provides its own custom options. Their
uses are not necessary, they're introduced there:

* `-I`: header include path, including the `/include` suffix, colon-separated.
* `-L`: library linking path, including the `/lib` suffix, colon-separated.
* `-R`: library runtime path, including the `/lib` suffix, colon-separated.

The following example assumes the CSXCAD/openEMS installation prefix is
`$HOME/opt/physics/`, and some libraries have been installed to
`/usr/local`.

```bash
export OPENEMS_INSTALL_PATH_IGNORE=1

python3 setup.py build_ext \
  -I "$HOME/opt/physics/include:/usr/local/include" \
  -L "$HOME/opt/physics/lib:/usr/local/lib" \
  -R "$HOME/opt/physics/lib"
```

## Troubleshooting

### ModuleNotFoundError: No module named 'openEMS.openEMS'

Ensure you're not running `python` under `openEMS/python` of the
source code tree. Otherwise, Python attempts to import the incomplete
source code instead of the complied Python extension.

```bash
cd /

source $HOME/opt/physics/venv/bin/activate
python3 -m "import openEMS"
```

If the error persists, debug the installation by running `pip` in
verbose mode.

```bash
export CSXCAD_INSTALL_PATH="$HOME/opt/physics"
export OPENEMS_INSTALL_PATH="$HOME/opt/physics"

pip3 install . --verbose
```

The `setup.py` method can also be used for troubleshooting.

If you are unable to solve the problem, create a post in the
[discussion forum](https://github.com/thliebig/openEMS-Project/discussions).
Make sure to provide detailed information about your system
(operating systems name and version, any error messages and
debugging outputs).

### FileNotFoundError: [Errno 2] No such file or directory: '/CSXCAD/python'

If you see a similar error message during installation:

```
ERROR: Could not install packages due to an OSError.
FileNotFoundError: [Errno 2] No such file or directory: '/CSXCAD/python'
error: subprocess-exited-with-error
```

It means the installer detected an incorrect CSXCAD Python source code path,
and it's unable to install the CSXCAD Python extension as a dependency.

Rerun pip with `pip install . --no-build-isolation`. Without
build isolation, the CSXCAD Python extension which you've already installed
will be used. Just ensure to install `CSXCAD/python` before `openEMS/python`.

Alternatively, provide the path manually via `CSXCAD_PYSRC_PATH` instead, and
rerun pip. See the section *Override `CSXCAD_PYSRC_PATH`* in the documentation.

### error: invalid command `bdist_wheel`

Under `--no-build-isolation`, one may encounter the following error:

    creating '/tmp/pip-modern-metadata-laicgdq2/CSXCAD.dist-info'
    error: invalid command 'bdist_wheel'
    Preparing metadata (pyproject.toml) ... error

On some systems, the command `bdist_wheel` is not provided by `setuptools`
but an additional package named `wheel`, such as `python3-wheel`.

This shouldn't happen if all dependencies listed in the latest documentation
have been correctly installed. If you need to manually install the package
`wheel` not already mentioned by the documentation, it means the documentation
is outdated, please submit a bug report.

### Unable to detect CSXCAD's Python source code path

If you see the following error message during installation:

```
  File "/home/fdtd/openEMS/python/bootstrap/setuptools_build_meta_custom.py", line 44, in add_csxcad
    raise RuntimeError(
      RuntimeError: Unable to detect CSXCAD's Python source code path. You're likely using an old pip without 'in-tree build' support. You can pick one solution below: (1) Rerun pip with 'pip install . --no-build-isolation' if CSXCAD Python extension is already installed (recommended). (2) Provide the path via CSXCAD_PYSRC_PATH and rerun pip (e.g. 'export CSXCAD_PYSRC_PATH=/home/user/openEMS-Project/CSXCAD/python/ && pip install . '). (3) Upgrade to pip 21.3 or newer.
      [end of output]

  note: This error originates from a subprocess, and is likely not a problem with pip.
```

It means the auto-CSXCAD detection result is ambiguous. You can choose
one solution below.

1. Choice 1: rerun pip with `pip install . --no-build-isolation`.

   Without build isolation, the CSXCAD Python extension which
   you've already installed will be used. Just ensure to install
   `CSXCAD/python` before `openEMS/python`. This is recommended,
   as it's the easiest solution.

2. Choice 2: provide the path via `CSXCAD_PYSRC_PATH` and rerun
   pip. See the section *Override `CSXCAD_PYSRC_PATH`* in the
   documentation.

3. Choice 3: upgrade to pip 21.3.

   ```bash
   # Activate the Python venv first (important) if you didn't.
   source $HOME/opt/physics/venv/bin/activate

   # Upgrade pip
   pip3 install --upgrade pip

   # Check version
   pip3 --version
   ```

   If the system's Python interpreter is older than Python 3.6,
   pip 21.3 cannot be installed. Upgrade your system to Python 3.6,
   or try an alternative solution above. Note that the openEMS
   extension is not fully-functional under Python 3.5 and below,
   `SyntaxError` may occur.
