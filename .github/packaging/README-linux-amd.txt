openEMS for Linux (x86_64) with the GPU engine, for AMD GPUs
============================================================

Package: @PACKAGE_NAME@

openEMS runs directly from this folder, there is nothing to install. The
examples below assume the package was extracted to ~/openEMS, i.e. that
~/openEMS/bin/openEMS exists. Adjust the path if you chose another folder.

Requirements: x86_64 with glibc 2.39 or newer (e.g. Ubuntu 24.04, Debian 13,
Fedora 40 or newer). All other libraries are included.


Contents
--------

    bin/            openEMS, nf2ff and sar_calc
    lib/            the openEMS, CSXCAD, fparser and nf2ff libraries
    lib/deps/       the libraries they need (HDF5, VTK, Boost, ...)
    share/          Octave/Matlab interface and tutorials
    python/         Python modules (wheels) and tutorials


Command line
------------

    ~/openEMS/bin/openEMS <simulation.xml>
    ~/openEMS/bin/openEMS <simulation.xml> --engine=gpu


Python
------

Install the modules (Python 3.13 or 3.14, e.g. in a virtual environment):

    python3 -m pip install numpy h5py matplotlib
    python3 -m pip install --no-index --find-links ~/openEMS/python openEMS

and tell them where the libraries are (e.g. in ~/.bashrc):

    export LD_LIBRARY_PATH=~/openEMS/lib

Check the installation:

    python3 -c "import openEMS; print(openEMS.__version__)"

lib/ holds only the openEMS libraries, so the setting does not affect other
programs.


GPU engine
----------

The GPU engine is built with HIP, which compiles for both GPU vendors from one
source; this package is the build for AMD, with code for the CDNA 2 and 3
(MI200, MI300) and RDNA 2 and 3 (RX 6000, RX 7000) architectures. It needs the
amdgpu kernel driver of the host, as the NVIDIA build needs NVIDIA's; the ROCm
libraries it loads are in the package. The build for NVIDIA GPUs is the package
named cuda. Select the engine with --engine=gpu, or in Python:

    FDTD.Run(sim_path, engine='gpu')

Without a supported GPU, the GPU engine falls back to its reference backend
on the CPU, which is slower than the default engine.


More
----

Project page:   https://openEMS.de
Documentation:  https://docs.openems.de
