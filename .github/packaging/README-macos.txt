openEMS for macOS (Apple silicon) with the Metal GPU engine
==========================================================

Package: @PACKAGE_NAME@

openEMS runs directly from this folder, there is nothing to install. The
examples below assume the package was extracted to ~/openEMS, i.e. that
~/openEMS/bin/openEMS exists. Adjust the path if you chose another folder.

Requirements: macOS 15 or newer on Apple silicon. All other libraries are
included.

The package is not notarized: if macOS refuses to run openEMS, remove the
quarantine flag of the downloaded folder once:

    xattr -dr com.apple.quarantine ~/openEMS


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

Install the modules (Python 3.13 or 3.14 from python.org or a virtual
environment):

    python3 -m pip install numpy h5py matplotlib
    python3 -m pip install --no-index --find-links ~/openEMS/python openEMS

and tell them where the libraries are (e.g. in ~/.zshrc):

    export DYLD_LIBRARY_PATH=~/openEMS/lib

Check the installation:

    python3 -c "import openEMS; print(openEMS.__version__)"

Start scripts as "python3 script.py". macOS removes DYLD_LIBRARY_PATH for
programs started through /usr/bin, so a script started as "./script.py" with
a "#!/usr/bin/env python3" line does not find the libraries. lib/ holds only
the openEMS libraries, so the setting does not affect other programs.


GPU engine
----------

The Metal GPU engine runs on the GPU of Apple silicon Macs. Select it with
--engine=gpu, or in Python:

    FDTD.Run(sim_path, engine='gpu')


More
----

Project page:   https://openEMS.de
Documentation:  https://docs.openems.de
