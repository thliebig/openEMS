.. _tutorial_horn_antenna:

Horn Antenna with Coaxial Pin Feed
===================================

* A pyramidal horn antenna fed by a coaxial probe (pin) inside the rectangular feed waveguide, terminated by a quarter-wave back-short instead of a waveguide port.

Introduction
-------------
**This tutorial covers:**

* Setup of a pyramidal horn antenna with a closed rectangular feed waveguide
* A coaxial pin feed modelled as a lumped port, exciting the TE10 waveguide mode
* A quarter guided-wavelength back-short instead of a waveguide port or PML termination
* Calculate the S-Parameter and input impedance
* Calculate the far-field pattern, directivity and aperture efficiency via NF2FF

Python Script
-------------
Get the latest version `from git <https://raw.githubusercontent.com/thliebig/openEMS/master/python/Tutorials/Horn_Antenna.py>`_.

.. include:: ./__Horn_Antenna.txt

Images
-------------
.. figure:: images/Horn_Ant.png
    :width: 80%
    :alt: 3D view of the horn antenna

    3D view of the Horn Antenna with coaxial pin feed (AppCSXCAD)

.. figure:: images/Horn_Ant_SPara.png
    :width: 80%
    :alt: S-Parameter

    S-Parameter of the horn antenna

.. figure:: images/Horn_Ant_Pattern.png
    :width: 80%
    :alt: Farfield pattern

    Farfield pattern on an E- and H-plane

