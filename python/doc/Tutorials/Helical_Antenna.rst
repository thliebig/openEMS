Helical Antenna
===============

Introduction
-------------
**This tutorial covers:**

* setup of a helix using the wire primitive
* setup a lumped feeding port (R_in = 120 Ohms)
* adding a near-field to far-field (nf2ff) box using an efficient subsampling
* calculate the S-Parameter of the antenna
* calculate and plot the far-field pattern

Python Script
-------------
Get the latest version `from git <https://raw.githubusercontent.com/thliebig/openEMS/master/python/Tutorials/Helical_Antenna.py>`_.

.. include:: ./__Helical_Antenna.txt

Images
-------------
.. figure:: images/Helix_Ant.png
    :width: 49%
    :alt: alternate text
    
    3D view of the Helical Antenna (AppCSXCAD)

.. figure:: images/Helix_Ant_Pattern.png
    :width: 49%
    :alt: alternate text
    
    Far-Field pattern showing a right-handed circular polarization.
