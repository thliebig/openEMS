.. _tutorial_parallel_plate_wg:

Parallel Plate Waveguide
=========================

* The simplest possible openEMS simulation: a parallel-plate waveguide excited with a sinusoidal TEM mode, showing the core workflow of geometry setup, field dump and result visualization.

Introduction
-------------
**This tutorial covers:**

* FDTD setup with a sinusoidal excitation and mixed boundary conditions
* Geometry and mesh definition with CSXCAD
* A time-domain E-field dump written as VTK files
* Geometry inspection with AppCSXCAD and animation in Paraview

Python Script
-------------
Get the latest version `from git <https://raw.githubusercontent.com/thliebig/openEMS/master/python/Tutorials/Parallel_Plate_Waveguide.py>`_.

.. include:: ./__Parallel_Plate_Waveguide.txt

Visualizing the Results
------------------------

The simulation writes the E-field dump to ``Et_*.vtr`` in the simulation
directory. To animate the propagating wave in Paraview:

1. **File → Open** and select the ``Et_..vtr`` group.
2. Click **Apply** in the Properties panel.
3. Set **Color by** to ``E-Field``.
4. Press **Play** in the Animation toolbar.
5. Use **Rescale to Data Range** occasionally to tune the colour mapping.

For a clearer view of the wave propagation, apply a **Warp By Vector** filter
(Filters → Alphabetical → Warp By Vector, then Apply).

.. seealso::
    The same tutorial for the :ref:`Octave/Matlab interface <octave_tutorial_parallel_plate_wg>`.
