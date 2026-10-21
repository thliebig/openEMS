.. _python_viewer:

Field Viewer
------------

Field results shown in 3D in the CSXCAD viewer (:ref:`csxcad_viewer`):
cutting planes through field and SAR dumps, iso-surfaces and far-field
patterns. Each is a viewer item: add it with ``Viewer.AddItem()``, change
its settings in the GUI or with ``Viewer.SetItemParameter()``; *File → Copy
view as Python* records them in the view state. See :ref:`concept_viewer`
for an introduction.

.. code-block:: python

    from CSXCAD.Viewer import Viewer
    from openEMS.Viewer import FieldSlices, FarField

    v = Viewer(CSX)
    v.AddItem(FieldSlices(os.path.join(Sim_Path, 'Ef.h5'), planes=[('z', 1.5), ('x', 0)],
                          quantity='E_z', freq=2.4e9))
    v.AddItem(FarField(nf2ff_res, quantity='directivity', scale='log'))
    v.Show()

The data comes from openEMS HDF5 dumps (``file_type=1``), from
``CalcNF2FF`` results on a full theta/phi grid, or from numpy arrays
(:meth:`FieldSource.FromArrays <openEMS.Viewer.Source.FieldSource.FromArrays>`).
Frequency-domain fields animate over the phase, time-domain fields over the
recorded steps, far fields over frequency (``sweep=True``). Dumps hold the
fields of the simulation's excitation as is: frequency-domain values are
spectra of the excitation pulse, so their absolute numbers depend on it.
``norm`` divides them by a factor, e.g. the incident port voltage at the
frequency, with ``unit`` naming the result:

.. code-block:: python

    u_inc = port.uf_inc[np.argmin(np.abs(freq - f0))]
    v.AddItem(FieldSlices(Ef_file, quantity='E_y', norm=u_inc, unit='V/m/V'))

.. autoclass:: openEMS.Viewer.Slices.FieldSlices
    :members: GetPlanes, SetPlanes, AddPlane, RemovePlane

.. autoclass:: openEMS.Viewer.Pattern.FarField
    :members: Values, MainLobe

.. autoclass:: openEMS.Viewer.IsoSurface.FieldIsoSurface

.. autoclass:: openEMS.Viewer.Source.FieldSource
    :members: FromArrays, Read, ReadAtFrequency, Lines, Domains
