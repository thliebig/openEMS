.. _python_results:

Setup File and Results
----------------------

Post-processing a simulation without rerunning its script: the XML holds
the structure, a JSON file next to it (``model.json`` for ``model.xml``)
the ports and NF2FF boxes, and the results to compute and plot.

.. code-block:: python

    FDTD.Write2XML(os.path.join(Sim_Path, 'model.xml'))   # writes model.json too
    FDTD.Run(Sim_Path)

    # later, in a new session
    FDTD = openEMS()
    FDTD.ReadFromXML(os.path.join(Sim_Path, 'model.xml'))
    FDTD.ReadSetup(os.path.join(Sim_Path, 'model.json'))
    port = FDTD.GetPort(1)

``ReadSetup`` adds nothing to the structure: each port is created again
against the loaded mesh, and every probe and excitation it would add has to
exist in the loaded structure already, unchanged. Ports created directly
from :mod:`openEMS.ports` are added with ``FDTD.AddPort()``;
:class:`~openEMS.setup_file.Setup` does the same without an openEMS object.

:class:`~openEMS.results.Results` describes what to compute and plot. A
quantity (``res.S(2, 1)``, ``res.Zin(1)``, ``res.FarFieldCut(...)``,
``res.DumpSlice(...)``) computes its data on first use and caches it;
figures collect quantities as curves or maps:

.. code-block:: python

    from openEMS.results import Results

    res = Results(FDTD, Sim_Path)            # frequencies: the excitation band
    fig = res.AddFigure('S-Parameter')
    fig.AddCurve(res.S(1, 1), fmt='dB')
    fig.AddCurve(res.S(2, 1), fmt='dB')

    pat = res.AddFigure('Pattern', kind='polar')
    pat.AddCurve(res.FarFieldCut('nf2ff', f0, phi=0), fmt='dBi')

    sar = res.AddFigure('SAR', grid=(1, 2), shared_clim=True)
    sar.AddMap(res.DumpSlice('SAR.h5', ('z', 0), freq=f0), norm=res.P_acc(1))

    f, s21 = res.S(2, 1).GetXY()             # the numbers
    print('Dmax: {:.1f} dBi'.format(10*np.log10(res.Dmax('nf2ff', f0).GetValue())))
    res.Write(os.path.join(Sim_Path, 'model.json'))
    res.Show()                               # matplotlib windows

With the CSXCAD viewer (``vtk``, and ``PySide6`` for its GUI) the figures
become tabs next to the structure, in one window with field displays:

.. code-block:: python

    from CSXCAD.Viewer import Viewer
    from openEMS.Viewer import FieldSlices

    v = Viewer(CSX)
    v.AddItem(FieldSlices(os.path.join(Sim_Path, 'Et.h5'), quantity='E_z', freq=f0))
    res.AddToViewer(v)
    v.Show()

.. automodule:: openEMS.results

    .. autoclass:: Results
        :members:

    .. autoclass:: Figure
        :members: AddCurve, AddMap, Draw

    .. autoclass:: Quantity
        :members: GetXY, GetMap, GetValue, At

    .. autofunction:: Format

.. automodule:: openEMS.setup_file

    .. autoclass:: Setup
        :members: AddPort, AddNF2FF, Write, Read
