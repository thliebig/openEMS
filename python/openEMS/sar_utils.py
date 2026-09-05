# -*- coding: utf-8 -*-

from openEMS.utilities import HDF5Dump

def readSAR(fn, f_idx=0):
    """Read a SAR result HDF5 file written by SAR_Calculation.

    Thin wrapper around `openEMS.utilities.HDF5Dump` for the common case of
    reading a complete SAR result.

    Parameters
    ----------
    fn : str
        Path to the SAR result HDF5 file.
    f_idx : int, optional
        Frequency index to read (default 0).

    Returns
    -------
    sar : ndarray, shape (nx, ny, nz)
        SAR values in W/kg.
    mesh : list of three 1-D ndarrays
        Mesh node coordinates [x, y, z] in metres.  This is ``mesh['lines']``
        from :meth:`HDF5Dump.GetMesh` -- a deliberate simplification for the
        common SAR use case (always Cartesian, coordinates already in metres).
        For full mesh metadata (type, scaling, names), for reading only a part
        of a large result, or for the ``/CellData`` and ``/CellWidth`` groups
        of a raw SAR dump, use :class:`HDF5Dump` directly.
    sar_data : dict
        Metadata from the file: 'mass' (kg), 'frequency' (Hz), 'power' (W),
        and any other dump attributes (e.g. 'maxSAR', 'dump_type').

    See Also
    --------
    openEMS.utilities.HDF5Dump
    """
    with HDF5Dump(fn) as dump:
        return (dump.GetFieldAtIndex(f_idx=f_idx),
                dump.GetMesh()['lines'],
                dump.GetAttributes(f_idx=f_idx))
