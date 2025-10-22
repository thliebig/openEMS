# -*- coding: utf-8 -*-

import numpy as np

def readSAR(fn, f_idx=0):
    sar_data = {}
    import h5py
    with h5py.File(fn, 'r') as h5:
        if 'openEMS_HDF5_version' in h5.attrs:
            sar = h5[f'/FieldData/FD/f{f_idx}']
            sar_data.update(sar.attrs)
            sar = np.array(sar)
            if h5.attrs['openEMS_HDF5_version'] <= 0.2:
                sar = sar.swapaxes(0,2)
            mesh = [None, None, None]
            for n, d in enumerate('xyz'):
                mesh[n] = np.array(h5['Mesh/'+d])
            return sar, mesh, sar_data

    return None, None, None
