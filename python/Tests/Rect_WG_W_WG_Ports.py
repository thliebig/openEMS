"""
Rectangular waveguide (WR-42) with generic WaveguidePort using HDF5 mode files.

This is the test counterpart of Tutorials/Rect_Waveguide.py.  Instead of
AddRectWaveGuidePort (which generates the analytic mode functions internally),
it writes the identical TE10 field distributions to HDF5 mode files and uses
the generic AddWaveGuidePort.  Results should be indistinguishable.

The test is repeated for all three propagation directions (x, y, z) to verify
that the HDF5-based WaveguidePort handles arbitrary orientations correctly.
The same mode files are reused for every direction: the file coordinates
(file-x, file-y) always map to (ny_P, ny_PP) of the port.

Assertions (per direction):
  - S11 < -40 dB across the passband
  - S21 > -0.1  dB across the passband
  - mode purity > 99 % at every time step where the recorded signal exceeds
    1 % of its peak amplitude (checked for both U and I probes of both ports)

(c) 2026 Thorsten Liebig <thorsten.liebig@gmx.de>
"""

import os
import tempfile
import numpy as np
import h5py

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0
from openEMS.utilities import check_mode_purity


def make_rect_wg_te10_hdf5(filename, field_type, a_draw, b_draw, N=60):
    """Write a TE10 rectangular waveguide mode profile to an HDF5 file.

    Parameters
    ----------
    filename : str
        Output HDF5 file path.
    field_type : 'E' or 'H'
        Which field to store.
    a_draw, b_draw : float
        Waveguide width / height in drawing units.
    N : int
        Number of grid points along the width (height scaled proportionally).
    """
    x = np.linspace(0, a_draw, N)
    y = np.linspace(0, b_draw, max(int(N * b_draw / a_draw), 2))
    XX, _ = np.meshgrid(x, y, indexing='ij')

    if field_type == 'E':
        # TE10: E_x = 0,  E_y = -sin(pi*x/a) / a
        Vx = np.zeros_like(XX)
        Vy = -np.sin(np.pi * XX / a_draw) / a_draw
    else:
        # TE10: H_x = sin(pi*x/a) / a,  H_y = 0
        Vx = np.sin(np.pi * XX / a_draw) / a_draw
        Vy = np.zeros_like(XX)

    with h5py.File(filename, 'w') as f:
        f.attrs['Version'] = 1.0
        f.create_dataset('x',  data=x)
        f.create_dataset('y',  data=y)
        f.create_dataset('Vx', data=Vx)
        f.create_dataset('Vy', data=Vy)


# ---------------------------------------------------------------------------
# Geometry & simulation parameters  (WR-42, identical to tutorial)
# ---------------------------------------------------------------------------
unit   = 1e-6       # drawing unit: 1 µm
a      = 10700      # waveguide width  (µm) — along ny_P axis
b      = 4300       # waveguide height (µm) — along ny_PP axis
length = 50000      # waveguide length (µm) — along propagation axis

f_start = 20e9
f_0     = 24e9
f_stop  = 26e9
lambda0 = C0 / f_0 / unit

mesh_res = lambda0 / 30

kc = np.pi / (a * unit)     # TE10 cutoff wavenumber (rad/m)

# ---------------------------------------------------------------------------
# Mode files — written once, reused for all directions
# (file-x maps to ny_P, file-y maps to ny_PP regardless of exc_dir)
# ---------------------------------------------------------------------------
mode_dir = os.path.join(tempfile.gettempdir(), 'Rect_WG_HDF5_modes')
os.makedirs(mode_dir, exist_ok=True)

E_file = os.path.join(mode_dir, 'TE10_E.h5')
H_file = os.path.join(mode_dir, 'TE10_H.h5')
make_rect_wg_te10_hdf5(E_file, 'E', a, b)
make_rect_wg_te10_hdf5(H_file, 'H', a, b)


def run_direction(exc_dir):
    exc_ny = 'xyz'.index(exc_dir)
    ny_P   = (exc_ny + 1) % 3
    ny_PP  = (exc_ny + 2) % 3

    Sim_Path = os.path.join(tempfile.gettempdir(), 'Rect_WG_HDF5_{}'.format(exc_dir))
    os.makedirs(Sim_Path, exist_ok=True)

    # -------------------------------------------------------------------
    # FDTD & CSX setup
    # -------------------------------------------------------------------
    FDTD = openEMS()
    FDTD.SetGaussExcite(0.5*(f_start+f_stop), 0.5*(f_stop-f_start))

    # PML on the two propagation-axis faces, PEC on the four transverse faces
    bc = [0] * 6
    bc[2*exc_ny]   = 3   # PML_8
    bc[2*exc_ny+1] = 3
    FDTD.SetBoundaryCond(bc)

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)

    mesh.AddLine('xyz'[exc_ny], [0, length])
    mesh.AddLine('xyz'[ny_P],   [0, a])
    mesh.AddLine('xyz'[ny_PP],  [0, b])

    # -------------------------------------------------------------------
    # Ports — local_origin='corner' so mode-file (x,y) aligns with
    # the (ny_P=0, ny_PP=0) corner of the port box
    # -------------------------------------------------------------------
    def make_coords(x_near, x_far):
        start = [0, 0, 0]; stop = [0, 0, 0]
        start[exc_ny] = x_near;  stop[exc_ny] = x_far
        start[ny_P]   = 0;       stop[ny_P]   = a
        start[ny_PP]  = 0;       stop[ny_PP]  = b
        return start, stop

    ports = []

    start, stop = make_coords(10*mesh_res, 15*mesh_res)
    mesh.AddLine('xyz'[exc_ny], [start[exc_ny], stop[exc_ny]])
    ports.append(FDTD.AddWaveGuidePort(1, start, stop, exc_dir,
                                        E_file=E_file, H_file=H_file,
                                        kc=kc, excite=1, local_origin='corner'))

    start, stop = make_coords(length - 10*mesh_res, length - 15*mesh_res)
    mesh.AddLine('xyz'[exc_ny], [start[exc_ny], stop[exc_ny]])
    ports.append(FDTD.AddWaveGuidePort(2, start, stop, exc_dir,
                                        E_file=E_file, H_file=H_file,
                                        kc=kc, excite=0, local_origin='corner'))

    mesh.SmoothMeshLines('all', mesh_res, ratio=1.4)

    # -------------------------------------------------------------------
    # Run
    # -------------------------------------------------------------------
    FDTD.Run(Sim_Path, cleanup=True, exact_endcriteria=True)

    # -------------------------------------------------------------------
    # Post-processing
    # -------------------------------------------------------------------
    freq = np.linspace(f_start, f_stop, 201)
    for port in ports:
        port.CalcPort(Sim_Path, freq)

    for port in ports:
        lbl = 'port {}'.format(port.number)
        check_mode_purity(lbl + ' U', port.u_data.ui_val[0], port.u_mode_purity[0], threshold=0.99, sig_frac=0.001)
        check_mode_purity(lbl + ' I', port.i_data.ui_val[0], port.i_mode_purity[0], threshold=0.99, sig_frac=0.001)

    s11 = ports[0].uf_ref / ports[0].uf_inc
    s21 = ports[1].uf_ref / ports[0].uf_inc

    s11_dB = 20 * np.log10(np.abs(s11))
    s21_dB = 20 * np.log10(np.abs(s21))

    print('  S11 max = {:.1f} dB'.format(np.max(s11_dB)))
    print('  S21 min = {:.1f} dB'.format(np.min(s21_dB)))

    assert np.max(s11_dB) < -40, \
        'FAIL [{}]: S11 too high: {:.1f} dB'.format(exc_dir, np.max(s11_dB))
    assert np.min(s21_dB) > -0.1, \
        'FAIL [{}]: S21 too low: {:.1f} dB'.format(exc_dir, np.min(s21_dB))


for direction in ('x', 'y', 'z'):
    print('Testing direction: {}'.format(direction))
    run_direction(direction)
    print('PASS [{}]'.format(direction))
