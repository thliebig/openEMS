# -*- coding: utf-8 -*-
"""
Rectangular resonant cavity

Simulates a PEC rectangular cavity excited by a Dirac delta pulse and
extracts the resonance frequencies via FFT.  Results are compared with the
analytical TM/TE mode formula.

Three orthogonally oriented lumped ports are placed in a single simulation so
that all resonance modes are excited regardless of polarisation.

Reference: https://en.wikipedia.org/wiki/Microwave_cavity#Rectangular_cavity

Tested with
  - Python 3.8+
  - openEMS v0.0.36+
"""

import os, tempfile
import numpy as np
import matplotlib.pyplot as plt

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
post_proc_only = False
# show the structure before simulating; set True here, or run with SHOW_STRUCTURE=1
show_structure = os.environ.get('SHOW_STRUCTURE', '0') == '1'

unit        = 1                          # metres
cavity_size = np.array([2.0, 3.0, 7.0]) # m

# Highest frequency to analyse: one wavelength across the shortest side
f_max = C0 / cavity_size.min()

# Frequency resolution: keep <1 % of the lowest resonance spacing
frequency_resolution = C0 / cavity_size.max() / 100

# Target ~30 cells per wavelength at f_max for uniform spatial resolution
cell_size = C0 / f_max / 30   # m

Sim_path = os.path.join(tempfile.gettempdir(), 'RectCavity')

# ---------------------------------------------------------------------------
# Analytical resonance frequencies
# ---------------------------------------------------------------------------

def _mode_frequency(cavity_size, mnp, E_dir):
    """Resonance frequency for mode mnp=(m,n,p) with given E-field direction."""
    m, n, p = mnp
    if m + n + p in (0, 1):
        return np.nan
    if E_dir == 'x' and (n == 0 or p == 0):
        return np.nan
    if E_dir == 'y' and (m == 0 or p == 0):
        return np.nan
    if E_dir == 'z' and (m == 0 or n == 0):
        return np.nan
    return C0 / 2 * np.sqrt(np.sum((mnp / cavity_size) ** 2))


def theoretical_resonances(cavity_size, f_max, E_dir, n_max=5):
    """Sorted array of analytical resonance frequencies below f_max."""
    freqs = []
    for m in range(n_max):
        for n in range(n_max):
            for p in range(n_max):
                f = _mode_frequency(cavity_size, np.array([m, n, p]), E_dir)
                if not np.isnan(f) and f <= f_max:
                    freqs.append(f)
    return np.unique(freqs)

# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

def run_simulation(sim_path, cavity_size, f_max, sim_time):
    os.makedirs(sim_path, exist_ok=True)

    CSX  = ContinuousStructure()
    FDTD = openEMS()
    FDTD.SetCSX(CSX)
    FDTD.SetBoundaryCond(['PEC'] * 6)
    FDTD.SetDiracExcite(f_max)
    FDTD.SetMaxTime(sim_time)
    #FDTD.SetOverSampling(2000)

    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    for i, ax in enumerate('xyz'):
        n_lines = max(int(np.ceil(cavity_size[i] / cell_size)) + 1, 5)
        mesh.AddLine(ax, np.linspace(0, cavity_size[i], n_lines).tolist())

    # Three ports, one per axis, placed off-centre to couple to many modes.
    # R≫cavity impedance → near-ideal voltage source with negligible loading.
    port_start = np.array([mesh.GetLine(n, np.argmin(abs(mesh.GetLines(n)-cavity_size[n]/3))) for n in range(3)])
    for n, E_dir in enumerate(['x', 'y', 'z']):
        port_size    = [0,0,0]
        port_size[n] = max(cavity_size)*0.1
        FDTD.AddLumpedPort(
            port_nr = n+1,
            R       = 1e99,
            start   = port_start,
            stop    = port_start + port_size,
            p_dir   = E_dir,
            excite  = 1,
        )

    CSX.Write2XML(os.path.join(sim_path, 'model.xml'))
    ### Show the structure (see show_structure above)
    if show_structure:
        view = dict(view='iso', grid=dict(outline=True))
        try:
            from CSXCAD.Viewer import Show
        except ImportError:
            print('Showing the structure needs vtk: pip install vtk')
        else:
            Show(CSX, **view)

    FDTD.Run(sim_path, cleanup=True)


def load_port_voltage(sim_path, port_nr):
    """Read the time-domain port voltage written by openEMS."""
    data = np.loadtxt(os.path.join(sim_path, f'port_ut_{port_nr}'), comments='%')
    return data[:, 0], data[:, 1]   # time (s), voltage (V)

# ---------------------------------------------------------------------------
# Post-processing and plotting
# ---------------------------------------------------------------------------

def plot_results(cavity_size, f_max, dirs, time_data, voltage_data):
    colors      = ['tab:blue', 'tab:orange', 'tab:green']
    line_styles = ['--', '-.', ':']

    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    size_str = ' × '.join(f'{v:g}' for v in cavity_size)
    fig.suptitle(f'Rectangular resonant cavity  {size_str} m')

    # Time domain
    ax = axes[0]
    for i, E_dir in enumerate(dirs):
        ax.plot(time_data[E_dir] * 1e9, voltage_data[E_dir],
                color=colors[i], label=f'E_{E_dir}', alpha=0.8)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Voltage (V)')
    ax.legend()
    ax.set_title('Time domain')

    # Frequency domain
    ax = axes[1]
    for i, E_dir in enumerate(dirs):
        t, v = time_data[E_dir], voltage_data[E_dir]
        dt    = t[1] - t[0]
        freqs = np.fft.rfftfreq(len(v), dt)
        spec  = np.abs(np.fft.rfft(v))
        mask  = (freqs <= f_max) & (freqs>0)
        ax.semilogy(freqs[mask] / 1e6, spec[mask], color=colors[i], label=f'E_{E_dir}')

        for f_res in theoretical_resonances(cavity_size, f_max, E_dir):
            ax.axvline(f_res / 1e6, color=colors[i],
                       linestyle=line_styles[i], linewidth=0.8, alpha=0.5)

    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('|FFT| (arb.)')
    ax.legend()
    ax.set_title('Frequency domain  (dashed lines = analytical modes)')

    plt.tight_layout()
    plt.show()

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
dirs     = ['x', 'y', 'z']
sim_time = 1.0 / frequency_resolution

if not post_proc_only:
    print('Running simulation ...')
    run_simulation(Sim_path, cavity_size, f_max, sim_time)

time_data, voltage_data = {}, {}
for port_nr, E_dir in enumerate(dirs, start=1):
    time_data[E_dir], voltage_data[E_dir] = load_port_voltage(Sim_path, port_nr)

plot_results(cavity_size, f_max, dirs, time_data, voltage_data)
