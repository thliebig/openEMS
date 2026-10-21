# -*- coding: utf-8 -*-
#
# Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""openEMS.Viewer.CLI
The ``openems-view`` command: a finished simulation, without its script.

It reads the setup :meth:`openEMS.openEMS.openEMS.Run` wrote into the
simulation path (``openEMS_setup.xml`` and ``.json``), recreates the ports,
NF2FF boxes and the figures of :class:`openEMS.results.Results`, adds the
HDF5 field dumps found there as hidden :class:`FieldSlices` (shown with the
view key ``visible`` or in the tree) and the metal edges of ``--debug-PEC``
(``PEC_dump.vtp``), also hidden, and then opens the GUI or writes the
outputs, with the view and output options of ``csxcad-view``.
"""

import argparse
import glob
import json
import os
import sys
import warnings

DESCRIPTION = 'View a finished openEMS simulation: structure, field dumps and figures.'
EPILOG = """\
SIM_PATH is the simulation path (Run() writes openEMS_setup.xml and .json
there) or an XML file written by Write2XML. Field dumps are items below
results/fields, hidden at first; --describe lists them with their settings.
The metal edges of a run with debug_PEC=True (PEC_dump.vtp) are the hidden
item results/debug/PEC_dump.vtp.
The view options are those of csxcad-view (--list-view-keys explains them).

examples:
  openems-view /tmp/Simp_Patch --describe
  openems-view /tmp/Simp_Patch --figures figures/
  openems-view /tmp/Simp_Patch --png Et.png --set view=+z --set 'visible={"results/fields/Et": true}'
  openems-view /tmp/Simp_Patch --png Ez.png --set 'visible={"results/fields/Et": true}' \\
      --set 'items={"results/fields/Et": {"quantity": "z", "scale": "log"}}'
  openems-view /tmp/Simp_Patch --png pec.png --set isolate=type:Metal \\
      --set 'visible={"results/debug/PEC_dump.vtp": true}'
"""


def _Parser():
    from CSXCAD.Viewer.CLI import AddViewArguments
    p = argparse.ArgumentParser(prog='openems-view', description=DESCRIPTION, epilog=EPILOG,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('sim_path', nargs='?', help='simulation path, or an XML file written by Write2XML')
    AddViewArguments(p)
    p.add_argument('--figures', metavar='DIR', help='write every figure as PNG into DIR')
    p.add_argument('--no-fields', action='store_true', help='do not add the field dumps')
    return p


def Load(sim_path):
    """ Load(sim_path)

    The simulation in `sim_path` (or of an XML file written by
    ``Write2XML``): ``(FDTD, results, sim_dir)``, the openEMS object with
    its ports and NF2FF boxes, the :class:`~openEMS.results.Results` stored
    with it (None if none) and the directory of the result files.
    """
    from openEMS.openEMS import openEMS
    from openEMS import setup_file
    if os.path.isdir(sim_path):
        xml = os.path.join(sim_path, setup_file.RUN_SETUP)
        if not os.path.isfile(xml):
            raise FileNotFoundError('{}: no {}; Run() writes it (since this version), or give an XML file '
                                    'written by Write2XML'.format(sim_path, setup_file.RUN_SETUP))
    else:
        xml = sim_path
        if not os.path.isfile(xml):
            raise FileNotFoundError('{}: no such file or directory'.format(sim_path))
    sim_dir = os.path.dirname(os.path.abspath(xml))
    FDTD = openEMS()
    FDTD.ReadFromXML(xml)
    sidecar = setup_file.SetupFileName(xml)
    results = None
    if os.path.isfile(sidecar):
        FDTD.ReadSetup(sidecar)
        with open(sidecar) as f:
            stored = 'results' in json.load(f)
        if stored:
            from openEMS.results import Results
            results = Results.Load(FDTD, sim_dir, sidecar)
    return FDTD, results, sim_dir


def FieldDumps(CSX, sim_dir):
    """ The HDF5 files of the field (and SAR) dumps of `CSX` in `sim_dir`. """
    from CSXCAD.CSProperties import PropertyType
    out = []
    for prop in CSX.GetAllProperties():
        if prop.GetType() & PropertyType.DUMPBOX:
            fn = os.path.join(sim_dir, prop.GetName() + '.h5')
            if os.path.isfile(fn) and fn not in out:
                out.append(fn)
    # e.g. SAR_Calculation results next to them
    for fn in sorted(glob.glob(os.path.join(sim_dir, 'SAR*.h5'))):
        if fn not in out:
            out.append(fn)
    return out


def DebugFiles(sim_dir):
    """ The metal edges openEMS wrote into `sim_dir` with ``--debug-PEC``
    (one file per multigrid level). """
    return sorted(glob.glob(os.path.join(sim_dir, 'PEC_dump*.vtp')))


def Main(argv=None):
    """ Entry point of ``openems-view``; returns the exit code. """
    from CSXCAD.Viewer.CLI import (ViewFromArgs, HasOutputs, WriteOutputs, ErrorText, KeyList, ItemSettings,
                                   UnknownPaths, AddFiles)
    args = _Parser().parse_args(argv)
    if args.list_view_keys:
        print(KeyList())
        return 0
    if args.sim_path is None:
        print('openems-view: give the simulation path', file=sys.stderr)
        return 2
    outputs = HasOutputs(args) or bool(args.figures)
    want_window = (not outputs or args.gui) and not args.offscreen
    try:
        FDTD, results, sim_dir = Load(args.sim_path)
        from CSXCAD.Viewer import Viewer
        v = Viewer(FDTD.GetCSX(), offscreen=not want_window, name=os.path.basename(os.path.normpath(sim_dir)),
                   **ViewFromArgs(args))
        AddFiles(v, args)
    except (FileNotFoundError, ValueError, KeyError) as err:
        print('openems-view: {}'.format(ErrorText(err)), file=sys.stderr)
        return 1
    with ItemSettings() as wrong:   # the view's settings of the items added here
        if not args.no_fields:
            from openEMS.Viewer import FieldSlices
            for fn in FieldDumps(FDTD.GetCSX(), sim_dir):
                try:
                    v.AddItem(FieldSlices(fn), visible=False)
                except Exception as err:   # e.g. a dump of another layout
                    warnings.warn('{}: not shown: {}'.format(os.path.basename(fn), err))
        for fn in DebugFiles(sim_dir):
            try:
                v.AddFile(fn, group='results/debug', visible=False, render=False)
            except Exception as err:
                warnings.warn('{}: not shown: {}'.format(os.path.basename(fn), err))
    if results is not None:
        results.AddToViewer(v, store=False, skip_failed=True)
    errors = wrong.errors + UnknownPaths(v, ViewFromArgs(args))
    if errors:
        for e in errors:
            print('openems-view: {}'.format(e), file=sys.stderr)
        return 1
    try:
        if args.figures:
            os.makedirs(args.figures, exist_ok=True)
            for name, vf in v.GetFigures().items():
                fn = ''.join(c if c.isalnum() or c in '-_' else '_' for c in name)
                vf.GetFigure().savefig(os.path.join(args.figures, fn + '.png'))
        WriteOutputs(v, args, prog='openems-view')
    except (ValueError, KeyError) as err:
        print('openems-view: {}'.format(ErrorText(err)), file=sys.stderr)
        return 1
    if want_window:
        v.Show(gui=False if args.no_gui else None)
    return 0


if __name__ == '__main__':
    sys.exit(Main())
