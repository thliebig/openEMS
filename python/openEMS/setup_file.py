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
#
"""
The setup file: ports and NF2FF boxes of a simulation, written as JSON next
to its XML so that they can be recreated after the XML is read back.

The XML holds only the CSXCAD properties a port created, not the port. Each
entry here stores how a port was created (its type and arguments); reading
it runs the same constructor against a scratch structure with the loaded
mesh, checks that every property it creates exists unchanged in the loaded
structure, and then points the new port at those properties. Nothing is
added to the loaded structure.

Ports and NF2FF boxes remember how they were created (see
:func:`openEMS.utilities.SetupType`), however they were created. Normally
used through the openEMS object (:meth:`~openEMS.openEMS.openEMS.WriteSetup`,
:meth:`~openEMS.openEMS.openEMS.ReadSetup`), which knows what its `Add*Port`
methods and `CreateNF2FFBox` created. Without one, e.g. in a tool that only
builds a structure:

>>> setup = Setup()
>>> setup.AddPort(ports.MSLPort(CSX, 1, metal, start, stop, 'x', 'z', excite=1))
>>> setup.Write('model.json')            # next to CSX.Write2XML('model.xml')
...
>>> CSX = ContinuousStructure()
>>> CSX.ReadFromXML('model.xml')
>>> setup = Setup()
>>> setup.Read('model.json', CSX)
>>> port = setup.GetPort(1)
"""

import inspect
import json
import os
import tempfile
import xml.etree.ElementTree as ET

import numpy as np

from CSXCAD import ContinuousStructure
from CSXCAD.CSProperties import CSProperties
from openEMS import ports, nf2ff

FORMAT  = 'openEMS-setup'
VERSION = 1

#: the setup :meth:`openEMS.openEMS.openEMS.Run` writes into the simulation
#: path (the XML; the setup file next to it, see :func:`SetupFileName`)
RUN_SETUP = 'openEMS_setup.xml'

def _types(module):
    return {cls.SETUP_TYPE: cls for cls in vars(module).values()
            if isinstance(cls, type) and 'SETUP_TYPE' in vars(cls)}

# neutral type name --> class
PORT_TYPES = _types(ports)

# properties the user passes in and the port adds geometry to; post-processing
# never needs them, so they are neither stored nor checked
_USER_PROPS = {'metal_prop': 'AddMetal', 'pec_prop': 'AddMetal', 'mat_prop': 'AddMaterial'}
_PLACEHOLDER = '__openEMS_setup_{}'

# XML that is look or bookkeeping, not physics
_IGNORE_ATTRS = ('ID', 'Visible')
_IGNORE_ELEMS = ('FillColor', 'EdgeColor', 'Attributes')


def _to_json(val, where):
    if isinstance(val, CSProperties):
        raise TypeError('{}: a CSXCAD property cannot be stored'.format(where))
    if isinstance(val, np.ndarray):
        return [_to_json(v, where) for v in val.tolist()]
    if isinstance(val, (list, tuple)):
        return [_to_json(v, where) for v in val]
    if isinstance(val, dict):
        return {str(k): _to_json(v, where) for k, v in val.items()}
    if isinstance(val, np.generic):
        val = val.item()
    if isinstance(val, float) and not np.isfinite(val):
        return str(val)   # 'inf', '-inf', 'nan': strict JSON has no Infinity
    if val is None or isinstance(val, (bool, int, float, str)):
        return val
    raise TypeError('{}: cannot store a value of type {}'.format(where, type(val).__name__))


def _from_json(val):
    if isinstance(val, list):
        return [_from_json(v) for v in val]
    if isinstance(val, dict):
        return {k: _from_json(v) for k, v in val.items()}
    if val in ('inf', '-inf', 'nan'):
        return float(val)
    return val


def SetupFileName(xml_file):
    """ The setup file belonging to `xml_file`: `.xml` replaced by `.json`. """
    root, ext = os.path.splitext(str(xml_file))
    if ext.lower() == '.xml':
        return root + '.json'
    return str(xml_file) + '.json'


class Setup:
    """ Ports and NF2FF boxes of a simulation, to write them to a setup file
    and recreate them from it. """
    def __init__(self):
        self.ports = []   # [type, args, port]
        self.nf2ff = []   # [args, nf2ff]

    def Clear(self):
        self.ports = []
        self.nf2ff = []

    def AddPort(self, port):
        """ AddPort(port)

        Add `port` (any :class:`openEMS.ports.Port` subclass of openEMS),
        unless it is in already. A subclass of your own is stored as the
        openEMS class it is built on. Returns the port. """
        if any(p is port for _, _, p in self.ports):
            return port
        cls, args = _InitArgs(port)
        if cls.SETUP_TYPE not in PORT_TYPES:
            raise TypeError('{} is not a port'.format(cls.__name__))
        # the user's properties the port adds geometry to are not needed
        # afterwards, see _USER_PROPS
        args = {k: v for k, v in args.items() if k not in _USER_PROPS or (k == 'mat_prop' and v is None)}
        self.ports.append([cls.SETUP_TYPE, args, port])
        return port

    def AddNF2FF(self, box):
        """ AddNF2FF(box)

        Add the :class:`openEMS.nf2ff.nf2ff` box `box`, unless it is in
        already. Returns the box. """
        if any(b is box for _, b in self.nf2ff):
            return box
        cls, args = _InitArgs(box)
        if cls is not nf2ff.nf2ff:
            raise TypeError('{} is not an nf2ff box'.format(cls.__name__))
        self.nf2ff.append([dict(args), box])
        return box

    def GetPorts(self):
        return [port for _, _, port in self.ports]

    def GetPort(self, port_nr, prefix=''):
        for _, _, port in self.ports:
            if port.number == port_nr and port.prefix == prefix:
                return port
        raise KeyError('no port {} with prefix "{}"'.format(port_nr, prefix))

    def GetNF2FFBoxes(self):
        return [box for _, box in self.nf2ff]

    def ToDict(self):
        out = {'format': FORMAT, 'version': VERSION, 'ports': [], 'nf2ff': []}
        for port_type, args, port in self.ports:
            where = '{} port {}'.format(port_type, port.number)
            out['ports'].append({'type': port_type, 'args': _to_json(args, where)})
        for args, box in self.nf2ff:
            out['nf2ff'].append({'args': _to_json(args, 'nf2ff box "{}"'.format(box.name))})
        return out

    def Write(self, file):
        """ Write(file)

        Write the setup file `file` (JSON). """
        data = self.ToDict()   # first, so a failure leaves no half-written file
        with open(file, 'w') as f:
            json.dump(data, f, indent=1, allow_nan=False)

    def Read(self, file, CSX):
        """ Read(file, CSX)

        Recreate the ports and NF2FF boxes of `file` against `CSX`, the
        structure read from the matching XML, see
        :meth:`openEMS.openEMS.openEMS.ReadSetup`. Replaces what this object
        held, unless any of them does not match `CSX`. """
        with open(file) as f:
            data = json.load(f)
        if data.get('format') != FORMAT:
            raise ValueError('{}: not an openEMS setup file'.format(file))
        if data.get('version', 0) > VERSION:
            raise ValueError('{}: setup file version {} is newer than supported ({})'.format(
                file, data.get('version'), VERSION))

        loaded = _PropertyXML(CSX)
        port_list = []
        for entry in data.get('ports', []):
            port_type = entry['type']
            if port_type not in PORT_TYPES:
                raise ValueError('{}: unknown port type "{}"'.format(file, port_type))
            args = _from_json(entry['args'])
            port = _Rebuild(CSX, loaded, '{} port {}'.format(port_type, args.get('port_nr')),
                            PORT_TYPES[port_type], args)
            port_list.append([port_type, args, port])
        nf2ff_list = []
        for entry in data.get('nf2ff', []):
            args = _from_json(entry['args'])
            box = _Rebuild(CSX, loaded, 'nf2ff box "{}"'.format(args.get('name')), nf2ff.nf2ff, args)
            nf2ff_list.append([args, box])
        self.ports = port_list
        self.nf2ff = nf2ff_list


def _InitArgs(obj):
    init_args = getattr(obj, '_init_args', None)
    if init_args is None:
        raise TypeError('{} does not record how it was created'.format(type(obj).__name__))
    return init_args


def _Rebuild(CSX, loaded, what, cls, args):
    """ Create `cls(scratch, **args)` on a scratch structure with the mesh of
    `CSX`, check its properties against `CSX` and rebind the result to
    `CSX`. """
    grid = CSX.GetGrid()
    scratch = ContinuousStructure()
    s_grid = scratch.GetGrid()
    s_grid.SetMeshType(grid.GetMeshType())
    s_grid.SetDeltaUnit(grid.GetDeltaUnit())
    for n in range(3):
        s_grid.SetLines(n, grid.GetLines(n))

    kw = dict(args)
    for name, add in _USER_PROPS.items():   # placeholders, see _USER_PROPS
        if name in inspect.signature(cls.__init__).parameters and not (name == 'mat_prop' and kw.get(name, 0) is None):
            kw[name] = getattr(scratch, add)(_PLACEHOLDER.format(name))
    obj = cls(scratch, **kw)

    created = _PropertyXML(scratch)
    mapping = {}
    for key, elem in created.items():
        if key[1].startswith(_PLACEHOLDER.format('')):
            mapping[key] = None
            continue
        found = loaded.get(key)
        if found is None:
            raise ValueError('{}: property "{}" ({}) not found in the structure'.format(what, key[1], key[0]))
        if found is _AMBIGUOUS:
            raise ValueError('{}: property "{}" ({}) exists more than once in the structure'.format(
                what, key[1], key[0]))
        diff = _Diff(elem, found, key[0])
        if diff:
            raise ValueError('{}: property "{}" differs from the structure at {}'.format(what, key[1], diff))
        mapping[key] = loaded.props[key]

    def rebind(val):
        if isinstance(val, CSProperties):
            return mapping[(created.tags[id(val)], val.GetName())]
        if isinstance(val, list):
            # placeholders of user properties are dropped from lists
            return [rebind(v) for v in val if not (isinstance(v, CSProperties) and rebind(v) is None)]
        return val
    for k, v in vars(obj).items():
        if v is scratch:
            setattr(obj, k, CSX)
        else:
            setattr(obj, k, rebind(v))
    obj._init_args = (cls, dict(args))   # not the placeholders
    return obj


_AMBIGUOUS = object()


class _PropertyXML(dict):
    """ (XML tag, name) --> XML element of every property in `CSX`. """
    def __init__(self, CSX):
        super().__init__()
        fd, fn = tempfile.mkstemp(suffix='.xml')
        os.close(fd)
        try:
            CSX.Write2XML(fn)
            root = ET.parse(fn).getroot()
        finally:
            os.remove(fn)
        elems = root.find('Properties')
        if elems is None:
            elems = []
        # XML order is the order of GetAllProperties(); CSXCAD keeps one
        # wrapper per property, so the wrapper identifies it
        self.props = {}
        self.tags = {}
        for elem, prop in zip(elems, CSX.GetAllProperties()):
            key = (elem.tag, elem.get('Name'))
            self.tags[id(prop)] = elem.tag
            if key in self:
                self[key] = _AMBIGUOUS
                continue
            self[key] = elem
            self.props[key] = prop


def _Diff(a, b, path):
    """ First difference between the XML elements `a` and `b`, ignoring look
    and bookkeeping, as a readable path, or None. """
    if a.tag != b.tag:
        return '{}: {} vs. {}'.format(path, a.tag, b.tag)
    attr_a = {k: v for k, v in a.attrib.items() if k not in _IGNORE_ATTRS}
    attr_b = {k: v for k, v in b.attrib.items() if k not in _IGNORE_ATTRS}
    for k in sorted(set(attr_a) | set(attr_b)):
        if attr_a.get(k) != attr_b.get(k):
            return '{}/@{}: {} vs. {}'.format(path, k, attr_a.get(k), attr_b.get(k))
    if (a.text or '').strip() != (b.text or '').strip():
        return '{}: text "{}" vs. "{}"'.format(path, (a.text or '').strip(), (b.text or '').strip())
    kids_a = [c for c in a if c.tag not in _IGNORE_ELEMS]
    kids_b = [c for c in b if c.tag not in _IGNORE_ELEMS]
    for n, (ca, cb) in enumerate(zip(kids_a, kids_b)):
        d = _Diff(ca, cb, '{}/{}[{}]'.format(path, ca.tag, n))
        if d:
            return d
    if len(kids_a) != len(kids_b):
        return '{}: {} vs. {} child elements'.format(path, len(kids_a), len(kids_b))
    return None
