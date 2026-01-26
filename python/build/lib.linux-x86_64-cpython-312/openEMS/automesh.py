# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 20:29:25 2017

@author: thorsten
"""

import sys
import numpy as np

from CSXCAD import CSPrimitives
from CSXCAD.Utilities import CheckNyDir, GetMultiDirs
from openEMS.physical_constants import C0

def mesh_combine(mesh1, mesh2, sort=True):
    mesh = [None, None, None]
    for ny in range(3):
        if mesh1[ny] is None and mesh2[ny] is None:
            continue
        elif mesh1[ny] is None:
            mesh[ny] = mesh2[ny]
        elif mesh2[ny] is None:
            mesh[ny] = mesh1[ny]
        else:
            mesh[ny] = list(sorted(mesh1[ny] + mesh2[ny]))
    return mesh

def mesh_hint_from_primitive(primitive, dirs, **kw):
    if primitive.GetType() == CSPrimitives.PrimitiveType.POINT:
        return mesh_hint_from_point(primitive, dirs, **kw)
    if primitive.GetType() == CSPrimitives.PrimitiveType.BOX:
        return mesh_hint_from_box(primitive, dirs, **kw)
    else:
        return None

def mesh_hint_from_point(point, dirs, **kw):
    """ mesh_hint_from_point(point, dirs)

    Get a grid hint for the coordinates of the point.

    :param dirs: str -- 'x','y','z' or 'xy', 'yz' or 'xyz' or 'all'
    :param mesh: combine mesh hint to existing mesh
    :returns: (3,) list of mesh hints
    """
    hint = [None, None, None]
    coord = point.GetCoord()
    for ny in GetMultiDirs(dirs):
        hint[ny] = [coord[ny],]
    if 'mesh' in kw:
        return mesh_combine(hint, kw['mesh'])
    return hint

def mesh_hint_from_box(box, dirs, **kw):
    """ mesh_hint_from_box(box, dirs, metal_edge_res=None, **kw)

    Get a grid hint for the edges of the given box with an an optional 2D metal
    edge resolution.

    :param dirs: str -- 'x','y','z' or 'xy', 'yz' or 'xyz' or 'all'
    :param metal_edge_res: float -- 2D flat edge resolution
    :param up_dir: bool -- Enable upper edge
    :param down_dir: bool -- Enable lower edge
    :param mesh: combine mesh hint to existing mesh
    :returns: (3,) list of mesh hints
    """
    metal_edge_res = kw.get('metal_edge_res', None)
    up_dir   = kw.get('up_dir'  , True)
    down_dir = kw.get('down_dir', True)

    if metal_edge_res is None:
        mer = 0
    else:
        mer = np.array([-1.0, 2.0])/3 * metal_edge_res
    if box.HasTransform():
        sys.stderr.write('FDTD::automesh: Warning, cannot add edges to grid with transformations enabled\n')
        return
    hint = [None, None, None]
    start = np.fmin(box.GetStart(), box.GetStop())
    stop  = np.fmax(box.GetStart(), box.GetStop())
    for ny in GetMultiDirs(dirs):
        hint[ny] = []
        if metal_edge_res is not None and stop[ny]-start[ny]>metal_edge_res:
            if down_dir:
                hint[ny].append(start[ny]-mer[0])
                hint[ny].append(start[ny]-mer[1])
            if up_dir:
                hint[ny].append(stop[ny]+mer[0])
                hint[ny].append(stop[ny]+mer[1])
        elif stop[ny]-start[ny]:
            if down_dir:
                hint[ny].append(start[ny])
            if up_dir:
                hint[ny].append(stop[ny])
        else:
            hint[ny].append(start[ny])
    if 'mesh' in kw:
        return mesh_combine(hint, kw['mesh'])
    return hint

def mesh_estimate_cfl_timestep(mesh):
    """ mesh_estimate_cfl_timestep(mesh)

    Estimate the maximum CFL time step of the given mesh needed to ensure numerical stability,
    assuming propagation in pure vacuum.

    :returns: the maximum CFL time step, in seconds
    """
    invMinDiff = [None, None, None]
    for ny in range(3):
        invMinDiff[ny] = np.min(np.diff(mesh.GetLines(ny))) ** -2

    delta_t = mesh.GetDeltaUnit() / (C0 * np.sqrt(np.sum(invMinDiff)))

    return delta_t

