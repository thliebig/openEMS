# -*- coding: utf-8 -*-
#
# Copyright (C) 2025 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

cimport openEMS.sar_calculation
import os

cdef class SAR_Calculation:
    """SAR averaging and calculation.

    Keyword arguments accepted by the constructor (all optional):

    mass : float
        Averaging mass in grams (0 = local SAR, default 0).
    method : str
        Spatial averaging method: 'SIMPLE' (default), 'IEEE_C95_3', or
        'IEEE_62704'. All three average over a cubical mass; see
        :meth:`SetAveragingMethod`. Has no effect when `mass` is 0.
    verbose : int
        Debug verbosity level.
    autoRange : float
        Restrict the calculation to the cells within this many dB of the
        peak local SAR.
    EnableCubeStats : bool
        Record per-cube averaging statistics in the output file.
    """
    def __cinit__(self, **kw):
        self.thisptr = new _SAR_Calculation()
        for k, v in kw.items():
            if k == 'mass':
                self.SetAveragingMass(v)
            elif k=='method':
                self.SetAveragingMethod(v)
            elif k=='verbose':
                self.SetDebugLevel(int(v))
            elif k=='debug':
                self.SetDebugLevel(int(v))
            elif k=='autoRange':
                self.EnableAutoRange(float(v))
            elif k=='EnableCubeStats':
                if v:
                    self.EnableCubeStats()
            else:
                raise Exception('Unknown keyword argument: "{}"'.format(k))

    def __dealloc__(self):
        del self.thisptr

    def SetDebugLevel(self, level):
        """Set the verbosity of the calculation.

        Parameters
        ----------
        level : int
            0 is silent, higher values print progressively more detail about
            the averaging (e.g. the auto range and per frequency summaries).
        """
        self.thisptr.SetDebugLevel(level)

    def EnableProgress(self, enable):
        """Enable or disable the progress indicator during averaging.

        Parameters
        ----------
        enable : bool
            Averaging a large mesh can take a while; this prints how far the
            calculation has come.
        """
        self.thisptr.EnableProgress(enable)

    def SetAveragingMass(self, mass):
        """Set averaging mass in grams (0 = local SAR)."""
        self.thisptr.SetAveragingMass(float(mass)/1000)

    def SetAveragingMethod(self, method, silent=True):
        """Set the spatial averaging method.

        All three methods perform a real cubical mass averaging: a cube is
        grown around each tissue cell until it encloses the averaging mass,
        and the SAR is the absorbed power in that cube divided by its mass.
        They differ only in how strictly the cube has to qualify as valid,
        and therefore in how many cells end up being filled in from a
        neighbouring cell's cube instead:

        * 'SIMPLE' (default) accepts every cube that reaches the target
          mass, even one clipped by the edge of the dump box, so a cube is
          built for essentially every tissue cell. This avoids the surface
          artifacts that the fill-in step can introduce, but does not
          follow the validity rules of the standards.
        * 'IEEE_C95_3' uses the same 5% mass tolerance, but additionally
          requires the cube to fit inside the dump box.
        * 'IEEE_62704' is the strictest: the mass has to match to within
          1e-6, and no more than 10% of the cube volume may be background
          (air). The latter is what rejects cubes sitting on a tissue
          surface.

        Note that this is unrelated to local SAR: an averaging mass of 0
        skips the averaging entirely and the method has no effect.

        Parameters
        ----------
        method : str
            'SIMPLE', 'IEEE_C95_3' or 'IEEE_62704'.
        silent : bool, optional
            Suppress the confirmation message printed by the C++ code.

        Returns
        -------
        bool
            True on success, False if the method name is unknown.

        Notes
        -----
        None of the methods is validated according to IEC/IEEE-62704-1.

        See Also
        --------
        SetAveragingMass
        """
        return self.thisptr.SetAveragingMethod(method.encode('UTF-8'), silent)

    def EnableAutoRange(self, dBmax):
        """Restrict the calculation to the region around the peak.

        Only the cells whose local SAR is within `dBmax` dB of the peak local
        SAR are averaged, which speeds up large meshes with a localised hot
        spot. The result is written on the reduced mesh, so the output covers a
        smaller region than the input.

        The averaged SAR of a cube is the mass weighted mean of the local SAR
        of its cells and can never exceed the largest local SAR inside that
        cube. Everything that is dropped here is therefore below the threshold
        after averaging as well, and the region is padded by roughly one
        averaging cube so that cubes centred just outside it are covered too.
        That padding is an estimate, so this remains a speedup and not a
        guarantee to find the global peak. A warning is printed if the peak
        that was found is itself below the threshold. Do not use the auto range
        for standard compliance work.

        Parameters
        ----------
        dBmax : float
            Range below the peak local SAR, in dB. Values <= 0 disable the
            auto range.
        """
        self.thisptr.EnableAutoRange(float(dBmax))

    def EnableCubeStats(self):
        """Record per-cell averaging cube statistics in the output file.

        This exposes how the averaging cube was found for every cell and is
        the main handle for validating the averaging against the conformance
        requirements of IEC/IEEE 62704-1 -- it shows which cells got a proper
        centred cube, which fell back to a neighbouring one, and what mass and
        volume each cube actually enclosed.

        Three datasets are added next to the SAR result:

        ``f{n}_CubeType`` (unsigned byte)
            How the cube for this cell was obtained:

            * 0 -- no averaging cube; background/air, or no cube could be built
            * 1 to 6 -- second pass: the cube was built with one face pinned to
              the cell, 1/2 = lower/upper x face, 3/4 = y, 5/6 = z
            * 7 -- first pass: a valid cube enclosing the target mass was found
              centred on the cell itself
            * 8 -- the cell was covered by another cell's cube but never got a
              valid cube of its own

        ``f{n}_CubeMass``
            Mass actually enclosed by the cube, in kg.  Compare against the
            requested averaging mass to check the convergence tolerance.

        ``f{n}_CubeVol``
            Volume of the cube, in m^3.

        .. note::
            Cube statistics can only be recorded for a **single** frequency.
            With more than one frequency of interest the calculation prints a
            warning and writes no statistics.

        See :ref:`concept_sar` for the averaging methods and the output format.
        """
        self.thisptr.EnableCubeStats()

    def CalcFromHDF5(self, h5_fn, out_name, export_cube_stats=False, numThreads=0):
        """Read raw field data from h5_fn, run the SAR calculation, and write
        results to out_name. Returns True on success.
        numThreads=0 uses all available hardware threads."""
        if not os.path.exists(h5_fn):
            raise Exception('File "{}" does not exist'.format(h5_fn))
        cdef string in_fn = h5_fn.encode('UTF-8')
        cdef string out_fn = out_name.encode('UTF-8')
        cdef unsigned int c_numThreads = numThreads
        if export_cube_stats:
            self.EnableCubeStats()
        with nogil:
            ok = self.thisptr.CalcFromHDF5(in_fn, out_fn, False, c_numThreads)
        return ok

