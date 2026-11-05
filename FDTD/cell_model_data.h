/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CELL_MODEL_DATA_H
#define CELL_MODEL_DATA_H

// The plain data of the cell model (CellModel) that the update kernels of the
// engines read. Plain data only: the CPU kernels are compiled once per
// instruction set and chosen at run time, and an inline function here,
// compiled with AVX2 in one file and without in another, could be merged by
// the linker into the AVX2 version for all callers.

#include <cstddef>
#include <cstdint>

//! The update class and the coefficients a run of cells shares
struct CellModel_Set
{
	enum Class : uint8_t {ZERO, LOSSLESS, LOSSY, UPML, DISP, UPML_PROFILE};

	double m;              //!< the material part of vi (iv): vi = m * g, see CellModel
	float vv;              //!< vv (ii), 1 for LOSSLESS
	float uvv, ufo, ufn;   //!< UPML: Operator_Ext_UPML::vv, vvfo, vvfn (ii, iifo, iifn); vv and m are the main update on the flux
	uint32_t offset;       //!< DISP: its first pole in CellModel_Component::poles; UPML_PROFILE: its first cell in profile
	uint8_t nl, nd;        //!< DISP: the number of Drude/Lorentz poles, then of Debye poles
	float solve;           //!< DISP: the implicit solve of the Debye poles (Operator_Ext_DebyeMaterial::v_solve_ADE), 1 without
	Class cls;
};

//! A pole of a dispersive material
/*!
  Drude/Lorentz (Operator_Ext_LorentzMaterial): v_int_ADE, v_ext_ADE and
  v_Lor_ADE (i_*), k_lor 0 for a Drude pole. Debye
  (Operator_Ext_DebyeMaterial): v_relax_ADE in k_int, v_drive_ADE in k_ext.
  */
struct CellModel_Pole
{
	float k_int, k_ext;
	float k_lor;
};

//! A run of cells along the fast axis with one coefficient set
struct CellModel_Run
{
	uint32_t start;   //!< first line along the fast axis
	uint32_t stop;    //!< one past the last
	uint32_t set;     //!< index in CellModel_Component::sets
	uint32_t u;       //!< UPML, UPML_PROFILE, DISP: the index of the run's state (CellModel_Component::state)
};

#endif // CELL_MODEL_DATA_H
