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

#ifndef CELL_MODEL_H
#define CELL_MODEL_H

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

#include "cell_model_data.h"

class Operator;
class Operator_Extension;
class Operator_Ext_UPML;
class Operator_Ext_Dispersive;
class Operator_Ext_LorentzMaterial;
class Operator_Ext_DebyeMaterial;

//! A set of mesh nodes by linear index (x*Ny + y)*Nz + z: the rank of a node among them in O(1)
/*!
  Init(), Add() each node, Finish(), then Find(). Added in ascending order, the rank is
  the node's index in that order (an extension's cell list, tests/test_cell_model.cpp).
  */
struct CellModel_NodeSet
{
	std::vector<uint64_t> bits;
	std::vector<uint32_t> before;   //!< per word the nodes in the words before it

	//! an empty set of \a nodes nodes
	void Init(size_t nodes);
	//! add \a lin, false if it is in the set already
	bool Add(size_t lin);
	//! the ranks, after the last Add()
	void Finish();
	//! the rank of \a lin, or -1 if it is not in the set
	int64_t Find(size_t lin) const;
};

//! The runs of one field component, per line of the fast axis
struct CellModel_Component
{
	std::vector<CellModel_Run> runs;
	std::vector<uint64_t> first;         //!< per line (split axis, then mid axis) the index of its first run, plus one entry at the end
	std::vector<CellModel_Set> sets;     //!< the distinct coefficient sets of the runs
	std::vector<CellModel_Pole> poles;   //!< the poles of the DISP sets
	std::vector<float> profile[5];       //!< the profiles of the UPML_PROFILE sets: m, vv, uvv, ufo, ufn per cell
	size_t state = 0;                    //!< the floats of state of the UPML and DISP runs (CellModel_Run::u)
};

//! The cells of an operator as runs of cell classes, for the update kernels of the engines
/*!
  Built by an engine from its operator (Analyse(), RecommendFastAxis(),
  Build(), FreeFused()); the engine keeps the runs in the form its kernels
  want and frees the model's analysis. Engine-neutral: the CPU engine
  (Engine_CPU) updates the runs directly, a GPU engine would expand them to an
  index per cell (notes: cell-model-SPEC).

  Per field component and line of the fast axis, the cells are runs of one
  class (CellModel_Set::Class): ZERO (both coefficients 0: PEC, the outer
  lines; the field is set to 0 as the update would), LOSSLESS (vv == 1:
  V += vi*curl), LOSSY, and the classes of the extensions the model takes
  over (UPML, UPML_PROFILE, DISP, below). The coefficients are factorised:
  vi = m * g, with g a product of per-axis mesh factors (edge length over
  area, G()) and m the material part, 1/eps up to constants. A run is a
  maximal set of neighbours with the same class, vv and m (equal within a
  relative 2^-20); so graded meshes do not break runs.

  The cells of a UPML (Operator_Ext_UPML) are the class UPML, runs with equal
  main and UPML coefficients and a flux per cell, which the kernel updates in
  one pass instead of the extension's pre and post passes around the main
  update:
    Fn = vv*F + vi*curl,  V = uvv*V - ufo*F + ufn*Fn,  F = Fn
  with F the flux. This is the extension's sequence, which leaves the real
  field in V outside the main update: its pre hook is the last before it, its
  post hook the first after it (its priority), so the other extensions see the
  same values. In a UPML normal to the fast axis the coefficients change from
  cell to cell along it: neighbouring UPML runs shorter than 4 cells on
  average are one run of the class UPML_PROFILE, with the coefficients per
  cell from a profile along the fast axis, which the lines through the same
  layers share.

  The cells of a dispersive material, Drude/Lorentz (Operator_Ext_LorentzMaterial,
  also the conducting sheets) or Debye (Operator_Ext_DebyeMaterial), are the
  class DISP, with per run the main coefficients and the poles of their ADEs,
  and per cell and pole their state. The kernel runs the extensions' pre
  passes, the main update, the Debye post pass and the Drude/Lorentz apply
  pass (X -= J) in one, so the latter two before the other extensions' post
  and apply hooks instead of in the order of their priority. The same up to
  rounding for additive hooks (excitation, TF/SF); the same for Mur, which
  overwrites its plane after the apply pass in the other engines too (equal
  priority runs in reverse order of creation); not the same for the
  cylinder's r=0 (so no fusion on a cylindrical mesh) and for a lumped RLC
  element with an engine extension (an inductor, or in series) on a
  Drude/Lorentz edge, which here sees the corrected field, in the other
  engines the uncorrected one: a modelling conflict either way, warned about
  (a lumped element is no material, so priority does not separate them).
  With a cylindrical mesh, a dispersive cell in a UPML, more than MAX_POLES
  poles of a kind or a second Debye extension, the extension stays.

  A lossless current inside metal, whose four voltages of the curl are ZERO
  cells, never changes: it is in no run, so the kernels skip it. Only if no
  extension writes such a voltage: none on the outer lines (Mur), in a UPML
  or excited, and only the extensions known not to; the engine decides
  whether to skip them (Analyse()).

  The engine updates the fused extensions' cells itself and drops their
  engine extensions (IsFused()); FreeFused() frees their operator data and
  the operator's own coefficients are freed by the operator after creating
  the engine (Operator::CreateEngine()), so an operator serves one engine.
  */
class CellModel
{
public:
	static const int MAX_POLES = 4;   //!< per kind, Drude/Lorentz and Debye

	explicit CellModel(const Operator* op);

	//! The extensions the model takes over (UPML, dispersive) and, if \a skipInert, the inert currents; before the engine's extensions exist
	void Analyse(bool skipInert = true);
	//! Whether the model updates the cells of the operator extension \a ext: the engine drops its engine extension
	bool IsFused(const Operator_Extension* ext) const;

	//! The axis along which the cells form the fewest runs: the fast axis the model is most compact on
	/*!
	  The run statistics only: an engine applies its own constraints on top
	  (Engine_CPU: the box plans need a closed cylinder's alpha or the local
	  absorbers' normal as the fast axis).
	  */
	int RecommendFastAxis() const;
	//! The runs along \a fast, the lines ordered by \a split, then \a mid
	void Build(int fast, int split, int mid);
	//! Free the operator data of the fused extensions and the analysis, after the engine took the runs
	void FreeFused();

	//! Write the model to the HDF5 file \a filename (no extension), after Build()
	/*!
	  For a script or a reader to find what the model made of a structure:
	  runs broken where the material does not change, an unexpected class,
	  inert regions missing, an extension not fused. `--debug-operator` writes
	  it as cellmodel_dump.h5 with the CPU engine. The layout:

	  - /Mesh: the mesh lines, as in the field dumps
	  - /CellModel, attributes: the axes of the runs (fast, split, mid: 0..2 for
	    x, y, z), skip_inert, the class names, the fused extensions and why
	    others were not fused
	  - /CellModel/volt, /CellModel/curr: per edge (vector fields [3][x][y][z])
	    `class` (CellModel_Set::Class; 6 an inert current in no run; -1 not
	    updated, the last line of the currents), `set` (its coefficient set,
	    -1 none) and `run` (the length of its run); attribute `counts`: the
	    edges per class 0..6
	  - /CellModel/volt/x .. z (and curr): the model of the component as the
	    kernels read it: `sets` [sets][10] (columns in its attribute), `runs`
	    [runs][4] (start, stop along the fast axis, set, state offset),
	    `first` [lines+1] (the first run of line split*N_mid + mid), `poles`
	    [poles][3] (k_int, k_ext, k_lor), `profile` [5][cells] (m, vv, uvv,
	    ufo, ufn of the UPML_PROFILE sets)
	  */
	void Dump(const std::string& filename) const;

	const CellModel_Component& Component(int f, int n) const {return m_Comp[f][n];}
	//! The factor of each line along \a axis of field \a f (0 voltages, 1 currents), component \a n: g = product over the axes
	const std::vector<double>& G(int f, int n, int axis) const {return m_G[f][n][axis];}

protected:
	const Operator* Op;
	unsigned int numLines[3];
	int m_F = 2, m_S = 0, m_M = 1;   //!< the axes of the runs (Build())

	//! the per-axis geometry factors g of field \a volt, component n
	void CalcGeometry();

	//! The class and coefficients of a cell: neighbours with equal keys form a run
	struct CellKey
	{
		CellModel_Set::Class cls;
		float vv;
		int64_t q;          //!< m quantised on a log scale, see GetCellKey(); its sign in sign
		int sign;
		float uvv, ufo, ufn;   //!< UPML only, else 0
		uint8_t nl, nd;        //!< DISP: the Drude/Lorentz poles, then the Debye poles (CellModel_Set); the others 0
		float solve;           //!< DISP: CellModel_Set::solve, 1 without Debye poles
		CellModel_Pole poles[2*MAX_POLES];   //!< k_ext, k_lor and solve rounded to 20 bits of mantissa
		int NumPoles() const {return nl+nd;}
		bool operator==(const CellKey& o) const
		{
			if (!(cls==o.cls && vv==o.vv && q==o.q && sign==o.sign && uvv==o.uvv && ufo==o.ufo && ufn==o.ufn
			      && nl==o.nl && nd==o.nd && solve==o.solve))
				return false;
			for (int p=0; p<NumPoles(); ++p)
				if (!(poles[p].k_int==o.poles[p].k_int && poles[p].k_ext==o.poles[p].k_ext && poles[p].k_lor==o.poles[p].k_lor))
					return false;
			return true;
		}
		double M() const;   //!< the m this key stands for
		uint32_t Hash() const;   //!< for counting runs only (RecommendFastAxis())
	};
	CellKey GetCellKey(int f, int n, const unsigned int pos[3]) const;

	//! the UPML region of a cell and the cell's position in it, or NULL
	const Operator_Ext_UPML* FindUPML(const unsigned int pos[3], unsigned int loc[3]) const;
	std::vector<const Operator_Ext_UPML*> m_UPML;

	//! the dispersive extensions the model takes over, with their cells
	void FindDispersive();
	//! a warning about lumped RLC edges with Drude/Lorentz poles, where the order of the two differs from the other engines
	void WarnLumpedInDispersive() const;
	struct Dispersive
	{
		const Operator_Ext_Dispersive* op;
		const Operator_Ext_LorentzMaterial* lorentz;   //!< or NULL
		const Operator_Ext_DebyeMaterial* debye;       //!< or NULL
		std::vector<CellModel_NodeSet> cells;   //!< per order (Debye: one list for all poles)
	};
	std::vector<Dispersive> m_Disp;
	//! the operator extensions IsFused() answers for, kept after FreeFused()
	std::vector<const Operator_Extension*> m_Fused;
	//! the dispersive extensions not fused, and why (Dump())
	std::vector<std::string> m_NotFused;
	//! free the cell lists and coefficients of the operator extensions of m_Disp
	void FreeDispersive();
	//! the poles of a cell into \a key
	void GetCellPoles(int f, int n, const unsigned int pos[3], CellKey& key) const;

	//! the excited voltages, and whether the inert currents may be skipped: only with known extensions (m_SkipInert from Analyse())
	void FindInert();
	bool m_SkipInert;
	std::unordered_set<size_t> m_ExcitedVolt;   //!< 3*lin + n
	//! a voltage that is always 0: ZERO coefficients, and no extension writes it
	bool IsZeroVolt(int n, const unsigned int pos[3]) const;
	//! a lossless current whose curl voltages are all IsZeroVolt()
	bool IsInertCurr(int n, const unsigned int pos[3]) const;

	//! [field: 0 volt, 1 curr][component][axis] the factor of each line, g = product over the axes
	std::vector<double> m_G[2][3][3];
	CellModel_Component m_Comp[2][3];
};

#endif // CELL_MODEL_H
