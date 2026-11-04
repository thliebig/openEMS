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

#ifndef ENGINE_CPU_H
#define ENGINE_CPU_H

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

#include "engine.h"
#include "engine_cpu_kernel.h"

//! A set of mesh nodes by linear index (x*Ny + y)*Nz + z: the rank of a node among them in O(1)
struct Engine_CPU_NodeSet
{
	std::vector<uint64_t> bits;
	std::vector<uint32_t> before;   //!< per word the nodes in the words before it
	//! the rank of \a lin, or -1 if it is not in the set
	int64_t Find(size_t lin) const;
};

//! The run lists of one field component, per line of the fast axis
struct Engine_CPU_Runs
{
	std::vector<Engine_CPU_Run> runs;
	std::vector<uint64_t> first;   //!< per line the index of its first run, plus one entry at the end
	std::vector<Engine_CPU_Set> sets;   //!< the distinct coefficient sets of the runs
};

//! CPU engine on plain linear field arrays, the coefficients as runs of cell classes
/*!
  The fields are two ArrayNIJK in a layout chosen per model: N slowest, then
  the two slow axes, the fast axis last (stride 1). The fast axis is the one
  along which the coefficients form the longest runs; the threads split the
  slowest axis (GetSplitAxis()).

  Per field component and line of the fast axis, the cells are runs of one
  class: ZERO (both coefficients 0: PEC, the outer lines; the field is set to
  0 as the update would), LOSSLESS (vv == 1: V += vi*curl) or LOSSY. The
  coefficients are factorised: vi = m * g, with g a product of per-axis mesh
  factors (edge length over area) and m the material part, 1/eps up to
  constants. A run is a maximal set of neighbours with the same class, vv
  and m (equal within a relative tolerance), and stores vv and a = m times
  the factors of its two slow lines; per cell the kernel multiplies with
  the fast factor. So graded meshes do not break runs, and the coefficient
  traffic is the fast factor array (in cache) plus a run per run.

  The cells of a UPML (Operator_Ext_UPML) are the class UPML, runs with equal
  main and UPML coefficients and a flux per cell, and the engine updates them
  itself in one pass instead of
  the extension's pre and post passes around the main update:
    Fn = vv*F + vi*curl,  V = uvv*V - ufo*F + ufn*Fn,  F = Fn
  with F the flux. This is the extension's sequence, which leaves the real
  field in V outside the main update: its pre hook is the last before it, its
  post hook the first after it (its priority), so the other extensions see the
  same values. The UPML engine extensions are dropped, the coefficients of
  the UPML operator extensions freed, and the operator frees its full
  coefficients after creating the engine (Operator::CreateEngine()).

  The cells of a dispersive material, Drude/Lorentz (Operator_Ext_LorentzMaterial,
  also the conducting sheets) or Debye (Operator_Ext_DebyeMaterial), are the
  class DISP, with per run the main coefficients and the poles of their ADEs,
  and per cell and pole their state. The kernel runs the extensions' pre
  passes, the main update, the Debye post pass and the Drude/Lorentz apply
  pass (X -= J) in one, so the latter two before the other extensions' post
  and apply hooks instead of in the order of their priority: the same for
  additive hooks up to rounding, but not for the cylinder's at r=0. With a
  cylindrical mesh, a dispersive cell in a UPML, more than MAX_POLES poles of
  a kind or a second Debye extension, the extension stays.

  A lossless current inside metal, whose four voltages of the curl are ZERO
  cells, never changes: it is in no run, so the kernels skip it. Only if no
  extension writes such a voltage: none on the outer lines (Mur), in a UPML
  or excited, and only the extensions known not to (FindInert()).

  In a UPML normal to the fast axis the coefficients change from cell to
  cell along it: neighbouring UPML runs shorter than 4 cells on average are
  one run of the class UPML_PROFILE, with the coefficients per cell from a
  profile along the fast axis, which the lines through the same layers
  share.

  The extensions see the fields through Engine::GetVolt() and friends, which
  follow the layout, and run in the dispatch of the basic engine.
  */
class Engine_CPU : public Engine
{
public:
	static Engine_CPU* New(const Operator* op, unsigned int numThreads = 0);
	virtual ~Engine_CPU();

	virtual void Init();
	virtual void Reset();

	virtual int GetSplitAxis() const {return m_SplitAxis;}

	//! Trapezoid tiles along the split axis sized for the L3, if every extension supports boxes and there is no sub-grid
	/*!
	  OPENEMS_CPU_TBLOCK=0 keeps the flat sweep, =k:W sets the timesteps per block
	  and the tile width in lines.
	  */
	virtual bool BuildBlockedPlans(unsigned int numThreads, std::vector<Engine_BlockedPlan>& plans);

	//! The fast axis (stride 1 in memory)
	int GetFastAxis() const {return m_FastAxis;}

protected:
	Engine_CPU(const Operator* op);

	using Engine::UpdateVoltages;
	using Engine::UpdateCurrents;
	virtual void UpdateVoltages(const Engine_Box& box);
	virtual void UpdateCurrents(const Engine_Box& box);
	//! the range updates of the old thread split: the box of x lines [startX, startX+numX)
	virtual void UpdateVoltages(unsigned int startX, unsigned int numX);
	virtual void UpdateCurrents(unsigned int startX, unsigned int numX);

	virtual std::string GetFieldLayout() const {return m_Layout;}

	//! the fast axis from the run statistics of all components (OPENEMS_CPU_FAST_AXIS=x|y|z overrides)
	void ChooseLayout();
	//! the per-axis geometry factors g of field \a volt, component n
	void CalcGeometry();
	void BuildRuns();

	//! The class and coefficients of a cell: neighbours with equal keys form a run
	struct CellKey
	{
		Engine_CPU_Set::Class cls;
		float vv;
		int64_t q;          //!< m quantised on a log scale, see GetCellKey(); its sign in sign
		int sign;
		float uvv, ufo, ufn;   //!< UPML only, else 0
		uint8_t nl, nd;        //!< DISP: the Drude/Lorentz poles, then the Debye poles (Engine_CPU_Set); the others 0
		float solve;           //!< DISP: Engine_CPU_Set::solve, 1 without Debye poles
		Engine_CPU_Pole poles[8];   //!< k_ext, k_lor and solve rounded to 20 bits of mantissa
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
		uint32_t Hash() const;   //!< for counting runs only (ChooseLayout())
	};
	CellKey GetCellKey(int f, int n, const unsigned int pos[3]) const;

	//! the UPML regions: take over their coefficients, drop their engine extensions
	void TakeOverUPML();
	//! the UPML region of a cell and the cell's position in it, or NULL
	const class Operator_Ext_UPML* FindUPML(const unsigned int pos[3], unsigned int loc[3]) const;
	std::vector<const class Operator_Ext_UPML*> m_UPML;

	static const int MAX_POLES = 4;   //!< per kind, Drude/Lorentz and Debye
	//! the dispersive extensions the engine updates itself, with their cells
	void FindDispersive();
	//! drop the engine extensions of m_Disp
	void TakeOverDispersive();
	//! free the cell lists and coefficients of the operator extensions of m_Disp, after BuildRuns()
	void FreeDispersive();
	struct Dispersive
	{
		const class Operator_Ext_Dispersive* op;
		const class Operator_Ext_LorentzMaterial* lorentz;   //!< or NULL
		const class Operator_Ext_DebyeMaterial* debye;       //!< or NULL
		std::vector<Engine_CPU_NodeSet> cells;   //!< per order (Debye: one list for all poles)
	};
	std::vector<Dispersive> m_Disp;
	//! whether the inert currents may be skipped: only known extensions, and the excited voltages
	void FindInert();
	bool m_SkipInert;
	std::unordered_set<size_t> m_ExcitedVolt;   //!< 3*lin + n
	//! a voltage that is always 0: ZERO coefficients, and no extension writes it
	bool IsZeroVolt(int n, const unsigned int pos[3]) const;
	//! a lossless current whose curl voltages are all IsZeroVolt()
	bool IsInertCurr(int n, const unsigned int pos[3]) const;

	//! the poles of a cell into \a key
	void GetCellPoles(int f, int n, const unsigned int pos[3], CellKey& key) const;

	//! the kernel data and the kernels of the instruction set of this CPU, after BuildRuns()
	void InitKernels();
	Engine_CPU_Kernel m_Kernel;
	Engine_CPU_Update m_UpdateVoltages;
	Engine_CPU_Update m_UpdateCurrents;
	std::string m_ISA;
	//! [field][component] the state of the UPML and DISP runs (Engine_CPU_Run::u)
	std::vector<float> m_State[2][3];
	//! [field][component] the poles of the DISP sets
	std::vector<Engine_CPU_Pole> m_Poles[2][3];
	//! [field][component] the profiles of the UPML_PROFILE sets: m, vv, uvv, ufo, ufn per cell
	std::vector<float> m_Profile[2][3][5];

	int m_FastAxis;
	int m_SplitAxis;     //!< the slowest axis
	int m_MidAxis;       //!< the other slow axis
	std::string m_Layout;

	//! [field: 0 volt, 1 curr][component][axis] the factor of each line, g = product over the axes
	std::vector<double> m_G[2][3][3];
	//! [field][component] the fast factor as float, for the kernels
	std::vector<float> m_Gf[2][3];
	Engine_CPU_Runs m_Runs[2][3];
};

#endif // ENGINE_CPU_H
