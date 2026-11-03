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
#include <vector>

#include "engine.h"

//! A run of cells along the fast axis with one update class and one coefficient set
struct Engine_CPU_Run
{
	enum Class : uint8_t {ZERO, LOSSLESS, LOSSY, UPML};

	uint32_t start;   //!< first line along the fast axis
	uint32_t stop;    //!< one past the last
	float vv;         //!< vv (ii) of the run, 1 for LOSSLESS
	float a;          //!< vi (iv) of a cell is a * g(fast line), see Engine_CPU
	uint64_t u;       //!< UPML: the index of the first cell in Engine_CPU_UPML
	Class cls;
};

//! The UPML cells of one field component, in the order of their runs: the main coefficients, the UPML ones and the flux
struct Engine_CPU_UPML
{
	std::vector<float> vv, vi;          //!< the main update, applied to the flux (Operator_Ext_UPML rewrites them in the PML)
	std::vector<float> uvv, ufo, ufn;   //!< Operator_Ext_UPML::vv, vvfo, vvfn (ii, iifo, iifn)
	std::vector<float> flux;
	void Clear() {vv.clear(); vi.clear(); uvv.clear(); ufo.clear(); ufn.clear(); flux.clear();}
};

//! The run lists of one field component, per line of the fast axis
struct Engine_CPU_Runs
{
	std::vector<Engine_CPU_Run> runs;
	std::vector<uint64_t> first;   //!< per line the index of its first run, plus one entry at the end
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

  The cells of a UPML (Operator_Ext_UPML) are the class UPML, with exact per
  cell coefficients, and the engine updates them itself in one pass instead of
  the extension's pre and post passes around the main update:
    Fn = vv*F + vi*curl,  V = uvv*V - ufo*F + ufn*Fn,  F = Fn
  with F the flux. This is the extension's sequence, which leaves the real
  field in V outside the main update: its pre hook is the last before it, its
  post hook the first after it (its priority), so the other extensions see the
  same values. The UPML engine extensions are dropped.

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
		Engine_CPU_Run::Class cls;
		float vv;
		int64_t q;          //!< m quantised on a log scale, see GetCellKey(); its sign in sign
		int sign;
		bool operator==(const CellKey& o) const {return cls==o.cls && vv==o.vv && q==o.q && sign==o.sign;}
		double M() const;   //!< the m this key stands for
	};
	CellKey GetCellKey(int f, int n, const unsigned int pos[3]) const;

	//! the UPML regions: take over their coefficients, drop their engine extensions
	void TakeOverUPML();
	//! the UPML region of a cell and the cell's position in it, or NULL
	const class Operator_Ext_UPML* FindUPML(const unsigned int pos[3], unsigned int loc[3]) const;
	std::vector<const class Operator_Ext_UPML*> m_UPML;
	Engine_CPU_UPML m_UPMLCells[2][3];

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
