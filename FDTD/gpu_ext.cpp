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

#include "gpu_ext.h"

#include <algorithm>
#include <cstring>
#include <numeric>
#include <set>

#include "engine_gpu.h"
#include "operator_cylindermultigrid.h"
#include "excitation.h"
#include "kernels/ext_ops.h"
#include "extensions/engine_ext_excitation.h"
#include "extensions/operator_ext_excitation.h"
#include "extensions/engine_ext_mur_abc.h"
#include "extensions/operator_ext_mur_abc.h"
#include "extensions/engine_ext_absorbing_bc.h"
#include "extensions/operator_ext_absorbing_bc.h"
#include "extensions/engine_ext_lumpedRLC.h"
#include "extensions/operator_ext_lumpedRLC.h"
#include "extensions/engine_ext_tfsf.h"
#include "extensions/operator_ext_tfsf.h"
#include "extensions/engine_ext_steadystate.h"
#include "extensions/operator_ext_steadystate.h"
#include "extensions/engine_ext_cylinder.h"
#include "extensions/operator_ext_cylinder.h"
#include "extensions/engine_ext_lorentzmaterial.h"
#include "extensions/operator_ext_lorentzmaterial.h"
#include "extensions/engine_ext_debyematerial.h"
#include "extensions/operator_ext_debyematerial.h"

#include <functional>
#include <map>
#include <stdexcept>

void GPU_PointList::Init(const std::vector<unsigned int>& x)
{
	m_Order.resize(x.size());
	std::iota(m_Order.begin(), m_Order.end(), 0);
	std::stable_sort(m_Order.begin(), m_Order.end(), [&](size_t a, size_t b) {return x[a]<x[b];});
	m_X.resize(x.size());
	for (size_t n=0; n<m_Order.size(); ++n)
		m_X[n] = x[m_Order[n]];
}

void GPU_PointList::Range(const Engine_Box& box, uint64_t& begin, uint64_t& end) const
{
	begin = std::lower_bound(m_X.begin(), m_X.end(), box.start[0]) - m_X.begin();
	end = std::lower_bound(m_X.begin(), m_X.end(), box.stop[0]) - m_X.begin();
}

namespace
{
//! the flat ArrayNIJK index of component \a n at \a pos
uint64_t FlatIndex(const unsigned int numLines[3], unsigned int n, unsigned int x, unsigned int y, unsigned int z)
{
	return (((uint64_t)n*numLines[0] + x)*numLines[1] + y)*numLines[2] + z;
}

//! the nodes (x*ny + y)*nz + z of flat ArrayNIJK indices
void AddNodes(std::vector<uint64_t>& nodes, const std::vector<uint64_t>& flat, const unsigned int numLines[3])
{
	const uint64_t N = (uint64_t)numLines[0]*numLines[1]*numLines[2];
	for (uint64_t f : flat)
		nodes.push_back(f%N);
}

template <typename T>
T* Upload(GPU_Backend* b, const std::vector<T>& v)
{
	return (T*)b->DeviceAlloc(std::max<size_t>(v.size(), 1)*sizeof(T), v.empty() ? NULL : v.data());
}
}

/***************************************************************************************************/

//! Engine_Ext_Excitation on the device
class GPU_Ext_Excitation : public GPU_Extension
{
public:
	GPU_Ext_Excitation(Operator_Ext_Excitation* op_ext, GPU_Backend* b, const unsigned int lines[3]) : m_Op(op_ext), m_B(b)
	{
		for (int a=0; a<3; ++a)
			numLines[a] = lines[a];
		Setup(m_Src[0], op_ext->Volt_Count, op_ext->Volt_index, op_ext->Volt_dir, op_ext->Volt_amp, op_ext->Volt_delay, op_ext->m_Exc->GetVoltageSignal());
		Setup(m_Src[1], op_ext->Curr_Count, op_ext->Curr_index, op_ext->Curr_dir, op_ext->Curr_amp, op_ext->Curr_delay, op_ext->m_Exc->GetCurrentSignal());
	}

	virtual void Hook(bool currents, GPU_Hook h, const Engine_Box& box, unsigned int ts) {if (h==GPU_APPLY) Apply(currents, box, ts);}
	virtual bool MidStepNodes(std::vector<uint64_t>& volt, std::vector<uint64_t>&) const {AddNodes(volt, m_Src[0].flat, numLines); return true;}

protected:
	struct Sources
	{
		GPU_PointList list;
		bool serial = false;   //!< points share an edge: in order
		std::vector<uint64_t> flat;   //!< the flat index of each point
		Op_Excitation op;
	};

	void Setup(Sources& s, size_t count, unsigned int* const index[3], const unsigned short* dir, const FDTD_FLOAT* amp,
	           const unsigned int* delay, const FDTD_FLOAT* signal)
	{
		std::vector<unsigned int> x(count);
		std::vector<uint64_t> flat(count);
		std::set<uint64_t> edges;
		for (size_t n=0; n<count; ++n)
		{
			x[n] = index[0][n];
			flat[n] = FlatIndex(numLines, dir[n], index[0][n], index[1][n], index[2][n]);
			s.serial |= !edges.insert(flat[n]).second;
		}
		s.list.Init(x);
		s.op.index = Upload(m_B, s.list.Sorted(flat));
		s.flat.swap(flat);
		s.op.amp = Upload(m_B, s.list.Sorted(std::vector<float>(amp, amp+count)));
		s.op.delay = Upload(m_B, s.list.Sorted(std::vector<uint32_t>(delay, delay+count)));
		const unsigned int length = m_Op->m_Exc->GetLength();
		s.op.signal = Upload(m_B, std::vector<float>(signal, signal+length));
		s.op.length = length;
	}

	void Apply(int f, const Engine_Box& box, unsigned int ts)
	{
		Sources& s = m_Src[f];
		uint64_t begin, end;
		s.list.Range(box, begin, end);
		if (begin>=end)
			return;
		s.op.field = m_B->DeviceField(f==1);
		s.op.numTS = (int)ts;
		s.op.period = (int)ts+1;
		if (m_Op->m_Exc->GetSignalPeriod()>0)
			s.op.period = int(m_Op->m_Exc->GetSignalPeriod()/m_Op->m_Exc->GetTimestep());
		m_B->RunPoints(GPU_OP_EXCITATION, &s.op, begin, end, s.serial);
	}

	Operator_Ext_Excitation* m_Op;
	GPU_Backend* m_B;
	unsigned int numLines[3];
	Sources m_Src[2];
};

/***************************************************************************************************/

//! A boundary sheet of first order Mur: Engine_Ext_Mur_ABC, and Engine_Ext_Absorbing_BC (with super-absorption of the currents)
class GPU_Ext_Sheet : public GPU_Extension
{
public:
	//! one field (voltages or currents) of the sheet: its points (i, j) from the origin \a p0 along nyP and nyPP
	struct Sheet
	{
		bool active = false;
		bool superAbsorption = false;   //!< the apply hook of the currents of Engine_Ext_Absorbing_BC
		unsigned int line, shift;       //!< along ny
		std::vector<uint64_t> nodes;    //!< the nodes on the sheet and next to it
		GPU_PointList list;
		Op_Sheet op;
	};

	GPU_Ext_Sheet(GPU_Backend* b, const unsigned int lines[3], int ny, unsigned int start_TS) : m_B(b), m_ny(ny), m_StartTS(start_TS)
	{
		for (int a=0; a<3; ++a)
			numLines[a] = lines[a];
	}

	//! the points of field \a f: (i, j) in [0, ni) x [0, nj) at p0 + (i, j) along nyP and nyPP, K1 and K2 by \a K(c, i, j)
	typedef std::function<float(int, unsigned int, unsigned int)> KFunc;
	void Setup(int f, unsigned int line, unsigned int shift, unsigned int p0P, unsigned int p0PP, unsigned int ni, unsigned int nj,
	           const KFunc& K1, const KFunc& K2, bool superAbsorption)
	{
		Sheet& S = m_Sheet[f];
		S.active = true;
		S.superAbsorption = superAbsorption;
		S.line = line;
		S.shift = shift;
		const int nyP = (m_ny+1)%3, nyPP = (m_ny+2)%3;
		const uint64_t nodes = (uint64_t)numLines[0]*numLines[1]*numLines[2];
		std::vector<unsigned int> x;
		std::vector<uint64_t> il, is;
		std::vector<float> k1[2], k2[2];
		unsigned int pos[3];
		for (unsigned int i=0; i<ni; ++i)
			for (unsigned int j=0; j<nj; ++j)
			{
				pos[nyP] = p0P+i;
				pos[nyPP] = p0PP+j;
				pos[m_ny] = line;
				x.push_back(pos[0]);
				il.push_back(((uint64_t)pos[0]*numLines[1] + pos[1])*numLines[2] + pos[2]);
				pos[m_ny] = shift;
				is.push_back(((uint64_t)pos[0]*numLines[1] + pos[1])*numLines[2] + pos[2]);
				for (int c=0; c<2; ++c)
				{
					k1[c].push_back(K1(c, i, j));
					k2[c].push_back(superAbsorption ? K2(c, i, j) : 0.0f);
				}
			}
		S.list.Init(x);
		S.nodes = il;
		S.nodes.insert(S.nodes.end(), is.begin(), is.end());
		S.op.line = Upload(m_B, S.list.Sorted(il));
		S.op.shift = Upload(m_B, S.list.Sorted(is));
		for (int c=0; c<2; ++c)
		{
			S.op.K1[c] = Upload(m_B, S.list.Sorted(k1[c]));
			S.op.K2[c] = Upload(m_B, S.list.Sorted(k2[c]));
			S.op.s[c] = Upload(m_B, std::vector<float>(x.size(), 0.0f));
		}
		S.op.offset[0] = nyP*nodes;
		S.op.offset[1] = nyPP*nodes;
	}

	virtual void Hook(bool currents, GPU_Hook h, const Engine_Box& box, unsigned int ts) {Run(currents, h, box, ts);}
	//! the voltages of the sheet are written after the update, the currents of the sheet and the line next to it read before
	virtual bool MidStepNodes(std::vector<uint64_t>& volt, std::vector<uint64_t>& curr) const
	{
		if (m_Sheet[0].active)
			volt.insert(volt.end(), m_Sheet[0].nodes.begin(), m_Sheet[0].nodes.begin() + m_Sheet[0].nodes.size()/2);
		if (m_Sheet[1].active)
			curr.insert(curr.end(), m_Sheet[1].nodes.begin(), m_Sheet[1].nodes.end());
		return true;
	}

protected:
	void Run(int f, GPU_Hook hook, const Engine_Box& box, unsigned int ts)
	{
		Sheet& S = m_Sheet[f];
		if (!S.active || (ts<m_StartTS))
			return;
		// the box must hold the line next to the sheet as well: the post hook reads it after the update of the box
		if ((S.line<box.start[m_ny]) || (S.line>=box.stop[m_ny]))
			return;
		if ((S.shift<box.start[m_ny]) || (S.shift>=box.stop[m_ny]))
			throw std::logic_error("GPU_Ext_Sheet: a box holds the sheet but not the line next to it");
		uint64_t begin, end;
		S.list.Range(box, begin, end);
		S.op.field = m_B->DeviceField(f==1);
		// Op_Sheet numbers its hooks as GPU_Hook, plus the apply hook with super-absorption
		S.op.hook = ((hook==GPU_APPLY) && S.superAbsorption) ? 3 : (int)hook;
		m_B->RunPoints(GPU_OP_SHEET, &S.op, begin, end, false);
	}

	GPU_Backend* m_B;
	unsigned int numLines[3];
	int m_ny;
	unsigned int m_StartTS;
	Sheet m_Sheet[2];
};

class GPU_Ext_Mur_ABC
{
public:
	static GPU_Extension* New(Operator_Ext_Mur_ABC* op, Engine_Ext_Mur_ABC* eng_ext, GPU_Backend* b, const unsigned int numLines[3])
	{
		GPU_Ext_Sheet* s = new GPU_Ext_Sheet(b, numLines, op->m_ny, eng_ext->GetStartTimestep());
		auto K1 = [op](int c, unsigned int i, unsigned int j) {return (c==0) ? op->m_Mur_Coeff_nyP(i, j) : op->m_Mur_Coeff_nyPP(i, j);};
		s->Setup(0, op->m_LineNr, op->m_LineNr_Shift, 0, 0, op->m_numLines[0], op->m_numLines[1], K1, K1, false);
		return s;
	}
};

class GPU_Ext_Absorbing_BC
{
public:
	static GPU_Extension* New(Operator_Ext_Absorbing_BC* op, GPU_Backend* b, const unsigned int numLines[3])
	{
		GPU_Ext_Sheet* s = new GPU_Ext_Sheet(b, numLines, op->m_ny, 0);
		const unsigned int* p0 = op->m_sheetX0;
		const bool pos = op->m_normalSignPositive;
		auto K1 = [op](int c, unsigned int i, unsigned int j) {return (c==0) ? op->m_K1_nyP(i, j) : op->m_K1_nyPP(i, j);};
		auto K2 = [op](int c, unsigned int i, unsigned int j) {return (c==0) ? op->m_K2_nyP(i, j) : op->m_K2_nyPP(i, j);};
		s->Setup(0, p0[op->m_ny], p0[op->m_ny] + (pos ? 1 : -1), p0[op->m_nyP], p0[op->m_nyPP], op->m_numLines[0], op->m_numLines[1], K1, K2, false);
		// the currents: super-absorption only, a line less (dual grid)
		if (op->m_ABCtype==Operator_Ext_Absorbing_BC::MUR_1ST_SA)
			s->Setup(1, p0[op->m_ny] + (pos ? 0 : -1), p0[op->m_ny] + (pos ? 1 : -2), p0[op->m_nyP], p0[op->m_nyPP],
			         op->m_numLines[0]-1, op->m_numLines[1]-1, K1, K2, true);
		return s;
	}
};

/***************************************************************************************************/

//! Engine_Ext_LumpedRLC on the device
class GPU_Ext_LumpedRLC : public GPU_Extension
{
public:
	GPU_Ext_LumpedRLC(Operator_Ext_LumpedRLC* op, GPU_Backend* b, const unsigned int numLines[3]) : m_B(b)
	{
		const size_t count = op->RLC_count;
		std::vector<unsigned int> x(count);
		std::vector<uint64_t> flat(count);
		std::set<uint64_t> edges;
		for (size_t n=0; n<count; ++n)
		{
			x[n] = op->v_RLC_pos[0][n];
			flat[n] = FlatIndex(numLines, op->v_RLC_dir[n], op->v_RLC_pos[0][n], op->v_RLC_pos[1][n], op->v_RLC_pos[2][n]);
			m_Serial |= !edges.insert(flat[n]).second;
		}
		m_List.Init(x);
		AddNodes(m_Nodes, flat, numLines);
		auto up = [&](const FDTD_FLOAT* v) {return Upload(m_B, m_List.Sorted(std::vector<float>(v, v+count)));};
		m_Op.index = Upload(m_B, m_List.Sorted(flat));
		m_Op.ilv = up(op->v_RLC_ilv);
		m_Op.i2v = up(op->v_RLC_i2v);
		m_Op.vv2 = up(op->v_RLC_vv2);
		m_Op.vj1 = up(op->v_RLC_vj1);
		m_Op.vj2 = up(op->v_RLC_vj2);
		m_Op.vvd = up(op->v_RLC_vvd);
		m_Op.ib0 = up(op->v_RLC_ib0);
		m_Op.b1 = up(op->v_RLC_b1);
		m_Op.b2 = up(op->v_RLC_b2);
		const std::vector<float> zero(std::max<size_t>(count, 1), 0.0f);
		m_Op.Il = Upload(m_B, zero);
		for (int r=0; r<3; ++r)
		{
			m_Vd[r] = Upload(m_B, zero);
			m_J[r] = Upload(m_B, zero);
		}
	}

	virtual void Hook(bool currents, GPU_Hook h, const Engine_Box& box, unsigned int ts)
	{
		if (!currents && (h!=GPU_POST))
			Run(box, ts, h==GPU_PRE);
	}
	virtual bool MidStepNodes(std::vector<uint64_t>& volt, std::vector<uint64_t>&) const {volt.insert(volt.end(), m_Nodes.begin(), m_Nodes.end()); return true;}

protected:
	// the rings hold timesteps ts, ts-1, ts-2 at [ts%3], [(ts+2)%3], [(ts+1)%3], see Engine_Ext_LumpedRLC
	void Run(const Engine_Box& box, unsigned int ts, bool pre)
	{
		uint64_t begin, end;
		m_List.Range(box, begin, end);
		if (begin>=end)
			return;
		m_Op.volt = m_B->DeviceField(false);
		m_Op.Vd0 = m_Vd[ts%3];
		m_Op.Vd1 = m_Vd[(ts+2)%3];
		m_Op.Vd2 = m_Vd[(ts+1)%3];
		m_Op.J0 = m_J[ts%3];
		m_Op.J1 = m_J[(ts+2)%3];
		m_Op.J2 = m_J[(ts+1)%3];
		if (pre)
		{
			m_Op.mode = 0;
			m_B->RunPoints(GPU_OP_LUMPED_RLC, &m_Op, begin, end, false);
			return;
		}
		// the passes of the engine extension: elements on one edge see each other's values as there
		for (int mode=1; mode<=3; ++mode)
		{
			m_Op.mode = mode;
			m_B->RunPoints(GPU_OP_LUMPED_RLC, &m_Op, begin, end, m_Serial && (mode==3));
		}
	}

	GPU_Backend* m_B;
	GPU_PointList m_List;
	std::vector<uint64_t> m_Nodes;
	bool m_Serial = false;
	Op_LumpedRLC m_Op;
	float* m_Vd[3];
	float* m_J[3];
};

/***************************************************************************************************/

//! Engine_Ext_TFSF on the device: the points of the planes in passes, a field value at most once per pass, in the order of the engine extension
class GPU_Ext_TFSF : public GPU_Extension
{
public:
	GPU_Ext_TFSF(Operator_Ext_TFSF* op, GPU_Backend* b, const unsigned int numLines[3]) : m_Op(op), m_B(b)
	{
		for (int f=0; f<2; ++f)
		{
			const bool volt = (f==0);
			// the points in the order of Engine_Ext_TFSF::PlanesImpl()
			std::vector<unsigned int> x;
			std::vector<uint64_t> flat;
			std::vector<uint32_t> delay;
			std::vector<float> dd, amp;
			unsigned int pos[3];
			for (int n=0; n<3; ++n)
			{
				const int nP = (n+1)%3, nPP = (n+2)%3;
				for (int side=0; side<2; ++side)
				{
					if (!op->m_ActiveDir[n][side])
						continue;
					pos[n] = (side==1) ? op->m_Stop[n] : (volt ? op->m_Start[n] : op->m_Start[n]-1);
					for (pos[nP]=op->m_Start[nP]; pos[nP]<op->m_Start[nP]+op->m_numLines[nP]; ++pos[nP])
						for (pos[nPP]=op->m_Start[nPP]; pos[nPP]<op->m_Start[nPP]+op->m_numLines[nPP]; ++pos[nPP])
						{
							if ((pos[0]>=numLines[0]) || (pos[1]>=numLines[1]) || (pos[2]>=numLines[2]))
								continue;
							const unsigned int ui_pos = (pos[nP]-op->m_Start[nP])*op->m_numLines[nPP] + (pos[nPP]-op->m_Start[nPP]);
							for (int c=0; c<2; ++c)
							{
								x.push_back(pos[0]);
								flat.push_back(FlatIndex(numLines, (c==0) ? nP : nPP, pos[0], pos[1], pos[2]));
								delay.push_back(volt ? op->m_VoltDelay[n][side][c][ui_pos] : op->m_CurrDelay[n][side][c][ui_pos]);
								dd.push_back(volt ? op->m_VoltDelayDelta[n][side][c][ui_pos] : op->m_CurrDelayDelta[n][side][c][ui_pos]);
								amp.push_back(volt ? op->m_VoltAmp[n][side][c][ui_pos] : op->m_CurrAmp[n][side][c][ui_pos]);
							}
						}
				}
			}
			if (volt)
				AddNodes(m_VoltNodes, flat, numLines);
			// pass k: the k-th addition to a field value
			std::map<uint64_t, unsigned int> seen;
			std::vector<unsigned int> pass(flat.size());
			for (size_t p=0; p<flat.size(); ++p)
				pass[p] = seen[flat[p]]++;
			unsigned int passes = 0;
			for (unsigned int k : pass)
				passes = std::max(passes, k+1);
			const unsigned int length = op->m_Exc->GetLength();
			// a voltage gets the H-field signal and the other way round
			const FDTD_FLOAT* signal = volt ? op->m_Exc->GetCurrentSignal() : op->m_Exc->GetVoltageSignal();
			float* sig = Upload(m_B, std::vector<float>(signal, signal+length));
			m_Pass[f].resize(passes);
			for (unsigned int k=0; k<passes; ++k)
			{
				std::vector<size_t> sel;
				for (size_t p=0; p<flat.size(); ++p)
					if (pass[p]==k)
						sel.push_back(p);
				auto pick = [&](const auto& v) {std::vector<typename std::decay<decltype(v)>::type::value_type> r; for (size_t p : sel) r.push_back(v[p]); return r;};
				Pass& P = m_Pass[f][k];
				P.list.Init(pick(x));
				P.op.index = Upload(m_B, P.list.Sorted(pick(flat)));
				P.op.delay = Upload(m_B, P.list.Sorted(pick(delay)));
				P.op.dd = Upload(m_B, P.list.Sorted(pick(dd)));
				P.op.amp = Upload(m_B, P.list.Sorted(pick(amp)));
				P.op.signal = sig;
				P.op.length = length;
				P.op.period = int(op->m_Exc->GetSignalPeriod()/op->m_Exc->GetTimestep());
			}
		}
	}

	virtual void Hook(bool currents, GPU_Hook h, const Engine_Box& box, unsigned int ts) {if (h==GPU_POST) Run(currents, box, ts);}
	virtual bool MidStepNodes(std::vector<uint64_t>& volt, std::vector<uint64_t>&) const {volt.insert(volt.end(), m_VoltNodes.begin(), m_VoltNodes.end()); return true;}

protected:
	struct Pass
	{
		GPU_PointList list;
		Op_TFSF op;
	};

	void Run(int f, const Engine_Box& box, unsigned int ts)
	{
		for (Pass& P : m_Pass[f])
		{
			uint64_t begin, end;
			P.list.Range(box, begin, end);
			P.op.field = m_B->DeviceField(f==1);
			P.op.ts = ts;
			m_B->RunPoints(GPU_OP_TFSF, &P.op, begin, end, false);
		}
	}

	Operator_Ext_TFSF* m_Op;
	GPU_Backend* m_B;
	std::vector<Pass> m_Pass[2];
	std::vector<uint64_t> m_VoltNodes;
};

/***************************************************************************************************/

//! Engine_Ext_SteadyState on the device: the probe voltages recorded per timestep, its check once per period on the host
class GPU_Ext_SteadyState : public GPU_Extension
{
public:
	GPU_Ext_SteadyState(Operator_Ext_SteadyState* op, Engine_Ext_SteadyState* eng_ext, GPU_Backend* b, const unsigned int numLines[3])
		: m_EngExt(eng_ext), m_B(b), m_Period(op->m_TS_period)
	{
		const size_t count = op->m_E_probe_dir.size();
		std::vector<unsigned int> x(count);
		std::vector<uint64_t> flat(count);
		for (size_t n=0; n<count; ++n)
		{
			x[n] = op->m_E_probe_pos[0].at(n);
			flat[n] = FlatIndex(numLines, op->m_E_probe_dir.at(n), op->m_E_probe_pos[0].at(n), op->m_E_probe_pos[1].at(n), op->m_E_probe_pos[2].at(n));
		}
		m_List.Init(x);
		m_Op.index = Upload(m_B, m_List.Sorted(flat));
		m_Op.length = 2*m_Period;
		m_Host.assign(std::max<size_t>(count, 1)*2*(size_t)m_Period, 0.0f);
		m_Op.records = Upload(m_B, m_Host);
	}

	virtual void Hook(bool currents, GPU_Hook h, const Engine_Box& box, unsigned int ts)
	{
		if (currents || (h!=GPU_APPLY))
			return;
		Record(box, ts);
		// the check reads the records of all probes: a timestep of the whole domain (NeedsFlatStep())
		if (NeedsFlatStep(ts))
		{
			// the records of the sorted probes, in the order of the engine extension
			m_B->DeviceDownload(m_Host.data(), m_Op.records, m_Host.size()*sizeof(float));
			std::vector<size_t> sorted(m_List.Size());
			for (size_t n=0; n<sorted.size(); ++n)
				sorted[m_List.Order()[n]] = n;
			const float* rec = m_Host.data();
			const unsigned int len = 2*m_Period;
			m_EngExt->CheckPeriod(ts%len, [&](size_t n, unsigned int t) {return (double)rec[sorted[n]*len + t];});
		}
	}
	virtual bool NeedsFlatStep(unsigned int ts) const {return m_EngExt->NeedsFlatStep(ts);}
	virtual bool MidStepNodes(std::vector<uint64_t>&, std::vector<uint64_t>&) const {return true;}   // reads voltages only

protected:
	void Record(const Engine_Box& box, unsigned int ts)
	{
		uint64_t begin, end;
		m_List.Range(box, begin, end);
		m_Op.volt = m_B->DeviceField(false);
		m_Op.pos = ts%(2*m_Period);
		m_B->RunPoints(GPU_OP_STEADY_STATE, &m_Op, begin, end, false);
	}

	Engine_Ext_SteadyState* m_EngExt;
	GPU_Backend* m_B;
	unsigned int m_Period;
	GPU_PointList m_List;
	Op_SteadyState m_Op;
	std::vector<float> m_Host;
};

/***************************************************************************************************/

//! Engine_Ext_Cylinder on the device (a closed alpha mesh)
class GPU_Ext_Cylinder : public GPU_Extension
{
public:
	GPU_Ext_Cylinder(Operator_Ext_Cylinder* op, GPU_Backend* b, const unsigned int lines[3]) : m_B(b), m_Closed(op->CC_closedAlpha), m_R0(op->CC_R0_included)
	{
		for (int a=0; a<3; ++a)
			numLines[a] = lines[a];
		m_Op.ny = numLines[1];
		m_Op.nz = numLines[2];
		m_Op.N = (uint64_t)numLines[0]*numLines[1]*numLines[2];
		m_Op.vv_R0 = Upload(m_B, m_R0 ? std::vector<float>(op->vv_R0, op->vv_R0+numLines[2]) : std::vector<float>());
		m_Op.vi_R0 = Upload(m_B, m_R0 ? std::vector<float>(op->vi_R0, op->vi_R0+numLines[2]) : std::vector<float>());
	}

	virtual void Hook(bool currents, GPU_Hook h, const Engine_Box& box, unsigned int) {if (h==GPU_POST) Run(!currents, box);}
	//! the voltages of the first alpha plane, and on the r=0 axis the voltages and the alpha currents it reads
	virtual bool MidStepNodes(std::vector<uint64_t>& volt, std::vector<uint64_t>& curr) const
	{
		if (!m_Closed)
			return true;
		for (unsigned int r=0; r<numLines[0]; ++r)
			for (unsigned int z=0; z<numLines[2]; ++z)
				volt.push_back((uint64_t)r*numLines[1]*numLines[2] + z);
		if (m_R0)
			for (uint64_t n=0; n<(uint64_t)numLines[1]*numLines[2]; ++n)
			{
				volt.push_back(n);
				curr.push_back(n);
			}
		return true;
	}

protected:
	// as Engine_Ext_Cylinder::BoxImpl(): the box spans the whole alpha axis
	void Run(bool volt, const Engine_Box& box)
	{
		if (!m_Closed)
			return;
		if ((box.start[1]!=0) || (box.stop[1]<numLines[1]))
			throw std::logic_error("GPU_Ext_Cylinder: a box does not span the whole alpha axis");
		m_Op.volt = m_B->DeviceField(false);
		m_Op.curr = m_B->DeviceField(true);
		const unsigned int r1 = std::min(box.stop[0], volt ? numLines[0] : numLines[0]-1);
		const unsigned int z1 = std::min(box.stop[2], volt ? numLines[2] : numLines[2]-1);
		const unsigned int r0 = box.start[0], z0 = box.start[2];
		if ((r0>=r1) || (z0>=z1))
			return;
		if (volt && m_R0 && (box.start[0]==0))
		{
			// the r=0 axis: per z its voltage, then per (alpha, z) its copies
			m_Op.mode = 0;
			m_B->RunPoints(GPU_OP_CYLINDER, &m_Op, z0, z1, false);
			m_Op.mode = 1;
			m_B->RunPoints(GPU_OP_CYLINDER, &m_Op, 0, (uint64_t)numLines[1]*numLines[2], false);
		}
		m_Op.mode = volt ? 2 : 3;
		m_Op.r0 = r0;
		m_Op.z0 = z0;
		m_Op.nzr = z1-z0;
		m_B->RunPoints(GPU_OP_CYLINDER, &m_Op, 0, (uint64_t)(r1-r0)*(z1-z0), false);
	}

	GPU_Backend* m_B;
	unsigned int numLines[3];
	bool m_Closed, m_R0;
	Op_Cylinder m_Op;
};

/***************************************************************************************************/

//! Engine_Ext_LorentzMaterial (also the conducting sheets) and Engine_Ext_DebyeMaterial on the device, for the dispersive
//! materials the cell model does not fuse (CellModel::Analyse(): a cylindrical mesh, more than four poles of a kind, ...)
class GPU_Ext_Dispersive : public GPU_Extension
{
public:
	GPU_Ext_Dispersive(Operator_Ext_LorentzMaterial* op, GPU_Backend* b, const unsigned int lines[3]) : m_B(b), m_Debye(false)
	{
		Init(lines);
		for (int o=0; o<op->m_Order; ++o)
			for (int f=0; f<2; ++f)
			{
				if (!(f ? op->m_curr_ADE_On[o] : op->m_volt_ADE_On[o]))
					continue;
				const bool lorentz = f ? op->m_curr_Lor_ADE_On[o] : op->m_volt_Lor_ADE_On[o];
				Order& ord = Add(f, op->m_LM_Count.at(o), op->m_LM_pos[o]);
				ord.op.lorentz = lorentz;
				ord.op.k_int = Upload(m_B, Gather(ord, 1, f ? &op->i_int_ADE[o] : &op->v_int_ADE[o]));
				ord.op.k_ext = Upload(m_B, Gather(ord, 1, f ? &op->i_ext_ADE[o] : &op->v_ext_ADE[o]));
				if (lorentz)
				{
					ord.op.k_lor = Upload(m_B, Gather(ord, 1, f ? &op->i_Lor_ADE[o] : &op->v_Lor_ADE[o]));
					ord.op.lor = Upload(m_B, std::vector<float>(3*ord.op.count, 0.0f));
				}
			}
	}

	GPU_Ext_Dispersive(Operator_Ext_DebyeMaterial* op, GPU_Backend* b, const unsigned int lines[3]) : m_B(b), m_Debye(true)
	{
		Init(lines);
		// one cell list for all poles, the voltages only (Operator_Ext_DebyeMaterial)
		Order& ord = Add(0, op->m_LM_Count.at(0), op->m_LM_pos[0]);
		ord.op.poles = op->m_PoleCount;
		ord.op.k_int = Upload(m_B, Gather(ord, op->m_PoleCount, op->v_relax_ADE));
		ord.op.k_ext = Upload(m_B, Gather(ord, op->m_PoleCount, op->v_drive_ADE));
		ord.op.solve = Upload(m_B, Gather(ord, 1, &op->v_solve_ADE));
		ord.op.pole = Upload(m_B, std::vector<float>(3*ord.op.count*op->m_PoleCount, 0.0f));
		ord.op.pre = Upload(m_B, std::vector<float>(3*ord.op.count, 0.0f));
	}

	//! the modes of Op_Dispersive per hook; the Debye material has voltages only
	virtual void Hook(bool currents, GPU_Hook h, const Engine_Box& box, unsigned int)
	{
		if (h==GPU_APPLY)
			Run(currents, 3, box);
		else if (h==GPU_PRE)
			Run(currents, m_Debye ? 1 : 0, box);
		else if (m_Debye)
			Run(currents, 2, box);
	}
	//! the voltages its apply pass writes; the currents its pre pass reads
	virtual bool MidStepNodes(std::vector<uint64_t>& volt, std::vector<uint64_t>& curr) const
	{
		for (int f=0; f<2; ++f)
			for (const Order& o : m_Orders[f])
				(f ? curr : volt).insert((f ? curr : volt).end(), o.nodes.begin(), o.nodes.end());
		return true;
	}

protected:
	//! the cells of one order of the voltages or currents, sorted along x
	struct Order
	{
		GPU_PointList list;
		bool serial = false;           //!< a cell twice: the apply passes in order
		std::vector<uint64_t> nodes;   //!< the node of each cell, in the order of the operator extension
		Op_Dispersive op;
	};

	void Init(const unsigned int lines[3])
	{
		for (int a=0; a<3; ++a)
			numLines[a] = lines[a];
	}

	//! the order of field \a f with \a count cells at \a pos, with its state of the apply pass
	Order& Add(int f, size_t count, unsigned int* const pos[3])
	{
		m_Orders[f].emplace_back();
		Order& o = m_Orders[f].back();
		std::memset(&o.op, 0, sizeof(o.op));
		std::vector<unsigned int> x(count);
		o.nodes.resize(count);
		for (size_t i=0; i<count; ++i)
		{
			x[i] = pos[0][i];
			o.nodes[i] = ((uint64_t)pos[0][i]*numLines[1] + pos[1][i])*numLines[2] + pos[2][i];
		}
		std::vector<uint64_t> sorted(o.nodes);
		std::sort(sorted.begin(), sorted.end());
		o.serial = std::adjacent_find(sorted.begin(), sorted.end())!=sorted.end();
		o.list.Init(x);
		o.op.node = Upload(m_B, o.list.Sorted(o.nodes));
		o.op.sN = (uint64_t)numLines[0]*numLines[1]*numLines[2];
		o.op.count = count;
		o.op.ade = Upload(m_B, std::vector<float>(3*count, 0.0f));
		return o;
	}

	//! the coefficients c[k][n][i] of \a poles poles at [(k*3 + n)*count + p] of the sorted points p
	std::vector<float> Gather(const Order& o, int poles, FDTD_FLOAT** const* c) const
	{
		const size_t count = o.op.count;
		std::vector<float> v(3*count*poles);
		for (int k=0; k<poles; ++k)
			for (int n=0; n<3; ++n)
				for (size_t p=0; p<count; ++p)
					v[(k*3 + n)*count + p] = c[k][n][o.list.Order()[p]];
		return v;
	}

	void Run(int f, int mode, const Engine_Box& box)
	{
		for (Order& o : m_Orders[f])
		{
			uint64_t begin, end;
			o.list.Range(box, begin, end);
			o.op.field = m_B->DeviceField(f==1);
			o.op.mode = mode;
			m_B->RunPoints(GPU_OP_DISPERSIVE, &o.op, begin, end, o.serial && (mode==3));
		}
	}

	GPU_Backend* m_B;
	bool m_Debye;
	unsigned int numLines[3];
	std::vector<Order> m_Orders[2];   //!< of the voltages, of the currents
};

GPU_Extension* GPU_CreateExtension(Engine_Extension* eng_ext, GPU_Backend* backend, const unsigned int numLines[3])
{
	if (dynamic_cast<Engine_Ext_Excitation*>(eng_ext))
		if (Operator_Ext_Excitation* op = dynamic_cast<Operator_Ext_Excitation*>(eng_ext->GetOperatorExtension()))
			return new GPU_Ext_Excitation(op, backend, numLines);
	if (Engine_Ext_Mur_ABC* e = dynamic_cast<Engine_Ext_Mur_ABC*>(eng_ext))
		if (Operator_Ext_Mur_ABC* op = dynamic_cast<Operator_Ext_Mur_ABC*>(eng_ext->GetOperatorExtension()))
			return GPU_Ext_Mur_ABC::New(op, e, backend, numLines);
	if (dynamic_cast<Engine_Ext_Absorbing_BC*>(eng_ext))
		if (Operator_Ext_Absorbing_BC* op = dynamic_cast<Operator_Ext_Absorbing_BC*>(eng_ext->GetOperatorExtension()))
			return GPU_Ext_Absorbing_BC::New(op, backend, numLines);
	if (dynamic_cast<Engine_Ext_LumpedRLC*>(eng_ext))
		if (Operator_Ext_LumpedRLC* op = dynamic_cast<Operator_Ext_LumpedRLC*>(eng_ext->GetOperatorExtension()))
			return new GPU_Ext_LumpedRLC(op, backend, numLines);
	if (dynamic_cast<Engine_Ext_TFSF*>(eng_ext))
		if (Operator_Ext_TFSF* op = dynamic_cast<Operator_Ext_TFSF*>(eng_ext->GetOperatorExtension()))
			return new GPU_Ext_TFSF(op, backend, numLines);
	if (Engine_Ext_SteadyState* e = dynamic_cast<Engine_Ext_SteadyState*>(eng_ext))
		if (Operator_Ext_SteadyState* op = dynamic_cast<Operator_Ext_SteadyState*>(eng_ext->GetOperatorExtension()))
			return new GPU_Ext_SteadyState(op, e, backend, numLines);
	if (dynamic_cast<Engine_Ext_Cylinder*>(eng_ext))
		if (Operator_Ext_Cylinder* op = dynamic_cast<Operator_Ext_Cylinder*>(eng_ext->GetOperatorExtension()))
			return new GPU_Ext_Cylinder(op, backend, numLines);
	if (dynamic_cast<Engine_Ext_LorentzMaterial*>(eng_ext))
		if (Operator_Ext_LorentzMaterial* op = dynamic_cast<Operator_Ext_LorentzMaterial*>(eng_ext->GetOperatorExtension()))
			return new GPU_Ext_Dispersive(op, backend, numLines);
	if (dynamic_cast<Engine_Ext_DebyeMaterial*>(eng_ext))
		if (Operator_Ext_DebyeMaterial* op = dynamic_cast<Operator_Ext_DebyeMaterial*>(eng_ext->GetOperatorExtension()))
			return new GPU_Ext_Dispersive(op, backend, numLines);
	return NULL;
}

/***************************************************************************************************/

namespace
{
//! the operation of the coupling of grids with \a nb and \a ns lines
Op_MultiGrid MultiGridOp(const unsigned int nb[3], const unsigned int ns[3])
{
	Op_MultiGrid op;
	std::memset(&op, 0, sizeof(op));
	for (int a=0; a<3; ++a)
	{
		op.nb[a] = nb[a];
		op.ns[a] = ns[a];
	}
	return op;
}
}

SubGrid_Coupling_GPU::SubGrid_Coupling_GPU(Engine_GPU* base, Engine_GPU* sub, const Operator_CylinderMultiGrid* op)
	: m_Base(base->GetBackend()), m_Sub(sub->GetBackend())
{
	const MultiGridInterpolation& interpol = op->GetInterpolation();
	m_SplitPos = interpol.split_pos;
	for (int a=0; a<3; ++a)
	{
		m_BaseLines[a] = op->GetNumberOfLines(a, true);
		m_SubLines[a] = op->GetInnerOperator()->GetNumberOfLines(a, true);
	}
	SetWeights(m_Volt, interpol.pos_v_2p, interpol.pos_v_2pp, interpol.w_v_2p, interpol.w_v_2pp);
	SetWeights(m_Curr, interpol.pos_i_2p, interpol.pos_i_2pp, interpol.w_i_2p, interpol.w_i_2pp);
}

void SubGrid_Coupling_GPU::SyncVoltages()
{
	if (m_SplitPos<1)
		return;
	Op_MultiGrid op = MultiGridOp(m_BaseLines, m_SubLines);
	op.base = m_Base->DeviceField(false);
	op.sub = m_Sub->DeviceField(false);
	op.r0 = m_SplitPos-1;
	op.mode = 0;
	m_Base->RunPoints(GPU_OP_MULTIGRID, &op, 0, (uint64_t)(m_BaseLines[1]/2)*m_BaseLines[2], false);
}

void SubGrid_Coupling_GPU::SyncCurrents()
{
	if (m_SplitPos>=2)
		Interpolate(true, m_SplitPos-2, 1);
}

void SubGrid_Coupling_GPU::Interpolate()
{
	if (m_SplitPos>=1)
		Interpolate(false, 0, m_SplitPos-1);
	if (m_SplitPos>=2)
		Interpolate(true, 0, m_SplitPos-2);
}

void SubGrid_Coupling_GPU::SetWeights(Weights& w, const std::vector<unsigned int> p2p[2], const std::vector<unsigned int> p2pp[2], const std::vector<float> w2p[2], const std::vector<float> w2pp[2])
{
	// in the memory of the base grid, which outlives the link
	std::vector<uint32_t> a, b;
	std::vector<float> c, d;
	for (int n=0; n<2; ++n)
	{
		a.insert(a.end(), p2p[n].begin(), p2p[n].end());
		b.insert(b.end(), p2pp[n].begin(), p2pp[n].end());
		c.insert(c.end(), w2p[n].begin(), w2p[n].end());
		d.insert(d.end(), w2pp[n].begin(), w2pp[n].end());
	}
	w.p2p = Upload(m_Base, a);
	w.p2pp = Upload(m_Base, b);
	w.w2p = Upload(m_Base, c);
	w.w2pp = Upload(m_Base, d);
}

void SubGrid_Coupling_GPU::Interpolate(bool currents, unsigned int r0, unsigned int count)
{
	const Weights& w = currents ? m_Curr : m_Volt;
	Op_MultiGrid op = MultiGridOp(m_BaseLines, m_SubLines);
	op.base = m_Base->DeviceField(currents);
	op.sub = m_Sub->DeviceField(currents);
	op.pos_2p = w.p2p;
	op.pos_2pp = w.p2pp;
	op.w_2p = w.w2p;
	op.w_2pp = w.w2pp;
	op.r0 = r0;
	op.mode = 1;
	m_Base->RunPoints(GPU_OP_MULTIGRID, &op, 0, (uint64_t)count*m_BaseLines[1]*m_BaseLines[2], false);
}
