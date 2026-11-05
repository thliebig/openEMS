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

#include "cell_model.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <thread>
#include <tuple>

#include "tools/global.h"
#include "tools/useful.h"
#include "operator.h"
#include "extensions/operator_ext_upml.h"
#include "extensions/operator_ext_cylinder.h"
#include "extensions/operator_ext_lorentzmaterial.h"
#include "extensions/operator_ext_debyematerial.h"
#include "extensions/operator_ext_excitation.h"
#include "extensions/operator_ext_mur_abc.h"
#include "extensions/operator_ext_lumpedRLC.h"
#include "extensions/operator_ext_steadystate.h"

using std::cout;
using std::endl;

// cells whose material part m agrees in the upper 20 bits of the mantissa (2^-20, about
// 1e-6 relative) share a run: the double rounded to them, as an integer
static const int M_DROP_BITS = 32;

static inline int64_t QuantiseM(double m)
{
	uint64_t b;
	std::memcpy(&b, &m, sizeof(b));
	return (int64_t)((b + (1ull<<(M_DROP_BITS-1))) >> M_DROP_BITS);
}

static inline double DequantiseM(int64_t q)
{
	const uint64_t b = (uint64_t)q << M_DROP_BITS;
	double m;
	std::memcpy(&m, &b, sizeof(m));
	return m;
}

// a float rounded to the upper 20 bits of its mantissa
static inline float RoundMantissa(float x)
{
	uint32_t b;
	std::memcpy(&b, &x, sizeof(b));
	b = (b + 4u) & ~7u;
	std::memcpy(&x, &b, sizeof(x));
	return x;
}

void CellModel_NodeSet::Init(size_t nodes)
{
	bits.assign((nodes+63)/64, 0);
	before.clear();
}

bool CellModel_NodeSet::Add(size_t lin)
{
	const uint64_t bit = 1ull<<(lin&63);
	if (bits[lin>>6] & bit)
		return false;
	bits[lin>>6] |= bit;
	return true;
}

void CellModel_NodeSet::Finish()
{
	before.resize(bits.size());
	uint32_t count = 0;
	for (size_t w=0; w<bits.size(); ++w)
	{
		before[w] = count;
		count += std::bitset<64>(bits[w]).count();
	}
}

int64_t CellModel_NodeSet::Find(size_t lin) const
{
	const uint64_t w = bits[lin>>6];
	const unsigned int b = lin & 63;
	if (((w>>b) & 1)==0)
		return -1;
	return before[lin>>6] + std::bitset<64>(w & ((1ull<<b)-1)).count();
}

CellModel::CellModel(const Operator* op) : Op(op), m_SkipInert(false)
{
	for (int n=0; n<3; ++n)
		numLines[n] = op->GetNumberOfLines(n, true);
}

void CellModel::Analyse(bool skipInert)
{
	m_SkipInert = skipInert;
	m_UPML.clear();
	m_Fused.clear();
	for (size_t n=0; n<Op->GetNumberOfExtentions(); ++n)
	{
		const Operator_Ext_UPML* upml = dynamic_cast<const Operator_Ext_UPML*>(Op->GetExtension(n));
		if (upml==NULL)
			continue;
		// FreeFused() frees the coefficients: an operator serves one engine (as its own, Operator::CreateEngine())
		if (upml->vv.data()==NULL)
			throw std::logic_error("CellModel: the UPML coefficients are freed, by an engine created on this operator before");
		m_UPML.push_back(upml);
		m_Fused.push_back(upml);
	}
	FindDispersive();
	for (const Dispersive& d : m_Disp)
		m_Fused.push_back(d.op);
	FindInert();
	CalcGeometry();
}

bool CellModel::IsFused(const Operator_Extension* ext) const
{
	return std::find(m_Fused.begin(), m_Fused.end(), ext)!=m_Fused.end();
}

void CellModel::FreeFused()
{
	FreeDispersive();
	// the UPML cells hold their coefficients now: free the extension's (24 bytes per node of a component pair)
	for (const Operator_Ext_UPML* upml : m_UPML)
	{
		Operator_Ext_UPML* w = const_cast<Operator_Ext_UPML*>(upml);
		w->vv.Reset();
		w->vvfo.Reset();
		w->vvfn.Reset();
		w->ii.Reset();
		w->iifo.Reset();
		w->iifn.Reset();
	}
	// the analysis: the cell sets of the dispersive extensions, the excited voltages
	m_Disp.clear();
	m_ExcitedVolt.clear();
}

void CellModel::FindDispersive()
{
	m_Disp.clear();
	bool cylinder = false;
	for (size_t n=0; n<Op->GetNumberOfExtentions(); ++n)
		if (dynamic_cast<const Operator_Ext_Cylinder*>(Op->GetExtension(n)))
			cylinder = true;
	const size_t N1 = numLines[1], N2 = numLines[2];
	const size_t nodes = (size_t)numLines[0]*N1*N2;
	int poles[2] = {0, 0};   // the Drude/Lorentz and the Debye poles a cell may have
	for (size_t n=0; n<Op->GetNumberOfExtentions(); ++n)
	{
		Dispersive d;
		d.op = dynamic_cast<const Operator_Ext_Dispersive*>(Op->GetExtension(n));
		d.lorentz = dynamic_cast<const Operator_Ext_LorentzMaterial*>(Op->GetExtension(n));
		d.debye = dynamic_cast<const Operator_Ext_DebyeMaterial*>(Op->GetExtension(n));
		if ((d.lorentz==NULL) && (d.debye==NULL))
			continue;
		const int kind = d.lorentz ? 0 : 1;
		const int np = d.lorentz ? d.lorentz->m_Order : d.debye->m_PoleCount;
		std::string why;
		if (cylinder)
			why = "a cylindrical mesh";
		else if (poles[kind]+np>MAX_POLES)
			why = "more than 4 poles of a kind";
		else if (d.debye && (poles[1]>0))
			why = "a second Debye extension";
		for (int o=0; (o<d.op->m_Order) && why.empty(); ++o)
			if ((d.op->m_LM_Count[o]>0) && (d.op->m_LM_pos[o][0]==NULL))
				why = "its cells freed by an engine before";
		d.cells.resize(d.op->m_Order);
		for (int o=0; (o<d.op->m_Order) && why.empty(); ++o)
		{
			CellModel_NodeSet& set = d.cells[o];
			set.Init(nodes);
			unsigned int* const* pos = d.op->m_LM_pos[o];
			for (size_t i=0; (i<d.op->m_LM_Count[o]) && why.empty(); ++i)
			{
				const unsigned int p[3] = {pos[0][i], pos[1][i], pos[2][i]};
				const size_t lin = (p[0]*N1 + p[1])*N2 + p[2];
				unsigned int loc[3];
				// the rank of a cell is its index in the list only if the list ascends
				if ((i>0) && (lin < (((size_t)pos[0][i-1]*N1 + pos[1][i-1])*N2 + pos[2][i-1])))
					why = "its cells out of order";
				else if (!set.Add(lin))
					why = "a cell twice";
				else if (FindUPML(p, loc))
					why = "a dispersive cell in a UPML";
			}
			set.Finish();
		}
		if (why.empty())
		{
			poles[kind] += np;
			m_Disp.push_back(d);
		}
		else if (g_settings.GetVerboseLevel()>0)
			cout << "CellModel: the " << d.op->GetExtensionName() << " stays an extension: " << why << endl;
	}
}

void CellModel::FindInert()
{
	m_ExcitedVolt.clear();
	for (size_t n=0; n<Op->GetNumberOfExtentions() && m_SkipInert; ++n)
	{
		const Operator_Extension* e = Op->GetExtension(n);
		if (const Operator_Ext_Excitation* exc = dynamic_cast<const Operator_Ext_Excitation*>(e))
		{
			for (size_t i=0; i<exc->Volt_Count; ++i)
			{
				const size_t lin = ((size_t)exc->Volt_index[0][i]*numLines[1] + exc->Volt_index[1][i])*numLines[2] + exc->Volt_index[2][i];
				m_ExcitedVolt.insert(3*lin + exc->Volt_dir[i]);
			}
			continue;
		}
		// extensions that write voltages only where the coefficients are not zero, or on the outer lines
		const bool known = dynamic_cast<const Operator_Ext_UPML*>(e) || dynamic_cast<const Operator_Ext_Mur_ABC*>(e)
		                || dynamic_cast<const Operator_Ext_LumpedRLC*>(e) || dynamic_cast<const Operator_Ext_SteadyState*>(e)
		                || dynamic_cast<const Operator_Ext_LorentzMaterial*>(e) || dynamic_cast<const Operator_Ext_DebyeMaterial*>(e);
		if (!known)
			m_SkipInert = false;
	}
}

bool CellModel::IsZeroVolt(int n, const unsigned int pos[3]) const
{
	for (int a=0; a<3; ++a)
		if ((pos[a]==0) || (pos[a]+1>=numLines[a]))
			return false;
	unsigned int p[3] = {pos[0], pos[1], pos[2]}, loc[3];
	if ((Op->GetVV(n, p)!=0) || (Op->GetVI(n, p)!=0) || FindUPML(pos, loc))
		return false;
	const size_t lin = ((size_t)pos[0]*numLines[1] + pos[1])*numLines[2] + pos[2];
	return m_ExcitedVolt.count(3*lin + n)==0;
}

bool CellModel::IsInertCurr(int n, const unsigned int pos[3]) const
{
	// I_n += iv * (V_nPP(pos) - V_nPP(pos+e_nP) - V_nP(pos) + V_nP(pos+e_nPP))
	const int nP = (n+1)%3, nPP = (n+2)%3;
	unsigned int q[3] = {pos[0], pos[1], pos[2]};
	if (!IsZeroVolt(nPP, q) || !IsZeroVolt(nP, q))
		return false;
	++q[nP];
	if (!IsZeroVolt(nPP, q))
		return false;
	--q[nP];
	++q[nPP];
	return IsZeroVolt(nP, q);
}

void CellModel::FreeDispersive()
{
	// the sets and the pole states hold it all now: up to 6 floats per cell, pole and component, and the
	// cell list; the counts and the flags stay (ShowStat()), the destructors take NULL
	auto free3 = [](FDTD_FLOAT** a) {
		if (a)
			for (int n=0; n<3; ++n)
			{
				delete[] a[n];
				a[n] = NULL;
			}
	};
	for (const Dispersive& d : m_Disp)
	{
		Operator_Ext_Dispersive* op = const_cast<Operator_Ext_Dispersive*>(d.op);
		for (int o=0; o<op->m_Order; ++o)
			for (int n=0; n<3; ++n)
			{
				delete[] op->m_LM_pos[o][n];
				op->m_LM_pos[o][n] = NULL;
			}
		if (d.lorentz)
		{
			Operator_Ext_LorentzMaterial* L = const_cast<Operator_Ext_LorentzMaterial*>(d.lorentz);
			for (int o=0; o<L->m_Order; ++o)
			{
				if (L->m_volt_ADE_On[o])
				{
					free3(L->v_int_ADE[o]);
					free3(L->v_ext_ADE[o]);
				}
				if (L->m_curr_ADE_On[o])
				{
					free3(L->i_int_ADE[o]);
					free3(L->i_ext_ADE[o]);
				}
				if (L->m_volt_Lor_ADE_On[o])
					free3(L->v_Lor_ADE[o]);
				if (L->m_curr_Lor_ADE_On[o])
					free3(L->i_Lor_ADE[o]);
			}
		}
		if (d.debye)
		{
			Operator_Ext_DebyeMaterial* D = const_cast<Operator_Ext_DebyeMaterial*>(d.debye);
			for (int o=0; o<D->m_PoleCount; ++o)
			{
				if (D->v_relax_ADE)
					free3(D->v_relax_ADE[o]);
				if (D->v_drive_ADE)
					free3(D->v_drive_ADE[o]);
			}
			free3(D->v_solve_ADE);
		}
	}
}

void CellModel::GetCellPoles(int f, int n, const unsigned int pos[3], CellKey& key) const
{
	const size_t lin = ((size_t)pos[0]*numLines[1] + pos[1])*numLines[2] + pos[2];
	// the Drude/Lorentz poles first, in the order of their extensions and orders: their apply passes subtract in it
	for (const Dispersive& d : m_Disp)
	{
		const Operator_Ext_LorentzMaterial* L = d.lorentz;
		if (L==NULL)
			continue;
		for (int o=0; o<L->m_Order; ++o)
		{
			if (!((f==0) ? L->m_volt_ADE_On[o] : L->m_curr_ADE_On[o]))
				continue;
			const int64_t i = d.cells[o].Find(lin);
			if (i<0)
				continue;
			const float ext = (f==0) ? L->v_ext_ADE[o][n][i] : L->i_ext_ADE[o][n][i];
			if (ext==0)
				continue;   // a component without the material: its J stays 0
			const bool lorentz = (f==0) ? L->m_volt_Lor_ADE_On[o] : L->m_curr_Lor_ADE_On[o];
			CellModel_Pole& p = key.poles[key.nl++];
			p.k_int = (f==0) ? L->v_int_ADE[o][n][i] : L->i_int_ADE[o][n][i];
			p.k_ext = RoundMantissa(ext);
			p.k_lor = lorentz ? RoundMantissa((f==0) ? L->v_Lor_ADE[o][n][i] : L->i_Lor_ADE[o][n][i]) : 0;
		}
	}
	// then the Debye poles (voltages only), of at most one extension (FindDispersive())
	for (const Dispersive& d : m_Disp)
	{
		const Operator_Ext_DebyeMaterial* D = d.debye;
		if ((D==NULL) || (f!=0) || !D->m_volt_ADE_On[0])
			continue;
		const int64_t i = d.cells[0].Find(lin);
		if (i<0)
			continue;
		for (int o=0; o<D->m_PoleCount; ++o)
		{
			const float drive = D->v_drive_ADE[o][n][i];
			if (drive==0)
				continue;   // a component without the material: its W stays 0
			CellModel_Pole& p = key.poles[key.nl + key.nd++];
			p.k_int = D->v_relax_ADE[o][n][i];
			p.k_ext = RoundMantissa(drive);
			p.k_lor = 0;
		}
		if (key.nd>0)
			key.solve = RoundMantissa(D->v_solve_ADE[n][i]);
	}
}

const Operator_Ext_UPML* CellModel::FindUPML(const unsigned int pos[3], unsigned int loc[3]) const
{
	for (const Operator_Ext_UPML* upml : m_UPML)
	{
		bool inside = true;
		for (int a=0; a<3 && inside; ++a)
			inside = (pos[a]>=upml->m_StartPos[a]) && (pos[a]<upml->m_StartPos[a]+upml->m_numLines[a]);
		if (!inside)
			continue;
		for (int a=0; a<3; ++a)
			loc[a] = pos[a] - upml->m_StartPos[a];
		return upml;
	}
	return NULL;
}

// The factors g = L/A of a cell, per axis: for a voltage the primary edge
// length along n over the dual area across it, for a current the dual edge
// length over the primary area. Any factorisation is exact (m = vi/g takes
// the rest); this one makes m depend on the material only, on a Cartesian mesh.
void CellModel::CalcGeometry()
{
	for (int a=0; a<3; ++a)
	{
		const unsigned int N = numLines[a];
		std::vector<double> prim(N), dual(N);
		for (unsigned int i=0; i+1<N; ++i)
			prim[i] = Op->GetDiscLine(a, i+1) - Op->GetDiscLine(a, i);
		prim[N-1] = N>1 ? prim[N-2] : 1;
		for (unsigned int i=0; i<N; ++i)
		{
			double lo = (i>0) ? prim[i-1] : prim[0];
			double hi = (i+1<N) ? prim[i] : prim[N-1];
			dual[i] = 0.5*(lo + hi);
		}
		for (int n=0; n<3; ++n)
		{
			std::vector<double>& gv = m_G[0][n][a];
			std::vector<double>& gi = m_G[1][n][a];
			gv.resize(N);
			gi.resize(N);
			for (unsigned int i=0; i<N; ++i)
			{
				gv[i] = (a==n) ? prim[i] : 1.0/dual[i];
				gi[i] = (a==n) ? dual[i] : 1.0/prim[i];
			}
		}
	}
}

double CellModel::CellKey::M() const
{
	if (cls==CellModel_Set::ZERO || sign==0)
		return 0;
	return sign*DequantiseM(q);
}

uint32_t CellModel::CellKey::Hash() const
{
	// FNV-1a over the fields; a collision only merges two runs in the statistics
	uint32_t h = 2166136261u;
	auto mix = [&h](const void* p, size_t n) {
		for (size_t i=0; i<n; ++i)
			h = (h ^ ((const unsigned char*)p)[i]) * 16777619u;
	};
	mix(&cls, sizeof(cls));
	mix(&vv, sizeof(vv));
	mix(&q, sizeof(q));
	mix(&sign, sizeof(sign));
	mix(&uvv, sizeof(uvv));
	mix(&ufo, sizeof(ufo));
	mix(&ufn, sizeof(ufn));
	mix(&nl, sizeof(nl));
	mix(&nd, sizeof(nd));
	mix(&solve, sizeof(solve));
	mix(poles, NumPoles()*sizeof(CellModel_Pole));
	return h;
}

CellModel::CellKey CellModel::GetCellKey(int f, int n, const unsigned int pos[3]) const
{
	CellKey key;
	unsigned int p[3] = {pos[0], pos[1], pos[2]};
	const float s = (f==0) ? Op->GetVV(n, p) : Op->GetII(n, p);
	const float c = (f==0) ? Op->GetVI(n, p) : Op->GetIV(n, p);
	key.vv = s;
	key.q = 0;
	key.sign = 0;
	key.uvv = key.ufo = key.ufn = 0;
	key.nl = key.nd = 0;
	key.solve = 1;
	std::memset(key.poles, 0, sizeof(key.poles));
	unsigned int loc[3];
	const Operator_Ext_UPML* upml = FindUPML(pos, loc);
	if (upml)
	{
		// the main coefficients (rewritten by the UPML for its flux) and the UPML ones, graded with the depth
		Operator_Ext_UPML* w = const_cast<Operator_Ext_UPML*>(upml);   // its Get*() are not const
		key.cls = CellModel_Set::UPML;
		key.uvv = (f==0) ? w->GetVV(n, loc)   : w->GetII(n, loc);
		key.ufo = (f==0) ? w->GetVVFO(n, loc) : w->GetIIFO(n, loc);
		key.ufn = (f==0) ? w->GetVVFN(n, loc) : w->GetIIFN(n, loc);
	}
	else if ((s==0) && (c==0))
	{
		key.cls = CellModel_Set::ZERO;
		return key;
	}
	else
	{
		key.cls = (s==1.0f) ? CellModel_Set::LOSSLESS : CellModel_Set::LOSSY;
		if (!m_Disp.empty())
			GetCellPoles(f, n, pos, key);
		if (key.NumPoles()>0)
			key.cls = CellModel_Set::DISP;
	}
	const double g = m_G[f][n][0][pos[0]] * m_G[f][n][1][pos[1]] * m_G[f][n][2][pos[2]];
	const double m = c/g;
	if (m!=0)
	{
		key.sign = (m>0) ? 1 : -1;
		key.q = QuantiseM(std::fabs(m));
	}
	return key;
}

int CellModel::RecommendFastAxis() const
{
	// the number of runs along each axis, over all components of both fields
	double runs[3] = {0, 0, 0};
	const size_t N0 = numLines[0], N1 = numLines[1], N2 = numLines[2];
	std::vector<uint32_t> keys(N0*N1*N2);   // 4 bytes per node, the peak of the setup
	// the threads fill and count slices of x lines
	const unsigned int numThreads = std::max(1u, std::min(AvailableThreads(), numLines[0]));
	std::vector<unsigned int> jpt = AssignJobs2Threads(numLines[0], numThreads, true);
	std::vector<unsigned int> x0(jpt.size()+1, 0);
	for (size_t t=0; t<jpt.size(); ++t)
		x0[t+1] = x0[t] + jpt[t];
	auto parallel = [&](const std::function<void(unsigned int, unsigned int, size_t)>& work) {
		std::vector<std::thread> threads;
		for (size_t t=0; t<jpt.size(); ++t)
			threads.push_back(std::thread(work, x0[t], x0[t+1], t));
		for (std::thread& th : threads)
			th.join();
	};
	std::vector<double> part(3*jpt.size());
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
		{
			parallel([&](unsigned int xa, unsigned int xb, size_t) {
				unsigned int pos[3];
				for (pos[0]=xa; pos[0]<xb; ++pos[0])
					for (pos[1]=0; pos[1]<N1; ++pos[1])
						for (pos[2]=0; pos[2]<N2; ++pos[2])
							keys[(pos[0]*N1 + pos[1])*N2 + pos[2]] = GetCellKey(f, n, pos).Hash();
			});
			// a run starts at the first line of the axis or where the key changes
			std::fill(part.begin(), part.end(), 0.0);
			parallel([&](unsigned int xa, unsigned int xb, size_t t) {
				for (size_t x=xa; x<xb; ++x)
					for (size_t y=0; y<N1; ++y)
					{
						const uint32_t* k = &keys[(x*N1 + y)*N2];
						for (size_t z=0; z<N2; ++z)
						{
							part[3*t+0] += (x==0) || (k[z]!=k[z-N1*N2]);
							part[3*t+1] += (y==0) || (k[z]!=k[z-N2]);
							part[3*t+2] += (z==0) || (k[z]!=k[z-1]);
						}
					}
			});
			for (size_t t=0; t<jpt.size(); ++t)
				for (int a=0; a<3; ++a)
					runs[a] += part[3*t+a];
		}

	int fast = 2;
	for (int a=0; a<2; ++a)
		if ((runs[a]<runs[fast]) || ((runs[a]==runs[fast]) && (numLines[a]>numLines[fast])))
			fast = a;
	if (g_settings.GetVerboseLevel()>0)
	{
		const double cells = 6.0*N0*N1*N2;
		cout << "CellModel: mean run length along x, y, z: " << cells/runs[0] << ", " << cells/runs[1]
		     << ", " << cells/runs[2] << "; the fewest runs along " << (char)('x'+fast) << endl;
	}
	return fast;
}

void CellModel::Build(int F, int S, int M)
{
	// one thread per component, each with its own cell counts
	size_t counts[2][3][6] = {};
	auto build = [&](int f, int n)
		{
			size_t* count = counts[f][n];
			CellModel_Component& r = m_Comp[f][n];
			std::vector<CellModel_Pole>& poles = r.poles;
			std::vector<float>* profile = r.profile;
			std::map<std::vector<float>, uint32_t> profileIndex;   // the start of each distinct profile
			for (int q=0; q<5; ++q)
				profile[q].clear();
			r.runs.clear();
			r.sets.clear();
			poles.clear();
			// the index of each distinct key in r.sets
			typedef std::tuple<int, float, int64_t, int, float, float, float, std::vector<float> > SetKey;
			std::map<SetKey, uint32_t> setIndex;
			r.first.assign((size_t)numLines[S]*numLines[M] + 1, 0);
			// the currents are updated on [0, numLines-1) per axis
			const unsigned int last = (f==0) ? 0 : 1;
			unsigned int pos[3];
			for (pos[S]=0; pos[S]<numLines[S]; ++pos[S])
				for (pos[M]=0; pos[M]<numLines[M]; ++pos[M])
				{
					const size_t line = (size_t)pos[S]*numLines[M] + pos[M];
					r.first[line] = r.runs.size();
					if ((pos[S]+last>=numLines[S]) || (pos[M]+last>=numLines[M]))
						continue;
					CellKey prev = CellKey();
					bool open = false;   // the last run ends at pos[F]
					for (pos[F]=0; pos[F]+last<numLines[F]; ++pos[F])
					{
						CellKey key = GetCellKey(f, n, pos);
						if ((f==1) && m_SkipInert && (key.cls==CellModel_Set::LOSSLESS) && IsInertCurr(n, pos))
						{
							++count[5];   // in no run
							open = false;
							continue;
						}
						++count[key.cls];
						if (open && (key==prev))
						{
							++r.runs.back().stop;
							continue;
						}
						std::vector<float> pk;
						pk.push_back(key.nl);
						pk.push_back(key.solve);
						for (int p=0; p<key.NumPoles(); ++p)
						{
							pk.push_back(key.poles[p].k_int);
							pk.push_back(key.poles[p].k_ext);
							pk.push_back(key.poles[p].k_lor);
						}
						const SetKey sk(key.cls, key.vv, key.q, key.sign, key.uvv, key.ufo, key.ufn, pk);
						auto it = setIndex.find(sk);
						if (it==setIndex.end())
						{
							CellModel_Set set;
							set.cls = key.cls;
							set.m = key.M();
							set.vv = (key.cls==CellModel_Set::LOSSLESS) ? 1.0f : key.vv;
							set.uvv = key.uvv;
							set.ufo = key.ufo;
							set.ufn = key.ufn;
							set.offset = (uint32_t)poles.size();
							set.nl = key.nl;
							set.nd = key.nd;
							set.solve = key.solve;
							poles.insert(poles.end(), key.poles, key.poles+key.NumPoles());
							it = setIndex.insert(std::make_pair(sk, (uint32_t)r.sets.size())).first;
							r.sets.push_back(set);
						}
						CellModel_Run run;
						run.start = pos[F];
						run.stop = pos[F]+1;
						run.set = it->second;
						run.u = 0;
						r.runs.push_back(run);
						prev = key;
						open = true;
					}
					// neighbouring short UPML runs: one run with a profile (a PML normal to the fast axis)
					std::vector<CellModel_Run> merged;
					const size_t first = r.first[line];
					for (size_t i=first; i<r.runs.size(); )
					{
						size_t j = i+1;
						if (r.sets[r.runs[i].set].cls==CellModel_Set::UPML)
							while ((j<r.runs.size()) && (r.sets[r.runs[j].set].cls==CellModel_Set::UPML) && (r.runs[j].start==r.runs[j-1].stop))
								++j;
						const unsigned int len = r.runs[j-1].stop - r.runs[i].start;
						if ((j-i<2) || (len >= 4*(j-i)))
						{
							merged.insert(merged.end(), r.runs.begin()+i, r.runs.begin()+j);
							i = j;
							continue;
						}
						std::vector<float> values(5*len);
						for (size_t q=i; q<j; ++q)
						{
							const CellModel_Set& s = r.sets[r.runs[q].set];
							for (unsigned int k=r.runs[q].start; k<r.runs[q].stop; ++k)
							{
								const unsigned int c = k - r.runs[i].start;
								values[0*len+c] = (float)s.m;
								values[1*len+c] = s.vv;
								values[2*len+c] = s.uvv;
								values[3*len+c] = s.ufo;
								values[4*len+c] = s.ufn;
							}
						}
						auto pit = profileIndex.find(values);
						if (pit==profileIndex.end())
						{
							if (profile[0].size()+len > 0xffffffffu)
								throw std::length_error("CellModel: more than 2^32 cells of UPML profiles in a component");
							pit = profileIndex.insert(std::make_pair(values, (uint32_t)profile[0].size())).first;
							for (int q=0; q<5; ++q)
								profile[q].insert(profile[q].end(), values.begin()+q*len, values.begin()+(q+1)*len);
						}
						const SetKey sk(CellModel_Set::UPML_PROFILE, 0, pit->second, 0, 0, 0, 0, std::vector<float>());
						auto it = setIndex.find(sk);
						if (it==setIndex.end())
						{
							CellModel_Set set = CellModel_Set();
							set.cls = CellModel_Set::UPML_PROFILE;
							set.m = 1;   // the kernel's a is then the slow factor, the profile has m
							set.solve = 1;
							set.offset = pit->second;
							it = setIndex.insert(std::make_pair(sk, (uint32_t)r.sets.size())).first;
							r.sets.push_back(set);
						}
						CellModel_Run run;
						run.start = r.runs[i].start;
						run.stop = r.runs[j-1].stop;
						run.set = it->second;
						run.u = 0;
						merged.push_back(run);
						i = j;
					}
					r.runs.resize(first);
					r.runs.insert(r.runs.end(), merged.begin(), merged.end());
				}
			r.first.back() = r.runs.size();
			// the state of the UPML and DISP runs, now that their lengths are known
			size_t state = 0;
			for (CellModel_Run& run : r.runs)
			{
				const CellModel_Set& set = r.sets[run.set];
				const size_t len = run.stop - run.start;
				if ((set.cls!=CellModel_Set::UPML) && (set.cls!=CellModel_Set::UPML_PROFILE) && (set.cls!=CellModel_Set::DISP))
					continue;
				if (state > 0xffffffffu)
					throw std::length_error("CellModel: more than 2^32 floats of UPML and ADE state in a component");
				run.u = (uint32_t)state;
				if ((set.cls==CellModel_Set::UPML) || (set.cls==CellModel_Set::UPML_PROFILE))
					state += len;
				else
				{
					for (unsigned int p=0; p<set.nl; ++p)
						state += ((poles[set.offset+p].k_lor!=0) ? 2 : 1)*len;
					state += set.nd*len;
				}
			}
			r.state = state;
		};
	// an exception in a thread would end the program: carried out and rethrown
	std::exception_ptr error[2][3];
	std::vector<std::thread> threads;
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
			threads.push_back(std::thread([&, f, n]() {
				try {build(f, n);}
				catch (...) {error[f][n] = std::current_exception();}
			}));
	for (std::thread& th : threads)
		th.join();
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
			if (error[f][n])
				std::rethrow_exception(error[f][n]);
	size_t count[6] = {0, 0, 0, 0, 0, 0};
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
			for (int c=0; c<6; ++c)
				count[c] += counts[f][n][c];

	if (g_settings.GetVerboseLevel()>0)
	{
		size_t runs = 0, sets = 0;
		for (int f=0; f<2; ++f)
			for (int n=0; n<3; ++n)
			{
				runs += m_Comp[f][n].runs.size();
				sets += m_Comp[f][n].sets.size();
			}
		const double cells = count[0] + count[1] + count[2] + count[3] + count[4] + count[5];
		cout << "CellModel: " << runs << " runs (" << runs*sizeof(CellModel_Run)/1048576.0 << " MiB), " << sets
		     << " coefficient sets, mean length " << cells/runs << "; cells zero "
		     << 100.0*count[0]/cells << " %, lossless " << 100.0*count[1]/cells << " %, lossy "
		     << 100.0*count[2]/cells << " %, upml " << 100.0*count[3]/cells << " %, dispersive "
		     << 100.0*count[4]/cells << " %, inert (skipped) " << 100.0*count[5]/cells << " %" << endl;
	}
}

