/*
*	Copyright (C) 2026 Roel Budé (roel.bude@proton.me)
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

#include "operator_ext_lossymetal.h"
#include "tools/arraylib/array_nijk.h"
#include "tools/constants.h"

#include "CSPropLossyMetal.h"

#include <array>
#include <cmath>

using std::cerr;
using std::cout;
using std::endl;

// number of R-L branches of the rational fit, see SurfaceImpedance::FitFace
#define LOSSYMETAL_BRANCHES 8

Operator_Ext_LossyMetal::Operator_Ext_LossyMetal(Operator* op, double f_max) : Operator_Ext_LorentzMaterial(op)
{
	m_f_max = f_max;
	m_NumBranches = LOSSYMETAL_BRANCHES;
	m_numSurface = m_numFinite = m_numPEC = m_numH = m_numNoThickness = 0;
}

Operator_Ext_LossyMetal::Operator_Ext_LossyMetal(Operator* op, Operator_Ext_LossyMetal* op_ext) : Operator_Ext_LorentzMaterial(op, op_ext)
{
	m_f_max = op_ext->m_f_max;
	m_NumBranches = op_ext->m_NumBranches;
	m_numSurface = m_numFinite = m_numPEC = m_numH = m_numNoThickness = 0;
}

Operator_Extension* Operator_Ext_LossyMetal::Clone(Operator* op)
{
	return new Operator_Ext_LossyMetal(op, this);
}

bool Operator_Ext_LossyMetal::IsLossy(const double* coord) const
{
	CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord, (CSProperties::PropertyType)(CSProperties::MATERIAL | CSProperties::METAL), false);
	return (prop!=NULL) && (prop->GetType() & CSProperties::LOSSYMETAL);
}

bool Operator_Ext_LossyMetal::IsMetal(const double* coord) const
{
	CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord, (CSProperties::PropertyType)(CSProperties::MATERIAL | CSProperties::METAL), false);
	return (prop!=NULL) && (prop->GetType() & CSProperties::METAL);
}

double Operator_Ext_LossyMetal::InsideLength(const double* coord, int dir, int s, double max) const
{
	double p[3] = {coord[0], coord[1], coord[2]};
	p[dir] = coord[dir] + s*max;
	if (IsLossy(p))
		return max;
	double lo = 0, hi = max;
	for (int i=0; i<30; ++i)
	{
		double mid = (lo+hi)/2;
		p[dir] = coord[dir] + s*mid;
		if (IsLossy(p))
			lo = mid;
		else
			hi = mid;
	}
	return lo;
}

double Operator_Ext_LossyMetal::Extent(const double* coord, const unsigned int* pos, int dir, int s) const
{
	int numLines = m_Op->GetNumberOfLines(dir,true);
	int idx = pos[dir];
	double p[3] = {coord[0], coord[1], coord[2]};
	while (true)
	{
		int next = idx + s;
		if ((next<0) || (next>=numLines))
			return std::fabs(m_Op->GetDiscLine(dir,idx) - coord[dir]);
		p[dir] = m_Op->GetDiscLine(dir,next);
		if (!IsLossy(p))
		{
			double start[3] = {coord[0], coord[1], coord[2]};
			start[dir] = m_Op->GetDiscLine(dir,idx);
			double len = std::fabs(m_Op->GetDiscLine(dir,next) - start[dir]);
			return std::fabs(start[dir] - coord[dir]) + InsideLength(start, dir, s, len);
		}
		idx = next;
	}
}

int Operator_Ext_LossyMetal::GetFit(CSPropLossyMetal* prop, double thickness)
{
	double sigma = prop->GetConductivity();
	double rq = (prop->GetRoughnessModel()==CSPropLossyMetal::ROUGHNESS_GRADIENT) ? prop->GetRq() : 0;
	if (rq<0)
		rq = 0;
	// beyond this thickness the face impedance does not depend on it in the fitted band
	double delta_lo = sqrt(2/(2*PI*m_f_max/1000*MUE0*sigma));
	FitKey key;
	key.prop = prop;
	key.t_index = -1;
	double t = INFINITY;
	if (thickness/2 < 6*rq + 15*delta_lo)
	{
		key.t_index = lround(log(thickness/1e-9)/log(1.005)); // 0.5% steps
		t = 1e-9*pow(1.005, key.t_index);
	}
	std::map<FitKey,int>::iterator it = m_FitIndex.find(key);
	if (it!=m_FitIndex.end())
		return it->second;

	SurfaceImpedance::Fit fit = SurfaceImpedance::FitFace(sigma, t, rq, m_f_max, m_w_k);
	int idx = m_Fits.size();
	m_Fits.push_back(fit);
	m_FitThickness.push_back(t);
	m_FitIndex[key] = idx;
	return idx;
}

bool Operator_Ext_LossyMetal::AddContributions(int n, const unsigned int* pos, const double* coord, CSPropLossyMetal* prop, std::vector<Contribution>& contrib)
{
	int dirs[2] = {(n+1)%3, (n+2)%3};
	double len = m_Op->GetEdgeLength(n,pos)/m_Op->GetGridDelta(); // drawing units, as all lengths here
	if (len<=0)
		return false;
	double delta = m_Op->GetGridDelta();
	double tol = m_Op->GetGeometryCSX()->GetDrawingTolerance();

	// half cell sizes towards the neighbouring cell centres: h[d][0] backwards, h[d][1] forwards,
	// zero on a (PMC) boundary plane
	double h[3][2];
	// probe offsets towards the quarter cells, on a boundary plane mirrored to the inside
	double off[3][2];
	for (int m=0; m<2; ++m)
	{
		int d = dirs[m];
		unsigned int numLines = m_Op->GetNumberOfLines(d,true);
		h[d][0] = (pos[d]>0) ? (m_Op->GetDiscLine(d,pos[d]) - m_Op->GetDiscLine(d,pos[d]-1))/2 : 0;
		h[d][1] = (pos[d]+1<numLines) ? (m_Op->GetDiscLine(d,pos[d]+1) - m_Op->GetDiscLine(d,pos[d]))/2 : 0;
		off[d][0] = (h[d][0]>0) ? -h[d][0] : h[d][1];
		off[d][1] = (h[d][1]>0) ? h[d][1] : -h[d][0];
	}

	// metal in the four quarter cells, metal[s0][s1] with s=0: backwards, s=1: forwards along dirs[0], dirs[1]
	bool metal[2][2];
	double p[3];
	for (int s0=0; s0<2; ++s0)
		for (int s1=0; s1<2; ++s1)
		{
			p[0]=coord[0]; p[1]=coord[1]; p[2]=coord[2];
			p[dirs[0]] += off[dirs[0]][s0];
			p[dirs[1]] += off[dirs[1]][s1];
			metal[s0][s1] = IsMetal(p);
		}

	// a body that ends before the neighbouring cell centres along dirs[m] is thinner than the mesh
	bool thin[2];
	double thick[2] = {-1, -1};
	for (int m=0; m<2; ++m)
	{
		int d = dirs[m];
		p[0]=coord[0]; p[1]=coord[1]; p[2]=coord[2];
		p[d] = coord[d] + off[d][0];
		thin[m] = !IsLossy(p);
		p[d] = coord[d] + off[d][1];
		thin[m] = thin[m] && !IsLossy(p);
	}

	bool missingThickness = false;
	// walk the four half segments of the dual face boundary, from the edge to the neighbouring cell centres
	for (int m=0; m<2; ++m)
	{
		int d = dirs[m];   // segment direction
		int e = dirs[1-m]; // the other direction, normal of a face lying along the segment
		for (int s=0; s<2; ++s)
		{
			double Lin = InsideLength(coord, d, 2*s-1, h[d][s]);
			if (Lin<=0)
				continue;
			if (thin[1-m])
			{
				// a sheet: both faces open, thickness of the body along e
				// the drawing tolerance makes a body appear larger by tol on either side
				if (thick[1-m]<0)
					thick[1-m] = std::max(0.0, InsideLength(coord, e, -1, h[e][0]) + InsideLength(coord, e, +1, h[e][1]) - 2*tol)*delta;
				double t = thick[1-m];
				if (t <= tol*delta)
					t = prop->GetThickness(); // a 2D primitive
				if (t<=0)
				{
					missingThickness = true;
					continue;
				}
				Contribution c;
				c.fit = GetFit(prop, t);
				c.weight = 2*Lin/len;
				contrib.push_back(c);
				continue;
			}
			// the two quarter cells on either side of this segment
			bool q0 = (m==0) ? metal[s][0] : metal[0][s];
			bool q1 = (m==0) ? metal[s][1] : metal[1][s];
			if (q0==q1)
				continue; // inside the metal or in the open, no face here
			// a face with normal e, the metal is on the side of the metal quarter cell
			double t = std::max(Extent(coord, pos, e, q0 ? -1 : +1) - tol, tol)*delta;
			Contribution c;
			c.fit = GetFit(prop, t);
			c.weight = Lin/len;
			contrib.push_back(c);
		}
	}
	if (missingThickness && contrib.empty())
		++m_numNoThickness;
	return !contrib.empty();
}

bool Operator_Ext_LossyMetal::BuildExtension()
{
	double dT = m_Op->GetTimestep();
	unsigned int pos[] = {0,0,0};
	double coord[3];
	unsigned int numLines[3] = {m_Op->GetNumberOfLines(0,true),m_Op->GetNumberOfLines(1,true),m_Op->GetNumberOfLines(2,true)};

	if (m_f_max<=0)
	{
		cerr << "Operator_Ext_LossyMetal::BuildExtension: Error, no max. frequency, lossy metal is PEC!" << endl;
		m_Order = 0;
		return false;
	}
	m_w_k = SurfaceImpedance::CornerFrequencies(m_f_max, m_NumBranches);
	m_Fits.clear();
	m_FitThickness.clear();
	m_FitIndex.clear();
	m_numSurface = m_numFinite = m_numPEC = m_numH = m_numNoThickness = 0;

	m_Order = 0;
	std::vector<unsigned int> v_pos[3];
	std::vector< std::array<std::vector<Contribution>,3> > v_contrib;
	// 0: no lossy metal, 1: inside a lossy metal body (PEC), 2: surface edge
	ArrayLib::ArrayNIJK<int8_t> state("state", numLines);

	CSPrimitives* prim = NULL;
	bool disable_pos;
	for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
	{
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			std::vector<CSPrimitives*> vPrims = m_Op->GetPrimitivesBoundBox(
				pos[0], pos[1], -1,
				(CSProperties::PropertyType)(CSProperties::MATERIAL | CSProperties::METAL)
			);

			for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
			{
				// disable the surface impedance inside a pml
				disable_pos = false;
				for (int m=0;m<3;++m)
				{
					int lo = m_Op->GetBCSize(2*m), hi = m_Op->GetBCSize(2*m+1);
					if (((lo>0) && (pos[m]<=(unsigned int)lo)) || ((hi>0) && (pos[m]>=numLines[m]-hi-1)))
						disable_pos = true;
				}

				std::array<std::vector<Contribution>,3> contrib;
				bool b_pos_on = false;
				for (int n=0; n<3; ++n)
				{
					// and on the boundary planes, unless it is a PMC
					bool disable_edge = disable_pos;
					for (int m=0;m<3;++m)
					{
						if (m==n)
							continue;
						if (((pos[m]==0) && (m_Op->m_BC[2*m]!=1)) || ((pos[m]==numLines[m]-1) && (m_Op->m_BC[2*m+1]!=1)))
							disable_edge = true;
					}
					state(n, pos[0], pos[1], pos[2]) = 0;
					if (m_Op->GetYeeCoords(n,pos,coord,false)==false)
						continue;

					CSProperties* prop = m_Op->GetGeometryCSX()->GetPropertyByCoordPriority(coord, vPrims, false, &prim);
					CSPropLossyMetal* lm_prop = dynamic_cast<CSPropLossyMetal*>(prop);
					if (lm_prop==NULL)
						continue;
					if (prim)
						prim->SetPrimitiveUsed(true);

					state(n, pos[0], pos[1], pos[2]) = 1;
					if (!disable_edge && (lm_prop->GetConductivity()>0) && AddContributions(n, pos, coord, lm_prop, contrib[n]))
					{
						state(n, pos[0], pos[1], pos[2]) = 2;
						b_pos_on = true;
						continue;
					}
					// inside the body, in a boundary condition or no conductivity: PEC
					m_Op->SetVV(n,pos[0],pos[1],pos[2], 0 );
					m_Op->SetVI(n,pos[0],pos[1],pos[2], 0 );
					++m_Op->m_Nr_PEC[n];
					++m_numPEC;
				}
				if (b_pos_on)
				{
					for (int n=0; n<3; ++n)
						v_pos[n].push_back(pos[n]);
					v_contrib.push_back(contrib);
				}
			}
		}
	}

	if (m_numNoThickness>0)
		cerr << "Operator_Ext_LossyMetal::BuildExtension: Warning: " << m_numNoThickness << " edges of zero-thickness LossyMetal primitives without a Thickness, fallback to PEC!" << endl;

	size_t numLM = v_pos[0].size();
	if (numLM==0)
		return false;

	m_Order = m_NumBranches;
	m_volt_ADE_On = new bool[m_Order];
	m_curr_ADE_On = new bool[m_Order];
	m_volt_Lor_ADE_On = new bool[m_Order];
	m_curr_Lor_ADE_On = new bool[m_Order];
	m_LM_pos = new unsigned int**[m_Order];
	v_int_ADE = new FDTD_FLOAT**[m_Order];
	v_ext_ADE = new FDTD_FLOAT**[m_Order];
	for (int o=0; o<m_Order; ++o)
	{
		m_LM_Count.push_back(numLM);
		m_volt_ADE_On[o] = true;
		m_curr_ADE_On[o] = false;
		m_volt_Lor_ADE_On[o] = false;
		m_curr_Lor_ADE_On[o] = false;
		m_LM_pos[o] = new unsigned int*[3];
		v_int_ADE[o] = new FDTD_FLOAT*[3];
		v_ext_ADE[o] = new FDTD_FLOAT*[3];
		for (int n=0; n<3; ++n)
		{
			m_LM_pos[o][n] = new unsigned int[numLM];
			for (size_t i=0; i<numLM; ++i)
				m_LM_pos[o][n][i] = v_pos[n].at(i);
			v_int_ADE[o][n] = new FDTD_FLOAT[numLM];
			v_ext_ADE[o][n] = new FDTD_FLOAT[numLM];
		}
	}

	std::vector<double> a(m_Order);
	for (size_t i=0; i<numLM; ++i)
	{
		pos[0]=v_pos[0].at(i); pos[1]=v_pos[1].at(i); pos[2]=v_pos[2].at(i);
		unsigned int index = m_Op->MainOp->SetPos(pos[0],pos[1],pos[2]);
		for (int n=0; n<3; ++n)
		{
			for (int o=0; o<m_Order; ++o)
			{
				v_int_ADE[o][n][i] = 0;
				v_ext_ADE[o][n][i] = 0;
			}
			const std::vector<Contribution>& contrib = v_contrib.at(i)[n];
			if (contrib.empty())
				continue;
			++m_numSurface;

			// edge admittance: sum of the face admittances times width/length, all fits share the corner frequencies
			double G = 0;
			bool finite = false;
			for (int o=0; o<m_Order; ++o)
				a[o] = 0;
			for (size_t c=0; c<contrib.size(); ++c)
			{
				const SurfaceImpedance::Fit& fit = m_Fits.at(contrib[c].fit);
				G += contrib[c].weight*fit.G;
				for (int o=0; o<m_Order; ++o)
					a[o] += contrib[c].weight*fit.a[o];
				finite = finite || (m_FitThickness.at(contrib[c].fit)<INFINITY);
			}
			if (finite)
				++m_numFinite;

			// branch o: R = 1/a, L = 1/(a w_k), stabilized by an additional capacitance as in Operator_Ext_ConductingSheet
			double sumInvL = 0, maxInvL = 0;
			for (int o=0; o<m_Order; ++o)
			{
				sumInvL += a[o]*m_w_k[o];
				maxInvL = std::max(maxInvL, a[o]*m_w_k[o]);
			}
			m_Op->EC_C[n][index] += dT*dT/4.0*(16.0*maxInvL + sumInvL);
			m_Op->EC_G[n][index] += G;
			m_Op->Calc_ECOperatorPos(n,pos);

			for (int o=0; o<m_Order; ++o)
			{
				double wdt = m_w_k[o]*dT;
				v_int_ADE[o][n][i] = (2.0-wdt)/(2.0+wdt);
				v_ext_ADE[o][n][i] = a[o]*wdt/(1.0+wdt/2.0)*m_Op->GetVI(n,pos[0],pos[1],pos[2]);
			}
		}
	}

	// the tangential H inside the body, between surface and PEC edges, would add the inductance of
	// a half cell in parallel to the surface admittance: switch it off
	for (pos[0]=0; pos[0]<numLines[0]-1; ++pos[0])
		for (pos[1]=0; pos[1]<numLines[1]-1; ++pos[1])
			for (pos[2]=0; pos[2]<numLines[2]-1; ++pos[2])
				for (int n=0; n<3; ++n)
				{
					int nP = (n+1)%3;
					int nPP = (n+2)%3;
					unsigned int pP[3] = {pos[0], pos[1], pos[2]};
					unsigned int pPP[3] = {pos[0], pos[1], pos[2]};
					++pP[nP];
					++pPP[nPP];
					int8_t s[4] = {state(nPP,pos[0],pos[1],pos[2]), state(nPP,pP[0],pP[1],pP[2]),
								   state(nP,pos[0],pos[1],pos[2]), state(nP,pPP[0],pPP[1],pPP[2])};
					if ((s[0]==0) || (s[1]==0) || (s[2]==0) || (s[3]==0))
						continue;
					if ((s[0]==2) && (s[1]==2) && (s[2]==2) && (s[3]==2))
						continue;
					m_Op->SetII(n,pos[0],pos[1],pos[2], 0 );
					m_Op->SetIV(n,pos[0],pos[1],pos[2], 0 );
					++m_numH;
				}

	cout << "Operator_Ext_LossyMetal: " << m_Fits.size() << " surface impedance fit(s), " << m_Order << " branches, f_max " << m_f_max/1e9 << " GHz:" << endl;
	for (size_t f=0; f<m_Fits.size(); ++f)
	{
		const SurfaceImpedance::Fit& fit = m_Fits.at(f);
		cout << "  face thickness ";
		if (m_FitThickness.at(f)<INFINITY)
			cout << m_FitThickness.at(f)*1e6 << " um";
		else
			cout << "inf";
		cout << ": max. error Re(Zs) " << fit.err_re*100 << "%, Im(Zs) " << fit.err_im*100 << "%";
		if (fit.d_ref>0)
			cout << ", referenced " << fit.d_ref*1e6 << " um outside the mean surface";
		cout << endl;
	}
	return true;
}

void Operator_Ext_LossyMetal::ShowStat(std::ostream &ostr) const
{
	Operator_Extension::ShowStat(ostr);
	ostr << " Surface impedance branches\t: " << m_Order << endl;
	ostr << " Surface edges\t\t\t: " << m_numSurface << " (" << m_numFinite << " of finite thickness)" << endl;
	ostr << " PEC edges inside\t\t: " << m_numPEC << endl;
	ostr << " Disabled H edges inside\t: " << m_numH << endl;
}
