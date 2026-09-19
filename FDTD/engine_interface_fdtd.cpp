/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "engine_interface_fdtd.h"
#include "engine_gpu.h"
#include <stdexcept>

using std::cerr;
using std::endl;

Engine_Interface_FDTD::Engine_Interface_FDTD(Operator* op) : Engine_Interface_Base(op)
{
	if (op==NULL)
		throw std::runtime_error("Engine_Interface_FDTD::Engine_Interface_FDTD: Error: Operator is not set!");
	m_Op = op;
	m_Eng = m_Op->GetEngine();
	if (m_Eng==NULL)
		throw std::runtime_error("Engine_Interface_FDTD::Engine_Interface_FDTD: Error: Engine is not set!");
}

Engine_Interface_FDTD::~Engine_Interface_FDTD()
{
}

double* Engine_Interface_FDTD::GetEField(const unsigned int* pos, double* out) const
{
	return GetRawInterpolatedField(pos, out, 0);
}

double* Engine_Interface_FDTD::GetJField(const unsigned int* pos, double* out) const
{
	return GetRawInterpolatedField(pos, out, 1);
}

double* Engine_Interface_FDTD::GetDField(const unsigned int* pos, double* out) const
{
	return GetRawInterpolatedField(pos, out, 3);
}

double* Engine_Interface_FDTD::GetRotHField(const unsigned int* pos, double* out) const
{
	return GetRawInterpolatedField(pos, out, 2);
}

double* Engine_Interface_FDTD::GetRawInterpolatedField(const unsigned int* pos, double* out, int type) const
{
	unsigned int iPos[] = {pos[0],pos[1],pos[2]};
	int nP,nPP;
	double delta;
	switch (m_InterpolType)
	{
	default:
	case NO_INTERPOLATION:
		for (int n=0; n<3; ++n)
			out[n] = GetRawField(n,pos,type);
		break;
	case NODE_INTERPOLATE:
		for (int n=0; n<3; ++n)
		{
			if (pos[n]==m_Op->GetNumberOfLines(n, true)-1)  // use only the "lower value" at the upper bound
			{
				--iPos[n];
				out[n] = (double)GetRawField(n,iPos,type);
				++iPos[n];
				continue;
			}
			delta = m_Op->GetEdgeLength(n,iPos);
			out[n] = GetRawField(n,iPos,type);
			if (delta==0)
			{
				out[n]=0;
				continue;
			}
			if (pos[n]==0) // use only the "upper value" at the lower bound
				continue;
			--iPos[n];
			double deltaDown = m_Op->GetEdgeLength(n,iPos);
			double deltaRel = delta / (delta+deltaDown);
			out[n] = out[n]*(1.0-deltaRel) + (double)GetRawField(n,iPos,type)*deltaRel;
			++iPos[n];
		}
		break;
	case CELL_INTERPOLATE:
		for (int n=0; n<3; ++n)
		{
			nP = (n+1)%3;
			nPP = (n+2)%3;
			if ((pos[0]==m_Op->GetNumberOfLines(0,true)-1) || (pos[1]==m_Op->GetNumberOfLines(1,true)-1) || (pos[2]==m_Op->GetNumberOfLines(2,true)-1))
			{
				out[n] = 0; //electric field outside the field domain is always zero
				continue;
			}
			out[n]=GetRawField(n,iPos,type);
			++iPos[nP];
			out[n]+=GetRawField(n,iPos,type);
			++iPos[nPP];
			out[n]+=GetRawField(n,iPos,type);
			--iPos[nP];
			out[n]+=GetRawField(n,iPos,type);
			--iPos[nPP];
			out[n]/=4;
		}
		break;
	}
	return out;
}

double* Engine_Interface_FDTD::GetHField(const unsigned int* pos, double* out) const
{
		return GetRawInterpolatedDualField(pos, out, 0);
}

double* Engine_Interface_FDTD::GetBField(const unsigned int* pos, double* out) const
{
		return GetRawInterpolatedDualField(pos, out, 1);
}

double Engine_Interface_FDTD::GetRawDualField(unsigned int n, const unsigned int* pos, int type) const
{
	double value = m_Eng->GetCurr(n,pos[0],pos[1],pos[2]);
	double delta = m_Op->GetEdgeLength(n,pos,true);
	if ((type==0) && (delta))
		return value/delta;
	if ((type==1) && (m_Op->m_mueR_ptr) && (delta))
	{
		ArrayLib::ArrayNIJK<float>& m_mueR = *m_Op->m_mueR_ptr;
		return value * m_mueR(n, pos[0], pos[1], pos[2]) / delta;
	}
	return 0.0;
}

double* Engine_Interface_FDTD::GetRawInterpolatedDualField(const unsigned int* pos, double* out, int type) const
{
	unsigned int iPos[] = {pos[0],pos[1],pos[2]};
	int nP,nPP;
	double delta;
	switch (m_InterpolType)
	{
	default:
	case NO_INTERPOLATION:
		out[0] = GetRawDualField(0, pos, type);
		out[1] = GetRawDualField(1, pos, type);
		out[2] = GetRawDualField(2, pos, type);
		break;
	case NODE_INTERPOLATE:
		for (int n=0; n<3; ++n)
		{
			nP = (n+1)%3;
			nPP = (n+2)%3;
			if ((pos[0]==m_Op->GetNumberOfLines(0,true)-1) || (pos[1]==m_Op->GetNumberOfLines(1,true)-1) || (pos[2]==m_Op->GetNumberOfLines(2,true)-1) || (pos[nP]==0) || (pos[nPP]==0))
			{
				out[n] = 0;
				continue;
			}
			out[n] = GetRawDualField(n, iPos, type);
			--iPos[nP];
			out[n]+= GetRawDualField(n, iPos, type);
			--iPos[nPP];
			out[n]+= GetRawDualField(n, iPos, type);
			++iPos[nP];
			out[n]+= GetRawDualField(n, iPos, type);
			++iPos[nPP];
			out[n]/=4;
		}
		break;
	case CELL_INTERPOLATE:
		for (int n=0; n<3; ++n)
		{
			delta = m_Op->GetEdgeLength(n,iPos,true);
			out[n] = GetRawDualField(n, iPos, type);
			if ((pos[n]>=m_Op->GetNumberOfLines(n,true)-1))
			{
				out[n] = 0; //magnetic field on the outer boundaries is always zero
				continue;
			}
			++iPos[n];
			double deltaUp = m_Op->GetEdgeLength(n,iPos,true);
			double deltaRel = delta / (delta+deltaUp);
			out[n] = out[n]*(1.0-deltaRel) + (double)GetRawDualField(n, iPos, type)*deltaRel;
			--iPos[n];
		}
		break;
	}

	return out;
}

double Engine_Interface_FDTD::CalcVoltageIntegral(const unsigned int* start, const unsigned int* stop) const
{
	if (((start[0]!=stop[0]) + (start[1]!=stop[1]) + (start[2]!=stop[2]))!=1)
	{
		cerr << "Engine_Interface_FDTD::CalcVoltageIntegral: Error, only a 1D/line integration is allowed" << endl;
		return 0;
	}
	//cerr << "CalcVoltageIntegral" << start[0] << ", " << start[1] << ", " << start[2] << " -> " << stop[0] << ", " << stop[1] << ", " << stop[2] << ", " << endl;
	double result=0;
	for (int n=0; n<3; ++n)
	{
		if (start[n]<stop[n])
		{
			unsigned int pos[3]={start[0],start[1],start[2]};
			for (; pos[n]<stop[n]; ++pos[n])
				result += m_Eng->GetVolt(n,pos[0],pos[1],pos[2]);

		}
		else
		{
			unsigned int pos[3]={stop[0],stop[1],stop[2]};
			for (; pos[n]<start[n]; ++pos[n])
				result -= m_Eng->GetVolt(n,pos[0],pos[1],pos[2]);
		}
	}
	return result;
}

double Engine_Interface_FDTD::GetRawField(unsigned int n, const unsigned int* pos, int type) const
{
	double value = m_Eng->GetVolt(n,pos[0],pos[1],pos[2]);
	double delta = m_Op->GetEdgeLength(n,pos);
	if ((type==0) && (delta))
		return value/delta;
	if ((type==1) && (m_Op->m_kappa_ptr) && (delta)) {
		ArrayLib::ArrayNIJK<float>& kappa = *m_Op->m_kappa_ptr;
		return value * kappa(n, pos[0], pos[1], pos[2]) / delta;
	}
	if ((type==3) && (m_Op->m_epsR_ptr) && (delta)) {
		ArrayLib::ArrayNIJK<float>& epsR = *m_Op->m_epsR_ptr;
		return value * epsR(n, pos[0], pos[1], pos[2]) / delta;
	}
	if (type==2) //calc rot(H)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;
		unsigned int locPos[] = {pos[0],pos[1],pos[2]};
		double area = m_Op->GetEdgeArea(n,pos);
		value  = m_Eng->GetCurr(nPP,pos);
		value -= m_Eng->GetCurr(nP,pos);
		if (pos[nPP]>0)
		{
			--locPos[nPP];
			value += m_Eng->GetCurr(nP,locPos);
			++locPos[nPP];
		}
		if (pos[nP]>0)
		{
			--locPos[nP];
			value -= m_Eng->GetCurr(nPP,locPos);
		}
		return value/area;
	}

	return 0.0;
}

namespace
{
//! Precomputed GetRawInterpolatedField()/GetRawInterpolatedDualField() (type 0) of the
//! nodes of a dump, for engines with the fields in the basic engine layout (GPU host mirror).
//! Every node component is one of a few forms, evaluated with the same operations as there.
class Field_Gather_FDTD : public Engine_Field_Gather
{
public:
	enum Form { ZERO, RAW, LERP, AVG4 };
	//! one component of a node: raw(k) = value(idx[k])/delta[k] (0 if delta is 0), see GetRawField()
	struct Entry
	{
		Form form;
		unsigned int idx[4];
		double delta[4];
		double rel;      //!< LERP: raw(0)*(1-rel) + raw(1)*rel
	};

	Field_Gather_FDTD(const Engine_GPU* eng, bool h_field, unsigned int nj, unsigned int nk)
		: m_Eng(eng), m_H(h_field), m_nj(nj), m_nk(nk) {}

	std::vector<Entry> entries;   //!< [line][k][n]

	virtual void Evaluate(size_t line_start, size_t line_stop, ArrayLib::ArrayNIJK<float> &field, const float* src=NULL) const
	{
		const float* f = src ? src : (m_H ? m_Eng->HostCurrents().data() : m_Eng->HostVoltages().data());
		for (size_t l=line_start; l<line_stop; ++l)
		{
			const unsigned int i = l/m_nj;
			const unsigned int j = l%m_nj;
			for (unsigned int k=0; k<m_nk; ++k)
				for (int n=0; n<3; ++n)
				{
					const Entry& e = entries[(l*m_nk + k)*3 + n];
					double out = 0;
					switch (e.form)
					{
					case ZERO:
						break;
					case RAW:
						out = Raw(f, e, 0);
						break;
					case LERP:
						out = Raw(f, e, 0)*(1.0-e.rel) + Raw(f, e, 1)*e.rel;
						break;
					case AVG4:
						out = Raw(f, e, 0);
						out+= Raw(f, e, 1);
						out+= Raw(f, e, 2);
						out+= Raw(f, e, 3);
						out/=4;
						break;
					}
					field(n, i, j, k) = out;
				}
		}
	}

protected:
	static inline double Raw(const float* f, const Entry& e, int k)
	{
		double value = f[e.idx[k]];
		if (e.delta[k])
			return value/e.delta[k];
		return 0.0;
	}

	const Engine_GPU* m_Eng;
	bool m_H;
	unsigned int m_nj, m_nk;
};
}

// The same values as GetEField()/GetHField() of this class: subclasses that change the
// field evaluation must override this (see Engine_Interface_Cylindrical_FDTD).
Engine_Field_Gather* Engine_Interface_FDTD::CreateFieldGather(bool h_field, const unsigned int numLines[3], unsigned int* const posLines[3]) const
{
	// only for the GPU engine, whose host mirror has the basic engine layout
	const Engine_GPU* eng_gpu = dynamic_cast<const Engine_GPU*>(m_Eng);
	if (!eng_gpu)
		return NULL;

	unsigned int N[3];
	for (int n=0; n<3; ++n)
		N[n] = m_Op->GetNumberOfLines(n, true);
	auto index = [&](int n, const unsigned int* p) -> unsigned int {return ((n*N[0] + p[0])*N[1] + p[1])*N[2] + p[2];};
	auto delta = [&](int n, const unsigned int* p) -> double {return m_Op->GetEdgeLength(n, p, h_field);};

	typedef Field_Gather_FDTD G;
	G* gather = new G(eng_gpu, h_field, numLines[1], numLines[2]);
	gather->entries.resize((size_t)numLines[0]*numLines[1]*numLines[2]*3);
	size_t e_idx = 0;
	for (unsigned int i=0; i<numLines[0]; ++i)
		for (unsigned int j=0; j<numLines[1]; ++j)
			for (unsigned int k=0; k<numLines[2]; ++k)
			{
				const unsigned int pos[3] = {posLines[0][i], posLines[1][j], posLines[2][k]};
				for (int n=0; n<3; ++n)
				{
					G::Entry& e = gather->entries[e_idx++];
					e.form = G::ZERO;
					e.rel = 0;
					for (int m=0; m<4; ++m) {e.idx[m] = 0; e.delta[m] = 0;}
					unsigned int p[3] = {pos[0], pos[1], pos[2]};
					const int nP = (n+1)%3;
					const int nPP = (n+2)%3;
					auto set = [&](int m) {e.idx[m] = index(n, p); e.delta[m] = delta(n, p);};
					const bool upper = (pos[0]==N[0]-1) || (pos[1]==N[1]-1) || (pos[2]==N[2]-1);
					switch (m_InterpolType)
					{
					default:
					case NO_INTERPOLATION:
						e.form = G::RAW;
						set(0);
						break;
					case NODE_INTERPOLATE:
						if (!h_field)
						{
							// see GetRawInterpolatedField()
							if (pos[n]==N[n]-1)
							{
								--p[n];
								e.form = G::RAW;
								set(0);
								break;
							}
							const double d = delta(n, p);
							if (d==0)
								break;   // ZERO
							e.form = G::RAW;
							set(0);
							if (pos[n]==0)
								break;
							--p[n];
							const double d_down = delta(n, p);
							e.form = G::LERP;
							e.rel = d / (d+d_down);
							set(1);
						}
						else
						{
							// see GetRawInterpolatedDualField()
							if (upper || (pos[nP]==0) || (pos[nPP]==0))
								break;   // ZERO
							e.form = G::AVG4;
							set(0);
							--p[nP];
							set(1);
							--p[nPP];
							set(2);
							++p[nP];
							set(3);
						}
						break;
					case CELL_INTERPOLATE:
						if (!h_field)
						{
							// see GetRawInterpolatedField()
							if (upper)
								break;   // ZERO
							e.form = G::AVG4;
							set(0);
							++p[nP];
							set(1);
							++p[nPP];
							set(2);
							--p[nP];
							set(3);
						}
						else
						{
							// see GetRawInterpolatedDualField()
							if (pos[n]>=N[n]-1)
								break;   // ZERO
							const double d = delta(n, p);
							e.form = G::LERP;
							set(0);
							++p[n];
							const double d_up = delta(n, p);
							e.rel = d / (d+d_up);
							set(1);
						}
						break;
					}
				}
			}
	return gather;
}

bool Engine_Interface_FDTD::TakeFieldSnapshot(unsigned int slot, const float* &volt, const float* &curr)
{
	Engine_GPU* eng_gpu = dynamic_cast<Engine_GPU*>(m_Eng);
	return eng_gpu && eng_gpu->SnapshotFields(slot, volt, curr);
}

void Engine_Interface_FDTD::WaitFieldSnapshot(unsigned int slot) const
{
	const Engine_GPU* eng_gpu = dynamic_cast<const Engine_GPU*>(m_Eng);
	if (eng_gpu)
		eng_gpu->WaitSnapshot(slot);
}

void Engine_Interface_FDTD::PrepareFieldAccess()
{
	// the GPU engine reads out-of-date values from the device, which is not thread-safe
	Engine_GPU* eng_gpu = dynamic_cast<Engine_GPU*>(m_Eng);
	if (eng_gpu)
		eng_gpu->UpdateHostMirror();
}

double Engine_Interface_FDTD::CalcFastEnergy() const
{
	double E_energy=0.0;
	double H_energy=0.0;

	unsigned int pos[3];
	if (m_Eng->GetType()==Engine::GPU)
	{
		// on the device, if the backend can
		Engine_GPU* eng_gpu = dynamic_cast<Engine_GPU*>(m_Eng);
		unsigned int numNodes[3];
		for (int n=0; n<3; ++n)
			numNodes[n] = m_Op->GetNumberOfLines(n)-1;
		if (eng_gpu && eng_gpu->CalcFastEnergy(numNodes, E_energy, H_energy))
			return EPS0*E_energy + MUE0*H_energy;
		// else on the host mirror below
		if (eng_gpu)
			eng_gpu->UpdateHostMirror();
	}
	// the GPU engine keeps a host mirror in the basic engine layout
	if ((m_Eng->GetType()==Engine::BASIC) || (m_Eng->GetType()==Engine::GPU))
	{
		for (pos[0]=0; pos[0]<m_Op->GetNumberOfLines(0)-1; ++pos[0])
		{
			for (pos[1]=0; pos[1]<m_Op->GetNumberOfLines(1)-1; ++pos[1])
			{
				for (pos[2]=0; pos[2]<m_Op->GetNumberOfLines(2)-1; ++pos[2])
				{
					E_energy+=m_Eng->Engine::GetVolt(0,pos[0],pos[1],pos[2]) * m_Eng->Engine::GetVolt(0,pos[0],pos[1],pos[2]);
					E_energy+=m_Eng->Engine::GetVolt(1,pos[0],pos[1],pos[2]) * m_Eng->Engine::GetVolt(1,pos[0],pos[1],pos[2]);
					E_energy+=m_Eng->Engine::GetVolt(2,pos[0],pos[1],pos[2]) * m_Eng->Engine::GetVolt(2,pos[0],pos[1],pos[2]);

					H_energy+=m_Eng->Engine::GetCurr(0,pos[0],pos[1],pos[2]) * m_Eng->Engine::GetCurr(0,pos[0],pos[1],pos[2]);
					H_energy+=m_Eng->Engine::GetCurr(1,pos[0],pos[1],pos[2]) * m_Eng->Engine::GetCurr(1,pos[0],pos[1],pos[2]);
					H_energy+=m_Eng->Engine::GetCurr(2,pos[0],pos[1],pos[2]) * m_Eng->Engine::GetCurr(2,pos[0],pos[1],pos[2]);
				}
			}
		}
	}
	else
	{
		for (pos[0]=0; pos[0]<m_Op->GetNumberOfLines(0)-1; ++pos[0])
		{
			for (pos[1]=0; pos[1]<m_Op->GetNumberOfLines(1)-1; ++pos[1])
			{
				for (pos[2]=0; pos[2]<m_Op->GetNumberOfLines(2)-1; ++pos[2])
				{
					E_energy+=m_Eng->GetVolt(0,pos[0],pos[1],pos[2]) * m_Eng->GetVolt(0,pos[0],pos[1],pos[2]);
					E_energy+=m_Eng->GetVolt(1,pos[0],pos[1],pos[2]) * m_Eng->GetVolt(1,pos[0],pos[1],pos[2]);
					E_energy+=m_Eng->GetVolt(2,pos[0],pos[1],pos[2]) * m_Eng->GetVolt(2,pos[0],pos[1],pos[2]);

					H_energy+=m_Eng->GetCurr(0,pos[0],pos[1],pos[2]) * m_Eng->GetCurr(0,pos[0],pos[1],pos[2]);
					H_energy+=m_Eng->GetCurr(1,pos[0],pos[1],pos[2]) * m_Eng->GetCurr(1,pos[0],pos[1],pos[2]);
					H_energy+=m_Eng->GetCurr(2,pos[0],pos[1],pos[2]) * m_Eng->GetCurr(2,pos[0],pos[1],pos[2]);
				}
			}
		}
	}
	return EPS0*E_energy + MUE0*H_energy;
}
