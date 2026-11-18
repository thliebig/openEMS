/*
*	Copyright (C) 2026 Roel Budé
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

#include "gpu_split.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "gpu_cell_model.h"
#include "tools/global.h"

namespace
{
//! The device extension of an engine extension on every part of a split grid, in order of the parts
class GPU_Ext_Split : public GPU_Extension
{
public:
	GPU_Ext_Split(GPU_Backend_Split* split, const std::vector<GPU_Extension*>& parts) : m_Split(split), m_Parts(parts) {}
	virtual ~GPU_Ext_Split()
	{
		for (GPU_Extension* e : m_Parts)
			delete e;
	}

	virtual void Hook(bool currents, GPU_Hook h, const Engine_Box& box, unsigned int ts)
	{
		// the halo lines the parts read: the line below part p at cut p-1, the one above at cut p
		for (size_t p=0; p<m_Parts.size(); ++p)
			if (const int halos = m_Parts[p]->HaloReads(currents, h))
			{
				if ((halos & 1) && (p>0))
					m_Split->ExchangeCut(p-1, currents, 1);
				if ((halos & 2) && (p+1<m_Parts.size()))
					m_Split->ExchangeCut(p, currents, 2);
			}
		Engine_Box local;
		for (size_t p=0; p<m_Parts.size(); ++p)
			if (m_Split->LocalBox(p, box, local))
				m_Parts[p]->Hook(currents, h, local, ts);
	}
	virtual bool NeedsFlatStep(unsigned int ts) const
	{
		for (const GPU_Extension* e : m_Parts)
			if (e->NeedsFlatStep(ts))
				return true;
		return false;
	}

protected:
	GPU_Backend_Split* m_Split;
	std::vector<GPU_Extension*> m_Parts;
};

}

std::vector<unsigned int> GPU_NumberList(const std::string& value, const std::string& option, const char* env, std::string& what, bool& from_env)
{
	std::string list = value;
	what = option + "=" + value;
	from_env = false;
	if (list.empty() && getenv(env))
	{
		list = getenv(env);
		what = std::string(env) + "=" + list;
		from_env = true;
	}
	std::vector<unsigned int> v;
	std::istringstream items(list);
	std::string item;
	while (std::getline(items, item, ','))
	{
		char* end = NULL;
		const long n = strtol(item.c_str(), &end, 10);
		if ((end==item.c_str()) || (*end!='\0') || (n<0))
			throw std::runtime_error(what + ": not a comma separated list of numbers");
		v.push_back((unsigned int)n);
	}
	return v;
}

std::vector<unsigned int> GPU_SplitCuts(unsigned int nx, size_t parts)
{
	std::string what;
	bool from_env;
	std::vector<unsigned int> cuts = GPU_NumberList(g_settings.GetGPUDomainCuts(), "--gpu-domain-cuts", "OPENEMS_GPU_DOMAIN_CUTS", what, from_env);
	if (cuts.empty())
	{
		// equal parts: the cross-section is the same along x
		for (size_t p=0; p<=parts; ++p)
			cuts.push_back((unsigned int)((uint64_t)nx*p/parts));
	}
	else
	{
		std::ostringstream err;
		if (cuts.size()!=parts+1)
			err << what << ": " << parts+1 << " x lines for " << parts << " devices";
		else if ((cuts.front()!=0) || (cuts.back()!=nx))
			err << what << ": from 0 to the number of x lines, " << nx;
		if (!err.str().empty())
			throw std::runtime_error(err.str());
	}
	for (size_t p=0; p<parts; ++p)
		if ((cuts[p+1]<cuts[p]) || (cuts[p+1]-cuts[p]<3))
		{
			std::ostringstream err;
			err << "GPU_Backend_Split: a part of fewer than 3 x lines (lines " << cuts[p] << " to " << cuts[p+1]
			    << " of " << nx << ") in the split over " << parts << " devices";
			throw std::runtime_error(err.str());
		}
	return cuts;
}

GPU_Backend_Split::GPU_Backend_Split(const std::vector<GPU_Backend*>& parts) : m_Parts(parts)
{
}

GPU_Backend_Split::~GPU_Backend_Split()
{
	for (GPU_Backend* b : m_Parts)
		delete b;
}

std::string GPU_Backend_Split::GetName() const
{
	std::string name;
	for (size_t p=0; p<m_Parts.size(); ++p)
		name += (p ? ", " : "") + m_Parts[p]->GetName();
	return std::to_string(m_Parts.size()) + " parts: " + name;
}

bool GPU_Backend_Split::Init(const Operator* op, const GPU_CellModel& model)
{
	const size_t P = m_Parts.size();
	if (!FitsOperator(op, model))
		return false;
	// before anything is on the devices: whether the parts can copy to each other at all
	for (size_t p=0; p+1<P; ++p)
		if (!m_Parts[p]->CanSendPlane(m_Parts[p+1]) || !m_Parts[p+1]->CanSendPlane(m_Parts[p]))
			throw std::runtime_error("GPU_Backend_Split: the backends " + m_Parts[p]->GetName() + " and " + m_Parts[p+1]->GetName()
			                         + " cannot split a grid between them");
	const std::vector<unsigned int> cuts = GPU_SplitCuts(model.numLines[0], P);
	m_Slabs.resize(P);
	for (size_t p=0; p<P; ++p)
	{
		GPU_Slab& s = m_Slabs[p];
		for (int a=0; a<3; ++a)
			s.lines[a] = model.numLines[a];
		s.x0 = cuts[p];
		s.x1 = cuts[p+1];
		s.off = s.x0 - (p>0 ? 1 : 0);
		s.nx = s.x1 + (p+1<P ? 1 : 0) - s.off;
		if (g_settings.GetVerboseLevel()>0)
			std::cout << "GPU_Backend_Split: part " << p << ", x lines " << s.x0 << " to " << s.x1-1 << ": " << m_Parts[p]->GetName() << std::endl;
		// the model of the part, freed once it is on its device; its grid is the model's (op NULL)
		if (!m_Parts[p]->Init(NULL, GPU_CellModel(model, s.off, s.off+s.nx)))
			return false;
	}
	return true;
}

void* GPU_Backend_Split::DeviceAlloc(size_t bytes, const void* host)
{
	UNUSED(bytes); UNUSED(host);
	throw std::logic_error("GPU_Backend_Split: no device memory of its own");
}

void GPU_Backend_Split::DeviceDownload(void* host, const void* dev, size_t bytes)
{
	UNUSED(host); UNUSED(dev); UNUSED(bytes);
	throw std::logic_error("GPU_Backend_Split: no device memory of its own");
}

float* GPU_Backend_Split::DeviceField(bool currents)
{
	UNUSED(currents);
	throw std::logic_error("GPU_Backend_Split: no device fields of its own");
}

void GPU_Backend_Split::RunPoints(GPU_PointOp kind, const void* op, uint64_t begin, uint64_t end, bool serial, const CM_PointFilter& filter)
{
	UNUSED(kind); UNUSED(op); UNUSED(begin); UNUSED(end); UNUSED(serial); UNUSED(filter);
	throw std::logic_error("GPU_Backend_Split: the device extensions run on the parts");
}

// the voltage update at x reads the currents at x-1: the line below each part; the current
// update at x the voltages at x+1: the line above
void GPU_Backend_Split::UpdateVoltages(const Engine_Box& box)
{
	Exchange(true, 1);
	Engine_Box local;
	for (size_t p=0; p<m_Parts.size(); ++p)
		if (LocalBox(p, box, local))
			m_Parts[p]->UpdateVoltages(local);
}

void GPU_Backend_Split::UpdateCurrents(const Engine_Box& box)
{
	Exchange(false, 2);
	Engine_Box local;
	for (size_t p=0; p<m_Parts.size(); ++p)
		if (LocalBox(p, box, local))
			m_Parts[p]->UpdateCurrents(local);
}

void GPU_Backend_Split::Exchange(bool currents, int halos)
{
	for (size_t c=0; c+1<m_Parts.size(); ++c)
		ExchangeCut(c, currents, halos);
}

void GPU_Backend_Split::ExchangeCut(size_t c, bool currents, int halos)
{
	const GPU_Slab& lo = m_Slabs.at(c);
	const GPU_Slab& hi = m_Slabs.at(c+1);
	const unsigned int cut = lo.x1;   // the first line of the part above
	if (halos & 1)
		m_Parts[c]->SendPlane(currents, cut-1-lo.off, m_Parts[c+1], cut-1-hi.off);
	if (halos & 2)
		m_Parts[c+1]->SendPlane(currents, cut-hi.off, m_Parts[c], cut-lo.off);
}

bool GPU_Backend_Split::LocalBox(size_t p, const Engine_Box& box, Engine_Box& local) const
{
	const GPU_Slab& s = m_Slabs.at(p);
	local = box;
	const unsigned int x0 = std::max(box.start[0], s.x0), x1 = std::min(box.stop[0], s.x1);
	if (x0>=x1)
		return false;
	local.start[0] = x0-s.off;
	local.stop[0] = x1-s.off;
	return true;
}

void GPU_Backend_Split::Download(bool currents, ArrayLib::ArrayNIJK<FDTD_FLOAT>& field)
{
	const GPU_Slab& s0 = m_Slabs.at(0);
	const size_t plane = (size_t)s0.lines[1]*s0.lines[2];
	if (field.size()!=3*plane*s0.lines[0])
		throw std::runtime_error("GPU_Backend_Split: host and device field size mismatch");
	// per part and component its owned lines, contiguous in both layouts; the parts copy at once, into
	// the mirror page-locked once (for all devices: the address space is unified)
	m_Parts.front()->PinHost(field.data(), field.size()*sizeof(FDTD_FLOAT));
	for (size_t p=0; p<m_Parts.size(); ++p)
	{
		const GPU_Slab& s = m_Slabs[p];
		const float* dev = m_Parts[p]->DeviceField(currents);
		for (unsigned int n=0; n<3; ++n)
			m_Parts[p]->DeviceDownloadAsync(field.data() + (n*(size_t)s.lines[0] + s.x0)*plane, dev + (n*(size_t)s.nx + s.x0-s.off)*plane,
			                                (s.x1-s.x0)*plane*sizeof(FDTD_FLOAT));
	}
	Synchronize();
}

void GPU_Backend_Split::ReadLine(bool currents, unsigned int n, unsigned int x, unsigned int y, const unsigned int numLines[3], FDTD_FLOAT* line)
{
	UNUSED(numLines);   // the parts' slabs
	for (size_t p=0; p<m_Parts.size(); ++p)
	{
		const GPU_Slab& s = m_Slabs[p];
		if (!s.Owns(x))
			continue;
		const float* dev = m_Parts[p]->DeviceField(currents);
		m_Parts[p]->DeviceDownload(line, dev + ((n*(size_t)s.nx + x-s.off)*s.lines[1] + y)*s.lines[2], s.lines[2]*sizeof(FDTD_FLOAT));
		return;
	}
	throw std::out_of_range("GPU_Backend_Split::ReadLine: no such x line");
}

void GPU_Backend_Split::Synchronize()
{
	for (GPU_Backend* b : m_Parts)
		b->Synchronize();
}

bool GPU_Backend_Split::CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy)
{
	E_energy = H_energy = 0;
	for (size_t p=0; p<m_Parts.size(); ++p)
	{
		const GPU_Slab& s = m_Slabs[p];
		const unsigned int x1 = std::min(numNodes[0], s.x1);
		if (x1<=s.x0)
			continue;
		const unsigned int local[3] = {x1-s.off, numNodes[1], numNodes[2]};
		double E = 0, H = 0;
		if (!m_Parts[p]->CalcFastEnergyLines(s.x0-s.off, local, E, H))
			return false;
		E_energy += E;
		H_energy += H;
	}
	return true;
}

GPU_Backend* GPU_Backend_Split::NewSubGridBackend()
{
	throw std::runtime_error("GPU_Backend_Split: a cylindrical multi-grid cannot be split over several devices");
}

GPU_Extension* GPU_Backend_Split::CreateExtension(Engine_Extension* eng_ext, const unsigned int numLines[3])
{
	UNUSED(numLines);   // the parts' slabs
	std::vector<GPU_Extension*> parts;
	for (size_t p=0; p<m_Parts.size(); ++p)
	{
		GPU_Extension* e = GPU_CreatePartExtension(eng_ext, m_Parts[p], m_Slabs[p], parts.empty() ? NULL : parts.front());
		if (e==NULL)
		{
			for (GPU_Extension* q : parts)
				delete q;
			return NULL;
		}
		parts.push_back(e);
	}
	return new GPU_Ext_Split(this, parts);
}
