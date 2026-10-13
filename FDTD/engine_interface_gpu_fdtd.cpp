/*
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
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

#include <stdexcept>
#include <vector>

#include "engine_interface_gpu_fdtd.h"
#include "engine_field_gather.h"
#include "engine_gpu.h"

Engine_Interface_GPU_FDTD::Engine_Interface_GPU_FDTD(Operator* op) : Engine_Interface_FDTD(op)
{
	m_Eng_GPU = dynamic_cast<Engine_GPU*>(m_Eng);
	if (m_Eng_GPU==NULL)
		throw std::runtime_error("Engine_Interface_GPU_FDTD::Engine_Interface_GPU_FDTD: Error: GPU-Engine is not set!");
}

Engine_Interface_GPU_FDTD::~Engine_Interface_GPU_FDTD()
{
	m_Eng_GPU = NULL;
}

double Engine_Interface_GPU_FDTD::CalcFastEnergy() const
{
	// the end criteria asks for this every Nyquist period, so it is worth the
	// device implementation; without one the host mirror has to be brought up to
	// date, it is not written at the end of a batch
	unsigned int numNodes[3];
	for (int n=0; n<3; ++n)
		numNodes[n] = m_Op->GetNumberOfLines(n)-1;
	double E_energy=0.0, H_energy=0.0;
	if (m_Eng_GPU->CalcFastEnergy(numNodes, E_energy, H_energy))
		return EPS0*E_energy + MUE0*H_energy;
	m_Eng_GPU->UpdateHostMirror();
	return Engine_Interface_FDTD::CalcFastEnergy();
}

void Engine_Interface_GPU_FDTD::PrepareFieldAccess()
{
	// several threads read the fields next; reading them from the device value by value is not thread-safe
	m_Eng_GPU->UpdateHostMirror();
}

bool Engine_Interface_GPU_FDTD::TakeFieldSnapshot(unsigned int slot, const float* &volt, const float* &curr)
{
	return m_Eng_GPU->SnapshotFields(slot, volt, curr);
}

void Engine_Interface_GPU_FDTD::WaitFieldSnapshot(unsigned int slot) const
{
	m_Eng_GPU->WaitSnapshot(slot);
}

bool Engine_Interface_GPU_FDTD::PrepareSnapshotGather(Engine_Field_Gather* gather)
{
	if (!m_Eng_GPU->CanSnapshotGather())
		return true;   // plain field snapshots, if the engine has any: evaluated on the host
	// the engine only snapshots the dumped nodes, so it has to know them beforehand
	Field_Gather_FDTD* g = dynamic_cast<Field_Gather_FDTD*>(gather);
	size_t offset = 0;
	if (!g || !m_Eng_GPU->AddSnapshotGather(g->HField(), g->entries, offset))
		return false;
	g->SetSnapshotEvaluated(offset);
	return true;
}

Engine_Field_Gather* Engine_Interface_GPU_FDTD::CreateFieldGather(bool h_field, const unsigned int numLines[3], unsigned int* const posLines[3]) const
{
	// the gather evaluates the nodes like GetEField()/GetHField(), which move the nodes
	// on the seam of a closed mesh to its other side (see GetSeamLine())
	std::vector<unsigned int> seam[3];
	unsigned int* lines[3];
	for (int n=0; n<3; ++n)
	{
		seam[n].assign(posLines[n], posLines[n]+numLines[n]);
		for (size_t l=0; l<seam[n].size(); ++l)
			seam[n][l] = GetSeamLine(n, seam[n][l]);
		lines[n] = seam[n].data();
	}
	return Field_Gather_FDTD::New(m_Op, m_InterpolType, &m_Eng_GPU->HostVoltages(), &m_Eng_GPU->HostCurrents(), h_field, numLines, lines);
}

int Engine_Interface_GPU_FDTD::CreateFieldDFT(const Engine_Field_Gather* gather, unsigned int count)
{
	const Field_Gather_FDTD* g = dynamic_cast<const Field_Gather_FDTD*>(gather);
	if (!g || g->entries.empty())
		return -1;
	const int id = m_Eng_GPU->AddFieldDFT(g->HField(), g->entries, count);
	if (id>=0)
	{
		const unsigned int nj = g->NumLinesJ(), nk = g->NumLinesK();
		m_FieldDFT[id] = {(unsigned int)(g->entries.size()/(3*(size_t)nj*nk)), nj, nk};
	}
	return id;
}

void Engine_Interface_GPU_FDTD::AccumulateFieldDFT(int id, const std::vector<std::complex<float>>& weights)
{
	m_Eng_GPU->AccumulateFieldDFT(id, weights);
}

bool Engine_Interface_GPU_FDTD::ReadFieldDFT(int id, std::vector<ArrayLib::ArrayNIJK<std::complex<float>>*>& fields)
{
	std::map<int, std::array<unsigned int, 3>>::const_iterator it = m_FieldDFT.find(id);
	std::vector<std::complex<float>> sums;
	if ((it==m_FieldDFT.end()) || !m_Eng_GPU->ReadFieldDFT(id, sums))
		return false;
	// the entries in the order of Field_Gather_FDTD: [line][k][n]
	const unsigned int ni = it->second[0], nj = it->second[1], nk = it->second[2];
	const size_t count = (size_t)ni*nj*nk*3;
	if (sums.size()!=count*fields.size())
		return false;
	for (size_t f=0; f<fields.size(); ++f)
	{
		const std::complex<float>* s = sums.data() + f*count;
		ArrayLib::ArrayNIJK<std::complex<float>>& field = *fields[f];
		for (unsigned int i=0; i<ni; ++i)
			for (unsigned int j=0; j<nj; ++j)
				for (unsigned int k=0; k<nk; ++k)
					for (int n=0; n<3; ++n)
						field(n, i, j, k) = *s++;
	}
	return true;
}
