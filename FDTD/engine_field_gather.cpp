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


#include <climits>

#include "engine_field_gather.h"
#include "engine_interface_fdtd.h"
#include "operator.h"

Field_Gather_FDTD::Field_Gather_FDTD(const Engine_Interface_FDTD* eng_if, const Operator* op,
                                     const ArrayLib::ArrayNIJK<FDTD_FLOAT>* volt, const ArrayLib::ArrayNIJK<FDTD_FLOAT>* curr,
                                     bool h_field, const unsigned int numLines[3], unsigned int* const posLines[3])
	: m_EngIf(eng_if), m_Op(op), m_Volt(volt), m_Curr(curr), m_H(h_field), m_RegionSnapshot(false), m_Offset(0)
{
	for (int n=0; n<3; ++n)
	{
		m_PosLines[n].assign(posLines[n], posLines[n]+numLines[n]);
		m_N[n] = op->GetNumberOfLines(n, true);
	}
	// the layout of the fields: the basic engine's (NIJK, z fastest) without arrays
	m_S[3] = 1;
	m_S[2] = m_N[2];
	m_S[1] = (size_t)m_N[1]*m_N[2];
	m_S[0] = (size_t)m_N[0]*m_N[1]*m_N[2];
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>* layout = volt ? volt : curr;
	if (layout)
		for (int d=0; d<4; ++d)
			m_S[d] = layout->stride(d);
	for (int n=0; n<3; ++n)
		for (int d=0; d<3; ++d)
			m_RegionPos[n][d].assign(m_N[d], UINT_MAX);

	// the region: the lines of the edges the interpolation of the dumped nodes reads
	auto mark = [this](unsigned int n, const unsigned int* p) -> double
	{
		for (int d=0; d<3; ++d)
			m_RegionPos[n][d][p[d]] = 0;
		return 0.0;
	};
	double out[3];
	for (unsigned int i=0; i<numLines[0]; ++i)
		for (unsigned int j=0; j<numLines[1]; ++j)
			for (unsigned int k=0; k<numLines[2]; ++k)
				Interpolate(i, j, k, out, mark);
	size_t start = 0;
	for (int n=0; n<3; ++n)
	{
		for (int d=0; d<3; ++d)
			for (unsigned int l=0; l<m_N[d]; ++l)
				if (m_RegionPos[n][d][l]!=UINT_MAX)
				{
					m_RegionPos[n][d][l] = region.lines[n][d].size();
					region.lines[n][d].push_back(l);
				}
		m_RegionStart[n] = start;
		start += region.Size(n);
	}
}

template <typename T, typename Raw>
void Field_Gather_FDTD::Interpolate(unsigned int i, unsigned int j, unsigned int k, T out[3], Raw raw) const
{
	const unsigned int pos[3] = {m_PosLines[0][i], m_PosLines[1][j], m_PosLines[2][k]};
	if (m_H)
		m_EngIf->InterpolateDualField(pos, out, raw);
	else
		m_EngIf->InterpolateField(pos, out, raw);
}

void Field_Gather_FDTD::BuildTable(size_t budget)
{
	m_TabStart.clear();
	m_TabIndex.clear();
	m_TabWeight.clear();
	const size_t ni = m_PosLines[0].size(), nj = m_PosLines[1].size(), nk = m_PosLines[2].size();
	const size_t nodes = ni*nj*nk;
	// estimated with 4 weights per component (the cell interpolation), 12 bytes each
	if ((nodes==0) || (nodes*3*(4*12 + 4) > budget))
		return;
	std::vector<uint32_t> start;
	std::vector<size_t> index;
	std::vector<float> weight;
	start.reserve(3*nodes+1);
	struct Edge {unsigned int n, p[3];};
	std::vector<Edge> edges;
	std::vector<std::pair<size_t, float> > comp[3];
	double out[3];
	for (unsigned int i=0; i<ni; ++i)
		for (unsigned int j=0; j<nj; ++j)
			for (unsigned int k=0; k<nk; ++k)
			{
				// the raw values the node reads, each once
				edges.clear();
				Interpolate(i, j, k, out, [&](unsigned int n, const unsigned int* p) -> double {
					for (const Edge& e : edges)
						if ((e.n==n) && (e.p[0]==p[0]) && (e.p[1]==p[1]) && (e.p[2]==p[2]))
							return 0.0;
					edges.push_back(Edge{n, {p[0], p[1], p[2]}});
					return 0.0;
				});
				for (int c=0; c<3; ++c)
					comp[c].clear();
				// the weight of each: the node for this raw value (divided by its edge length) 1, the others 0
				for (const Edge& e : edges)
				{
					const double delta = m_Op->GetEdgeLength(e.n, e.p, m_H);
					if (delta==0)
						continue;   // read as 0
					Interpolate(i, j, k, out, [&](unsigned int n, const unsigned int* p) -> double {
						return ((n==e.n) && (p[0]==e.p[0]) && (p[1]==e.p[1]) && (p[2]==e.p[2])) ? 1.0 : 0.0;
					});
					const size_t idx = e.n*m_S[0] + e.p[0]*m_S[1] + e.p[1]*m_S[2] + e.p[2]*m_S[3];
					for (int c=0; c<3; ++c)
						if (out[c]!=0)
							comp[c].push_back(std::make_pair(idx, (float)(out[c]/delta)));
				}
				for (int c=0; c<3; ++c)
				{
					if (index.size() > 0xffffffffu)
						return;   // the old path
					start.push_back((uint32_t)index.size());
					for (const std::pair<size_t, float>& w : comp[c])
					{
						index.push_back(w.first);
						weight.push_back(w.second);
					}
				}
				if (start.size()*4 + index.size()*12 > budget)
					return;   // more weights than estimated: the old path
			}
	start.push_back((uint32_t)index.size());
	m_TabStart.swap(start);
	m_TabIndex.swap(index);
	m_TabWeight.swap(weight);
}

void Field_Gather_FDTD::Evaluate(size_t line_start, size_t line_stop, ArrayLib::ArrayNIJK<float> &field, const float* src) const
{
	if (!src && !m_TabStart.empty())
	{
		const float* f = m_H ? m_Curr->data() : m_Volt->data();
		const unsigned int nj = m_PosLines[1].size(), nk = m_PosLines[2].size();
		for (size_t l=line_start; l<line_stop; ++l)
			for (unsigned int k=0; k<nk; ++k)
			{
				const size_t node = l*nk + k;
				for (int c=0; c<3; ++c)
				{
					double sum = 0;
					for (uint32_t w=m_TabStart[3*node+c]; w<m_TabStart[3*node+c+1]; ++w)
						sum += (double)m_TabWeight[w]*f[m_TabIndex[w]];
					field(c, l/nj, l%nj, k) = sum;
				}
			}
		return;
	}
	// the raw values of the region in a snapshot (see SetSnapshotRegion()), else the whole field
	const bool in_region = src && m_RegionSnapshot;
	const float* f = in_region ? src + m_Offset : (src ? src : (m_H ? m_Curr->data() : m_Volt->data()));
	// as GetRawField()/GetRawDualField() (type 0) of Engine_Interface_FDTD
	auto raw = [&](unsigned int n, const unsigned int* p) -> double
	{
		double value = Value(f, in_region, n, p);
		double delta = m_Op->GetEdgeLength(n, p, m_H);
		if (delta)
			return value/delta;
		return 0.0;
	};
	const unsigned int nj = m_PosLines[1].size();
	double out[3];
	for (size_t l=line_start; l<line_stop; ++l)
		for (unsigned int k=0; k<m_PosLines[2].size(); ++k)
		{
			Interpolate(l/nj, l%nj, k, out, raw);
			for (int n=0; n<3; ++n)
				field(n, l/nj, l%nj, k) = out[n];
		}
}

void Field_Gather_FDTD::EvaluateSums(const std::complex<float>* raw_sums, ArrayLib::ArrayNIJK<std::complex<float>> &field) const
{
	auto raw = [&](unsigned int n, const unsigned int* p) -> std::complex<double>
	{
		std::complex<double> value = Value(raw_sums, true, n, p);
		double delta = m_Op->GetEdgeLength(n, p, m_H);
		if (delta)
			return value/delta;
		return 0.0;
	};
	std::complex<double> out[3];
	for (unsigned int i=0; i<m_PosLines[0].size(); ++i)
		for (unsigned int j=0; j<m_PosLines[1].size(); ++j)
			for (unsigned int k=0; k<m_PosLines[2].size(); ++k)
			{
				Interpolate(i, j, k, out, raw);
				for (int n=0; n<3; ++n)
					field(n, i, j, k) = out[n];
			}
}
