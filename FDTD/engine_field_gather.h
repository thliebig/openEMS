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

#ifndef ENGINE_FIELD_GATHER_H
#define ENGINE_FIELD_GATHER_H

#include <climits>
#include <complex>
#include <vector>
#include "Common/engine_interface_base.h"
#include "tools/constants.h"
#include "tools/arraylib/array_nijk.h"

class Operator;

//! One component of a dumped node: the (interpolated) E or H field from up to four voltages or currents
/*!
  idx[k] indexes the voltages or currents in the flat ArrayNIJK layout,
  ((n*numLines[0] + x)*numLines[1] + y)*numLines[2] + z.
  raw(k) = field[idx[k]] / delta[k] (0 if delta[k] is 0), and the value is
  - ZERO: 0
  - RAW: raw(0)
  - LERP: raw(0)*(1-rel) + raw(1)*rel
  - AVG4: (raw(0) + raw(1) + raw(2) + raw(3)) / 4
  in double precision, rounded to float. See Field_Gather_FDTD.
  */
struct Engine_GatherEntry
{
	enum Form { ZERO, RAW, LERP, AVG4 };
	unsigned int form;
	unsigned int idx[4];
	double delta[4];
	double rel;
};

//! The raw voltages or currents a gather reads
/*!
  For each component n the edges on the mesh lines lines[n][0] x lines[n][1] x
  lines[n][2], component after component, in the order of the flat ArrayNIJK layout.
  A backend copies or sums these raw values (a few line tables on the device instead of
  an entry per dumped value); the host applies the entries to them, see
  Field_Gather_FDTD::RegionIndex().
  */
struct Engine_GatherRegion
{
	std::vector<unsigned int> lines[3][3];   //!< [n][direction], ascending

	size_t Size(int n) const {return lines[n][0].size()*lines[n][1].size()*lines[n][2].size();}
	size_t Size() const {return Size(0)+Size(1)+Size(2);}
};

//! Precomputed GetRawInterpolatedField()/GetRawInterpolatedDualField() (type 0) of the
//! nodes of a dump, for fields in the basic engine layout (e.g. the host mirror of Engine_GPU).
//! Every node component is one of a few forms, evaluated with the same operations as there.
class Field_Gather_FDTD : public Engine_Field_Gather
{
public:
	//! one component of a node: raw(k) = value(idx[k])/delta[k] (0 if delta is 0), see GetRawField().
	//! Evaluated on the host; a backend only copies or sums the raw values of \a region.
	typedef Engine_GatherEntry Entry;
	typedef Engine_GatherEntry G;

	Field_Gather_FDTD(const ArrayLib::ArrayNIJK<FDTD_FLOAT>* volt, const ArrayLib::ArrayNIJK<FDTD_FLOAT>* curr,
	                  bool h_field, unsigned int nj, unsigned int nk)
		: m_Volt(volt), m_Curr(curr), m_H(h_field), m_nj(nj), m_nk(nk), m_Evaluated(false), m_Offset(0) {}

	//! The evaluation of GetEField() (\a h_field false) or GetHField() of Engine_Interface_FDTD with the interpolation \a interpol at the nodes posLines, reading \a volt or \a curr
	static Field_Gather_FDTD* New(const Operator* op, Engine_Interface_Base::InterpolationType interpol,
	                              const ArrayLib::ArrayNIJK<FDTD_FLOAT>* volt, const ArrayLib::ArrayNIJK<FDTD_FLOAT>* curr,
	                              bool h_field, const unsigned int numLines[3], unsigned int* const posLines[3]);

	std::vector<Entry> entries;   //!< [line][k][n]
	Engine_GatherRegion region;   //!< the raw values the entries read

	//! The values of the entries from the raw values of \a region (snapshot \a src, or the sums of a frequency domain dump)
	template <typename R, typename T>
	R Value(const T* raw, const Entry& e) const
	{
		return Combine<R>(e, [&](int k) -> R {return Raw<R>(raw[RegionIndex(e.idx[k])], e.delta[k]);});
	}
	//! Position of the field value \a idx (flat ArrayNIJK index) in \a region
	size_t RegionIndex(unsigned int idx) const
	{
		const unsigned int z = idx%m_N[2];
		idx /= m_N[2];
		const unsigned int y = idx%m_N[1];
		idx /= m_N[1];
		const unsigned int x = idx%m_N[0];
		const unsigned int n = idx/m_N[0];
		return m_RegionStart[n] + ((size_t)m_RegionPos[n][0][x]*region.lines[n][1].size() + m_RegionPos[n][1][y])*region.lines[n][2].size() + m_RegionPos[n][2][z];
	}
	//! The frequency domain values of the dump from the sums \a raw of the raw values of \a region
	void EvaluateSums(const std::complex<float>* raw, ArrayLib::ArrayNIJK<std::complex<float>> &field) const
	{
		const size_t lines = entries.size()/(3*(size_t)m_nk);   // (i, j) lines
		for (size_t l=0; l<lines; ++l)
			for (unsigned int k=0; k<m_nk; ++k)
				for (int n=0; n<3; ++n)
					field(n, l/m_nj, l%m_nj, k) = Value<std::complex<double>>(raw, entries[(l*m_nk + k)*3 + n]);
	}

	bool HField() const {return m_H;}
	unsigned int NumLinesJ() const {return m_nj;}
	unsigned int NumLinesK() const {return m_nk;}
	//! The snapshots hold the raw values of \a region from \a offset on, see PrepareSnapshotGather()
	void SetSnapshotEvaluated(size_t offset) {m_Evaluated = true; m_Offset = offset;}

	virtual void Evaluate(size_t line_start, size_t line_stop, ArrayLib::ArrayNIJK<float> &field, const float* src=NULL) const
	{
		// a snapshot of the raw values of the region (see SetSnapshotEvaluated())
		const bool raw = src && m_Evaluated;
		const float* f = raw ? src + m_Offset : (src ? src : (m_H ? m_Curr->data() : m_Volt->data()));
		for (size_t l=line_start; l<line_stop; ++l)
		{
			const unsigned int i = l/m_nj;
			const unsigned int j = l%m_nj;
			for (unsigned int k=0; k<m_nk; ++k)
				for (int n=0; n<3; ++n)
				{
					const Entry& e = entries[(l*m_nk + k)*3 + n];
					if (raw)
						field(n, i, j, k) = Value<double>(f, e);
					else
						field(n, i, j, k) = Combine<double>(e, [&](int c) -> double {return Raw<double>(f[e.idx[c]], e.delta[c]);});
				}
		}
	}

	//! Find \a region and its index tables, after the entries are set
	void SetRegion(const unsigned int N[3]);

protected:
	template <typename R, typename T>
	static inline R Raw(T value, double delta)
	{
		if (delta)
			return R(value)/delta;
		return R(0);
	}

	//! The value of entry \a e from raw(k), see Engine_GatherEntry
	template <typename R, typename F>
	static inline R Combine(const Entry& e, F raw)
	{
		R out = 0;
		switch (e.form)
		{
		case G::ZERO:
			break;
		case G::RAW:
			out = raw(0);
			break;
		case G::LERP:
			out = raw(0)*(1.0-e.rel) + raw(1)*e.rel;
			break;
		case G::AVG4:
			out = raw(0);
			out+= raw(1);
			out+= raw(2);
			out+= raw(3);
			out/=4;
			break;
		}
		return out;
	}

	unsigned int m_N[3];                          //!< mesh lines of the fields
	size_t m_RegionStart[3];                      //!< first raw value of component n in the region
	std::vector<unsigned int> m_RegionPos[3][3];  //!< [n][direction][line]: position of the line in region.lines[n][direction]

	const ArrayLib::ArrayNIJK<FDTD_FLOAT>* m_Volt;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>* m_Curr;
	bool m_H;
	unsigned int m_nj, m_nk;
	bool m_Evaluated;   //!< a snapshot holds the evaluated entries, see SetSnapshotEvaluated()
	size_t m_Offset;
};

#endif // ENGINE_FIELD_GATHER_H
