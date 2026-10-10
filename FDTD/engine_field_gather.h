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

//! Precomputed GetRawInterpolatedField()/GetRawInterpolatedDualField() (type 0) of the
//! nodes of a dump, for fields in the basic engine layout (see Engine::GetBasicFieldArrays()).
//! Every node component is one of a few forms, evaluated with the same operations as there.
class Field_Gather_FDTD : public Engine_Field_Gather
{
public:
	//! one component of a node: raw(k) = value(idx[k])/delta[k] (0 if delta is 0), see GetRawField().
	//! The layout a backend evaluates on the device, see Engine::AddSnapshotGather().
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

	bool HField() const {return m_H;}
	unsigned int NumLinesJ() const {return m_nj;}
	unsigned int NumLinesK() const {return m_nk;}
	//! The snapshots hold the values of these entries from \a offset on, see PrepareSnapshotGather()
	void SetSnapshotEvaluated(size_t offset) {m_Evaluated = true; m_Offset = offset;}

	virtual void Evaluate(size_t line_start, size_t line_stop, ArrayLib::ArrayNIJK<float> &field, const float* src=NULL) const
	{
		// the device evaluated the entries into the snapshot, only unpack them
		if (src && m_Evaluated)
		{
			for (size_t l=line_start; l<line_stop; ++l)
			{
				const float* v = src + m_Offset + l*m_nk*3;
				const unsigned int i = l/m_nj;
				const unsigned int j = l%m_nj;
				for (unsigned int k=0; k<m_nk; ++k)
					for (int n=0; n<3; ++n)
						field(n, i, j, k) = v[k*3 + n];
			}
			return;
		}
		const float* f = src ? src : (m_H ? m_Curr->data() : m_Volt->data());
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
					case G::ZERO:
						break;
					case G::RAW:
						out = Raw(f, e, 0);
						break;
					case G::LERP:
						out = Raw(f, e, 0)*(1.0-e.rel) + Raw(f, e, 1)*e.rel;
						break;
					case G::AVG4:
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

	const ArrayLib::ArrayNIJK<FDTD_FLOAT>* m_Volt;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>* m_Curr;
	bool m_H;
	unsigned int m_nj, m_nk;
	bool m_Evaluated;   //!< a snapshot holds the evaluated entries, see SetSnapshotEvaluated()
	size_t m_Offset;
};

#endif // ENGINE_FIELD_GATHER_H
