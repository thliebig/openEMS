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

#include <complex>
#include <vector>
#include "Common/engine_interface_base.h"
#include "tools/constants.h"
#include "tools/arraylib/array_nijk.h"

class Engine_Interface_FDTD;
class Operator;

//! The raw voltages or currents a dump reads
/*!
  For each component n the edges on the mesh lines lines[n][0] x lines[n][1] x
  lines[n][2], component after component, in the order of the flat ArrayNIJK layout.
  A backend copies these raw values into a snapshot or sums them for a frequency domain
  dump (a few line tables on the device); the host interpolates them, see Field_Gather_FDTD.
  */
struct Engine_GatherRegion
{
	std::vector<unsigned int> lines[3][3];   //!< [n][direction], ascending

	size_t Size(int n) const {return lines[n][0].size()*lines[n][1].size()*lines[n][2].size();}
	size_t Size() const {return Size(0)+Size(1)+Size(2);}
};

//! GetEField()/GetHField() of Engine_Interface_FDTD at the nodes of a dump, from other sources of the raw values
/*!
  The fields of a dump of the GPU engine come from the host mirror, from a snapshot of the
  raw values the dump reads (its region), or, for a frequency domain dump, from the sums
  of these raw values. All are interpolated with Engine_Interface_FDTD::InterpolateField()
  and InterpolateDualField(), the code of the node by node evaluation, so the dumps from the
  host mirror and the snapshots are bit-identical to it. Nothing is stored per node.
  */
class Field_Gather_FDTD : public Engine_Field_Gather
{
public:
	//! The dump of the E field (\a h_field false) or H field of \a eng_if at the nodes posLines, reading \a volt or \a curr
	Field_Gather_FDTD(const Engine_Interface_FDTD* eng_if, const Operator* op,
	                  const ArrayLib::ArrayNIJK<FDTD_FLOAT>* volt, const ArrayLib::ArrayNIJK<FDTD_FLOAT>* curr,
	                  bool h_field, const unsigned int numLines[3], unsigned int* const posLines[3]);

	Engine_GatherRegion region;   //!< the raw values the dump reads

	bool HField() const {return m_H;}
	//! The snapshots hold the raw values of \a region from \a offset on, see PrepareSnapshotGather()
	void SetSnapshotRegion(size_t offset) {m_RegionSnapshot = true; m_Offset = offset;}

	virtual void Evaluate(size_t line_start, size_t line_stop, ArrayLib::ArrayNIJK<float> &field, const float* src=NULL) const;
	//! The frequency domain values of the dump from the sums \a raw of the raw values of \a region
	void EvaluateSums(const std::complex<float>* raw, ArrayLib::ArrayNIJK<std::complex<float>> &field) const;

protected:
	//! Interpolate node (i, j, k) of the dump from raw(n, pos), the value of edge n at pos divided by its (dual) length
	template <typename T, typename Raw>
	void Interpolate(unsigned int i, unsigned int j, unsigned int k, T out[3], Raw raw) const;
	//! the field value of edge n at pos (raw voltage or current) in \a f, which holds the whole field or \a region
	template <typename T>
	T Value(const T* f, bool in_region, unsigned int n, const unsigned int* pos) const
	{
		if (in_region)
			return f[m_RegionStart[n] + ((size_t)m_RegionPos[n][0][pos[0]]*region.lines[n][1].size() + m_RegionPos[n][1][pos[1]])*region.lines[n][2].size() + m_RegionPos[n][2][pos[2]]];
		return f[(((size_t)n*m_N[0] + pos[0])*m_N[1] + pos[1])*m_N[2] + pos[2]];
	}

	const Engine_Interface_FDTD* m_EngIf;
	const Operator* m_Op;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>* m_Volt;
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>* m_Curr;
	bool m_H;
	std::vector<unsigned int> m_PosLines[3];      //!< the dumped mesh lines
	unsigned int m_N[3];                          //!< mesh lines of the fields
	size_t m_RegionStart[3];                      //!< first raw value of component n in the region
	std::vector<unsigned int> m_RegionPos[3][3];  //!< [n][direction][line]: position of the line in region.lines[n][direction]
	bool m_RegionSnapshot;   //!< a snapshot holds the raw values of the region, see SetSnapshotRegion()
	size_t m_Offset;
};

#endif // ENGINE_FIELD_GATHER_H
