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
#include "operator.h"

// The same values as GetEField()/GetHField() of Engine_Interface_FDTD for the given nodes;
// the caller has moved the nodes on the seam of a closed mesh as those do
// (see Engine_Interface_GPU_FDTD::CreateFieldGather()).
Field_Gather_FDTD* Field_Gather_FDTD::New(const Operator* op, Engine_Interface_Base::InterpolationType interpol,
                                          const ArrayLib::ArrayNIJK<FDTD_FLOAT>* volt, const ArrayLib::ArrayNIJK<FDTD_FLOAT>* curr,
                                          bool h_field, const unsigned int numLines[3], unsigned int* const posLines[3])
{
	unsigned int N[3];
	for (int n=0; n<3; ++n)
		N[n] = op->GetNumberOfLines(n, true);
	// the entries index the fields with 32 bit: larger grids dump without a gather
	if (3*(size_t)N[0]*N[1]*N[2] > UINT_MAX)
		return NULL;
	auto index = [&](int n, const unsigned int* p) -> unsigned int {return ((n*N[0] + p[0])*N[1] + p[1])*N[2] + p[2];};
	auto delta = [&](int n, const unsigned int* p) -> double {return op->GetEdgeLength(n, p, h_field);};

	typedef Engine_GatherEntry G;
	Field_Gather_FDTD* gather = new Field_Gather_FDTD(volt, curr, h_field, numLines[1], numLines[2]);
	gather->entries.resize((size_t)numLines[0]*numLines[1]*numLines[2]*3);
	size_t e_idx = 0;
	for (unsigned int i=0; i<numLines[0]; ++i)
		for (unsigned int j=0; j<numLines[1]; ++j)
			for (unsigned int k=0; k<numLines[2]; ++k)
			{
				const unsigned int pos[3] = {posLines[0][i], posLines[1][j], posLines[2][k]};
				for (int n=0; n<3; ++n)
				{
					G& e = gather->entries[e_idx++];
					e.form = G::ZERO;
					e.rel = 0;
					for (int m=0; m<4; ++m) {e.idx[m] = 0; e.delta[m] = 0;}
					unsigned int p[3] = {pos[0], pos[1], pos[2]};
					const int nP = (n+1)%3;
					const int nPP = (n+2)%3;
					auto set = [&](int m) {e.idx[m] = index(n, p); e.delta[m] = delta(n, p);};
					const bool upper = (pos[0]==N[0]-1) || (pos[1]==N[1]-1) || (pos[2]==N[2]-1);
					switch (interpol)
					{
					default:
					case Engine_Interface_Base::NO_INTERPOLATION:
						e.form = G::RAW;
						set(0);
						break;
					case Engine_Interface_Base::NODE_INTERPOLATE:
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
					case Engine_Interface_Base::CELL_INTERPOLATE:
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
