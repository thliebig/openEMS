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

#include <algorithm>
#include <stdexcept>

#include "gpu_backend_reference.h"
#include "operator.h"

GPU_Backend_Reference::GPU_Backend_Reference()
{
	for (int n=0; n<3; ++n)
		numLines[n] = 0;
}

bool GPU_Backend_Reference::Init(const Operator* op)
{
	for (int n=0; n<3; ++n)
		numLines[n] = op->GetNumberOfLines(n, true);

	m_volt.Init("gpu_volt", numLines);
	m_curr.Init("gpu_curr", numLines);
	m_vv.Init("gpu_vv", numLines);
	m_vi.Init("gpu_vi", numLines);
	m_ii.Init("gpu_ii", numLines);
	m_iv.Init("gpu_iv", numLines);

	// upload the final operator coefficients, including all changes by operator extensions
	unsigned int pos[3];
	for (int n=0; n<3; ++n)
		for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
				for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
				{
					m_vv(n, pos[0], pos[1], pos[2]) = op->GetVV(n, pos[0], pos[1], pos[2]);
					m_vi(n, pos[0], pos[1], pos[2]) = op->GetVI(n, pos[0], pos[1], pos[2]);
					m_ii(n, pos[0], pos[1], pos[2]) = op->GetII(n, pos[0], pos[1], pos[2]);
					m_iv(n, pos[0], pos[1], pos[2]) = op->GetIV(n, pos[0], pos[1], pos[2]);
				}
	return true;
}

// operation by operation like Engine::UpdateVoltages(): reordering the float
// operations alone would break the bit-identical results
void GPU_Backend_Reference::UpdateVoltages()
{
	unsigned int pos[3];
	bool shift[3];

	// shift is 0 on the lower boundary, where the neighbour term cancels itself
	for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
	{
		shift[0]=pos[0];
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			shift[1]=pos[1];
			for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
			{
				shift[2]=pos[2];
				//for x
				m_volt(0, pos[0], pos[1], pos[2]) *= m_vv(0, pos[0], pos[1], pos[2]);
				m_volt(0, pos[0], pos[1], pos[2]) +=
				    m_vi(0, pos[0], pos[1], pos[2]) * (
				        m_curr(2, pos[0], pos[1]         , pos[2]         ) -
				        m_curr(2, pos[0], pos[1]-shift[1], pos[2]         ) -
				        m_curr(1, pos[0], pos[1]         , pos[2]         ) +
				        m_curr(1, pos[0], pos[1]         , pos[2]-shift[2])
				    );

				//for y
				m_volt(1, pos[0], pos[1], pos[2]) *= m_vv(1, pos[0], pos[1], pos[2]);
				m_volt(1, pos[0], pos[1], pos[2]) +=
				    m_vi(1, pos[0], pos[1], pos[2]) * (
				        m_curr(0, pos[0]         , pos[1], pos[2]         ) -
				        m_curr(0, pos[0]         , pos[1], pos[2]-shift[2]) -
				        m_curr(2, pos[0]         , pos[1], pos[2]         ) +
				        m_curr(2, pos[0]-shift[0], pos[1], pos[2]         )
				    );

				//for z
				m_volt(2, pos[0], pos[1], pos[2]) *= m_vv(2, pos[0], pos[1], pos[2]);
				m_volt(2, pos[0], pos[1], pos[2]) +=
				    m_vi(2, pos[0], pos[1], pos[2]) * (
				        m_curr(1, pos[0]         , pos[1]         , pos[2]) -
				        m_curr(1, pos[0]-shift[0], pos[1]         , pos[2]) -
				        m_curr(0, pos[0]         , pos[1]         , pos[2]) +
				        m_curr(0, pos[0]         , pos[1]-shift[1], pos[2])
				    );
			}
		}
	}
}

// see Engine::UpdateCurrents()
void GPU_Backend_Reference::UpdateCurrents()
{
	unsigned int pos[3];

	// the last line of each direction has no upper neighbour and stays untouched
	for (pos[0]=0; pos[0]<numLines[0]-1; ++pos[0])
	{
		for (pos[1]=0; pos[1]<numLines[1]-1; ++pos[1])
		{
			for (pos[2]=0; pos[2]<numLines[2]-1; ++pos[2])
			{
				//for x
				m_curr(0, pos[0], pos[1], pos[2]) *= m_ii(0, pos[0], pos[1], pos[2]);
				m_curr(0, pos[0], pos[1], pos[2]) +=
				    m_iv(0, pos[0], pos[1], pos[2]) * (
				        m_volt(2, pos[0], pos[1]  , pos[2]  ) -
				        m_volt(2, pos[0], pos[1]+1, pos[2]  ) -
				        m_volt(1, pos[0], pos[1]  , pos[2]  ) +
				        m_volt(1, pos[0], pos[1]  , pos[2]+1)
				    );

				//for y
				m_curr(1, pos[0], pos[1], pos[2]) *= m_ii(1, pos[0], pos[1], pos[2]);
				m_curr(1, pos[0], pos[1], pos[2]) +=
				    m_iv(1, pos[0], pos[1], pos[2]) * (
				        m_volt(0, pos[0]  , pos[1], pos[2]  ) -
				        m_volt(0, pos[0]  , pos[1], pos[2]+1) -
				        m_volt(2, pos[0]  , pos[1], pos[2]  ) +
				        m_volt(2, pos[0]+1, pos[1], pos[2]  )
				    );

				//for z
				m_curr(2, pos[0], pos[1], pos[2]) *= m_ii(2, pos[0], pos[1], pos[2]);
				m_curr(2, pos[0], pos[1], pos[2]) +=
				    m_iv(2, pos[0], pos[1], pos[2]) * (
				        m_volt(1, pos[0]  , pos[1]  , pos[2]) -
				        m_volt(1, pos[0]+1, pos[1]  , pos[2]) -
				        m_volt(0, pos[0]  , pos[1]  , pos[2]) +
				        m_volt(0, pos[0]  , pos[1]+1, pos[2])
				    );
			}
		}
	}
}

static void CopyField(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& src, ArrayLib::ArrayNIJK<FDTD_FLOAT>& dst)
{
	if (src.size()!=dst.size())
		throw std::runtime_error("GPU_Backend_Reference: host and device field size mismatch");
	std::copy_n(src.data(), src.size(), dst.data());
}

void GPU_Backend_Reference::DownloadVoltages(ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt)
{
	CopyField(m_volt, volt);
}

void GPU_Backend_Reference::DownloadCurrents(ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr)
{
	CopyField(m_curr, curr);
}

void GPU_Backend_Reference::UploadVoltages(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt)
{
	CopyField(volt, m_volt);
}

void GPU_Backend_Reference::UploadCurrents(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr)
{
	CopyField(curr, m_curr);
}
