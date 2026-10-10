/*
*	Copyright (C) 2011 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "engine_interface_cylindrical_fdtd.h"
#include <vector>

Engine_Interface_Cylindrical_FDTD::Engine_Interface_Cylindrical_FDTD(Operator_Cylinder* op) : Engine_Interface_FDTD(op)
{
	m_Op_Cyl = op;
}

Engine_Interface_Cylindrical_FDTD::~Engine_Interface_Cylindrical_FDTD()
{
}

double* Engine_Interface_Cylindrical_FDTD::GetHField(const unsigned int* pos, double* out) const
{
	if (m_Op_Cyl->GetClosedAlpha()==false)
		return Engine_Interface_FDTD::GetHField(pos, out);

	unsigned int iPos[] = {pos[0],pos[1],pos[2]};

	if ((m_InterpolType==CELL_INTERPOLATE) && (pos[1]==m_Op->GetNumberOfLines(1)))
		iPos[1]=0;

	if ((m_InterpolType==NODE_INTERPOLATE) && (iPos[1]==0))
		iPos[1]=m_Op->GetNumberOfLines(1);

	return Engine_Interface_FDTD::GetHField(iPos, out);
}

double* Engine_Interface_Cylindrical_FDTD::GetRawInterpolatedField(const unsigned int* pos, double* out, int type) const
{
	if (m_Op_Cyl->GetClosedAlpha()==false)
		return Engine_Interface_FDTD::GetRawInterpolatedField(pos,out,type);

	unsigned int iPos[] = {pos[0],pos[1],pos[2]};

	if ((m_InterpolType==NODE_INTERPOLATE) && (pos[1]==0))
		iPos[1]=m_Op->GetNumberOfLines(1);

	if ((m_InterpolType==CELL_INTERPOLATE) && (pos[1]==m_Op->GetNumberOfLines(1)))
		iPos[1]=0;

	return Engine_Interface_FDTD::GetRawInterpolatedField(iPos,out,type);
}

Engine_Field_Gather* Engine_Interface_Cylindrical_FDTD::CreateFieldGather(bool h_field, const unsigned int numLines[3], unsigned int* const posLines[3]) const
{
	if (m_Op_Cyl->GetClosedAlpha()==false)
		return Engine_Interface_FDTD::CreateFieldGather(h_field, numLines, posLines);

	// the nodes on the alpha seam are evaluated on its other side, see GetHField() and GetRawInterpolatedField()
	std::vector<unsigned int> alpha(posLines[1], posLines[1]+numLines[1]);
	for (size_t n=0; n<alpha.size(); ++n)
	{
		if ((m_InterpolType==NODE_INTERPOLATE) && (alpha[n]==0))
			alpha[n]=m_Op->GetNumberOfLines(1);
		else if ((m_InterpolType==CELL_INTERPOLATE) && (alpha[n]==m_Op->GetNumberOfLines(1)))
			alpha[n]=0;
	}
	unsigned int* const lines[3] = {posLines[0], alpha.data(), posLines[2]};
	return Engine_Interface_FDTD::CreateFieldGather(h_field, numLines, lines);
}
