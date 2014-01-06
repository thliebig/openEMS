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

Engine_Interface_Cylindrical_FDTD::Engine_Interface_Cylindrical_FDTD(Operator_sse* op) : Engine_Interface_SSE_FDTD(op)
{
	m_Op_Cyl = dynamic_cast<Operator_Cylinder*>(op);
	if (m_Op_Cyl==NULL)
	{
		cerr << "Engine_Interface_Cylindrical_FDTD::Engine_Interface_Cylindrical_FDTD: Error: Operator is not a cylindrical operator! Exit!" << endl;
		exit(1);
	}
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
