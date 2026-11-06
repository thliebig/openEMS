/*
*	Copyright (C) 2026 Ivan Iudice
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

#include "operator_ext_mpi.h"
#include "engine_ext_mpi.h"
#include "tools/mpi_domain.h"

Operator_Ext_MPI::Operator_Ext_MPI(Operator* op) : Operator_Extension(op)
{
	for (int n=0; n<3; ++n)
	{
		m_numLines[n] = 0;
		m_NeighborUp[n] = -1;
		m_NeighborDown[n] = -1;
	}
}

Operator_Ext_MPI::~Operator_Ext_MPI()
{
}

bool Operator_Ext_MPI::BuildExtension()
{
	for (int n=0; n<3; ++n)
	{
		m_numLines[n] = m_Op->GetNumberOfLines(n, true);
		m_NeighborUp[n] = g_mpi.GetNeighbor(n, true);
		m_NeighborDown[n] = g_mpi.GetNeighbor(n, false);
	}
	return true;
}

Engine_Extension* Operator_Ext_MPI::CreateEngineExtention()
{
	Engine_Ext_MPI* eng_ext = new Engine_Ext_MPI(this);
	return eng_ext;
}

void Operator_Ext_MPI::ShowStat(std::ostream &ostr) const
{
	Operator_Extension::ShowStat(ostr);
	ostr << " Rank " << g_mpi.Rank() << " of " << g_mpi.Size() << ", neighbours (down/up):";
	for (int n=0; n<3; ++n)
		ostr << " " << m_NeighborDown[n] << "/" << m_NeighborUp[n];
	ostr << std::endl;
}
