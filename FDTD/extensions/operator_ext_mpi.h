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

#ifndef OPERATOR_EXT_MPI_H
#define OPERATOR_EXT_MPI_H

#include "FDTD/operator.h"
#include "operator_extension.h"

//! The halo exchange of an MPI domain decomposition (see MPI_Domain)
class Operator_Ext_MPI : public Operator_Extension
{
	friend class Engine_Ext_MPI;
public:
	Operator_Ext_MPI(Operator* op);
	virtual ~Operator_Ext_MPI();

	virtual bool BuildExtension();
	virtual Engine_Extension* CreateEngineExtention();

	virtual std::string GetExtensionName() const {return std::string("MPI halo exchange");}
	virtual void ShowStat(std::ostream &ostr) const;

protected:
	unsigned int m_numLines[3];
	int m_NeighborUp[3];
	int m_NeighborDown[3];
};

#endif // OPERATOR_EXT_MPI_H
