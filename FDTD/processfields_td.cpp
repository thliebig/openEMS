/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "processfields_td.h"
#include <iomanip>
#include <sstream>
#include <string>

ProcessFieldsTD::ProcessFieldsTD(Operator* op, Engine* eng) : ProcessFields(op, eng)
{
	pad_length = 8;
}

ProcessFieldsTD::~ProcessFieldsTD()
{
}

void ProcessFieldsTD::DumpCellInterpol(ofstream &file)
{
	if (DumpType==0)
	{
		//create array
		FDTD_FLOAT**** E_T = Create_N_3DArray(numDLines);
		unsigned int pos[3] = {start[0],start[1],start[2]};
		unsigned int OpPos[3];
		double delta;
//		cerr << "processing e-fields... " << endl;
		for (pos[0]=0;pos[0]<numDLines[0];++pos[0])
		{
			OpPos[0]=start[0]+pos[0];
			for (pos[1]=0;pos[1]<numDLines[1];++pos[1])
			{
				OpPos[1]=start[1]+pos[1];
				for (pos[2]=0;pos[2]<numDLines[2];++pos[2])
				{
					OpPos[2]=start[2]+pos[2];
					//in x
					delta  = Op->discLines[0][OpPos[0]+1] - Op->discLines[0][OpPos[0]];
					E_T[0][pos[0]][pos[1]][pos[2]] = Eng->volt[0][OpPos[0]][OpPos[1]][OpPos[2]] + Eng->volt[0][OpPos[0]][OpPos[1]+1][OpPos[2]] + Eng->volt[0][OpPos[0]][OpPos[1]][OpPos[2]+1] + Eng->volt[0][OpPos[0]][OpPos[1]+1][OpPos[2]+1];
					E_T[0][pos[0]][pos[1]][pos[2]] /= (4*delta*Op->gridDelta);
					//in y
					delta  = Op->discLines[1][OpPos[1]+1] - Op->discLines[1][OpPos[1]];
					E_T[1][pos[0]][pos[1]][pos[2]] = Eng->volt[1][OpPos[0]][OpPos[1]][OpPos[2]] + Eng->volt[1][OpPos[0]+1][OpPos[1]][OpPos[2]] + Eng->volt[1][OpPos[0]][OpPos[1]][OpPos[2]+1] + Eng->volt[1][OpPos[0]+1][OpPos[1]][OpPos[2]+1];
					E_T[1][pos[0]][pos[1]][pos[2]] /= (4*delta*Op->gridDelta);
					//in z
					delta  = Op->discLines[2][OpPos[2]+1] - Op->discLines[2][OpPos[2]];
					E_T[2][pos[0]][pos[1]][pos[2]] = Eng->volt[2][OpPos[0]][OpPos[1]][OpPos[2]] + Eng->volt[2][OpPos[0]][OpPos[1]+1][OpPos[2]] + Eng->volt[2][OpPos[0]+1][OpPos[1]][OpPos[2]] + Eng->volt[2][OpPos[0]+1][OpPos[1]+1][OpPos[2]];
					E_T[2][pos[0]][pos[1]][pos[2]] /= (4*delta*Op->gridDelta);
				}
			}
		}
		DumpVectorArray2VTK(file,string("E-Field"),E_T,discDLines,numDLines);
		Delete_N_3DArray(E_T,numDLines);
		E_T = NULL;
	}

	if (DumpType==1)
	{
		//create array
		FDTD_FLOAT**** H_T = Create_N_3DArray(numDLines);
		unsigned int pos[3] = {start[0],start[1],start[2]};
		unsigned int OpPos[3];
		double delta;
//		cerr << "processing h-fields... " << endl;
		for (pos[0]=0;pos[0]<numDLines[0];++pos[0])
		{
			OpPos[0]=start[0]+pos[0];
			for (pos[1]=0;pos[1]<numDLines[1];++pos[1])
			{
				OpPos[1]=start[1]+pos[1];
				for (pos[2]=0;pos[2]<numDLines[2];++pos[2])
				{
					OpPos[2]=start[2]+pos[2];
					//in x
					if (OpPos[0]==0) delta  = Op->discLines[0][OpPos[0]+1] - Op->discLines[0][OpPos[0]];
					else delta = 0.5* (Op->discLines[0][OpPos[0]+1] - Op->discLines[0][OpPos[0]-1]);
					H_T[0][pos[0]][pos[1]][pos[2]] = Eng->curr[0][OpPos[0]][OpPos[1]][OpPos[2]] + Eng->curr[0][OpPos[0]+1][OpPos[1]][OpPos[2]];
					H_T[0][pos[0]][pos[1]][pos[2]] /= (2*delta*Op->gridDelta);
					//in y
					delta  = Op->discLines[1][OpPos[1]+1] - Op->discLines[1][OpPos[1]];
					H_T[1][pos[0]][pos[1]][pos[2]] = Eng->curr[1][OpPos[0]][OpPos[1]][OpPos[2]] + Eng->curr[1][OpPos[0]][OpPos[1]+1][OpPos[2]];
					H_T[1][pos[0]][pos[1]][pos[2]] /= (2*delta*Op->gridDelta);
					//in z
					delta  = Op->discLines[2][OpPos[2]+1] - Op->discLines[2][OpPos[2]];
					H_T[2][pos[0]][pos[1]][pos[2]] = Eng->curr[2][OpPos[0]][OpPos[1]][OpPos[2]] + Eng->curr[2][OpPos[0]][OpPos[1]][OpPos[2]+1];
					H_T[2][pos[0]][pos[1]][pos[2]] /= (2*delta*Op->gridDelta);
				}
			}
		}
		DumpVectorArray2VTK(file,string("H-Field"),H_T,discDLines,numDLines);
		Delete_N_3DArray(H_T,numDLines);
		H_T = NULL;
	}
}

void ProcessFieldsTD::DumpNoInterpol(ofstream &file)
{
	if (DumpType==0)
	{
		//create array
		FDTD_FLOAT**** E_T = Create_N_3DArray(numLines);
		unsigned int pos[3];
		for (pos[0]=0;pos[0]<numLines[0];++pos[0])
		{
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				for (pos[2]=0;pos[2]<numLines[2];++pos[2])
				{
					E_T[0][pos[0]][pos[1]][pos[2]] = Eng->volt[0][pos[0]+start[0]][pos[1]+start[1]][pos[2]+start[2]];
					E_T[1][pos[0]][pos[1]][pos[2]] = Eng->volt[1][pos[0]+start[0]][pos[1]+start[1]][pos[2]+start[2]];
					E_T[2][pos[0]][pos[1]][pos[2]] = Eng->volt[2][pos[0]+start[0]][pos[1]+start[1]][pos[2]+start[2]];
				}
			}
		}
		DumpVectorArray2VTK(file,string("E-Field"),E_T,discLines,numLines);
		Delete_N_3DArray(E_T,numLines);
		E_T = NULL;
	}

	if (DumpType==1)
	{
		//create array
		FDTD_FLOAT**** H_T = Create_N_3DArray(numLines);
		unsigned int pos[3] = {start[0],start[1],start[2]};
		for (pos[0]=0;pos[0]<numLines[0];++pos[0])
		{
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				for (pos[2]=0;pos[2]<numLines[2];++pos[2])
				{
					//in x
					H_T[0][pos[0]][pos[1]][pos[2]] = Eng->curr[0][pos[0]+start[0]][pos[1]+start[1]][pos[2]+start[2]];
					H_T[1][pos[0]][pos[1]][pos[2]] = Eng->curr[1][pos[0]+start[0]][pos[1]+start[1]][pos[2]+start[2]];
					H_T[2][pos[0]][pos[1]][pos[2]] = Eng->curr[2][pos[0]+start[0]][pos[1]+start[1]][pos[2]+start[2]];
				}
			}
		}
		DumpVectorArray2VTK(file,string("H-Field"),H_T,discLines,numLines);
		Delete_N_3DArray(H_T,numLines);
		H_T = NULL;
	}
}

int ProcessFieldsTD::Process()
{
	if (Enabled==false) return -1;
	if (filePattern.empty()) return -1;
	if (CheckTimestep()==false) return GetNextInterval();
	stringstream ss;
	ss << std::setw( pad_length ) << std::setfill( '0' ) << Eng->numTS;

	string filename = filePattern + ss.str() + ".vtk";
	ofstream file(filename.c_str());
	if (file.is_open()==false) { cerr << "ProcessFieldsTD::Process: can't open file '" << filename << "' for writing... abort! " << endl; return GetNextInterval();};

	if (DumpMode==0)
		DumpNoInterpol(file);
	if (DumpMode==2)
		DumpCellInterpol(file);
	file.close();
	return GetNextInterval();
}
