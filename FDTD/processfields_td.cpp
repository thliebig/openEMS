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

void ProcessFieldsTD::DumpCellInterpol(string filename)
{
	FDTD_FLOAT**** volt = Eng->GetVoltages();
	FDTD_FLOAT**** curr = Eng->GetCurrents();

	if (m_DumpType==E_FIELD_DUMP)
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
					E_T[0][pos[0]][pos[1]][pos[2]] = volt[0][OpPos[0]][OpPos[1]][OpPos[2]] + volt[0][OpPos[0]][OpPos[1]+1][OpPos[2]] + volt[0][OpPos[0]][OpPos[1]][OpPos[2]+1] + volt[0][OpPos[0]][OpPos[1]+1][OpPos[2]+1];
					E_T[0][pos[0]][pos[1]][pos[2]] /= (4*delta*Op->gridDelta);
					//in y
					delta  = Op->discLines[1][OpPos[1]+1] - Op->discLines[1][OpPos[1]];
					E_T[1][pos[0]][pos[1]][pos[2]] = volt[1][OpPos[0]][OpPos[1]][OpPos[2]] + volt[1][OpPos[0]+1][OpPos[1]][OpPos[2]] + volt[1][OpPos[0]][OpPos[1]][OpPos[2]+1] + volt[1][OpPos[0]+1][OpPos[1]][OpPos[2]+1];
					E_T[1][pos[0]][pos[1]][pos[2]] /= (4*delta*Op->gridDelta);
					//in z
					delta  = Op->discLines[2][OpPos[2]+1] - Op->discLines[2][OpPos[2]];
					E_T[2][pos[0]][pos[1]][pos[2]] = volt[2][OpPos[0]][OpPos[1]][OpPos[2]] + volt[2][OpPos[0]][OpPos[1]+1][OpPos[2]] + volt[2][OpPos[0]+1][OpPos[1]][OpPos[2]] + volt[2][OpPos[0]+1][OpPos[1]+1][OpPos[2]];
					E_T[2][pos[0]][pos[1]][pos[2]] /= (4*delta*Op->gridDelta);
				}
			}
		}

		if (m_fileType==VTK_FILETYPE)
		{
			ofstream file(filename.c_str());
			if (file.is_open()==false) { cerr << "ProcessFieldsTD::Process: can't open file '" << filename << "' for writing... abort! " << endl;};
			DumpVectorArray2VTK(file,string("E-Field"),E_T,discDLines,numDLines);
			file.close();
		}
		else if (m_fileType==HDF5_FILETYPE)
		{
			stringstream ss;
			ss << std::setw( pad_length ) << std::setfill( '0' ) << Eng->GetNumberOfTimesteps();
			DumpVectorArray2HDF5(filename.c_str(),string( ss.str() ),E_T,numDLines);
		}
		else
			cerr << "ProcessFieldsTD::DumpCellInterpol: unknown File-Type" << endl;
		Delete_N_3DArray(E_T,numDLines);
		E_T = NULL;
	}

	if (m_DumpType==1)
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
					H_T[0][pos[0]][pos[1]][pos[2]] = curr[0][OpPos[0]][OpPos[1]][OpPos[2]] + curr[0][OpPos[0]+1][OpPos[1]][OpPos[2]];
					H_T[0][pos[0]][pos[1]][pos[2]] /= (2*delta*Op->gridDelta);
					//in y
					delta  = Op->discLines[1][OpPos[1]+1] - Op->discLines[1][OpPos[1]];
					H_T[1][pos[0]][pos[1]][pos[2]] = curr[1][OpPos[0]][OpPos[1]][OpPos[2]] + curr[1][OpPos[0]][OpPos[1]+1][OpPos[2]];
					H_T[1][pos[0]][pos[1]][pos[2]] /= (2*delta*Op->gridDelta);
					//in z
					delta  = Op->discLines[2][OpPos[2]+1] - Op->discLines[2][OpPos[2]];
					H_T[2][pos[0]][pos[1]][pos[2]] = curr[2][OpPos[0]][OpPos[1]][OpPos[2]] + curr[2][OpPos[0]][OpPos[1]][OpPos[2]+1];
					H_T[2][pos[0]][pos[1]][pos[2]] /= (2*delta*Op->gridDelta);
				}
			}
		}
		if (m_fileType==VTK_FILETYPE)
		{
			ofstream file(filename.c_str());
			if (file.is_open()==false) { cerr << "ProcessFieldsTD::Process: can't open file '" << filename << "' for writing... abort! " << endl;};
			DumpVectorArray2VTK(file,string("H-Field"),H_T,discDLines,numDLines);
			file.close();
		}
		else if (m_fileType==HDF5_FILETYPE)
		{
			stringstream ss;
			ss << std::setw( pad_length ) << std::setfill( '0' ) << Eng->GetNumberOfTimesteps();
			DumpVectorArray2HDF5(filename.c_str(),string( ss.str() ),H_T,numDLines);
		}
		else
			cerr << "ProcessFieldsTD::DumpCellInterpol: unknown File-Type" << endl;
		Delete_N_3DArray(H_T,numDLines);
		H_T = NULL;
	}
}

void ProcessFieldsTD::DumpNoInterpol(string filename)
{
	FDTD_FLOAT**** volt = Eng->GetVoltages();
	FDTD_FLOAT**** curr = Eng->GetCurrents();

	unsigned int pos[3];
	double delta[3];
	if (m_DumpType==E_FIELD_DUMP)
	{
		//create array
		FDTD_FLOAT**** E_T = Create_N_3DArray(numLines);
		for (pos[0]=0;pos[0]<numLines[0];++pos[0])
		{
			delta[0]=fabs(Op->MainOp->GetIndexDelta(0,pos[0]+start[0]));
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				delta[1]=fabs(Op->MainOp->GetIndexDelta(1,pos[1]+start[1]));
				for (pos[2]=0;pos[2]<numLines[2];++pos[2])
				{
					delta[2]=fabs(Op->MainOp->GetIndexDelta(2,pos[2]+start[2]));
					E_T[0][pos[0]][pos[1]][pos[2]] = volt[0][pos[0]+start[0]][pos[1]+start[1]][pos[2]+start[2]]/delta[0]/Op->gridDelta;
					E_T[1][pos[0]][pos[1]][pos[2]] = volt[1][pos[0]+start[0]][pos[1]+start[1]][pos[2]+start[2]]/delta[1]/Op->gridDelta;
					E_T[2][pos[0]][pos[1]][pos[2]] = volt[2][pos[0]+start[0]][pos[1]+start[1]][pos[2]+start[2]]/delta[2]/Op->gridDelta;
				}
			}
		}
		if (m_fileType==VTK_FILETYPE)
		{
			ofstream file(filename.c_str());
			if (file.is_open()==false) { cerr << "ProcessFieldsTD::Process: can't open file '" << filename << "' for writing... abort! " << endl;};
			DumpVectorArray2VTK(file,string("E-Field"),E_T,discLines,numLines);
			file.close();
		}
		else if (m_fileType==HDF5_FILETYPE)
		{
			stringstream ss;
			ss << std::setw( pad_length ) << std::setfill( '0' ) << Eng->GetNumberOfTimesteps();
			DumpVectorArray2HDF5(filename.c_str(),string( ss.str() ),E_T,numLines);
		}
		else
			cerr << "ProcessFieldsTD::DumpCellInterpol: unknown File-Type" << endl;

		Delete_N_3DArray(E_T,numLines);
		E_T = NULL;
	}

	if (m_DumpType==H_FIELD_DUMP)
	{
		//create array
		FDTD_FLOAT**** H_T = Create_N_3DArray(numLines);
		for (pos[0]=0;pos[0]<numLines[0];++pos[0])
		{
			delta[0]=fabs(Op->MainOp->GetIndexWidth(0,pos[0]+start[0]));
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				delta[1]=fabs(Op->MainOp->GetIndexWidth(1,pos[1]+start[1]));
				for (pos[2]=0;pos[2]<numLines[2];++pos[2])
				{
					delta[2]=fabs(Op->MainOp->GetIndexWidth(2,pos[2]+start[2]));
					//in x
					H_T[0][pos[0]][pos[1]][pos[2]] = curr[0][pos[0]+start[0]][pos[1]+start[1]][pos[2]+start[2]]/delta[0]/Op->gridDelta;
					H_T[1][pos[0]][pos[1]][pos[2]] = curr[1][pos[0]+start[0]][pos[1]+start[1]][pos[2]+start[2]]/delta[1]/Op->gridDelta;
					H_T[2][pos[0]][pos[1]][pos[2]] = curr[2][pos[0]+start[0]][pos[1]+start[1]][pos[2]+start[2]]/delta[2]/Op->gridDelta;
				}
			}
		}
		if (m_fileType==VTK_FILETYPE)
		{
			ofstream file(filename.c_str());
			if (file.is_open()==false) { cerr << "ProcessFieldsTD::Process: can't open file '" << filename << "' for writing... abort! " << endl;};
			DumpVectorArray2VTK(file,string("H-Field"),H_T,discLines,numLines);
			file.close();
		}
		else if (m_fileType==HDF5_FILETYPE)
		{
			stringstream ss;
			ss << std::setw( pad_length ) << std::setfill( '0' ) << Eng->GetNumberOfTimesteps();
			DumpVectorArray2HDF5(filename.c_str(),string( ss.str() ),H_T,numLines);
		}
		else
			cerr << "ProcessFieldsTD::DumpCellInterpol: unknown File-Type" << endl;

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
	ss << filePattern << std::setw( pad_length ) << std::setfill( '0' ) << Eng->GetNumberOfTimesteps();

	if (m_fileType==VTK_FILETYPE)
	{
		ss << ".vtk";
		if (m_DumpMode==NO_INTERPOLATION)
			DumpNoInterpol(ss.str());
		if (m_DumpMode==CELL_INTERPOLATE)
			DumpCellInterpol(ss.str());
	}
	else if (m_fileType==HDF5_FILETYPE)
	{
		if (m_DumpMode==NO_INTERPOLATION)
			DumpNoInterpol(m_fileName);
		if (m_DumpMode==CELL_INTERPOLATE)
			DumpCellInterpol(m_fileName);
	}
	return GetNextInterval();
}
