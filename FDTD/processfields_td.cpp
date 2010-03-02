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

void ProcessFieldsTD::Process()
{
	stringstream ss;
	ss << std::setw( pad_length ) << std::setfill( '0' ) << Eng->numTS;

	string filename = filePattern + ss.str() + ".vtk";
	ofstream file(filename.c_str());
	if (file.is_open()==false) { cerr << "ProcessFieldsTD::Process: can't open file '" << filename << "' for writing... abort! " << endl; return;};

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
		DumpFieldArray2VTK(file,string("E-Field"),E_T,discDLines,numDLines);
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
					H_T[1][pos[0]][pos[1]][pos[2]] = Eng->curr[0][OpPos[0]][OpPos[1]][OpPos[2]] + Eng->curr[0][OpPos[0]][OpPos[1]+1][OpPos[2]];
					H_T[1][pos[0]][pos[1]][pos[2]] /= (2*delta*Op->gridDelta);
					//in z
					delta  = Op->discLines[2][OpPos[2]+1] - Op->discLines[2][OpPos[2]];
					H_T[2][pos[0]][pos[1]][pos[2]] = Eng->curr[0][OpPos[0]][OpPos[1]][OpPos[2]] + Eng->curr[0][OpPos[0]][OpPos[1]][OpPos[2]+1];
					H_T[2][pos[0]][pos[1]][pos[2]] /= (2*delta*Op->gridDelta);
				}
			}
		}
		DumpFieldArray2VTK(file,string("H-Field"),H_T,discDLines,numDLines);
		Delete_N_3DArray(H_T,numDLines);
		H_T = NULL;
	}
	file.close();
}
