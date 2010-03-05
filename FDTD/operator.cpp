#include <fstream>
#include "operator.h"
#include "tools/array_ops.h"

Operator::Operator()
{
	Operator::Init();
}

Operator::~Operator()
{
	Operator::Reset();
}

void Operator::Init()
{
	CSX = NULL;

	ExciteSignal = NULL;
	E_Ex_delay = NULL;
	vv=NULL;
	vi=NULL;
	iv=NULL;
	ii=NULL;
	for (int n=0;n<3;++n)
	{
		discLines[n]=NULL;
		E_Ex_amp[n]=NULL;
		E_Ex_index[n]=NULL;
	}

	MainOp=NULL;
	DualOp=NULL;

	for (int n=0;n<3;++n)
	{
		EC_C[n]=NULL;
		EC_G[n]=NULL;
		EC_L[n]=NULL;
		EC_R[n]=NULL;
	}
}

void Operator::Reset()
{
	delete[] ExciteSignal;
	delete[] E_Ex_delay;
	Delete_N_3DArray(vv,numLines);
	Delete_N_3DArray(vi,numLines);
	Delete_N_3DArray(iv,numLines);
	Delete_N_3DArray(ii,numLines);
	for (int n=0;n<3;++n)
	{
		delete[] discLines[n];
		delete[] E_Ex_amp[n];
		delete[] E_Ex_index[n];
	}
	delete MainOp;
	delete DualOp;
	for (int n=0;n<3;++n)
	{
		delete[] EC_C[n];
		delete[] EC_G[n];
		delete[] EC_L[n];
		delete[] EC_R[n];
	}

	Operator::Init();
}

unsigned int Operator::GetNyquistNum(double fmax)
{
	if (dT==0) return 1;
	double T0 = 1/fmax;
	return floor(T0/2/dT);
}

bool Operator::SnapToMesh(double* dcoord, unsigned int* uicoord, bool lower)
{
	bool ok=true;
	for (int n=0;n<3;++n)
	{
		uicoord[n]=0;
		if (dcoord[n]<discLines[n][0]) {ok=false;uicoord[n]=0;}
		else if (dcoord[n]==discLines[n][0]) {uicoord[n]=0;}
		else if (dcoord[n]>discLines[n][numLines[n]-1]) {ok=false;uicoord[n]=numLines[n]-1; if (lower) uicoord[n]=numLines[n]-2;}
		else if (dcoord[n]==discLines[n][numLines[n]-1]) {uicoord[n]=numLines[n]-1; if (lower) uicoord[n]=numLines[n]-2;}
		else
			for (unsigned int i=1;i<numLines[n]-1;++i)
			{
				if (dcoord[n]<discLines[n][i])
				{
					if (fabs(dcoord[n]-discLines[n][i])<(fabs(dcoord[n]-discLines[n][i-1])))
						uicoord[n]=i;
					else
						uicoord[n]=i-1;
					if (lower) uicoord[n]=i-1;
					i = numLines[n];
				}
			}
	}
//	cerr << "Operator::SnapToMesh Wish: " << dcoord[0] << " " << dcoord[1] << " " << dcoord[2] << endl;
//	cerr << "Operator::SnapToMesh Found: " << discLines[0][uicoord[0]] << " " << discLines[1][uicoord[1]] << " " << discLines[2][uicoord[2]] << endl;
//	cerr << "Operator::SnapToMesh Index: " << uicoord[0] << " " << uicoord[1] << " " << uicoord[2] << endl;
	return ok;
}

double Operator::GetNumberCells()
{
	if (numLines)
		return (numLines[0]-1)*(numLines[1]-1)*(numLines[2]-1);
	return 0;
}

void Operator::ShowSize()
{
	unsigned int OpSize = 12*numLines[0]*numLines[1]*numLines[2]*sizeof(FDTD_FLOAT);
	unsigned int FieldSize = 6*numLines[0]*numLines[1]*numLines[2]*sizeof(FDTD_FLOAT);
	double MBdiff = 1024*1024;

	cout << "---- Stat: FDTD Operator ----" << endl;
	cout << "Dimensions        : " << numLines[0] << "x" << numLines[1] << "x" << numLines[2] << " = " <<  numLines[0]*numLines[1]*numLines[2] << " Cells (" << numLines[0]*numLines[1]*numLines[2]/1e6 << " MCells)" << endl;
	cout << "Size of Operator  : " << OpSize << " Byte (" << (double)OpSize/MBdiff << " MB) " << endl;
	cout << "Size of Field-Data: " << FieldSize << " Byte (" << (double)FieldSize/MBdiff << " MB) " << endl;
	cout << "-----------------------------" << endl;
}

void Operator::CalcGaussianPulsExcitation(double f0, double fc)
{
	if (dT==0) return;

	ExciteLength = (unsigned int)(2.0 * 9.0/(2.0*PI*fc) / dT);
	cerr << "Operator::CalcGaussianPulsExcitation: Length of the excite signal: " << ExciteLength << " timesteps" << endl;
	delete[] ExciteSignal;
	ExciteSignal = new FDTD_FLOAT[ExciteLength+1];
	ExciteSignal[0]=0.0;
	for (unsigned int n=1;n<ExciteLength+1;++n)
	{
		ExciteSignal[n] = cos(2.0*PI*f0*(n*dT-9.0/(2.0*PI*fc)))*exp(-1*pow(2.0*PI*fc*n*dT/3.0-3,2));
//		cerr << ExciteSignal[n] << endl;
	}
}

void Operator::DumpOperator2File(string filename)
{
	ofstream file(filename.c_str(),ios_base::out);
//	file.open;
	if (file.is_open()==false)
	{
		cerr << "Operator::DumpOperator2File: Can't open file: " << filename << endl;
		return;
	}
	file << "########### Operator vv ###########" << endl;
	file << "ix\tiy\tiz\tvv_x\tvv_y\tvv_z" << endl;
	Dump_N_3DArray2File(file,vv,numLines);

	file << "########### Operator vi ###########" << endl;
	file << "ix\tiy\tiz\tvi_x\tvi_y\tvi_z" << endl;
	Dump_N_3DArray2File(file,vi,numLines);

	file << "########### Operator iv ###########" << endl;
	file << "ix\tiy\tiz\tiv_x\tiv_y\tiv_z" << endl;
	Dump_N_3DArray2File(file,iv,numLines);

	file << "########### Operator ii ###########" << endl;
	file << "ix\tiy\tiz\tii_x\tii_y\tii_z" << endl;
	Dump_N_3DArray2File(file,ii,numLines);

	file.close();
}

void Operator::SetGeometryCSX(ContinuousStructure* geo)
{
	if (geo==NULL) return;

	Reset();
	CSX = geo;

	CSRectGrid* grid=CSX->GetGrid();
	for (int n=0;n<3;++n)
	{
		discLines[n] = grid->GetLines(n,discLines[n],numLines[n],true);
		if (numLines[n]<3) {cerr << "CartOperator::SetGeometryCSX: you need at least 3 disc-lines in every direction (3D!)!!!" << endl; Reset(); return;}
	}
	MainOp = new AdrOp(numLines[0],numLines[1],numLines[2]);
	MainOp->SetGrid(discLines[0],discLines[1],discLines[2]);
	if (grid->GetDeltaUnit()<=0)  {cerr << "CartOperator::SetGeometryCSX: grid delta unit must not be <=0 !!!" << endl; Reset(); return;}
	else gridDelta=grid->GetDeltaUnit();
	MainOp->SetGridDelta(1);
	MainOp->AddCellAdrOp();
}

void Operator::InitOperator()
{
	Delete_N_3DArray(vv,numLines);
	Delete_N_3DArray(vi,numLines);
	Delete_N_3DArray(iv,numLines);
	Delete_N_3DArray(ii,numLines);
	vv = Create_N_3DArray(numLines);
	vi = Create_N_3DArray(numLines);
	iv = Create_N_3DArray(numLines);
	ii = Create_N_3DArray(numLines);
}

int Operator::CalcECOperator()
{
	if (Calc_EC()==0)
		return -1;

	CalcTimestep();

	InitOperator();

	unsigned int i=0;
	unsigned int pos[3];

	for (int n=0;n<3;++n)
	{
		for (pos[0]=0;pos[0]<numLines[0];++pos[0])
		{
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				for (pos[2]=0;pos[2]<numLines[2];++pos[2])
				{
					i = MainOp->SetPos(pos[0],pos[1],pos[2]);
					vv[n][pos[0]][pos[1]][pos[2]] = (1-dT*EC_G[n][i]/2/EC_C[n][i])/(1+dT*EC_G[n][i]/2/EC_C[n][i]);
					vi[n][pos[0]][pos[1]][pos[2]] = (dT/EC_C[n][i])/(1+dT*EC_G[n][i]/2/EC_C[n][i]);

					ii[n][pos[0]][pos[1]][pos[2]] = (1-dT*EC_R[n][i]/2/EC_L[n][i])/(1+dT*EC_R[n][i]/2/EC_L[n][i]);
					iv[n][pos[0]][pos[1]][pos[2]] = (dT/EC_L[n][i])/(1+dT*EC_R[n][i]/2/EC_L[n][i]);
				}
			}
		}
	}

	//cleanup
	for (int n=0;n<3;++n)
	{
		delete[] EC_C[n];EC_C[n]=NULL;
		delete[] EC_G[n];EC_G[n]=NULL;
		delete[] EC_L[n];EC_L[n]=NULL;
		delete[] EC_R[n];EC_R[n]=NULL;
	}

	//Always apply PEC to all boundary's
	bool PEC[6]={1,1,1,1,1,1};
	ApplyElectricBC(PEC);

	if (CalcEFieldExcitation()==false) return -1;
	CalcPEC();

	return 0;
}

void Operator::ApplyElectricBC(bool* dirs)
{
	if (dirs==NULL) return;
	unsigned int pos[3];
	unsigned int ipos;
	for (int n=0;n<3;++n)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;
		for (pos[nP]=0;pos[nP]<numLines[nP];++pos[nP])
		{
			for (pos[nPP]=0;pos[nPP]<numLines[nPP];++pos[nPP])
			{
				pos[n]=0;
				vv[nP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n];
				vi[nP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n];
				vv[nPP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n];
				vi[nPP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n];
				pos[n]=numLines[n]-1;
				vv[nP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n+1];
				vi[nP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n+1];
				vv[nPP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n+1];
				vi[nPP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n+1];
			}
		}
	}
}

void Operator::ApplyMagneticBC(bool* dirs)
{
	if (dirs==NULL) return;
	unsigned int pos[3];
	unsigned int ipos;
	for (int n=0;n<3;++n)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;
		for (pos[nP]=0;pos[nP]<numLines[nP];++pos[nP])
		{
			for (pos[nPP]=0;pos[nPP]<numLines[nPP];++pos[nPP])
			{
				pos[n]=0;
				ii[nP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n];
				iv[nP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n];
				ii[nPP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n];
				iv[nPP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n];

				pos[n]=numLines[n]-2;
				ii[nP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n+1];
				iv[nP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n+1];
				ii[nPP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n+1];
				iv[nPP][pos[0]][pos[1]][pos[2]] *= (FDTD_FLOAT)!dirs[2*n+1];
			}
		}
	}
}


bool Operator::Calc_ECPos(int n, unsigned int* pos, double* inEC)
{
	double coord[3];
	double shiftCoord[3];
	int nP = (n+1)%3;
	int nPP = (n+2)%3;
	coord[0] = discLines[0][pos[0]];
	coord[1] = discLines[1][pos[1]];
	coord[2] = discLines[2][pos[2]];
	unsigned int ipos = MainOp->SetPos(pos[0],pos[1],pos[2]);
	double delta=MainOp->GetIndexDelta(n,pos[n]);
	double deltaP=MainOp->GetIndexDelta(nP,pos[nP]);
	double deltaPP=MainOp->GetIndexDelta(nPP,pos[nPP]);
	double delta_M=MainOp->GetIndexDelta(n,pos[n]-1);
	double deltaP_M=MainOp->GetIndexDelta(nP,pos[nP]-1);
	double deltaPP_M=MainOp->GetIndexDelta(nPP,pos[nPP]-1);

	//******************************* epsilon,kappa averaging *****************************//
	//shift up-right
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]+deltaP*0.25;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.25;
	CSProperties* prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[0] = mat->GetEpsilonWeighted(n,shiftCoord)*fabs(deltaP*deltaPP);
		inEC[1] = mat->GetKappaWeighted(n,shiftCoord)*fabs(deltaP*deltaPP);
	}
	else
	{
		inEC[0] = 1*fabs(deltaP*deltaPP);
		inEC[1] = 0;
	}
	//shift up-left
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]-deltaP_M*0.25;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.25;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[0] += mat->GetEpsilonWeighted(n,shiftCoord)*fabs(deltaP*deltaPP);
		inEC[1] += mat->GetKappaWeighted(n,shiftCoord)*fabs(deltaP*deltaPP);
	}
	else
	{
		inEC[0] += 1*fabs(deltaP*deltaPP);
		inEC[1] += 0;
	}

	//shift down-right
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]+deltaP*0.25;
	shiftCoord[nPP] = coord[nPP]-deltaPP_M*0.25;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[0] += mat->GetEpsilonWeighted(n,shiftCoord)*fabs(deltaP*deltaPP);
		inEC[1] += mat->GetKappaWeighted(n,shiftCoord)*fabs(deltaP*deltaPP);
	}
	else
	{
		inEC[0] += 1*fabs(deltaP*deltaPP);
		inEC[1] += 0;
	}

	//shift down-left
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]-deltaP_M*0.25;
	shiftCoord[nPP] = coord[nPP]-deltaPP_M*0.25;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[0] += mat->GetEpsilonWeighted(n,shiftCoord)*fabs(deltaP*deltaPP);
		inEC[1] += mat->GetKappaWeighted(n,shiftCoord)*fabs(deltaP*deltaPP);
	}
	else
	{
		inEC[0] += 1*fabs(deltaP*deltaPP);
		inEC[1] += 0;
	}

	inEC[0]*=gridDelta/fabs(delta)/4*__EPS0__;
	inEC[1]*=gridDelta/fabs(delta)/4;

	//******************************* mu,sigma averaging *****************************//
	//shift down
	shiftCoord[n] = coord[n]-delta_M*0.25;
	shiftCoord[nP] = coord[nP]+deltaP*0.5;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.5;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[2] = fabs(delta_M) / mat->GetMueWeighted(n,shiftCoord);
		if (mat->GetSigma(n))
			inEC[3] = fabs(delta_M) / mat->GetSigmaWeighted(n,shiftCoord);
		else
			inEC[3] = 0;
	}
	else
	{
		inEC[2] = fabs(delta_M);
		inEC[3] = 0;
	}
	//shift up
	shiftCoord[n] = coord[n]+delta*0.25;
	shiftCoord[nP] = coord[nP]+deltaP*0.5;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.5;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[2] += mat->GetMue(n)*fabs(delta);
		if (mat->GetSigmaWeighted(n,shiftCoord))
			inEC[3] += fabs(delta)/mat->GetSigmaWeighted(n,shiftCoord);
		else
			inEC[3] = 0;
	}
	else
	{
		inEC[2] += 1*fabs(delta);
		inEC[3] = 0;
	}

	inEC[2] = gridDelta * fabs(deltaP*deltaPP) * 2 * __MUE0__ / inEC[2];
	if (inEC[3]) inEC[3]=gridDelta*fabs(deltaP*deltaPP) * 2 / inEC[3];

	return true;
}

bool Operator::Calc_EffMatPos(int n, unsigned int* pos, double* inMat)
{
	int nP = (n+1)%3;
	int nPP = (n+2)%3;

	unsigned int ipos = MainOp->SetPos(pos[0],pos[1],pos[2]);
	double delta=MainOp->GetIndexDelta(n,pos[n]);
	double deltaP=MainOp->GetIndexDelta(nP,pos[nP]);
	double deltaPP=MainOp->GetIndexDelta(nPP,pos[nPP]);

	double delta_M=MainOp->GetIndexDelta(n,pos[n]-1);
	double deltaP_M=MainOp->GetIndexDelta(nP,pos[nP]-1);
	double deltaPP_M=MainOp->GetIndexDelta(nPP,pos[nPP]-1);

	this->Calc_ECPos(n,pos,inMat);

	inMat[0] *= (delta*delta)/MainOp->GetNodeVolume(ipos)/gridDelta;
	inMat[1] *= (delta*delta)/MainOp->GetNodeVolume(ipos)/gridDelta;

	inMat[2] *= 0.5*(fabs(delta_M) + fabs(delta)) / fabs(deltaP*deltaPP) / gridDelta;
	inMat[3] *= 0.5*(fabs(delta_M) + fabs(delta)) / fabs(deltaP*deltaPP) / gridDelta;

	return true;
}


bool Operator::Calc_EC()
{
	if (CSX==NULL) {cerr << "CartOperator::Calc_EC: CSX not given or invalid!!!" << endl; return false;}

	unsigned int ipos;
	unsigned int pos[3];
	double inEC[4];
	for (int n=0;n<3;++n)
	{
		//init x-cell-array
		delete[] EC_C[n];
		delete[] EC_G[n];
		delete[] EC_L[n];
		delete[] EC_R[n];
		EC_C[n] = new double[MainOp->GetSize()];
		EC_G[n] = new double[MainOp->GetSize()];
		EC_L[n] = new double[MainOp->GetSize()];
		EC_R[n] = new double[MainOp->GetSize()];
		for (unsigned int i=0;i<MainOp->GetSize();i++) //init all
		{
			EC_C[n][i]=0;
			EC_G[n][i]=0;
			EC_L[n][i]=0;
			EC_R[n][i]=0;
		}
		for (pos[2]=0;pos[2]<numLines[2];++pos[2])
		{
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				for (pos[0]=0;pos[0]<numLines[0];++pos[0])
				{
					Calc_ECPos(n,pos,inEC);
					ipos = MainOp->SetPos(pos[0],pos[1],pos[2]);
					EC_C[n][ipos]=inEC[0];
					EC_G[n][ipos]=inEC[1];
					EC_L[n][ipos]=inEC[2];
					EC_R[n][ipos]=inEC[3];
				}
			}
		}
	}

	return true;
}

double Operator::CalcTimestep()
{
	dT=1e200;
	double newT;
	unsigned int pos[3];
	unsigned int ipos;
	unsigned int ipos_PM;
	unsigned int ipos_PPM;
	MainOp->SetReflection2Cell();
	for (int n=0;n<3;++n)
	{
		int nP = (n+1)%3;
		int nPP = (n+2)%3;

		for (pos[2]=0;pos[2]<numLines[2];++pos[2])
		{
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				for (pos[0]=0;pos[0]<numLines[0];++pos[0])
				{
					ipos = MainOp->SetPos(pos[0],pos[1],pos[2]);
					ipos_PM = MainOp->Shift(nP,-1);
					MainOp->ResetShift();
					ipos_PPM= MainOp->Shift(nPP,-1);
					MainOp->ResetShift();
					newT = 2/sqrt( ( 4/EC_L[nP][ipos] + 4/EC_L[nP][ipos_PPM] + 4/EC_L[nPP][ipos] + 4/EC_L[nPP][ipos_PM]) / EC_C[n][ipos] );
					if (newT<dT) dT=newT;
				}
			}
		}
	}
//	cerr << "Timestep: " << dT << endl;
	return 0;
}

bool Operator::CalcEFieldExcitation()
{
	if (dT==0) return false;
	vector<unsigned int> vIndex[3];
	vector<FDTD_FLOAT> vExcit[3];
	vector<unsigned int> vDelay;
	unsigned int ipos;
	unsigned int pos[3];
	double coord[3];

	for (pos[2]=0;pos[2]<numLines[2];++pos[2])
	{
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			for (pos[0]=0;pos[0]<numLines[0];++pos[0])
			{
				coord[0] = discLines[0][pos[0]];
				coord[1] = discLines[1][pos[1]];
				coord[2] = discLines[2][pos[2]];
				CSProperties* prop = CSX->GetPropertyByCoordPriority(coord,CSProperties::ELECTRODE);
				if (prop)
				{
					CSPropElectrode* elec = prop->ToElectrode();
					if ((elec->GetExcitType()==0) || (elec->GetExcitType()==1)) //soft or hard E-Field excite!
					{
						vDelay.push_back((unsigned int)(elec->GetDelay()/dT));
						for (int n=0;n<3;++n)
						{
							vIndex[n].push_back(pos[n]);
							double delta=MainOp->GetIndexDelta(n,pos[n])*gridDelta;
							if (elec->GetActiveDir(n))
								vExcit[n].push_back(elec->GetWeightedExcitation(n,coord)*delta);
							else
								vExcit[n].push_back(0);
							if ((elec->GetExcitType()==1) && (elec->GetActiveDir(n))) //hard excite
							{
								vv[(n+1)%3][pos[0]][pos[1]][pos[2]] = 0;
								vi[(n+1)%3][pos[0]][pos[1]][pos[2]] = 0;
								vv[(n+2)%3][pos[0]][pos[1]][pos[2]] = 0;
								vi[(n+2)%3][pos[0]][pos[1]][pos[2]] = 0;
							}
						}
					}
				}
			}
		}
	}
	E_Ex_Count = vIndex[0].size();
	for (int n=0;n<3;++n)
	{
		delete[] E_Ex_index[n];
		E_Ex_index[n] = new unsigned int[E_Ex_Count];
		for (unsigned int i=0;i<E_Ex_Count;++i)
			E_Ex_index[n][i]=vIndex[n].at(i);
	}
	delete[] E_Ex_delay;
	E_Ex_delay = new unsigned int[E_Ex_Count];
	for (unsigned int i=0;i<E_Ex_Count;++i)
		E_Ex_delay[i]=vDelay.at(i);
	for (int n=0;n<3;++n)
	{
		delete[] E_Ex_amp[n];
		E_Ex_amp[n] = new FDTD_FLOAT[E_Ex_Count];
		for (unsigned int i=0;i<E_Ex_Count;++i)
			E_Ex_amp[n][i]=vExcit[n].at(i);
	}
	return true;
}

bool Operator::CalcPEC()
{
	unsigned int pos[3];
	double coord[3];
	double delta;

	for (int n=0;n<3;++n)
	{
		for (pos[2]=0;pos[2]<numLines[2];++pos[2])
		{
			coord[2] = discLines[2][pos[2]];
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				coord[1] = discLines[1][pos[1]];
				for (pos[0]=0;pos[0]<numLines[0];++pos[0])
				{
					coord[0] = discLines[0][pos[0]];
					MainOp->SetPos(pos[0],pos[1],pos[2]);
					delta=MainOp->GetIndexDelta(n,pos[n]);
					coord[n]= discLines[n][pos[n]] + delta*0.5;
					CSProperties* prop = CSX->GetPropertyByCoordPriority(coord, (CSProperties::PropertyType)(CSProperties::MATERIAL | CSProperties::METAL));
					if (prop)
					{
						if (prop->GetType()==CSProperties::METAL) //set to PEC
						{
							vv[n][pos[0]][pos[1]][pos[2]] = 0;
							vi[n][pos[0]][pos[1]][pos[2]] = 0;
//							cerr << "CartOperator::CalcPEC: PEC found at " << pos[0] << " ; "  << pos[1] << " ; " << pos[2] << endl;
						}
					}
				}
			}
		}
	}
}

