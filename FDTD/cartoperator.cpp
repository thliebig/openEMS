#include "cartoperator.h"

CartOperator::CartOperator()
{
    Init();
}

CartOperator::~CartOperator()
{
	Reset();
}

void CartOperator::Init()
{
    CSX = NULL;
    MainOp=NULL;
    DualOp=NULL;
	E_Ex_index = NULL;
	E_Ex_delay = NULL;
    for (int n=0;n<3;++n)
    {
		discLines[n]=NULL;
        EC_C[n]=NULL;
        EC_G[n]=NULL;
        EC_L[n]=NULL;
        EC_R[n]=NULL;
		vv[n]=NULL;
		vi[n]=NULL;
		iv[n]=NULL;
		ii[n]=NULL;
		E_Ex_amp[n]=NULL;
	}
}

void CartOperator::Reset()
{
	delete[] E_Ex_index;
	delete[] E_Ex_delay;
	for (int n=0;n<3;++n)
	{
		delete[] EC_C[n];
		delete[] EC_G[n];
		delete[] EC_L[n];
		delete[] EC_R[n];
		delete[] vv[n];
		delete[] vi[n];
		delete[] iv[n];
		delete[] ii[n];
		delete[] E_Ex_amp[n];
	}
	delete MainOp;
	delete DualOp;
	Init();
}


void CartOperator::SetGeometryCSX(ContinuousStructure* geo)
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

int CartOperator::CalcECOperator()
{
    if (Calc_EC()==0)
        return -1;

	CalcTimestep();

	for (int n=0;n<3;++n)
	{
		delete[] vv[n];
		vv[n] = new FDTD_FLOAT[MainOp->GetSize()];
		delete[] vi[n];
		vi[n] = new FDTD_FLOAT[MainOp->GetSize()];
		delete[] iv[n];
		iv[n] = new FDTD_FLOAT[MainOp->GetSize()];
		delete[] ii[n];
		ii[n] = new FDTD_FLOAT[MainOp->GetSize()];

		for (unsigned int i=0;i<MainOp->GetSize();++i)
		{
			vv[n][i] = (1-dT*EC_G[n][i]/2/EC_C[n][i])/(1+dT*EC_G[n][i]/2/EC_C[n][i]);
			vi[n][i] = (dT/EC_C[n][i])/(1+dT*EC_G[n][i]/2/EC_C[n][i]);

			ii[n][i] = (1-dT*EC_R[n][i]/2/EC_L[n][i])/(1+dT*EC_R[n][i]/2/EC_L[n][i]);
			iv[n][i] = (dT/EC_L[n][i])/(1+dT*EC_R[n][i]/2/EC_L[n][i]);
//			cerr << iv[n][i] << endl;
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
	bool PEC[6]={0,0,0,0,0,0};
	ApplyElectricBC(PEC);

	if (CalcEFieldExcitation()==false) return -1;

    return 0;
}

void CartOperator::ApplyElectricBC(bool* dirs)
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
				ipos=MainOp->SetPos(pos[0],pos[1],pos[2]);
				vv[n][ipos] *= (FDTD_FLOAT)!dirs[2*n];
				vi[n][ipos] *= (FDTD_FLOAT)!dirs[2*n];

				pos[n]=numLines[n]-1;
				ipos=MainOp->SetPos(pos[0],pos[1],pos[2]);
				vv[n][ipos] *= (FDTD_FLOAT)!dirs[2*n+1];
				vi[n][ipos] *= (FDTD_FLOAT)!dirs[2*n+1];
			}
		}
	}
}

void CartOperator::ApplyMagneticBC(bool* dirs)
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
				ipos=MainOp->SetPos(pos[0],pos[1],pos[2]);
				ii[n][ipos] *= (FDTD_FLOAT)!dirs[2*n];
				iv[n][ipos] *= (FDTD_FLOAT)!dirs[2*n];

				pos[n]=numLines[n]-2;
				ipos=MainOp->SetPos(pos[0],pos[1],pos[2]);
				ii[n][ipos] *= (FDTD_FLOAT)!dirs[2*n+1];
				iv[n][ipos] *= (FDTD_FLOAT)!dirs[2*n+1];
			}
		}
	}
}


bool CartOperator::Calc_ECPos(int n, unsigned int* pos, double* inEC)
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
		inEC[0] = mat->GetEpsilon(n)*fabs(deltaP*deltaPP);
		inEC[1] = mat->GetKappa(n)*fabs(deltaP*deltaPP);
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
		inEC[0] += mat->GetEpsilon(n)*fabs(deltaP*deltaPP);
		inEC[1] += mat->GetKappa(n)*fabs(deltaP*deltaPP);
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
		inEC[0] += mat->GetEpsilon(n)*fabs(deltaP*deltaPP);
		inEC[1] += mat->GetKappa(n)*fabs(deltaP*deltaPP);
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
		inEC[0] += mat->GetEpsilon(n)*fabs(deltaP*deltaPP);
		inEC[1] += mat->GetKappa(n)*fabs(deltaP*deltaPP);
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
		inEC[2] = fabs(delta_M) / mat->GetMue(n);
		if (mat->GetSigma(n))
			inEC[3] = fabs(delta_M) / mat->GetSigma(n);
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
		if (mat->GetSigma(n))
			inEC[3] += fabs(delta)/mat->GetSigma(n);
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

bool CartOperator::Calc_EffMatPos(int n, unsigned int* pos, double* inMat)
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

}


bool CartOperator::Calc_EC()
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

double CartOperator::CalcTimestep()
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

bool CartOperator::CalcEFieldExcitation()
{
	vector<unsigned int> vIndex;
	vector<FDTD_FLOAT> vExcit[3];
	vector<FDTD_FLOAT> vDelay;
	unsigned int ipos;
	unsigned int pos[3];
	double coord[3];

	for (pos[2]=0;pos[2]<numLines[2];++pos[2])
	{
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			for (pos[0]=0;pos[0]<numLines[0];++pos[0])
			{
				ipos = MainOp->SetPos(pos[0],pos[1],pos[2]);
				coord[0] = discLines[0][pos[0]];
				coord[1] = discLines[1][pos[1]];
				coord[2] = discLines[2][pos[2]];
				CSProperties* prop = CSX->GetPropertyByCoordPriority(coord,CSProperties::ELECTRODE);
				if (prop)
				{
					CSPropElectrode* elec = prop->ToElectrode();
					if (elec->GetExcitType()==0)
					{
						vIndex.push_back(ipos);
						for (int n=0;n<3;++n)
						{
							double delta=MainOp->GetIndexDelta(n,pos[n])*gridDelta;
							vExcit[n].push_back(elec->GetWeightedExcitation(n,coord)*delta);
						}
					}
				}
			}
		}
	}
	E_Ex_Count = vIndex.size();
	delete[] E_Ex_index;
	E_Ex_index = new unsigned int[E_Ex_Count];
	delete[] E_Ex_delay;
	E_Ex_delay = new FDTD_FLOAT[E_Ex_Count];
	for (unsigned int i=0;i<E_Ex_Count;++i)
		E_Ex_delay[i]=vIndex.at(i);
	for (int n=0;n<3;++n)
	{
		delete[] E_Ex_amp[n];
		E_Ex_amp[n] = new FDTD_FLOAT[E_Ex_Count];
		for (unsigned int i=0;i<E_Ex_Count;++i)
			E_Ex_amp[n][i]=vExcit[n].at(i);
	}
	return true;
}
