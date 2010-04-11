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

#include "operator_cylinder.h"

Operator_Cylinder* Operator_Cylinder::New()
{
	Operator_Cylinder* op = new Operator_Cylinder();
	op->Init();
	return op;
}

Operator_Cylinder::Operator_Cylinder()
{
}

Operator_Cylinder::~Operator_Cylinder()
{
	Operator::Reset();
}

void Operator_Cylinder::Init()
{
	CC_closedAlpha = false;
	CC_R0_included = false;
	Operator::Init();
}

void Operator_Cylinder::Reset()
{
	Operator::Reset();
}

inline unsigned int Operator_Cylinder::GetNumberOfLines(int ny) const
{
	//this is necessary for a correct field processing... cylindrical engine has to reset this by adding +1
	if (CC_closedAlpha && ny==1)
		return numLines[1]-1;

	return numLines[ny];
}


bool Operator_Cylinder::SetGeometryCSX(ContinuousStructure* geo)
{
	if (Operator::SetGeometryCSX(geo)==false) return false;

	double minmaxA = fabs(discLines[1][numLines[1]-1]-discLines[1][0]);
//			cerr << minmaxA -2*PI << " < " << (2*PI)/10/numLines[1] << endl;
	if (fabs(minmaxA-2*PI) < (2*PI)/10/numLines[1]) //check minmaxA smaller then a tenth of average alpha-width
	{
		cout << "Operator_Cylinder::SetGeometryCSX: Alpha is a full 2*PI => closed Cylinder..." << endl;
		CC_closedAlpha = true;
		discLines[1][numLines[1]-1] = discLines[1][0] + 2*PI;
		cerr << "Operator_Cylinder::SetGeometryCSX: Warning, not handling the disc-line width and material averaging correctly yet for a closed cylinder..." << endl;
		if (MainOp->GetIndexDelta(1,0)-MainOp->GetIndexDelta(1,numLines[1]-2) > (2*PI)/10/numLines[1])
		{
			cerr << "Operator_Cylinder::SetGeometryCSX: first and last angle delta must be the same... deviation to large..." << MainOp->GetIndexDelta(1,0) - MainOp->GetIndexDelta(1,numLines[1]-2) << endl;
			exit(1);
		}
		if (MainOp->GetIndexDelta(1,0)-MainOp->GetIndexDelta(1,numLines[1]-2) > 0)
		{
			cerr << "Operator_Cylinder::SetGeometryCSX: first and last angle delta must be the same... auto correction of deviation: " << MainOp->GetIndexDelta(1,0) - MainOp->GetIndexDelta(1,numLines[1]-2) << endl;
			discLines[1][numLines[1]-2] = discLines[1][numLines[1]-1]-MainOp->GetIndexDelta(1,0);
		}
	}
	else if (minmaxA>2*PI)
		{cerr << "Operator_Cylinder::SetGeometryCSX: Alpha Max-Min must not be larger than 2*PI!!!" << endl; Reset(); return false;}
	else CC_closedAlpha=false;

	if (discLines[0][0]<0)
		{cerr << "Operator_Cylinder::SetGeometryCSX: r<0 not allowed in Cylinder Coordinates!!!" << endl; Reset(); return false;}
	else if (discLines[0][0]==0.0)
	{
		cout << "Operator_Cylinder::SetGeometryCSX: r=0 included..." << endl;
		cerr << "Operator_Cylinder::SetGeometryCSX: r=0 included not yet implemented... exit... " << endl; exit(1);
		CC_R0_included=true;
	}

	return true;
}

inline void Operator_Cylinder::Calc_ECOperatorPos(int n, unsigned int* pos)
{
	unsigned int i = MainOp->SetPos(pos[0],pos[1],pos[2]);
	if (EC_C[n][i]>0)
	{
		vv[n][pos[0]][pos[1]][pos[2]] = (1-dT*EC_G[n][i]/2/EC_C[n][i])/(1+dT*EC_G[n][i]/2/EC_C[n][i]);
		vi[n][pos[0]][pos[1]][pos[2]] = (dT/EC_C[n][i])/(1+dT*EC_G[n][i]/2/EC_C[n][i]);
	}
	else
	{
		vv[n][pos[0]][pos[1]][pos[2]] = 0;
		vi[n][pos[0]][pos[1]][pos[2]] = 0;
	}

	if (EC_L[n][i]>0)
	{
		ii[n][pos[0]][pos[1]][pos[2]] = (1-dT*EC_R[n][i]/2/EC_L[n][i])/(1+dT*EC_R[n][i]/2/EC_L[n][i]);
		iv[n][pos[0]][pos[1]][pos[2]] = (dT/EC_L[n][i])/(1+dT*EC_R[n][i]/2/EC_L[n][i]);
	}
	else
	{
		ii[n][pos[0]][pos[1]][pos[2]] = 0;
		iv[n][pos[0]][pos[1]][pos[2]] = 0;
	}
}

void Operator_Cylinder::ApplyElectricBC(bool* dirs)
{
	if (dirs==NULL) return;
	if (CC_closedAlpha)
	{
		dirs[2]=0;dirs[3]=0; //no PEC in alpha directions...
	}
	if (CC_R0_included)
	{
		dirs[2]=0;  //no PEC in r_min directions...
	}

	Operator::ApplyElectricBC(dirs);
}

void Operator_Cylinder::ApplyMagneticBC(bool* dirs)
{
	if (dirs==NULL) return;
	if (CC_closedAlpha)
	{
		dirs[2]=0;dirs[3]=0; //no PMC in alpha directions...
	}
	if (CC_R0_included)
	{
		dirs[2]=0;  //no PMC in r_min directions...
	}
	Operator::ApplyMagneticBC(dirs);
}

bool Operator_Cylinder::Calc_ECPos(int n, unsigned int* pos, double* inEC)
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
	double geom_factor,A_n;

	//******************************* epsilon,kappa averaging *****************************//
	//shift up-right
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]+deltaP*0.25;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.25;
	CSProperties* prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	switch (n)
	{
	case 0:
		geom_factor = fabs(deltaPP*deltaP/delta)*(coord[0]+fabs(delta)/2)*0.25;
		break;
	case 1:
		geom_factor = fabs(deltaP*deltaPP/(delta*coord[0]))*0.25;
		break;
	case 2:
		geom_factor = fabs(deltaPP/delta) * (pow(coord[0]+fabs(deltaP)/2.0,2.0) - pow(coord[0],2.0))*0.25;
		break;
	}
	geom_factor*=gridDelta;
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[0] = mat->GetEpsilonWeighted(n,shiftCoord)*geom_factor*__EPS0__;
		inEC[1] = mat->GetKappaWeighted(n,shiftCoord)*geom_factor;
	}
	else
	{
		inEC[0] = 1*geom_factor*__EPS0__;
		inEC[1] = 0;
	}

	//shift up-left
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]-deltaP_M*0.25;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.25;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	switch (n)
	{
	case 0:
		geom_factor = fabs(deltaPP*deltaP_M/delta)*(coord[0]+fabs(delta)/2)*0.25;
		break;
	case 1:
		geom_factor = fabs(deltaP_M*deltaPP/(delta*coord[0]))*0.25;
		break;
	case 2:
		geom_factor = fabs(deltaPP/delta) * (pow(coord[0],2.0) - pow(coord[0]-fabs(deltaP_M)/2.0,2.0))*0.25;
		break;
	}
	geom_factor*=gridDelta;
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[0] += mat->GetEpsilonWeighted(n,shiftCoord)*geom_factor*__EPS0__;
		inEC[1] += mat->GetKappaWeighted(n,shiftCoord)*geom_factor;
	}
	else
	{
		inEC[0] += 1*geom_factor*__EPS0__;
		inEC[1] += 0;
	}

	//shift down-right
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]+deltaP*0.25;
	shiftCoord[nPP] = coord[nPP]-deltaPP_M*0.25;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	switch (n)
	{
	case 0:
		geom_factor = fabs(deltaPP_M*deltaP/delta)*(coord[0]+fabs(delta)/2)*0.25;
		break;
	case 1:
		geom_factor = fabs(deltaP*deltaPP_M/(delta*coord[0]))*0.25;
		break;
	case 2:
		geom_factor = fabs(deltaPP_M/delta) * (pow(coord[0]+fabs(deltaP)/2.0,2.0) - pow(coord[0],2.0))*0.25;
		break;
	}
	geom_factor*=gridDelta;
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[0] += mat->GetEpsilonWeighted(n,shiftCoord)*geom_factor*__EPS0__;
		inEC[1] += mat->GetKappaWeighted(n,shiftCoord)*geom_factor;
	}
	else
	{
		inEC[0] += 1*geom_factor*__EPS0__;
		inEC[1] += 0;
	}

	//shift down-left
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]-deltaP_M*0.25;
	shiftCoord[nPP] = coord[nPP]-deltaPP_M*0.25;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	switch (n)
	{
	case 0:
		geom_factor = fabs(deltaPP_M*deltaP_M/delta)*(coord[0]+fabs(delta)/2)*0.25;
		break;
	case 1:
		geom_factor = fabs(deltaP_M*deltaPP_M/(delta*coord[0]))*0.25;
		break;
	case 2:
		geom_factor = fabs(deltaPP_M/delta) * (pow(coord[0],2.0) - pow(coord[0]-fabs(deltaP_M)/2.0,2.0))*0.25;
		break;
	}
	geom_factor*=gridDelta;
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[0] += mat->GetEpsilonWeighted(n,shiftCoord)*geom_factor*__EPS0__;
		inEC[1] += mat->GetKappaWeighted(n,shiftCoord)*geom_factor;
	}
	else
	{
		inEC[0] += 1*geom_factor*__EPS0__;
		inEC[1] += 0;
	}

	if (CC_R0_included && (n>0) && (coord[0]==0))
	{
		inEC[0]=0;
		inEC[1]=0;
	}

//	if ((n==2) && (pos[1]==0) && (pos[2]==0))
//	cerr << n << " -> " <<  pos[0] << " " << pos[1] << " " << pos[2] << " " << inEC[0] << endl;

	//******************************* mu,sigma averaging *****************************//
	//shift down
	shiftCoord[n] = coord[n]-delta_M*0.25;
	shiftCoord[nP] = coord[nP]+deltaP*0.5;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.5;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	double delta_n = fabs(delta_M);
	if (n==1)
	{
		delta_n = delta_n * (coord[0]+0.5*fabs(deltaPP)); //multiply delta-angle by radius
	}
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[2] = delta_n / mat->GetMueWeighted(n,shiftCoord);
		if (mat->GetSigma(n))
			inEC[3] = delta_n / mat->GetSigmaWeighted(n,shiftCoord);
		else
			inEC[3] = 0;
	}
	else
	{
		inEC[2] = delta_n;
		inEC[3] = 0;
	}
	//shift up
	shiftCoord[n] = coord[n]+delta*0.25;
	shiftCoord[nP] = coord[nP]+deltaP*0.5;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.5;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	delta_n = fabs(delta);
	if (n==1)
	{
		delta_n = delta_n * (coord[0]+0.5*fabs(deltaPP)); //multiply delta-angle by radius
	}
	if (prop)
	{
		CSPropMaterial* mat = prop->ToMaterial();
		inEC[2] += mat->GetMue(n)*delta_n;
		if (mat->GetSigmaWeighted(n,shiftCoord))
			inEC[3] += delta_n/mat->GetSigmaWeighted(n,shiftCoord);
		else
			inEC[3] = 0;
	}
	else
	{
		inEC[2] += 1*delta_n;
		inEC[3] = 0;
	}

	A_n = fabs(deltaP*deltaPP);
	if (n==0) //z-direction n==0 -> r; nP==1 -> alpha; nPP==2 -> z
	{
		A_n = A_n * coord[0];
	}
	if (n==2) //z-direction n==2 -> z; nP==0 -> r; nPP==1 -> alpha
	{
		A_n = fabs(deltaPP) * (pow(coord[0]+fabs(deltaP),2.0) - pow(coord[0],2.0))*0.5;
	}

	inEC[2] = gridDelta * A_n * 2 * __MUE0__ / inEC[2];
	if (inEC[3]) inEC[3]=gridDelta * A_n * 2 / inEC[3];

//	if ((n==0) && (pos[1]==0) && (pos[2]==0))
//		cerr <<  inEC[2]/(coord[0]) <<  endl;
//	cerr << n << " -> " <<  pos[0] << " " << pos[1] << " " << pos[2] << " " << inEC[2] << endl;

	return true;
}

bool Operator_Cylinder::Calc_EffMatPos(int n, unsigned int* pos, double* inMat)
{
	return false; //fixme

//	int nP = (n+1)%3;
//	int nPP = (n+2)%3;
//
//	unsigned int ipos = MainOp->SetPos(pos[0],pos[1],pos[2]);
//	double delta=MainOp->GetIndexDelta(n,pos[n]);
//	double deltaP=MainOp->GetIndexDelta(nP,pos[nP]);
//	double deltaPP=MainOp->GetIndexDelta(nPP,pos[nPP]);
//
//	double delta_M=MainOp->GetIndexDelta(n,pos[n]-1);
//	double deltaP_M=MainOp->GetIndexDelta(nP,pos[nP]-1);
//	double deltaPP_M=MainOp->GetIndexDelta(nPP,pos[nPP]-1);
//
//	this->Calc_ECPos(n,pos,inMat);
//
//	inMat[0] *= (delta*delta)/MainOp->GetNodeVolume(ipos)/gridDelta;
//	inMat[1] *= (delta*delta)/MainOp->GetNodeVolume(ipos)/gridDelta;
//
//	inMat[2] *= 0.5*(fabs(delta_M) + fabs(delta)) / fabs(deltaP*deltaPP) / gridDelta;
//	inMat[3] *= 0.5*(fabs(delta_M) + fabs(delta)) / fabs(deltaP*deltaPP) / gridDelta;
//
//	return true;
}

