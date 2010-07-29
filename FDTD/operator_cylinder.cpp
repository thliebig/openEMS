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

#include "engine.h"
#include "processfields.h"
#include "operator_cylinder.h"
#include "operator_extension.h"
#include "operator_ext_cylinder.h"

Operator_Cylinder* Operator_Cylinder::New(unsigned int numThreads)
{
	cout << "Create cylindrical FDTD operator" << endl;
	Operator_Cylinder* op = new Operator_Cylinder();
	op->setNumThreads(numThreads);
	op->Init();
	return op;
}

Operator_Cylinder::Operator_Cylinder() : __OP_CYLINDER_BASE_CLASS__()
{
	m_MeshType = ProcessFields::CYLINDRICAL_MESH;
}

Operator_Cylinder::~Operator_Cylinder()
{
	__OP_CYLINDER_BASE_CLASS__::Reset();
}

void Operator_Cylinder::Init()
{
	CC_closedAlpha = false;
	CC_R0_included = false;
	__OP_CYLINDER_BASE_CLASS__::Init();
}

void Operator_Cylinder::Reset()
{
	__OP_CYLINDER_BASE_CLASS__::Reset();
}

void Operator_Cylinder::InitOperator()
{
	__OP_CYLINDER_BASE_CLASS__::InitOperator();

	if (CC_closedAlpha || CC_R0_included)
		this->AddExtension(new Operator_Ext_Cylinder(this));
}

inline unsigned int Operator_Cylinder::GetNumberOfLines(int ny) const
{
	//this is necessary for a correct field processing... cylindrical engine has to reset this by adding +1
	if (CC_closedAlpha && ny==1)
		return numLines[1]-1;

	return numLines[ny];
}

string Operator_Cylinder::GetDirName(int ny) const
{
	if (ny==0) return "rho";
	if (ny==1) return "alpha";
	if (ny==2) return "z";
	return "";
}

double Operator_Cylinder::GetMeshDelta(int n, const int* pos, bool dualMesh) const
{
	double delta = __OP_CYLINDER_BASE_CLASS__::GetMeshDelta(n,pos,dualMesh);
	if (delta==0) return delta;
	if (n==1)
	{
		return delta * GetDiscLine(0,pos[0],dualMesh);
	}
	return delta;
}

double Operator_Cylinder::GetNodeArea(int ny, const int pos[3], bool dualMesh) const
{
	if (ny==2)
	{
		double da = __OP_CYLINDER_BASE_CLASS__::GetMeshDelta(1,pos,dualMesh)/gridDelta;
		double r1,r2;
		if (!dualMesh)
		{
			r1 = (discLines[0][pos[0]] - fabs(MainOp->GetIndexDelta(0,pos[0]-1))/2.0)*gridDelta;
			r2 = (discLines[0][pos[0]] + fabs(MainOp->GetIndexDelta(0,pos[0]  ))/2.0)*gridDelta;
		}
		else
		{
			r1 =  discLines[0][pos[0]]*gridDelta;
			r2 = (discLines[0][pos[0]] + fabs(MainOp->GetIndexDelta(0,pos[0])))*gridDelta;
		}
		if (r1<0)
			return da * pow(r2,2);
		return da/2* (pow(r2,2) - pow(r1,2));
	}
	return __OP_CYLINDER_BASE_CLASS__::GetNodeArea(ny,pos,dualMesh);
}

bool Operator_Cylinder::SetGeometryCSX(ContinuousStructure* geo)
{
	if (__OP_CYLINDER_BASE_CLASS__::SetGeometryCSX(geo)==false) return false;

	double minmaxA = fabs(discLines[1][numLines[1]-1]-discLines[1][0]);
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
	else
	{
		CC_closedAlpha=false;
	}

	if (discLines[0][0]<0)
		{cerr << "Operator_Cylinder::SetGeometryCSX: r<0 not allowed in Cylinder Coordinates!!!" << endl; Reset(); return false;}
	else if (discLines[0][0]==0.0)
	{
		cout << "Operator_Cylinder::SetGeometryCSX: r=0 included..." << endl;
		CC_R0_included= true;  //also needed for correct ec-calculation
	}

	return true;
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
		// E in alpha direction ( aka volt[1][x][y][z] ) is not defined for r==0 --> always zero...
		unsigned int pos[3] = {0,0,0};
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			for (pos[2]=0;pos[2]<numLines[2];++pos[2])
			{
				GetVV(1,pos[0],pos[1],pos[2]) = 0;
				GetVI(1,pos[0],pos[1],pos[2]) = 0;
			}
		}
	}
	__OP_CYLINDER_BASE_CLASS__::ApplyElectricBC(dirs);
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
		dirs[0]=0;  //no PMC in r_min directions...
	}
	__OP_CYLINDER_BASE_CLASS__::ApplyMagneticBC(dirs);
}

bool Operator_Cylinder::Calc_ECPos(int n, const unsigned int* pos, double* inEC) const
{
	double coord[3];
	double shiftCoord[3];
	int nP = (n+1)%3;
	int nPP = (n+2)%3;
	coord[0] = discLines[0][pos[0]];
	coord[1] = discLines[1][pos[1]];
	coord[2] = discLines[2][pos[2]];
	double delta=MainOp->GetIndexDelta(n,pos[n]);
	double deltaP=MainOp->GetIndexDelta(nP,pos[nP]);
	double deltaPP=MainOp->GetIndexDelta(nPP,pos[nPP]);
	double delta_M=MainOp->GetIndexDelta(n,pos[n]-1);
	double deltaP_M=MainOp->GetIndexDelta(nP,pos[nP]-1);
	double deltaPP_M=MainOp->GetIndexDelta(nPP,pos[nPP]-1);
	double geom_factor=0,A_n=0;

	//******************************* epsilon,kappa averaging *****************************//
	//shift up-right
	shiftCoord[n] = coord[n]+delta*0.5;
	shiftCoord[nP] = coord[nP]+deltaP*0.25;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.25;
	CSProperties* prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	switch (n)
	{
	case 0:
		geom_factor = fabs((deltaPP*deltaP/delta)*(coord[0]+fabs(delta)/2))*0.25;
		break;
	case 1:
		geom_factor = fabs(deltaP*deltaPP/(delta*coord[0]))*0.25;
		break;
	case 2:
		geom_factor = fabs((deltaPP/delta) * (pow(coord[0]+fabs(deltaP)/2.0,2.0) - pow(coord[0],2.0)))*0.25;
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
		geom_factor = fabs((deltaPP*deltaP_M/delta)*(coord[0]+fabs(delta)/2))*0.25;
		break;
	case 1:
		geom_factor = fabs(deltaP_M*deltaPP/(delta*coord[0]))*0.25;
		break;
	case 2:
		geom_factor = fabs((deltaPP/delta) * (pow(coord[0],2.0) - pow(coord[0]-fabs(deltaP_M)/2.0,2.0)))*0.25;
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
		geom_factor = fabs((deltaPP_M*deltaP/delta)*(coord[0]+fabs(delta)/2))*0.25;
		break;
	case 1:
		geom_factor = fabs(deltaP*deltaPP_M/(delta*coord[0]))*0.25;
		break;
	case 2:
		geom_factor = fabs((deltaPP_M/delta) * (pow(coord[0]+fabs(deltaP)/2.0,2.0) - pow(coord[0],2.0)))*0.25;
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
		geom_factor = fabs((deltaPP_M*deltaP_M/delta)*(coord[0]+fabs(delta)/2))*0.25;
		break;
	case 1:
		geom_factor = fabs(deltaP_M*deltaPP_M/(delta*coord[0]))*0.25;
		break;
	case 2:
		geom_factor = fabs((deltaPP_M/delta) * (pow(coord[0],2.0) - pow(coord[0]-fabs(deltaP_M)/2.0,2.0)))*0.25;
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

	//******************************* mu,sigma averaging *****************************//
	//shift down
	shiftCoord[n] = coord[n]-delta_M*0.25;
	shiftCoord[nP] = coord[nP]+deltaP*0.5;
	shiftCoord[nPP] = coord[nPP]+deltaPP*0.5;
	prop = CSX->GetPropertyByCoordPriority(shiftCoord,CSProperties::MATERIAL);
	double delta_n = fabs(delta_M);
	if (n==1)
	{
		delta_n = delta_n * fabs(coord[0]+0.5*fabs(deltaPP)); //multiply delta-angle by radius
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
		delta_n = delta_n * fabs(coord[0]+0.5*fabs(deltaPP)); //multiply delta-angle by radius
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
	if (n==0) //x-direction n==0 -> r; nP==1 -> alpha; nPP==2 -> z
	{
		A_n = A_n * coord[0];
	}
	if (n==2) //z-direction n==2 -> z; nP==0 -> r; nPP==1 -> alpha
	{
		A_n = fabs(deltaPP) * (pow(coord[0]+fabs(deltaP),2.0) - pow(coord[0],2.0))*0.5;
	}

	inEC[2] = gridDelta * A_n * 2 * __MUE0__ / inEC[2];
	if (inEC[3]) inEC[3]=gridDelta * A_n * 2 / inEC[3];

//	if ((n==1) && (pos[1]==0) && (pos[2]==0))
//		cerr <<  inEC[2]/(coord[0]) <<  endl;
//	cerr << n << " -> " <<  pos[0] << " " << pos[1] << " " << pos[2] << " " << inEC[2] << endl;

	return true;
}

bool Operator_Cylinder::Calc_EffMatPos(int n, const unsigned int* pos, double* inMat) const
{
	__OP_CYLINDER_BASE_CLASS__::Calc_EffMatPos(n, pos, inMat);

	// H_rho is not defined at position r==0
	if (CC_R0_included && (n==0) && (pos[0]==0))
	{
		inMat[2] = 0;
		inMat[3] = 0;
	}

	// E_alpha is not defined at position r==0
	if (CC_R0_included && (n==1) && (pos[0]==0))
	{
		inMat[0]=0;
		inMat[1]=0;
	}

	return true;
}

void Operator_Cylinder::AddExtension(Operator_Extension* op_ext)
{
	if (op_ext->IsCylinderCoordsSave())
		m_Op_exts.push_back(op_ext);
	else
		cerr << "Operator_Cylinder::AddExtension: Warning: Operator extension \"" << op_ext->GetExtensionName() << "\" is not compatible with cylinder-coords!! skipping...!" << endl;
}

double Operator_Cylinder::CalcTimestep()
{
	return CalcTimestep_Var1();
}
