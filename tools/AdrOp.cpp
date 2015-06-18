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

#include "AdrOp.h"

using namespace std;

AdrOp::AdrOp(unsigned int muiImax, unsigned int muiJmax, unsigned int muiKmax, unsigned int muiLmax)
{
	//error-handling...
	error = new ErrorMsg(9);
	if (error==NULL)
	{
		fprintf(stderr,"Memory allocation failed!! exiting...");
		exit(1);
	}
	error->SetMsg(1,"Adress Operator: Memory allocation failed!! exiting...");
	error->SetMsg(2,"Adress Operator: Invalid Adress requested!! exiting...");
	error->SetMsg(3,"Adress Operator: Invalid Position set!! exiting...");
	error->SetMsg(4,"Adress Operator: Invalid jump or passing end of iteration!! exiting...");
	error->SetMsg(5,"Adress Operator: 4D not yet implemented!! exiting...");
	error->SetMsg(6,"Adress Operator: Position not set!! exiting...");
	error->SetMsg(7,"Adress Operator: Cells not added to Adress Operator!! exiting...");
	error->SetMsg(8,"Adress Operator: Invalid Node!! exiting...");
	error->SetMsg(9,"Adress Operator: Grid invalid!! exiting...");

	//if (muiImax<0) muiImax=0;
	//if (muiJmax<0) muiJmax=0;
	//if (muiKmax<0) muiKmax=0;
	//if (muiLmax<0) muiLmax=0;

	uiDimension=0;
	if (muiImax>0) uiDimension++;
	else exit(-1);
	if (muiJmax>0) uiDimension++;
	else exit(-2);
	if (muiKmax>0) uiDimension++;
	if ( (muiLmax>0) && (muiKmax>0) ) uiDimension++;
//	cout << "\n-----Adress Operator created: Dimension: " << uiDimension << "----" <<endl;
	uiImax=muiImax;
	uiJmax=muiJmax;
	uiKmax=muiKmax;
	uiLmax=muiLmax=0;
	if (uiDimension==2) uiSize=uiImax*uiJmax;
	else if (uiDimension==3) uiSize=uiImax*uiJmax*uiKmax;
	else if (uiDimension==4) uiSize=uiImax*uiJmax*uiKmax*uiLmax;
	else uiSize=0;
	bPosSet=false;
	if (uiDimension==4) error->Error(5);
	iIshift=iJshift=iKshift=0;
	reflect=false;
	uiTypeOffset=0;
	clCellAdr=NULL;
	dGrid[0]=NULL;
	dGrid[1]=NULL;
	dGrid[2]=NULL;
	dGrid[3]=NULL;
	dDeltaUnit=1;
	bDebug=false;
}

AdrOp::AdrOp(AdrOp* origOP)
{
	clCellAdr=NULL;
	error=NULL; // has to be done!!!

	uiDimension=origOP->uiDimension;
	uiSize=origOP->uiSize;
	uiImax=origOP->uiImax;
	uiJmax=origOP->uiJmax;
	uiKmax=origOP->uiKmax;
	uiLmax=origOP->uiLmax;
	uiIpos=origOP->uiIpos;
	uiJpos=origOP->uiJpos;
	uiKpos=origOP->uiKpos;
	uiLpos=origOP->uiLpos;
	for (int ii=0; ii<4; ++ii) dGrid[ii]=origOP->dGrid[ii];
	dDeltaUnit=origOP->dDeltaUnit;
	iIshift=origOP->iIshift;
	iJshift=origOP->iJshift;
	iKshift=origOP->iKshift;
	for (int ii=0; ii<3; ++ii) iCellShift[ii]=origOP->iCellShift[ii];
	i=origOP->i;
	j=origOP->j;
	k=origOP->k;
	l=origOP->l;
	reflect=origOP->reflect;
	uiTypeOffset=origOP->uiTypeOffset;

	bPosSet=origOP->bPosSet;
	bDebug=origOP->bDebug;
//	return;
	if (origOP->clCellAdr!=NULL) clCellAdr= new AdrOp(origOP->clCellAdr);
}

AdrOp::~AdrOp()
{
//	cerr << "\n------Adress Operator deconstructed-----\n" << endl;
	delete error;
	error=NULL;
	delete clCellAdr;
	clCellAdr=NULL;
}

unsigned int AdrOp::SetPos(unsigned int muiIpos, unsigned int muiJpos, unsigned int muiKpos, unsigned int muiLpos)
{
	if (bDebug) fprintf(stderr,"AdrOp Debug:: SetPos(%d,%d,%d,%d) Max(%d,%d,%d,%d) \n",muiIpos,muiJpos,muiKpos,muiLpos,uiImax,uiJmax,uiKmax,uiLmax);
	bPosSet=false;
	if (muiIpos<uiImax) uiIpos=muiIpos;
	else error->Error(3);
	if (muiJpos<uiJmax) uiJpos=muiJpos;
	else error->Error(3);
	if ((muiKpos>=uiKmax) && (uiDimension>2)) error->Error(3);
	else if (uiDimension>2) uiKpos=muiKpos;
	if ((muiLpos>=uiLmax) && (uiDimension>3)) error->Error(3);
	else if (uiDimension>3) uiLpos=muiLpos;
	bPosSet=true;
//	cerr << "Position i:" << uiIpos << " j: " << uiJpos << " k: " << uiKpos << " l: " << 0 << " MAX: i:" << uiImax << " j: "  << uiJmax << " k: " << uiKmax << endl; //debug
	ADRESSEXPENSE(0,0,0,0,uiDimension+1,18)
	return GetPos();
}

bool  AdrOp::SetPosChecked(unsigned int muiIpos, unsigned int muiJpos, unsigned int muiKpos, unsigned int muiLpos)
{
	bPosSet=true;
	if (muiIpos<uiImax) uiIpos=muiIpos;
	else bPosSet=false;
	if (muiJpos<uiJmax) uiJpos=muiJpos;
	else bPosSet=false;
	if ((muiKpos>=uiKmax) && (uiDimension>2)) bPosSet=false;
	else if (uiDimension>2) uiKpos=muiKpos;
	if ((muiLpos>=uiLmax) && (uiDimension>3)) bPosSet=false;
	else if (uiDimension>3) uiLpos=muiLpos;
	ADRESSEXPENSE(0,0,0,0,uiDimension+1,18)
	return bPosSet;
}

void AdrOp::SetGrid(double *gridI,double *gridJ,double *gridK,double *gridL)
{
	dGrid[0]=gridI;
	dGrid[1]=gridJ;
	dGrid[2]=gridK;
	dGrid[3]=gridL;
	ADRESSEXPENSE(0,0,0,0,4,0)
}

bool AdrOp::CheckPos(unsigned int muiIpos, unsigned int muiJpos, unsigned int muiKpos, unsigned int muiLpos)
{
	bPosSet=true;
	if ((muiIpos>=uiImax)) bPosSet=false;
	if ((muiJpos>=uiJmax)) bPosSet=false;
	if ((muiKpos>=uiKmax) && (uiDimension>2)) bPosSet=false;
	if ((muiLpos>=uiLmax) && (uiDimension>3)) bPosSet=false;
	ADRESSEXPENSE(0,0,0,0,uiDimension+1,18)
	return bPosSet;
}

bool AdrOp::CheckRelativePos(int muiIrel,int muiJrel,int muiKrel, int muiLrel)
{
	bPosSet=true;
	if ((muiIrel+(int)uiIpos<0) || (muiIrel+(int)uiIpos>=(int)uiImax)) bPosSet=false;
	if ((muiJrel+(int)uiJpos<0) || (muiJrel+(int)uiJpos>=(int)uiJmax)) bPosSet=false;
	if (((muiKrel+(int)uiKpos<0) || (muiKrel+(int)uiKpos>=(int)uiKmax)) && (uiDimension>2)) bPosSet=false;
	if (((muiLrel+(int)uiLpos<0) || (muiLrel+(int)uiLpos>=(int)uiLmax)) && (uiDimension>3)) bPosSet=false;
	ADRESSEXPENSE(2*uiDimension,0,0,0,uiDimension+1,18)
	return bPosSet;
}

unsigned int AdrOp::GetPos(int muiIrel, int muiJrel, int muiKrel, int /*muiLrel*/)
{
	if (bPosSet==false) error->Error(6);
	if (reflect)
	{
#if EXPENSE_LOG==1
		if (muiIrel+(int)uiIpos<0) ADRESSEXPENSE(2,1,0,0,1,0)
			if (muiIrel+(int)uiIpos>(int)uiImax-1) ADRESSEXPENSE(4,1,0,0,1,0)
				if (muiJrel+(int)uiJpos<0) ADRESSEXPENSE(2,1,0,0,1,0)
					if (muiJrel+(int)uiJpos>(int)uiJmax-1) ADRESSEXPENSE(4,1,0,0,1,0)
						if (muiKrel+(int)uiKpos<0) ADRESSEXPENSE(2,1,0,0,1,0)
							if (muiKrel+(int)uiKpos>(int)uiKmax-1) ADRESSEXPENSE(4,1,0,0,1,0)
#endif

								if (muiIrel+(int)uiIpos<0) muiIrel=-2*uiIpos-muiIrel-uiTypeOffset;
		if (muiIrel+(int)uiIpos>(int)uiImax-1) muiIrel=2*(uiImax-1-uiIpos)-muiIrel+uiTypeOffset;
		if (muiJrel+(int)uiJpos<0) muiJrel=-2*uiJpos-muiJrel-uiTypeOffset;
		if (muiJrel+(int)uiJpos>(int)uiJmax-1) muiJrel=2*(uiJmax-1-uiJpos)-muiJrel+uiTypeOffset;
		if (muiKrel+(int)uiKpos<0) muiKrel=-2*uiKpos-muiKrel-uiTypeOffset;
		if (muiKrel+(int)uiKpos>(int)uiKmax-1) muiKrel=2*(uiKmax-1-uiKpos)-muiKrel+uiTypeOffset;
	}
	if (uiDimension==2)
	{
		ADRESSEXPENSE(7,1,0,0,0,7)
		if ( (muiIrel+uiIpos<uiImax) && (muiJrel+uiJpos<uiJmax) )
			return (muiIrel+uiIpos)+(muiJrel+uiJpos)*uiImax;
		else error->Error(2);
		return 0;
	}
	else if (uiDimension==3)
	{
		ADRESSEXPENSE(9,3,0,0,0,11)
		if ( (muiIrel+uiIpos<uiImax) && (muiJrel+uiJpos<uiJmax) && (muiKrel+uiKpos<uiKmax) )
			return (muiIrel+uiIpos) + (muiJrel+uiJpos)*uiImax + (muiKrel+uiKpos)*uiJmax*uiImax;
		else error->Error(2);
		return 0;
	}
	else return 0;
}

unsigned int AdrOp::GetPosFromNode(int ny, unsigned int uiNode)
{
	while (ny<0) ny+=uiDimension;
	ny=ny%uiDimension;
	unsigned int help=uiNode,i=0,j=0,k=0,l=0;
	i=help%uiImax;
	help=help/uiImax;
	j=help%uiJmax;
	help=help/uiJmax;
	if (uiKmax>0)
	{
		k=help%uiKmax;
		help=help/uiKmax;
		l=help;
	}
	if (!CheckPos(i,j,k,l)) error->Error(8);
	ADRESSEXPENSE(0,7,0,0,13,4)
	switch (ny)
	{
	case 0:
		{
			return i;
			break;
		}
	case 1:
		{
			return j;
			break;
		}
	case 2:
		{
			return k;
			break;
		}
	case 3:
		{
			return l;
			break;
		}
	}
	return 0;
}

double AdrOp::GetNodeVolume(unsigned int uiNode)
{
	for (unsigned int n=0; n<uiDimension; n++) if (dGrid[n]==NULL) error->Error(9);
	double dVol=1;
	unsigned int uiMax[4]={uiImax,uiJmax,uiKmax,uiLmax};
	unsigned int uiPos[4]={GetPosFromNode(0,uiNode),GetPosFromNode(1,uiNode),GetPosFromNode(2,uiNode),GetPosFromNode(3,uiNode)};
	for (unsigned int n=0; n<uiDimension; n++)
	{
		if ((uiPos[n]>0) && (uiPos[n]<uiMax[n]-1))
		{
			dVol*=0.5*dDeltaUnit*(dGrid[n][uiPos[n]+1]-dGrid[n][uiPos[n]-1]);
			ADRESSEXPENSE(4,0,1,3,0,4)
		}
		else if ((uiPos[n]==0) && (uiPos[n]<uiMax[n]-1))
		{
			dVol*=dDeltaUnit*(dGrid[n][uiPos[n]+1]-dGrid[n][uiPos[n]]);
			ADRESSEXPENSE(3,0,1,2,0,4)
		}
		else if ((uiPos[n]>0) && (uiPos[n]==uiMax[n]-1))
		{
			dVol*=dDeltaUnit*(dGrid[n][uiPos[n]]-dGrid[n][uiPos[n]-1]);
			ADRESSEXPENSE(3,0,1,2,0,4)
		}
	}
	return dVol;
}

double AdrOp::GetIndexWidth(int ny, int index)
{
	for (unsigned int n=0; n<uiDimension; n++) if (dGrid[n]==NULL) error->Error(9);
	double width=0;
	while (ny<0) ny+=uiDimension;
	ny=ny%uiDimension;
	unsigned int uiMax[4]={uiImax,uiJmax,uiKmax,uiLmax};
	if ((index>0) && (index<(int)uiMax[ny]-1)) width= (dGrid[ny][index+1]-dGrid[ny][index-1])/2*dDeltaUnit;
	else if ((index==0) && (index<(int)uiMax[ny]-1)) width=(dGrid[ny][index+1]-dGrid[ny][index])*dDeltaUnit;
	else if ((index>0) && (index==(int)uiMax[ny]-1)) width=(dGrid[ny][index]-dGrid[ny][index-1])*dDeltaUnit;
	else width= 0;
	return width;
}

double AdrOp::GetIndexCoord(int ny, int index)
{
	for (unsigned int n=0; n<uiDimension; n++) if (dGrid[n]==NULL) error->Error(9);
	while (ny<0) ny+=uiDimension;
	ny=ny%uiDimension;
	unsigned int uiMax[4]={uiImax,uiJmax,uiKmax,uiLmax};
	if (index<0) index=0;
	if (index>=(int)uiMax[ny]) index=uiMax[ny]-1;
	return dGrid[ny][index]*dDeltaUnit;
}

double AdrOp::GetIndexDelta(int ny, int index)
{
	if (index<0) return GetIndexCoord(ny, 0) - GetIndexCoord(ny, 1);
	unsigned int uiMax[4]={uiImax,uiJmax,uiKmax,uiLmax};
	if (index>=(int)uiMax[ny]-1) return GetIndexCoord(ny, (int)uiMax[ny]-2) - GetIndexCoord(ny, (int)uiMax[ny]-1);
	return GetIndexCoord(ny, index+1) - GetIndexCoord(ny, index);
}


unsigned int AdrOp::Shift(int ny, int step)
{
	if (bPosSet==false) error->Error(6);
	while (ny<0) ny+=uiDimension;
	ny=ny%uiDimension;
	switch (ny)
	{
	case 0:
		{
			iIshift=step;
//			if ((int)uiIpos+step<0) iIshift=-2*uiIpos-iIshift;
//			else if ((int)uiIpos+step>=(int)uiImax) iIshift=-1*iIshift+2*(uiImax-1-uiIpos);
			break;
		}
	case 1:
		{
			iJshift=step;
//			if ((int)uiJpos+iJshift<0) iJshift=-2*uiJpos-iJshift;
//			else if ((int)uiJpos+iJshift>=(int)uiJmax) iJshift=-1*iJshift+2*(uiJmax-1-uiJpos);
			break;
		}
	case 2:
		{
			iKshift=step;
//			if ((int)uiKpos+iKshift<0) iKshift=-2*uiKpos-iKshift;
//			else if ((int)uiKpos+iKshift>=(int)uiKmax) iKshift=-1*iKshift+2*(uiKmax-1-uiKpos);
			break;
		}
	}
	ADRESSEXPENSE(1,1,0,0,2,3)
	return GetPos(iIshift,iJshift,iKshift);
}

bool AdrOp::CheckShift(int ny, int step)
{
	while (ny<0) ny+=uiDimension;
	ny=ny%uiDimension;
	int shift[3]={0,0,0};
	shift[ny]=step;
	if (this->CheckPos(uiIpos+shift[0],uiJpos+shift[1],uiKpos+shift[2]))
	{
		this->Shift(ny,step);
		return true;
	}
	else return false;
}

unsigned int AdrOp::GetShiftedPos(int ny, int step)
{
	if ((ny<0) || (ny>2))
		return GetPos(iIshift,iJshift,iKshift);
	int pos[3] = {iIshift, iJshift, iKshift};
	pos[ny]+=step;
	return GetPos(pos[0],pos[1],pos[2]);
}

void AdrOp::ResetShift()
{
	iIshift=iJshift=iKshift=0;
}

unsigned int AdrOp::Iterate(int jump)
{
	if (abs(jump)>=(int)uiImax) error->Error(4);
	i=uiIpos+jump;
	if (i>=uiImax)
	{
		i=i-uiImax;
		j=uiJpos+1;
		if (j>=uiJmax)
		{
			j=0;
			if (uiDimension==3)
			{
				k=uiKpos+1;
				if (k>=uiKmax) k=0;
				uiKpos=k;
			}
		}
		uiIpos=i;
		uiJpos=j;
		return GetPos();
	}
	else
	{
		uiIpos=i;
		return GetPos();
	}
}

unsigned int AdrOp::GetSize()
{
	return uiSize;
}


void AdrOp::SetReflection2Node()
{
	reflect=true;
	uiTypeOffset=0;
}

void AdrOp::SetReflection2Cell()
{
	reflect=true;
	uiTypeOffset=1;
}

void AdrOp::SetReflectionOff()
{
	reflect=false;
}

AdrOp* AdrOp::AddCellAdrOp()
{
	if (clCellAdr!=NULL) return clCellAdr;
	if (uiDimension==3) clCellAdr = new AdrOp(uiImax-1,uiJmax-1,uiKmax-1);
	else if (uiDimension==2) clCellAdr = new AdrOp(uiImax-1,uiJmax-1);
	else clCellAdr=NULL;
	if (clCellAdr!=NULL)
	{
		clCellAdr->SetPos(0,0,0);
		clCellAdr->SetReflection2Cell();
	}
	iCellShift[0]=iCellShift[1]=iCellShift[2]=0;
	return clCellAdr;
}

AdrOp* AdrOp::DeleteCellAdrOp()
{
	delete clCellAdr;
	clCellAdr=NULL;
	return NULL;
}

unsigned int AdrOp::ShiftCell(int ny, int step)
{
	if (clCellAdr==NULL) error->Error(7);
	while (ny<0) ny+=uiDimension;
	ny=ny%uiDimension;
	iCellShift[ny]=step;
	ADRESSEXPENSE(3,1,0,0,1,2)
	return clCellAdr->GetPos(uiIpos+iCellShift[0],uiJpos+iCellShift[1],uiKpos+iCellShift[2]);
}

bool AdrOp::ShiftCellCheck(int ny, int step)
{
	return clCellAdr->CheckShift(ny,step);
}

void AdrOp::ResetCellShift()
{
	if (clCellAdr==NULL) error->Error(7);
	iCellShift[0]=iCellShift[1]=iCellShift[2]=0;
}

unsigned int AdrOp::GetCellPos(bool incShift)
{
	if (bPosSet==false) error->Error(6);
	if (clCellAdr==NULL) error->Error(7);
#if EXPENSE_LOG==1
	if (incShift) ADRESSEXPENSE(3,0,0,0,0,2)
#endif
		if (incShift) return clCellAdr->GetPos(uiIpos+iCellShift[0],uiJpos+iCellShift[1],uiKpos+iCellShift[2]);
		else return clCellAdr->GetPos(uiIpos,uiJpos,uiKpos);
}

unsigned int AdrOp::GetCellPos(int i, int j, int k)
{
	if (bPosSet==false) error->Error(6);
	return clCellAdr->GetPos(uiIpos+i,uiJpos+j,uiKpos+k);
}


double AdrOp::GetShiftCellVolume(int ny, int step)
{
	for (unsigned int n=0; n<uiDimension; n++) if (dGrid[n]==NULL) error->Error(9);
	int uiMax[4]={(int)uiImax-1,(int)uiJmax-1,(int)uiKmax-1,(int)uiLmax-1};
	while (ny<0) ny+=uiDimension;
	ny=ny%uiDimension;
	iCellShift[ny]=step;
	int uiPos[4]={(int)uiIpos+iCellShift[0],(int)uiJpos+iCellShift[1],(int)uiKpos+iCellShift[2]};
	double dVol=1;
	for (unsigned int n=0; n<uiDimension; ++n)
	{
		if (uiMax[n]>0)
		{
			while ((uiPos[n]<0) || (uiPos[n]>=uiMax[n]))
			{
				if (uiPos[n]<0) uiPos[n]=-1*uiPos[n]-1;
				else if (uiPos[n]>=uiMax[n]) uiPos[n]=uiMax[n]-(uiPos[n]-uiMax[n]+1);
			}
			dVol*=(dGrid[n][uiPos[n]+1]-dGrid[n][uiPos[n]])*dDeltaUnit;
		}
	}
	return dVol;
}


deltaAdrOp::deltaAdrOp(unsigned int max)
{
	uiMax=max;
}

deltaAdrOp::~deltaAdrOp()
{
}

void deltaAdrOp::SetMax(unsigned int max)
{
	uiMax=max;
}

unsigned int deltaAdrOp::GetAdr(int pos)
{
	if (uiMax==1) return 0;
	if (pos<0)	pos=pos*-1;
	else if (pos>(int)uiMax-1) pos=2*(uiMax-1)-pos+1;
	if ((pos<0) || (pos>(int)uiMax-1))
	{
		fprintf(stderr," Error exiting... ");
		getchar();
		exit(-1);
	}
	return pos;
}


