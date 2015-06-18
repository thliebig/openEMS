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

/*!
\class AdrOp
\author Thorsten Liebig
\version $Revision: 1.10 $
\date $Date: 2006/10/29 18:50:44 $
*/

#ifndef ADROP_H
#define ADROP_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "ExpenseLog.h"
#include "ErrorMsg.h"

class AdrOp
{
public:
	///Constructor, define dimension/size here
	AdrOp(unsigned int muiImax, unsigned int muiYmax, unsigned int muiKmax=0, unsigned int muiLmax=0);
	///Copy-Constructor
	AdrOp(AdrOp* origOP);
	///Deconstructor
	virtual ~AdrOp();
	///Set the current n-dim position, get 1-dim array position as return value
	/*!A position has to be set or all other methodes will case error! \n The methode will exit with error message if invalid position is set! \sa ErrorMsg */
	unsigned int SetPos(unsigned int muiIpos, unsigned int muiJpos, unsigned int muiKpos=0, unsigned int muiLpos=0);

	bool SetPosChecked(unsigned int muiIpos, unsigned int muiJpos, unsigned int muiKpos=0, unsigned int muiLpos=0);

	void SetGrid(double *gridI,double *gridJ,double *gridK=NULL,double *gridL=NULL);
	void SetGridDelta(double delta) {this->dDeltaUnit=delta;};

	bool CheckPos(unsigned int muiIpos, unsigned int muiJpos, unsigned int muiKpos=0, unsigned int muiLpos=0);
	bool CheckRelativePos(int muiIrel=0,int muiJrel=0,int muiKrel=0, int muiLrel=0);
	///will return current 1-dim position, in addition to a relative n-dim shift
	/*!In case of crossing the boundaries, activate reflection or an error will be invoked\sa SetReflection2Node \sa SetReflection2Cell \sa SetReflectionOff */
	unsigned int GetPos(int muiIrel=0,int muiJrel=0,int muiKrel=0, int muiLrel=0);

	double GetNodeVolume(unsigned int uiNode);

	double GetIndexWidth(int ny, int index);
	double GetIndexCoord(int ny, int index);
	///Get the gird delta at the given  index of direction ny. (if index<0 return negative value as index=0 would give, if index>=max-1 returns negative value as index=max-2 would give)
	double GetIndexDelta(int ny, int index);

//	double GetCellVolume(unsigned int uiCell);

	unsigned int GetPosFromNode(int ny, unsigned int uiNode);
	///Set a shift in ny direction (e.g. 0 for i-direction)
	/*!Shift set by this methode will be ignored by methode GetPos*/
	unsigned int Shift(int ny, int step);
	///Set a checked shift in ny direction (e.g. 0 for i-direction)
	/*!Shift set by this methode will be ignored by methode GetPos*/
	bool CheckShift(int ny, int step);
	///Returns the current 1-dim position including shift by methode "Shift" + additional (transitory) shift
	unsigned int GetShiftedPos(int ny=-1, int step=0);
	///Reset shift set by "Shift"-methode
	void ResetShift();
	///Iterates through AdrOp; --- obsolete ---
	unsigned int Iterate(int jump=1);
	///Retruns size of array
	unsigned int GetSize();
	///Set mode to reflect by node
	/*!1D-example (6 nodes): \image html node_reflect.PNG order: 0,1,2,3,4,5,4,3,...*/
	void SetReflection2Node();
	///Set mode to reflect by cell
	/*!1D-example (5 cells): \image html cell_reflect.PNG order: 0,1,2,3,4,4,3,...*/
	void SetReflection2Cell();
	///Deactivate reflection (default)
	void SetReflectionOff();
	///Add a cell adress operator (dimensions: i-1 j-1 k-1 l-1)
	/*!\image html cells_nodes.png */
	AdrOp* AddCellAdrOp();

	AdrOp* GetCellAdrOp() {return clCellAdr;};

	///Deconstructed cell adress operator if no longer needed
	AdrOp* DeleteCellAdrOp();
	///Shift cell in ny dircetion; cell reflection is active
	unsigned int ShiftCell(int ny, int step);
	///Shift cell in ny dircetion; cell reflection is active; return check
	bool ShiftCellCheck(int ny, int step);
	///Reset cell shift
	void ResetCellShift();
	///Get current cell position from cell adress operator
	unsigned int GetCellPos(bool incShift=true);
	///Get cell position from cell adress operator
	unsigned int GetCellPos(int i, int j, int k=0);
	//get volume of cell; incl shift
	double GetShiftCellVolume(int ny, int step);


	void SetDebugOn() {this->bDebug=true;};
	void SetDebugOff() {this->bDebug=false;};

protected:
	AdrOp *clCellAdr;
	unsigned int uiDimension;
	unsigned int uiSize;
	unsigned int uiImax,uiJmax,uiKmax,uiLmax;
	unsigned int uiIpos, uiJpos, uiKpos, uiLpos;
	double *dGrid[4];
	double dDeltaUnit;
	int iIshift, iJshift, iKshift;
	int iCellShift[3];
	unsigned int i,j,k,l;
	bool reflect;
	unsigned int uiTypeOffset;

	bool bPosSet;
	bool bDebug;
	ErrorMsg *error;
};


class deltaAdrOp
{
public:
	deltaAdrOp(unsigned int max=0);
	virtual ~deltaAdrOp();
	void SetMax(unsigned int max);
	unsigned int GetAdr(int pos);

protected:
	unsigned int uiMax;
};


#endif // ADROP_H
