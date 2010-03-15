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

#ifndef OPENEMS_H
#define OPENEMS_H

class Operator;
class Engine;
class ProcessingArray;

class openEMS
{
public:
    openEMS();
	~openEMS();

	int SetupFDTD(const char* file);

	void RunFDTD();

	void Reset();

	void SetEnableDumps(bool val) {Enable_Dumps=val;}
	void SetEndCriteria(double val) {endCrit=val;}

	void DebugMaterial() {DebugMat=true;}

protected:
	//! Number of Timesteps
	int NrTS;
	bool Enable_Dumps;
	bool DebugMat;
	double endCrit;
	Operator* FDTD_Op;
	Engine* FDTD_Eng;
	ProcessingArray* PA;
};

#endif // OPENEMS_H
