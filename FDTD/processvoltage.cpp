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

#include "processvoltage.h"

ProcessVoltage::ProcessVoltage(Operator* op, Engine* eng) : Processing(op, eng)
{
}

ProcessVoltage::~ProcessVoltage()
{
	file.close();
}

void ProcessVoltage::OpenFile(string outfile)
{
	if (file.is_open()) file.close();

	file.open(outfile.c_str());
	if (file.is_open()==false)
	{
		cerr << "Can't open file: " << outfile << endl;
		return;
	}
}

int ProcessVoltage::Process()
{
	if (Enabled==false) return -1;
	if (CheckTimestep()==false) return GetNextInterval();
	FDTD_FLOAT voltage=CalcLineIntegral(start,stop,0);
//	cerr << voltage << endl;
	voltages.push_back(voltage);
	file << (double)Eng->numTS*Op->GetTimestep() << "\t" << voltage << endl;
	return GetNextInterval();
}
