/*
*	Copyright (C) 2011 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef PROCESSFIELDS_SAR_H
#define PROCESSFIELDS_SAR_H

#include "processfields_fd.h"

class ProcessFieldsSAR : public ProcessFieldsFD
{
public:
	ProcessFieldsSAR(Engine_Interface_Base* eng_if);
	virtual ~ProcessFieldsSAR();

	virtual string GetProcessingName() const {return "SAR dump";}

	virtual void InitProcess();

protected:
	virtual void DumpFDData();

	double GetKappaDensityRatio(const unsigned int* pos);
};

#endif // PROCESSFIELDS_SAR_H
