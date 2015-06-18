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

#ifndef PROCESSFIELDS_FD_H
#define PROCESSFIELDS_FD_H

#include "processfields.h"

class ProcessFieldsFD : public ProcessFields
{
public:
	ProcessFieldsFD(Engine_Interface_Base* eng_if);
	virtual ~ProcessFieldsFD();

	virtual std::string GetProcessingName() const {return "frequency domain field dump";}

	virtual void InitProcess();

	virtual int Process();
	virtual void PostProcess();

protected:
	virtual void DumpFDData();

	//! frequency domain field storage
	std::vector<std::complex<float>****> m_FD_Fields;
};

#endif // PROCESSFIELDS_FD_H
