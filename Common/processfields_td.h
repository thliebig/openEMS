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

#ifndef PROCESSFIELDS_TD_H
#define PROCESSFIELDS_TD_H

#include "processfields.h"

class ProcessFieldsTD : public ProcessFields
{
public:
	ProcessFieldsTD(Engine_Interface_Base* eng_if);
	virtual ~ProcessFieldsTD();

	virtual std::string GetProcessingName() const {return "time domain field dump";}

	virtual void InitProcess();

	virtual int Process();

	//! Set the length of the filename timestep pad filled with zeros (default is 8)
	void SetPadLength(int val) {pad_length=val;};

protected:
	int pad_length;
};

#endif // PROCESSFIELDS_TD_H
