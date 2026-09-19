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

#include <atomic>

#include "processfields.h"
#include "tools/arraylib/array_nijk.h"

class ProcessFieldsFD : public ProcessFields
{
public:
	ProcessFieldsFD(Engine_Interface_Base* eng_if);
	virtual ~ProcessFieldsFD();

	virtual std::string GetProcessingName() const {return "frequency domain field dump";}

	virtual void InitProcess();

	virtual int Process();
	virtual void PostProcess();
	//! Wait for the sums of the samples taken from snapshots
	virtual void FinishAsync();

protected:
	virtual void DumpFDData();
	//! Add \a weights[n] times \a field to the frequency domain field n
	void AddSample(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& field, const std::vector<std::complex<float>>& weights);

	//! frequency domain field storage
	std::vector<ArrayLib::ArrayNIJK<std::complex<float>>*> m_FD_Fields;
	//! the sums kept by the engine instead (see Engine_Interface_Base::CreateFieldDFT()): -1: none, else its id
	int m_FieldDFT;
	//! else the samples may be summed from field snapshots in the background (see AsyncDumps)
	bool m_Snapshots;
	bool m_AsyncUsed;
	std::atomic<bool> m_AsyncFailed;
};

#endif // PROCESSFIELDS_FD_H
