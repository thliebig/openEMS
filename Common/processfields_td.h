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

#include <atomic>

#include "processfields.h"

class ProcessFieldsTD : public ProcessFields
{
public:
	ProcessFieldsTD(Engine_Interface_Base* eng_if);
	virtual ~ProcessFieldsTD();

	virtual std::string GetProcessingName() const {return "time domain field dump";}

	virtual void InitProcess();

	virtual int Process();

	//! Close the dump file, so other readers see all data
	virtual void PostProcess();

	//! Wait for the dumps written in the background
	virtual void FinishAsync();

	//! Set the length of the filename timestep pad filled with zeros (default is 8)
	void SetPadLength(int val) {pad_length=val;};

protected:
	int pad_length;

	//! Write the dump of timestep \a ts, see Process()
	bool WriteHDF5(unsigned int ts, float time, ArrayLib::ArrayNIJK<float>* field, const float* src);
	std::atomic<bool> m_AsyncFailed;   //!< a background write failed
	bool m_Snapshots;   //!< the dumps may read field snapshots (see Engine_Interface_Base::PrepareSnapshotGather())
	bool m_AsyncUsed;
};

#endif // PROCESSFIELDS_TD_H
