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

#ifndef PROCESSING_H
#define PROCESSING_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "operator.h"
#include "engine.h"

class Processing
{
public:
	Processing(Operator* op, Engine* eng);
	virtual ~Processing();

	virtual void SetName(string val) {m_Name=val;}

	virtual void InitProcess() {};
	virtual void Reset();

	virtual void DefineStartStopCoord(double* dstart, double* dstop);

	void SetProcessInterval(unsigned int interval) {ProcessInterval=max((unsigned int)1,interval);}

	void AddStep(unsigned int step);
	void AddSteps(vector<unsigned int> steps);

	void AddFrequency(double freq);
	void AddFrequency(vector<double>  freqs);

	bool CheckTimestep();
	virtual int Process() {return GetNextInterval();}

	//! If disabled, Process() will do nothing...
	virtual void SetEnable(bool val) {Enabled=val;}
	//! If disabled, Process() will do nothing...
	virtual bool GetEnable() const {return Enabled;}

	virtual void SetWeight(double weight) {m_weight=weight;}
	virtual double GetWeight() {return m_weight;}

	//! Invoke this flag to flush all stored data to disk
	virtual void FlushNext() {m_Flush = true;}
	virtual void FlushData() {};

	//! Set the dump precision
	void SetPrecision(unsigned int val) {m_precision = val;}

	//! Dump probe geometry to file (will obay main or dual mesh property)
	virtual void DumpBox2File( string vtkfilenameprefix) const {DumpBox2File(vtkfilenameprefix,m_dualMesh);}

	//! Dump probe geometry to file
	virtual void DumpBox2File( string vtkfilenameprefix, bool dualMesh) const;

protected:
	Operator* Op;
	Engine* Eng;
	unsigned int m_precision;

	string m_Name;

	bool m_Flush;

	double m_weight;

	bool Enabled;

	int GetNextInterval() const;
	unsigned int ProcessInterval;

	size_t m_PS_pos; //! current position in list of processing steps
	vector<unsigned int> m_ProcessSteps; //! list of processing steps

	//! Vector of frequency samples
	vector<double> m_FD_Samples;
	//! Number of samples already processed
	unsigned int m_FD_SampleCount;
	//! Sampling interval needed for the FD_Samples
	unsigned int m_FD_Interval;

	void Dump_FD_Data(vector<_Complex double> value, double factor, string filename);

	//! define if given coords are on main or dualMesh (default is false)
	bool m_dualMesh;

	//! define/store snapped start/stop coords as mesh index
	unsigned int start[3];
	unsigned int stop[3];

	//! define/store if snapped start/stop coords are inside the field domain
	bool m_start_inside[3];
	bool m_stop_inside[3];

	ofstream file;
	string m_filename;

	virtual void OpenFile(string outfile);

	double CalcLineIntegral(unsigned int* start, unsigned int* stop, int field) const;
	double CalcLineIntegral_I(unsigned int* start, unsigned int* stop) const;
	double CalcLineIntegral_V(unsigned int* start, unsigned int* stop) const;
};

class ProcessingArray
{
public:
	ProcessingArray(unsigned int maximalInterval) {maxInterval=maximalInterval;}
	~ProcessingArray() {};

	void AddProcessing(Processing* proc);

	//! Invoke this flag to flush all stored data to disk for all processings on next Process()
	void FlushNext();

	void Reset();

	//! Deletes all given processing's, can be helpful, but use carefull!!!
	void DeleteAll();

	int Process();

	void DumpBoxes2File( string vtkfilenameprefix ) const;

protected:
	unsigned int maxInterval;
	vector<Processing*> ProcessArray;
};

#endif // PROCESSING_H
