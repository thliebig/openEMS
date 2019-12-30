/*
*	Copyright (C) 2010-2015 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef __GNUC__ // not GCC
#include <emmintrin.h>
#endif

#include <complex>
typedef std::complex<double> double_complex;
#define _I double_complex(0.0,1.0)

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>

#define _USE_MATH_DEFINES

#include "Common/engine_interface_base.h"

class Operator_Base;

class Processing
{
public:
	virtual ~Processing();

	enum MeshType { CARTESIAN_MESH, CYLINDRICAL_MESH};

	//! Set the interface to the engine. Each processing needs its own engine interface. This class will take ownership and cleanup the interface on deletion!
	void SetEngineInterface(Engine_Interface_Base* eng_if);

	virtual void SetName(std::string val) {m_Name=val;}
	virtual void SetName(std::string val, int number);
	virtual std::string GetName() const {return m_Name;}

	//! Get the name for this processing, will be used in file description.
	virtual std::string GetProcessingName() const = 0;

	virtual void InitProcess() {};
	virtual void Reset();

	virtual void DefineStartStopCoord(double* dstart, double* dstop);

	virtual void ShowSnappedCoords();

	void SetProcessInterval(unsigned int interval);
	void SetProcessStartStopTime(double start, double stop);

	void AddStep(unsigned int step);
	void AddSteps(std::vector<unsigned int> steps);

	void AddFrequency(double freq);
	void AddFrequency(std::vector<double> *freqs);

	bool CheckTimestep();

	//! Process data prior to the simulation run.
	virtual void PreProcess() {};

	//! Process data during simulation run.
	virtual int Process() {return GetNextInterval();}

	//! Process data after simulation has finished.
	virtual void PostProcess();

	//! If disabled, Process() will do nothing...
	virtual void SetEnable(bool val) {Enabled=val;}
	//! If disabled, Process() will do nothing...
	virtual bool GetEnable() const {return Enabled;}

	virtual void SetWeight(double weight) {m_weight=weight;}
	virtual double GetWeight() {return m_weight;}

	//! Invoke this flag to flush all stored data to disk
	virtual void FlushNext() {m_Flush = true;}
	virtual void FlushData() {};

	void SetMeshType(MeshType meshType) {m_Mesh_Type=meshType;}

	//! Set the dump precision
	void SetPrecision(unsigned int val) {m_precision = val;}

	//! Dump probe geometry to file (will obay main or dual mesh property)
	virtual void DumpBox2File(std::string vtkfilenameprefix) const {DumpBox2File(vtkfilenameprefix,m_dualMesh);}

	//! Dump probe geometry to file
	virtual void DumpBox2File(std::string vtkfilenameprefix, bool dualMesh) const;

	virtual void SetDualMesh(bool val) {m_dualMesh=val;}
	virtual void SetDualTime(bool val) {m_dualTime=val;}

protected:
	Processing(Engine_Interface_Base* eng_if);
	Engine_Interface_Base* m_Eng_Interface;
	const Operator_Base* Op;
	MeshType m_Mesh_Type;

	unsigned int m_precision;

	std::string m_Name;

	bool m_Flush;

	double m_weight;

	bool Enabled;

	int GetNextInterval() const;
	unsigned int ProcessInterval;

	size_t m_PS_pos; //! current position in list of processing steps
	std::vector<unsigned int> m_ProcessSteps; //! list of processing steps

	//! Vector of frequency samples
	std::vector<double> m_FD_Samples;
	//! Number of samples already processed
	unsigned int m_FD_SampleCount;
	//! Sampling interval needed for the FD_Samples
	unsigned int m_FD_Interval;

	//! define if given coords are on main or dualMesh (default is false)
	bool m_dualMesh;

	//! define if given processing uses the dual time concept (default is false);
	bool m_dualTime;

	//! define the snap method used for this processing
	int m_SnapMethod;

	//! dimension of the snapped box
	int m_Dimension;

	//! define/store snapped start/stop coords as mesh index
	unsigned int start[3];
	unsigned int stop[3];

	//! start/stop timestep
	unsigned int startTS, stopTS;

	//! define/store if snapped start/stop coords are inside the field domain
	bool m_start_inside[3];
	bool m_stop_inside[3];

	std::ofstream file;
	std::string m_filename;

	virtual void OpenFile(std::string outfile);
};

class ProcessingArray
{
public:
	ProcessingArray(unsigned int maximalInterval) {maxInterval=maximalInterval;}
	~ProcessingArray() {};

	void AddProcessing(Processing* proc);

	void InitAll();

	//! Invoke this flag to flush all stored data to disk for all processings on next Process()
	void FlushNext();

	void Reset();

	//! Deletes all given processing's, can be helpful, but use carefull!!!
	void DeleteAll();

	//! Invoke PreProcess() on all Processings.
	void PreProcess();

	//! Invoke Process() on all Processings. Will return the smallest next iteration interval.
	int Process();

	//! Invoke PostProcess() on all Processings.
	void PostProcess();

	void DumpBoxes2File(std::string vtkfilenameprefix ) const;

	size_t GetNumberOfProcessings() const {return ProcessArray.size();}

	Processing* GetProcessing(size_t number) {return ProcessArray.at(number);}

protected:
	unsigned int maxInterval;
	std::vector<Processing*> ProcessArray;
};

#endif // PROCESSING_H
