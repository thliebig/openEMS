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

#ifndef PROCESSFIELDS_H
#define PROCESSFIELDS_H

#include "processing.h"
#include "tools/array_ops.h"

#define __VTK_DATA_TYPE__ "double"

class VTK_File_Writer;
class HDF5_File_Writer;

class ProcessFields : public Processing
{
public:
	ProcessFields(Engine_Interface_Base* eng_if);
	virtual ~ProcessFields();

	//! File type definition.
	enum FileType { VTK_FILETYPE, HDF5_FILETYPE};

	//! Dump type definitions.
	/*!
	  Current dump types are electric field (E_FIELD_DUMP), magnetic field (H_FIELD_DUMP),
	  (conduction) electric current density (kappa*E) (J_FIELD_DUMP) and total current density (rotH)
	  */
	enum DumpType { E_FIELD_DUMP=0, H_FIELD_DUMP=1, J_FIELD_DUMP=2, ROTH_FIELD_DUMP=3, D_FIELD_DUMP=4, B_FIELD_DUMP=5, SAR_LOCAL_DUMP=20, SAR_1G_DUMP=21, SAR_10G_DUMP=22, SAR_RAW_DATA=29};

	virtual std::string GetProcessingName() const {return "common field processing";}

	virtual void InitProcess();

	virtual void DefineStartStopCoord(double* dstart, double* dstop);

	//! Define a field dump sub sampling rate for a given direction (default: \a dir = -1 means all directions)
	virtual void SetSubSampling(unsigned int subSampleRate, int dir=-1);

	//! Define a field dump optimal resolution for a given direction (default: \a dir = -1 means all directions)
	virtual void SetOptResolution(double optRes, int dir=-1);

	//! Set the filename for a hdf5 data group file (HDF5 FileType only) \sa SetFileType()
	void SetFileName(std::string fn) {m_filename=fn;}
	std::string SetFileName() const {return m_filename;}

	//! Define the Dump-Mode
	void SetDumpMode(Engine_Interface_Base::InterpolationType mode);
	//! This methode will dump all fields on a main cell node using 2 E-field and 4 H-fields per direction.
	void SetDumpMode2Node() {SetDumpMode(Engine_Interface_Base::NODE_INTERPOLATE);}
	//! This methode will dump all fields in the center of a main cell (dual-node) using 4 E-field and 2 H-fields per direction.
	void SetDumpMode2Cell() {SetDumpMode(Engine_Interface_Base::CELL_INTERPOLATE);}

	//! Set dump type: 0 for E-fields, 1 for H-fields, 2 for D-fields, 3 for B-fields, 4 for J-fields, etc...
	virtual void SetDumpType(DumpType type) {m_DumpType=type;}

	double CalcTotalEnergyEstimate() const;

	void SetFileType(FileType fileType) {m_fileType=fileType;}

	static std::string GetFieldNameByType(DumpType type);

	virtual bool NeedConductivity() const;
	virtual bool NeedPermittivity() const;
	virtual bool NeedPermeability() const;

protected:
	DumpType m_DumpType;
	FileType m_fileType;

	VTK_File_Writer* m_Vtk_Dump_File;
	HDF5_File_Writer* m_HDF5_Dump_File;

	enum SampleType {NONE, SUBSAMPLE, OPT_RESOLUTION} m_SampleType;
	virtual void CalcMeshPos();

	//! field dump sub-sampling (if enabled)
	unsigned int subSample[3];

	//! field dump optimal resolution (if enabled)
	double optResolution[3];

	//! dump mesh information
	unsigned int numLines[3];	//number of lines to dump
	unsigned int* posLines[3];	//grid positions to dump
	double* discLines[3];		//mesh disc lines to dump

	//! Calculate and return the defined field. Caller has to cleanup the array.
	FDTD_FLOAT**** CalcField();
};

#endif // PROCESSFIELDS_H
