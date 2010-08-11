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

class ProcessFields : public Processing
{
public:
	ProcessFields(Operator* op, Engine* eng);
	virtual ~ProcessFields();

	enum FileType { VTK_FILETYPE, HDF5_FILETYPE};
	enum DumpType { E_FIELD_DUMP, H_FIELD_DUMP};
	enum DumpMode { NO_INTERPOLATION, NODE_INTERPOLATE, CELL_INTERPOLATE};

	virtual void InitProcess();

	virtual void DefineStartStopCoord(double* dstart, double* dstop);

	//! Define a field dump sub sampling rate for a given direction (default: \a dir = -1 means all directions)
	virtual void SetSubSampling(unsigned int subSampleRate, int dir=-1);

	//! Used file pattern e.g. pattern="tmp/efield_" --> "tmp/efield_000045.vtk" for timestep 45 or "tmp/efield_2.40000e9.vtk" for 2.4GHz E-field dump. (VTK FileType only) \sa SetFileType()
	void SetFilePattern(string fp) {m_filename=filePattern=fp;}

	//! Set the filename for a hdf5 data group file (HDF5 FileType only) \sa SetFileType()
	void SetFileName(string fn) {m_filename=fn;}

	//! Define the Dump-Mode
	void SetDumpMode(DumpMode mode) {m_DumpMode=mode;}
	//! This methode will dump all fields on a main cell node using 2 E-field and 4 H-fields per direction.
	void SetDumpMode2Node() {m_DumpMode=NODE_INTERPOLATE;}
	//! This methode will dump all fields in the center of a main cell (dual-node) using 4 E-field and 2 H-fields per direction.
	void SetDumpMode2Cell() {m_DumpMode=CELL_INTERPOLATE;}

	//! Set dump type: 0 for E-fields, 1 for H-fields, 2 for D-fields, 3 for B-fields, 4 for J-fields, etc...
	void SetDumpType(DumpType type) {m_DumpType=type;}

	//! Write a vtk header to an already open file with given mesh-type
	static void WriteVTKHeader(ofstream &file, double const* const* discLines, unsigned int const* numLines, unsigned int precision=12, string header_info = string(), MeshType meshT = CARTESIAN_MESH, double discLines_scaling = 1);
	//! Write a vtk header to an already open file (cartesian grid)
	static void WriteVTKCartesianGridHeader(ofstream &file, double const* const* discLines, unsigned int const* numLines, unsigned int precision=12, string header_info = string(), double discLines_scaling = 1);
	//! Write a vtk header to an already open file (cylindrical grid)
	static void WriteVTKCylindricalGridHeader(ofstream &file, double const* const* discLines, unsigned int const* numLines, unsigned int precision=12, string header_info = string(), double discLines_scaling = 1);
	//! Append a vtk vector array to an already open vtk file, write a header first! \sa WriteVTKHeader()
	static void WriteVTKVectorArray(ofstream &file, string name, FDTD_FLOAT const* const* const* const* array, unsigned int const* numLines, unsigned int precision=12);
	//! Append a vtk scalar array to an already open vtk file, write a header first! \sa WriteVTKHeader()
	static void WriteVTKScalarArray(ofstream &file, string name, FDTD_FLOAT const* const* const* array, unsigned int const* numLines, unsigned int precision=12);

	static bool DumpVectorArray2VTK(ofstream &file, string name, FDTD_FLOAT const* const* const * const* array, double const* const* discLines, unsigned int const* numLines, unsigned int precision=12, string header_info = string(), MeshType meshT = CARTESIAN_MESH, double discLines_scaling = 1);
	static bool DumpMultiVectorArray2VTK(ofstream &file, string names[], FDTD_FLOAT const* const* const* const* const* array, unsigned int numFields, double const* const* discLines, unsigned int const* numLines, unsigned int precision=12, string header_info = string(), MeshType meshT = CARTESIAN_MESH, double discLines_scaling = 1);
	static bool DumpScalarArray2VTK(ofstream &file, string name, FDTD_FLOAT const* const* const* array, double const* const* discLines, unsigned int const* numLines, unsigned int precision=12, string header_info = string(), MeshType meshT = CARTESIAN_MESH, double discLines_scaling = 1);
	static bool DumpMultiScalarArray2VTK(ofstream &file, string names[], FDTD_FLOAT const* const* const* const* array, unsigned int numFields, double const* const* discLines, unsigned int const* numLines, unsigned int precision=12, string header_info = string(), MeshType meshT = CARTESIAN_MESH, double discLines_scaling = 1);

	static bool DumpVectorArray2HDF5(string filename, string name, FDTD_FLOAT const* const* const* const* array, unsigned int const* numLines, float time=0);

	double CalcTotalEnergy() const;

	void SetFileType(FileType fileType) {m_fileType=fileType;}

	static string GetFieldNameByType(DumpType type);
	static string GetInterpolationNameByType(DumpMode mode);

//	virtual void Process();
protected:
	DumpMode m_DumpMode;
	DumpType m_DumpType;
	string filePattern;
	FileType m_fileType;

	//! field dump sub-sampling
	unsigned int subSample[3];

	//! dump mesh
	unsigned int numLines[3];
	double* discLines[3];
	//! dual dump mesh
	unsigned int numDLines[3];
	double* discDLines[3];
};

#endif // PROCESSFIELDS_H
