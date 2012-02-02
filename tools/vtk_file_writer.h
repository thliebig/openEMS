/*
*	Copyright (C) 2011,2012 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef VTK_FILE_WRITER_H
#define VTK_FILE_WRITER_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <complex>

class vtkDataSet;

class VTK_File_Writer
{
public:
	VTK_File_Writer(std::string filename, int meshType=0);
	virtual ~VTK_File_Writer();

	//! Set the filename
	virtual void SetFilename(std::string filename) {m_filename=filename;}
	//! Set the header information. May not be supported by all file types or setting.
	virtual void SetHeader(std::string header) {m_header=header;}

	//! Tell write to append data. May fail if filename has changed or filetype doesn't support this.
	virtual void SetAppendMode(bool val) {m_AppendMode=val;}
	//! Set binary flag (if the file type supports it)
	virtual void SetBinary(bool val) {m_Binary=val;}
	//! Set compression flag (if the file type supports it)
	virtual void SetCompress(bool val) {m_Compress=val;}

	void SetNativeDump(bool val) {m_NativeDump=val;}

	virtual void SetMeshLines(double const* const* lines, unsigned int const* count, double scaling=1);

	virtual void AddScalarField(std::string fieldname, double const* const* const* field);
	virtual void AddScalarField(std::string fieldname, float const* const* const* field);
	virtual void AddVectorField(std::string fieldname, double const* const* const* const* field);
	virtual void AddVectorField(std::string fieldname, float const* const* const* const* field);

	virtual int GetNumberOfFields() const;
	virtual void ClearAllFields();

	//! Get if timestep file series is active. \sa SetTimestepActive
	virtual bool GetTimestepActive() {return m_ActiveTS;}
	//! Set the timestep file series flag. \sa GetTimestepActive \sa SetTimestep
	virtual void SetTimestepActive(bool val) {m_ActiveTS = val;}
	//! Set the current timestep, this will set the timestep flag to true. \sa SetTimestepActive
	virtual void SetTimestep(unsigned int ts) {m_timestep=ts;SetTimestepActive(true);}

	virtual bool Write();

	virtual bool WriteASCII();
	virtual bool WriteXML();

protected:
	std::string m_filename;
	std::string m_header;

	//timestep properties
	bool m_ActiveTS;
	unsigned int m_timestep;

	vtkDataSet* m_GridData;

	//mesh information
	int m_MeshType;
	std::vector<double> m_MeshLines[3];
	bool m_NativeDump;

	bool m_AppendMode;
	bool m_Binary;
	bool m_Compress;

	virtual std::string GetTimestepFilename(int pad_length=10) const;
};

#endif // VTK_FILE_Writer_H
