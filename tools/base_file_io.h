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

#ifndef BASE_FILE_IO_H
#define BASE_FILE_IO_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>

//! Abstract base class for dumping scalar or vector field data
class Base_File_IO
{
public:
	virtual ~Base_File_IO();

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

	//! Set the mesh lines for the given mesh type.
	virtual void SetMeshLines(double const* const* lines, unsigned int const* count, double scaling=1) = 0 ;

	void SetNativeDump(bool val) {m_NativeDump=val;};

	//! Add a scalar field. \sa GetNumberOfFields \sa ClearAllFields
	virtual void AddScalarField(std::string fieldname, double const* const* const* field, unsigned int const* size) = 0;
	//! Add a scalar field. \sa GetNumberOfFields \sa ClearAllFields
	virtual void AddScalarField(std::string fieldname, float const* const* const* field, unsigned int const* size) = 0;
	//! Add a vector field. \sa GetNumberOfFields \sa ClearAllFields
	virtual void AddVectorField(std::string fieldname, double const* const* const* const* field, unsigned int const* size) = 0;
	//! Add a vector field. \sa GetNumberOfFields \sa ClearAllFields
	virtual void AddVectorField(std::string fieldname, float const* const* const* const* field, unsigned int const* size) = 0;

	//! Get the number of fields. \sa ClearAllFields
	virtual int GetNumberOfFields() const = 0;
	//! Clear all included fields. \sa GetNumberOfFields
	virtual void ClearAllFields() = 0;

	//! Get if timestep file series is active. \sa SetTimestepActive
	virtual bool GetTimestepActive() {return m_ActiveTS;}
	//! Set the timestep file series flag. \sa GetTimestepActive \sa SetTimestep
	virtual void SetTimestepActive(bool val) {m_ActiveTS = val;}
	//! Set the current timestep, this will set the timestep flag to true. \sa SetTimestepActive
	virtual void SetTimestep(unsigned int ts) {m_timestep=ts;SetTimestepActive(true);}

	virtual bool Write() = 0;

protected:
	Base_File_IO(std::string filename, int meshType=0);
	std::string m_filename;
	std::string m_header;

	bool m_ActiveTS;
	unsigned int m_timestep;

	int m_MeshType;
	bool m_NativeDump;

	bool m_AppendMode;
	bool m_Binary;
	bool m_Compress;
};


#endif // BASE_FILE_IO_H
