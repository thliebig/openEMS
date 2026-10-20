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

#include <cstdint>
#include <string>
#include <vector>

#include "tools/arraylib/array_nijk.h"
#include "tools/arraylib/array_ijk.h"

//! Write fields on the mesh as VTK XML files (.vtr for a Cartesian, .vts for a cylindrical mesh),
//! the arrays zlib compressed, without the VTK library
class VTK_File_Writer
{
public:
	VTK_File_Writer(std::string filename, int meshType=0);
	virtual ~VTK_File_Writer();

	//! Set the filename, without extension
	virtual void SetFilename(std::string filename) {m_filename=filename;}
	//! Set the header information. Not written to the XML files.
	virtual void SetHeader(std::string header) {m_header=header;}

	void SetNativeDump(bool val) {m_NativeDump=val;}

	virtual void SetMeshLines(double const* const* lines, unsigned int const* count, double scaling=1);

	virtual void AddScalarField(std::string fieldname, ArrayLib::ArrayIJK<float> &field);
	virtual void AddVectorField(std::string fieldname, ArrayLib::ArrayNIJK<float>& field);

	virtual int GetNumberOfFields() const;
	virtual void ClearAllFields();

	//! Get if timestep file series is active. \sa SetTimestepActive
	virtual bool GetTimestepActive() {return m_ActiveTS;}
	//! Set the timestep file series flag. \sa GetTimestepActive \sa SetTimestep
	virtual void SetTimestepActive(bool val) {m_ActiveTS = val;}
	//! Set the current timestep, this will set the timestep flag to true. \sa SetTimestepActive
	virtual void SetTimestep(unsigned int ts) {m_timestep=ts;SetTimestepActive(true);}

	virtual bool Write();

	//! Write lines as a VTK XML poly data file (.vtp): points (x,y,z each) and lines (two point
	//! indices each); optional a value per line as cell data, and the metres per coordinate unit
	//! (if >0) as field data
	static bool WriteLines(std::string filename, const std::vector<float>& points, const std::vector<int64_t>& lines,
						   std::string cellDataName="", const std::vector<unsigned char>& cellData=std::vector<unsigned char>(),
						   double lengthUnit=0);

protected:
	std::string m_filename;
	std::string m_header;

	//timestep properties
	bool m_ActiveTS;
	unsigned int m_timestep;

	//mesh information
	int m_MeshType;
	std::vector<double> m_MeshLines[3];
	bool m_NativeDump;
	//! the points of a cylindrical mesh in Cartesian coordinates
	std::vector<float> m_Points;

	struct Field
	{
		std::string name;
		int components;
		std::vector<float> values;
	};
	std::vector<Field> m_Fields;
	//! remove a field of that name, as a new one replaces it
	void ClearField(std::string fieldname);

	virtual std::string GetTimestepFilename(int pad_length=10) const;
};

#endif // VTK_FILE_Writer_H
