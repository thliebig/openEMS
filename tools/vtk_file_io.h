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

#ifndef VTK_FILE_IO_H
#define VTK_FILE_IO_H

#include "base_file_io.h"

class vtkDataSet;

class VTK_File_IO : public Base_File_IO
{
public:
	VTK_File_IO(std::string filename, int meshType=0);
	virtual ~VTK_File_IO();

	virtual void SetMeshLines(double const* const* lines, unsigned int const* count, double scaling=1);

	virtual void AddScalarField(std::string fieldname, double const* const* const* field, unsigned int const* size);
	virtual void AddScalarField(std::string fieldname, float const* const* const* field, unsigned int const* size);
	virtual void AddVectorField(std::string fieldname, double const* const* const* const* field, unsigned int const* size);
	virtual void AddVectorField(std::string fieldname, float const* const* const* const* field, unsigned int const* size);

	virtual int GetNumberOfFields() const;
	virtual void ClearAllFields();

	virtual bool Write();

	virtual bool WriteASCII();
	virtual bool WriteXML();

protected:
	vtkDataSet* m_GridData;

	virtual std::string GetTimestepFilename(int pad_length=10) const;
};

#endif // VTK_FILE_IO_H
