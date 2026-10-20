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

using namespace std;

#include "vtk_file_writer.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <zlib.h>

namespace
{

const char* ByteOrder()
{
	const uint16_t one = 1;
	return *reinterpret_cast<const unsigned char*>(&one) ? "LittleEndian" : "BigEndian";
}

string Escape(const string& s)
{
	string out;
	for (char c : s)
	{
		switch (c)
		{
		case '&': out += "&amp;"; break;
		case '<': out += "&lt;"; break;
		case '>': out += "&gt;"; break;
		case '"': out += "&quot;"; break;
		default: out += c;
		}
	}
	return out;
}

void Base64(ostream& out, const unsigned char* data, size_t n)
{
	static const char table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
	string buf;
	buf.reserve(4*(n/3+1));
	size_t i = 0;
	for (; i+2<n; i+=3)
	{
		uint32_t v = (uint32_t(data[i])<<16) | (uint32_t(data[i+1])<<8) | data[i+2];
		buf += table[v>>18]; buf += table[(v>>12)&63]; buf += table[(v>>6)&63]; buf += table[v&63];
	}
	if (i<n)
	{
		uint32_t v = uint32_t(data[i])<<16;
		if (i+1<n)
			v |= uint32_t(data[i+1])<<8;
		buf += table[v>>18]; buf += table[(v>>12)&63];
		buf += (i+1<n) ? table[(v>>6)&63] : '=';
		buf += '=';
	}
	out << buf;
}

// the data of an array as VTK's zlib compressor writes it: a header of UInt32 (the number of
// blocks, the block size, the size of the last block, the compressed size of each block), then
// the blocks of 32 KiB, each zlib compressed; header and blocks base64 encoded separately
bool WriteCompressed(ostream& out, const void* data, size_t bytes)
{
	const size_t block = 32768;
	size_t nblocks = (bytes+block-1)/block;
	vector<uint32_t> header(3+nblocks);
	header[0] = (uint32_t)nblocks;
	header[1] = (uint32_t)block;
	header[2] = nblocks>0 ? (uint32_t)(bytes-(nblocks-1)*block) : 0;
	vector<unsigned char> compressed;
	vector<unsigned char> buf(compressBound(block));
	const unsigned char* src = static_cast<const unsigned char*>(data);
	for (size_t b=0;b<nblocks;++b)
	{
		uLongf len = (uLongf)buf.size();
		size_t size = min(block, bytes-b*block);
		if (compress2(buf.data(), &len, src+b*block, (uLong)size, Z_DEFAULT_COMPRESSION)!=Z_OK)
			return false;
		header[3+b] = (uint32_t)len;
		compressed.insert(compressed.end(), buf.begin(), buf.begin()+len);
	}
	Base64(out, reinterpret_cast<const unsigned char*>(header.data()), header.size()*sizeof(uint32_t));
	Base64(out, compressed.data(), compressed.size());
	return true;
}

bool WriteArray(ostream& out, const string& indent, const char* type, const string& name, int components,
				const void* data, size_t bytes, const string& attributes="")
{
	out << indent << "<DataArray type=\"" << type << "\"";
	if (!name.empty())
		out << " Name=\"" << Escape(name) << "\"";
	if (components>1)
		out << " NumberOfComponents=\"" << components << "\"";
	out << attributes << " format=\"binary\">\n" << indent << "  ";
	bool ok = WriteCompressed(out, data, bytes);
	out << "\n" << indent << "</DataArray>\n";
	return ok;
}

void WriteFileStart(ostream& out, const char* type)
{
	out << "<?xml version=\"1.0\"?>\n"
		<< "<VTKFile type=\"" << type << "\" version=\"0.1\" byte_order=\"" << ByteOrder()
		<< "\" header_type=\"UInt32\" compressor=\"vtkZLibDataCompressor\">\n";
}

}

VTK_File_Writer::VTK_File_Writer(string filename, int meshType)
{
	SetFilename(filename);
	m_MeshType = meshType;
	m_NativeDump = false;
	m_ActiveTS = false;
	m_timestep = 0;

	if ((m_MeshType!=0) && (m_MeshType!=1))
		cerr << "VTK_File_Writer::VTK_File_Writer: Error, unknown mesh type: " << m_MeshType << endl;
}

VTK_File_Writer::~VTK_File_Writer()
{
}

void VTK_File_Writer::SetMeshLines(double const* const* lines, unsigned int const* count, double scaling)
{
	if (m_MeshType==0) //cartesian mesh
	{
		for (int n=0;n<3;++n)
		{
			m_MeshLines[n].clear();
			m_MeshLines[n].reserve(count[n]);
			for (unsigned int i=0; i<count[n]; i++)
				m_MeshLines[n].push_back(lines[n][i]*scaling);
		}
	}
	else if (m_MeshType==1) //cylindrical mesh
	{
		for (int n=0;n<3;++n)
		{
			m_MeshLines[n].clear();
			m_MeshLines[n].reserve(count[n]);
			double scale=1;
			if (n!=1)
				scale*=scaling;
			for (unsigned int i=0; i<count[n]; i++)
				m_MeshLines[n].push_back(lines[n][i]*scale);
		}

		m_Points.resize((size_t)3*count[0]*count[1]*count[2]);
		size_t id=0;
		for (unsigned int k=0; k<count[2]; ++k)
			for (unsigned int j=0; j<count[1]; ++j)
				for (unsigned int i=0; i<count[0]; ++i)
				{
					m_Points[id++] = (float)(lines[0][i] * cos(lines[1][j]) * scaling);
					m_Points[id++] = (float)(lines[0][i] * sin(lines[1][j]) * scaling);
					m_Points[id++] = (float)(lines[2][k] * scaling);
				}
	}
	else
	{
		cerr << "VTK_File_Writer::SetMeshLines: Error, unknown mesh type: " << m_MeshType << endl;
	}
}

void VTK_File_Writer::AddScalarField(std::string fieldname, ArrayLib::ArrayIJK<float> &field)
{
	Field f;
	f.name = fieldname;
	f.components = 1;
	f.values.resize(m_MeshLines[0].size()*m_MeshLines[1].size()*m_MeshLines[2].size());
	size_t id=0;
	for (unsigned int k=0;k<m_MeshLines[2].size();++k)
		for (unsigned int j=0;j<m_MeshLines[1].size();++j)
			for (unsigned int i=0;i<m_MeshLines[0].size();++i)
				f.values[id++] = field(i, j, k);
	ClearField(fieldname);
	m_Fields.push_back(f);
}

void VTK_File_Writer::AddVectorField(string fieldname, ArrayLib::ArrayNIJK<float>& field)
{
	Field f;
	f.name = fieldname;
	f.components = 3;
	f.values.resize(3*m_MeshLines[0].size()*m_MeshLines[1].size()*m_MeshLines[2].size());
	size_t id=0;
	for (unsigned int k=0;k<m_MeshLines[2].size();++k)
	{
		for (unsigned int j=0;j<m_MeshLines[1].size();++j)
		{
			float cos_a = cos(m_MeshLines[1].at(j)); //needed only for m_MeshType==1 (cylindrical mesh)
			float sin_a = sin(m_MeshLines[1].at(j)); //needed only for m_MeshType==1 (cylindrical mesh)
			for (unsigned int i=0;i<m_MeshLines[0].size();++i)
			{
				if ((m_MeshType==0) || (m_NativeDump))
				{
					f.values[id++] = field(0, i, j, k);
					f.values[id++] = field(1, i, j, k);
					f.values[id++] = field(2, i, j, k);
				}
				else
				{
					f.values[id++] = field(0, i, j, k) * cos_a - field(1, i, j, k) * sin_a;
					f.values[id++] = field(0, i, j, k) * sin_a + field(1, i, j, k) * cos_a;
					f.values[id++] = field(2, i, j, k);
				}
			}
		}
	}
	ClearField(fieldname);
	m_Fields.push_back(f);
}

void VTK_File_Writer::ClearField(std::string fieldname)
{
	for (size_t n=0;n<m_Fields.size();)
	{
		if (m_Fields[n].name==fieldname)
			m_Fields.erase(m_Fields.begin()+n);
		else
			++n;
	}
}

int VTK_File_Writer::GetNumberOfFields() const
{
	return (int)m_Fields.size();
}

void VTK_File_Writer::ClearAllFields()
{
	m_Fields.clear();
}

string VTK_File_Writer::GetTimestepFilename(int pad_length) const
{
	if (m_ActiveTS==false)
		return m_filename;

	stringstream ss;
	ss << m_filename << "_" << std::setw( pad_length ) << std::setfill( '0' ) << m_timestep;

	return ss.str();
}

bool VTK_File_Writer::Write()
{
	const char* type;
	const char* ext;
	if (m_MeshType==0) //cartesian mesh
	{
		type = "RectilinearGrid";
		ext = ".vtr";
	}
	else if (m_MeshType==1) //cylindrical mesh
	{
		type = "StructuredGrid";
		ext = ".vts";
	}
	else
	{
		cerr << "VTK_File_Writer::Write: Error, unknown mesh type: " << m_MeshType << endl;
		return false;
	}

	string filename = GetTimestepFilename() + ext;
	ofstream file(filename.c_str(), ios::binary);
	if (!file.is_open())
	{
		cerr << "VTK_File_Writer::Write: Error, can't open file: " << filename << endl;
		return false;
	}

	stringstream extent;
	extent << "0 " << (long)m_MeshLines[0].size()-1 << " 0 " << (long)m_MeshLines[1].size()-1 << " 0 " << (long)m_MeshLines[2].size()-1;
	WriteFileStart(file, type);
	file << "  <" << type << " WholeExtent=\"" << extent.str() << "\">\n";
	file << "  <Piece Extent=\"" << extent.str() << "\">\n";
	file << "    <PointData>\n";
	bool ok = true;
	for (size_t n=0;n<m_Fields.size();++n)
		ok &= WriteArray(file, "      ", "Float32", m_Fields[n].name, m_Fields[n].components,
						 m_Fields[n].values.data(), m_Fields[n].values.size()*sizeof(float));
	file << "    </PointData>\n";
	file << "    <CellData>\n";
	file << "    </CellData>\n";
	if (m_MeshType==0)
	{
		file << "    <Coordinates>\n";
		for (int n=0;n<3;++n)
			ok &= WriteArray(file, "      ", "Float64", "", 1, m_MeshLines[n].data(), m_MeshLines[n].size()*sizeof(double));
		file << "    </Coordinates>\n";
	}
	else
	{
		file << "    <Points>\n";
		ok &= WriteArray(file, "      ", "Float32", "Points", 3, m_Points.data(), m_Points.size()*sizeof(float));
		file << "    </Points>\n";
	}
	file << "  </Piece>\n";
	file << "  </" << type << ">\n";
	file << "</VTKFile>\n";
	file.close();

	if (!ok || file.fail())
	{
		cerr << "VTK_File_Writer::Write: Error, writing file failed: " << filename << endl;
		return false;
	}
	return true;
}

bool VTK_File_Writer::WriteLines(std::string filename, const std::vector<float>& points, const std::vector<int64_t>& lines,
								 std::string cellDataName, const std::vector<unsigned char>& cellData, double lengthUnit)
{
	ofstream file(filename.c_str(), ios::binary);
	if (!file.is_open())
	{
		cerr << "VTK_File_Writer::WriteLines: Error, can't open file: " << filename << endl;
		return false;
	}

	// the end of each line in the connectivity
	vector<int64_t> offsets(lines.size()/2);
	for (size_t n=0;n<offsets.size();++n)
		offsets[n] = 2*(int64_t)(n+1);

	WriteFileStart(file, "PolyData");
	file << "  <PolyData>\n";
	bool ok = true;
	if (lengthUnit>0)
	{
		file << "    <FieldData>\n";
		ok &= WriteArray(file, "      ", "Float64", "LengthUnit", 1, &lengthUnit, sizeof(double), " NumberOfTuples=\"1\"");
		file << "    </FieldData>\n";
	}
	file << "    <Piece NumberOfPoints=\"" << points.size()/3 << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << offsets.size()
		 << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
	file << "      <PointData>\n";
	file << "      </PointData>\n";
	file << "      <CellData>\n";
	if (!cellDataName.empty())
		ok &= WriteArray(file, "        ", "UInt8", cellDataName, 1, cellData.data(), cellData.size());
	file << "      </CellData>\n";
	file << "      <Points>\n";
	ok &= WriteArray(file, "        ", "Float32", "Points", 3, points.data(), points.size()*sizeof(float));
	file << "      </Points>\n";
	file << "      <Lines>\n";
	ok &= WriteArray(file, "        ", "Int64", "connectivity", 1, lines.data(), lines.size()*sizeof(int64_t));
	ok &= WriteArray(file, "        ", "Int64", "offsets", 1, offsets.data(), offsets.size()*sizeof(int64_t));
	file << "      </Lines>\n";
	file << "    </Piece>\n";
	file << "  </PolyData>\n";
	file << "</VTKFile>\n";
	file.close();

	if (!ok || file.fail())
	{
		cerr << "VTK_File_Writer::WriteLines: Error, writing file failed: " << filename << endl;
		return false;
	}
	return true;
}
