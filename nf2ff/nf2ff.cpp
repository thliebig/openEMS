/*
*	Copyright (C) 2012 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "nf2ff.h"
#include "nf2ff_calc.h"
#include "../tools/array_ops.h"
#include "../tools/useful.h"
#include "../tools/hdf5_file_reader.h"
#include "../tools/hdf5_file_writer.h"
#include <hdf5.h>
#include <boost/algorithm/string.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <complex>
#include <iostream>
#include <sstream>

//external libs
#include "tinyxml.h"

nf2ff::nf2ff(vector<float> freq, vector<float> theta, vector<float> phi, vector<float> center, unsigned int numThreads)
{
	m_freq = freq;

	m_numTheta = theta.size();
	m_theta = new float[m_numTheta];
	for (size_t n=0;n<m_numTheta;++n)
		 m_theta[n]=theta.at(n);

	m_numPhi = phi.size();
	m_phi = new float[m_numPhi];
	for (size_t n=0;n<m_numPhi;++n)
		 m_phi[n]=phi.at(n);

	m_nf2ff.resize(freq.size(),NULL);
	for (size_t fn=0;fn<freq.size();++fn)
	{
		m_nf2ff.at(fn) = new nf2ff_calc(freq.at(fn),theta, phi, center);
		if (numThreads)
			m_nf2ff.at(fn)->SetNumThreads(numThreads);
	}
	m_radius = 1;
	m_Verbose = 0;
}

nf2ff::~nf2ff()
{
	m_freq.clear();
	for (size_t fn=0;fn<m_nf2ff.size();++fn)
		delete m_nf2ff.at(fn);
	m_nf2ff.clear();

	delete[] m_phi;
	m_phi = NULL;
	delete[] m_theta;
	m_theta = NULL;
}

bool nf2ff::AnalyseXMLNode(TiXmlElement* ti_nf2ff)
{
	if (ti_nf2ff==NULL)
		return false;

	unsigned int numThreads=0;
	int ihelp=0;
	if (ti_nf2ff->QueryIntAttribute("NumThreads",&ihelp) == TIXML_SUCCESS)
	{
		numThreads = ihelp;
		cerr << "nf2ff: Set number of threads to: " << numThreads << endl;
	}
	int Verbose=0;
	if (ti_nf2ff->QueryIntAttribute("Verbose",&Verbose) == TIXML_SUCCESS)
		cerr << "nf2ff: Set verbose level to " << Verbose << endl;
	else
		Verbose = 0;

	const char* attr = NULL;
	attr = ti_nf2ff->Attribute("freq");
	if (attr==NULL)
	{
		cerr << "nf2ff::AnalyseXMLNode: Can't read frequency inforamtions ... " << endl;
		return false;
	}
	vector<float> freq = SplitString2Float(attr);

	vector<float> center;
	attr = ti_nf2ff->Attribute("Center");
	if (attr!=NULL)
		center = SplitString2Float(attr);

	attr = ti_nf2ff->Attribute("Outfile");
	if (attr==NULL)
	{
		cerr << "nf2ff::AnalyseXMLNode: Can't read frequency inforamtions ... " << endl;
		return false;
	}
	string outfile = string(attr);
	if (outfile.empty())
	{
		cerr << "nf2ff::AnalyseXMLNode: outfile is empty, skipping nf2ff... " << endl;
		return false;
	}

	TiXmlElement* ti_theta = ti_nf2ff->FirstChildElement("theta");
	if (ti_theta==NULL)
	{
		cerr << "nf2ff::AnalyseXMLNode: Can't read theta values ... " << endl;
		return false;
	}
	TiXmlNode* ti_theta_node = ti_theta->FirstChild();
	if (ti_theta_node==NULL)
	{
		cerr << "nf2ff::AnalyseXMLNode: Can't read theta text child ... " << endl;
		return false;
	}
	TiXmlText* ti_theta_text = ti_theta_node->ToText();
	if (ti_theta_text==NULL)
	{
		cerr << "nf2ff::AnalyseXMLNode: Can't read theta text values ... " << endl;
		return false;
	}
	vector<float> theta = SplitString2Float(ti_theta_text->Value());

	TiXmlElement* ti_phi = ti_nf2ff->FirstChildElement("phi");
	if (ti_phi==NULL)
	{
		cerr << "nf2ff::AnalyseXMLNode: Can't read phi values ... " << endl;
		return false;
	}
	TiXmlNode* ti_phi_node = ti_phi->FirstChild();
	if (ti_phi_node==NULL)
	{
		cerr << "nf2ff::AnalyseXMLNode: Can't read phi text child ... " << endl;
		return false;
	}
	TiXmlText* ti_phi_text = ti_phi_node->ToText();
	if (ti_phi_text==NULL)
	{
		cerr << "nf2ff::AnalyseXMLNode: Can't read phi text values ... " << endl;
		return false;
	}
	vector<float> phi = SplitString2Float(ti_phi_text->Value());

	nf2ff* l_nf2ff = new nf2ff(freq,theta,phi,center,numThreads);
	l_nf2ff->SetVerboseLevel(Verbose);

	TiXmlElement* ti_Planes = ti_nf2ff->FirstChildElement();
	string E_name;
	string H_name;
	while (ti_Planes!=NULL)
	{
		E_name = string(ti_Planes->Attribute("E_Field"));
		H_name = string(ti_Planes->Attribute("H_Field"));
		if ((!E_name.empty()) && (!H_name.empty()))
		{
			if (l_nf2ff->AnalyseFile(E_name,H_name)==false)
			{
				cerr << "nf2ff::AnalyseXMLNode: Error, analysing Plane ... " << endl;
				return false;
			}
		}
		else
		{
			cerr << "nf2ff::AnalyseXMLNode: Error, invalid plane entry ... " << endl;
			return false;
		}
		ti_Planes = ti_Planes->NextSiblingElement("Planes");
	}
	l_nf2ff->Write2HDF5(outfile);
	delete l_nf2ff;
	return true;
}

bool nf2ff::AnalyseXMLFile(string filename)
{
	TiXmlDocument doc(filename.c_str());
	if (!doc.LoadFile())
	{
		cerr << "nf2ff::AnalyseXMLFile: Error loading xml-file failed!!! File: " << filename << endl;
		return false;
	}
	TiXmlElement* ti_nf2ff = doc.FirstChildElement("nf2ff");
	if (ti_nf2ff==NULL)
	{
		cerr << "nf2ff::AnalyseXMLFile: Can't read nf2ff ... " << endl;
		return false;
	}

	return AnalyseXMLNode(ti_nf2ff);
}

bool nf2ff::AnalyseFile(string E_Field_file, string H_Field_file)
{
	HDF5_File_Reader E_file(E_Field_file);
	HDF5_File_Reader H_file(H_Field_file);

	if (m_Verbose>0)
		cerr << "nf2ff: Reading planes: " << E_Field_file << " & " << E_Field_file << endl;

	// read E-mesh
	float* E_lines[3]={NULL,NULL,NULL};
	unsigned int E_numLines[3];
	int E_meshType;
	if (E_file.ReadMesh(E_lines, E_numLines, E_meshType) == false)
	{
		cerr << "nf2ff::AnalyseFile: Error reading  E-field mesh..." << endl;
		return false;
	}

	// read H-mesh
	float* H_lines[3]={NULL,NULL,NULL};
	unsigned int H_numLines[3];
	int H_meshType;
	if (H_file.ReadMesh(H_lines, H_numLines, H_meshType) == false)
	{
		cerr << "nf2ff::AnalyseFile: Error reading H-Field mesh..." << endl;
		return false;
	}

	// compare E/H meshs
	if (E_meshType!=H_meshType)
	{
		cerr << "nf2ff::AnalyseFile: Error mesh types don't agree" << endl;
		return false;
	}
	if ((E_numLines[0]!=H_numLines[0]) || (E_numLines[1]!=H_numLines[1]) || (E_numLines[2]!=H_numLines[2]))
	{
		cerr << "nf2ff::AnalyseFile: Error mesh dimensions don't agree" << endl;
		return false;
	}
	for (int n=0;n<3;++n)
		for (unsigned int m=0;m<E_numLines[n];++m)
			if (E_lines[n][m]!=H_lines[n][m])
			{
				cerr << "nf2ff::AnalyseFile: Error mesh lines don't agree" << endl;
				return false;
			}

	if (m_Verbose>1)
		cerr << "nf2ff: Data-Size: " << E_numLines[0] << "x" << E_numLines[1] << "x"  << E_numLines[2] << endl;
	if (m_Verbose>1)
		cerr << "nf2ff: calculate dft..." << endl;

	unsigned int data_size[4];
	vector<complex<float>****> E_fd_data;
	E_file.CalcFDVectorData(m_freq,E_fd_data,data_size);

	vector<complex<float>****> H_fd_data;
	H_file.CalcFDVectorData(m_freq,H_fd_data,data_size);

	if (m_Verbose>0)
		cerr << "nf2ff: Analysing far-field for " <<  m_nf2ff.size() << " frequencies.  " << endl;

	for (size_t fn=0;fn<m_nf2ff.size();++fn)
	{
		if (m_Verbose>1)
			cerr << "nf2ff: f = " << m_freq.at(fn) << "Hz (" << fn+1 << "/" << m_nf2ff.size() << ") ...";
		m_nf2ff.at(fn)->AddPlane(E_lines, E_numLines, E_fd_data.at(fn), H_fd_data.at(fn),E_meshType);
		if (m_Verbose>1)
			cerr << " done." << endl;
	}
	for (int n=0;n<3;++n)
	{
		delete[] H_lines[n];
		delete[] E_lines[n];
	}
	return true;
}

bool nf2ff::Write2HDF5(string filename)
{
	HDF5_File_Writer hdf_file(filename);

	//write mesh information
	hdf_file.SetCurrentGroup("/Mesh");
	size_t meshsize[1]={m_numTheta};
	if (hdf_file.WriteData(string("theta"),m_theta,1,meshsize)==false)
		return false;
	meshsize[0]=m_numPhi;
	if (hdf_file.WriteData(string("phi"),m_phi,1,meshsize)==false)
		return false;
	meshsize[0]=1;
	float rad[1]={m_radius};
	if (hdf_file.WriteData(string("r"),rad,1,meshsize)==false)
		return false;

	float attr_value = 2;
	hdf_file.WriteAtrribute("/Mesh", "MeshType", &attr_value, 1);

	//write field data
	size_t dim = 2;
	size_t pos = 0;
	size_t datasize[2]={m_numPhi,m_numTheta};
	size_t size = datasize[0]*datasize[1];
	float* buffer = new float[size];
	complex<float>** field_data;
	string field_names[2]={"E_theta", "E_phi"};
	for (int n=0;n<2;++n)
	{
		hdf_file.SetCurrentGroup("/nf2ff/" + field_names[n] + "/FD");
		for (size_t fn=0;fn<m_freq.size();++fn)
		{
			stringstream ss;
			ss << "f" << fn;
			pos = 0;
			if (n==0)
				field_data = GetETheta(fn);
			else
				field_data = GetEPhi(fn);
			for (size_t j=0;j<m_numPhi;++j)
				for (size_t i=0;i<m_numTheta;++i)
				{
					buffer[pos++]=real(field_data[i][j]);
				}
			if (hdf_file.WriteData(ss.str() + "_real",buffer,dim,datasize)==false)
			{
				delete[] buffer;
				cerr << "nf2ff::Write2HDF5: Error writing field data" << endl;
				return false;
			}

			pos = 0;
			for (size_t j=0;j<m_numPhi;++j)
				for (size_t i=0;i<m_numTheta;++i)
				{
					buffer[pos++]=imag(field_data[i][j]);
				}
			if (hdf_file.WriteData(ss.str() + "_imag",buffer,dim,datasize)==false)
			{
				delete[] buffer;
				cerr << "nf2ff::Write2HDF5: Error writing field data" << endl;
				return false;
			}
		}
	}
	delete[] buffer;

	//write frequency attribute
	hdf_file.WriteAtrribute("/nf2ff", "Frequency",m_freq);

	buffer = new float[m_freq.size()];
	//write radiated power attribute
	for (size_t fn=0;fn<m_freq.size();++fn)
		buffer[fn] = GetRadPower(fn);
	hdf_file.WriteAtrribute("/nf2ff", "Prad",buffer,m_freq.size());

	//write max directivity attribute
	for (size_t fn=0;fn<m_freq.size();++fn)
		buffer[fn] = GetMaxDirectivity(fn);
	hdf_file.WriteAtrribute("/nf2ff", "Dmax",buffer,m_freq.size());

	delete[] buffer;
	return true;
}
