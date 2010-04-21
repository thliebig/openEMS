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

#ifndef OPERATOR_H
#define OPERATOR_H

#include "ContinuousStructure.h"
#include "tools/AdrOp.h"
#include "tools/constants.h"

#define FDTD_FLOAT float


#if __SIZEOF_FLOAT__ != 4
	#error wrong size of float
#endif

typedef float v4sf __attribute__ ((vector_size (16))); // vector of four single floats

union f4vector
{
  v4sf v;
  float f[4];
};



//! Abstract base-class for the FDTD-operator
class Operator
{
public:
	//! Create a new operator
	static Operator* New();
	virtual ~Operator();

	virtual bool SetGeometryCSX(ContinuousStructure* geo);

	virtual int CalcECOperator();

	//! Calculate a custom signal \return number of Nyquist timesteps defined by f0
	virtual unsigned int CalcCustomExcitation(double f0, int nTS, string signal);
	//! Calculate an excitation with center of f0 and the half bandwidth fc \return number of Nyquist timesteps
	virtual unsigned int CalcGaussianPulsExcitation(double f0, double fc);
	//! Calculate a sinusoidal excitation with frequency f0 and a duration of nTS number of timesteps \return number of Nyquist timesteps
	virtual unsigned int CalcSinusExcitation(double f0, int nTS);
	//! Calculate a dirac impuls excitation \return number of Nyquist timesteps
	virtual unsigned int CalcDiracPulsExcitation();
	//! Calculate a step excitation \return number of Nyquist timesteps
	virtual unsigned int CalcStepExcitation();

	//! Get the excitation timestep with the (first) max amplitude
	virtual unsigned int GetMaxExcitationTimestep() const;

	virtual void ApplyElectricBC(bool* dirs); //applied by default to all boundaries
	virtual void ApplyMagneticBC(bool* dirs);

	double GetTimestep() const {return dT;};
	double GetNumberCells() const;

	virtual unsigned int GetNumberOfLines(int ny) const {return numLines[ny];}

	void SetNyquistNum(unsigned int nyquist) {m_nyquistTS=nyquist;}
	unsigned int GetNyquistNum() const {return m_nyquistTS;}
	unsigned int CalcNyquistNum(double fmax);

	void ShowStat() const;

	void DumpOperator2File(string filename);
	void DumpMaterial2File(string filename);

	virtual double GetGridDelta() const {return gridDelta;}
	//! Get the mesh delta times the grid delta for a 3D position
	virtual double GetMeshDelta(int n, const int* pos, bool dualMesh=false) const;
	virtual double GetMeshDelta(int n, const unsigned int* pos, bool dualMesh=false) const;

	//! Get the disc line in n direction
	virtual double GetDiscLine(int n, int pos, bool dualMesh=false) const;
	virtual double GetDiscLine(int n, unsigned int pos, bool dualMesh=false) const;
	virtual bool SnapToMesh(double* coord, unsigned int* uicoord, bool lower=false, bool* inside=NULL);

protected:
	//! use New() for creating a new Operator
	Operator();

	virtual void Init();
	virtual void Reset();
	virtual void InitOperator();

	struct Grid_Path
	{
		vector<unsigned int> posPath[3];
		vector<unsigned short> dir;
	};
	struct Grid_Path FindPath(double start[], double stop[]);

	ContinuousStructure* CSX;

	//E-Field Excitation
	//!	  Calc the electric field excitation.
	virtual bool CalcEFieldExcitation();

	virtual bool CalcPEC();

	//Calc timestep only internal use
	virtual double CalcTimestep();
	double dT; //FDTD timestep!
	unsigned int m_nyquistTS;

	//! Calc operator at certain pos
	virtual inline void Calc_ECOperatorPos(int n, unsigned int* pos);

	//EC elements, internal only!
	virtual bool Calc_EC();
	virtual bool Calc_ECPos(int n, unsigned int* pos, double* inEC);
	virtual bool Calc_EffMatPos(int n, unsigned int* pos, double* inMat);
	double* EC_C[3];
	double* EC_G[3];
	double* EC_L[3];
	double* EC_R[3];

	unsigned int numLines[3];
	double* discLines[3];
	double gridDelta;
	AdrOp* MainOp;
	AdrOp* DualOp;

	// engine/post-proc needs access
public:
	//EC operator
	FDTD_FLOAT**** vv; //calc new voltage from old voltage
	FDTD_FLOAT**** vi; //calc new voltage from old current
	FDTD_FLOAT**** ii; //calc new current from old current
	FDTD_FLOAT**** iv; //calc new current from old voltage
	f4vector**** vv_; //calc new voltage from old voltage
	f4vector**** vi_; //calc new voltage from old current
	f4vector**** iv_; //calc new current from old current
	f4vector**** ii_; //calc new current from old voltage

	//Excitation time-signal
	unsigned int ExciteLength;
	FDTD_FLOAT* ExciteSignal;

	//E-Field Excitation
	unsigned int E_Exc_Count;
	unsigned int* E_Exc_index[3];
	unsigned short* E_Exc_dir;
	FDTD_FLOAT* E_Exc_amp; //represented as edge-voltages!!
	unsigned int* E_Exc_delay;
};

#endif // OPERATOR_H
