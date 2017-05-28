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

#ifndef ENGINE_INTERFACE_BASE_H
#define ENGINE_INTERFACE_BASE_H

#include "tools/global.h"

class Operator_Base;

//! This is the abstact base for all Engine Interface classes.
/*!
	 This is the abstact base for all Engine Interface classes. It will provide unified access to the field information of the corresponding engine.
	 All processing methods should only access this base class.
*/
class Engine_Interface_Base
{
public:
	enum InterpolationType { NO_INTERPOLATION, NODE_INTERPOLATE, CELL_INTERPOLATE };

	virtual ~Engine_Interface_Base() {;} //!< provide a virtual destructor to correctly delete derived objects

	//! Set the operator used for this engine interface.
	virtual void SetOperator(Operator_Base* base_op) {m_Op_Base=base_op;}
	//! Get the operator used for this engine interface.
	virtual const Operator_Base* GetOperator() const {return m_Op_Base;}

	//! Set the current interpolation type \sa GetInterpolationType
	void SetInterpolationType(InterpolationType type) {m_InterpolType=type;}
	//! Set the current interpolation type \sa GetInterpolationType
	void SetInterpolationType(int type) {m_InterpolType=(InterpolationType)type;}
	//! Get the current interpolation type as string \sa SetInterpolationType GetInterpolationType GetInterpolationNameByType
	std::string GetInterpolationTypeString() {return GetInterpolationNameByType(m_InterpolType);}
	//! Get the current interpolation type \sa SetInterpolationType
	InterpolationType GetInterpolationType() {return m_InterpolType;}

	//! Get the (interpolated) electric field at \p pos. \sa SetInterpolationType
	virtual double* GetEField(const unsigned int* pos, double* out) const =0;
	//! Get the (interpolated) magnetic field at \p pos. \sa SetInterpolationType
	virtual double* GetHField(const unsigned int* pos, double* out) const =0;
	//! Get the (interpolated) electric current density field at \p pos. \sa SetInterpolationType
	virtual double* GetJField(const unsigned int* pos, double* out) const =0;
	//! Get the total current density field by rot(H) at \p pos. \sa SetInterpolationType
	virtual double* GetRotHField(const unsigned int* pos, double* out) const =0;
	//! Get the (interpolated) electric flux density field at \p pos. \sa SetInterpolationType
	virtual double* GetDField(const unsigned int* pos, double* out) const =0;
	//! Get the (interpolated) magnetic flux density field at \p pos. \sa SetInterpolationType
	virtual double* GetBField(const unsigned int* pos, double* out) const =0;

	//! Calculate the electric field integral along a given line
	virtual double CalcVoltageIntegral(const unsigned int* start, const unsigned int* stop) const =0;

	//! Convert the interpolation type into a string.
	static std::string GetInterpolationNameByType(InterpolationType mode);

	//! Get the current simulation time
	virtual double GetTime(bool dualTime=false) const =0;

	//! Get the current number of timesteps
	virtual unsigned int GetNumberOfTimesteps() const =0;

	//! Calc (roughly) the total energy
	/*!
	  This method only calculates a very rough estimate of the total energy in the simulation domain.
	  The result may even be roughly proportional to the real system energy only.
	  Primary goal is speed, not accuracy.
	  */
	virtual double CalcFastEnergy() const =0;

protected:
	Engine_Interface_Base(Operator_Base* base_op);

	Operator_Base* m_Op_Base;

	InterpolationType m_InterpolType;
};

#endif // ENGINE_INTERFACE_BASE_H
