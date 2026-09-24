/*
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
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

#ifndef OPERATOR_GPU_H
#define OPERATOR_GPU_H

#include "operator.h"

//! Operator for the GPU engine
/*!
  The operator is built on the host like the basic operator, with the material
  and PEC evaluation split over threads like Operator_Multithread. Engine_GPU
  uploads the final coefficients to the device.

  It derives from the basic Operator and not from Operator_Multithread; only the
  two threaded build steps are repeated here.

  This is the operator of a Cartesian mesh only. A cylindrical mesh keeps its
  Operator_Cylinder / Operator_CylinderMultiGrid, which create the GPU engine
  once Operator_Cylinder::SetGPUBackend() named a backend.
  */
class Operator_GPU : public Operator
{
public:
	//! Create a new operator for the GPU backend \a backend (see GPU_Backend::New)
	static Operator_GPU* New(const std::string& backend="auto");
	virtual ~Operator_GPU();

	//! Create the Engine_GPU on the backend of this operator
	virtual Engine* CreateEngine();

	const std::string& GetBackendName() const {return m_Backend;}

protected:
	//! use New() for creating a new Operator
	Operator_GPU(const std::string& backend);

	//! Calc_EC_Range() and CalcPEC_Range() over x line ranges in parallel
	virtual bool Calc_EC();
	virtual bool CalcPEC();

	//! x line ranges of the threads, at most one per available CPU, bounds inclusive
	void ThreadRanges(std::vector<unsigned int>& start, std::vector<unsigned int>& stop) const;

	std::string m_Backend;
};

#endif // OPERATOR_GPU_H
