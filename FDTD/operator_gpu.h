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
  The operator is built on the host like the basic operator. Engine_GPU uploads
  the final coefficients to the device.
  */
class Operator_GPU : public Operator
{
public:
	//! Create a new operator for the GPU backend \a backend (see GPU_Backend::New)
	static Operator_GPU* New(const std::string& backend="auto");
	virtual ~Operator_GPU();

	virtual Engine* CreateEngine();

	const std::string& GetBackendName() const {return m_Backend;}

protected:
	//! use New() for creating a new Operator
	Operator_GPU(const std::string& backend);

	std::string m_Backend;
};

#endif // OPERATOR_GPU_H
