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

#include "engine_interface_sse_fdtd.h"
#include <stdexcept>

Engine_Interface_SSE_FDTD::Engine_Interface_SSE_FDTD(Operator* op) : Engine_Interface_FDTD(op)
{
	m_Eng_SSE = dynamic_cast<Engine_sse*>(m_Eng);
	if (m_Eng_SSE==NULL)
		throw std::runtime_error("Engine_Interface_SSE_FDTD::Engine_Interface_SSE_FDTD: Error: SSE-Engine is not set!");
}

Engine_Interface_SSE_FDTD::~Engine_Interface_SSE_FDTD()
{
	m_Eng_SSE=NULL;
}

double Engine_Interface_SSE_FDTD::CalcFastEnergy() const
{
	f4vector E_energy;
	E_energy.f[0]=0;
	E_energy.f[1]=0;
	E_energy.f[2]=0;
	E_energy.f[3]=0;
	f4vector H_energy;
	H_energy = E_energy;

	unsigned int pos[3];
	for (pos[0]=m_EnergyStart[0]; pos[0]<m_Op->GetNumberOfLines(0)-1; ++pos[0])
	{
		for (pos[1]=m_EnergyStart[1]; pos[1]<m_Op->GetNumberOfLines(1)-1; ++pos[1])
		{
			for (pos[2]=0; pos[2]<m_Eng_SSE->GetNumberOfVectors(); ++pos[2])
			{
				ArrayLib::ArrayNIJK<f4vector>& f4_volt = *m_Eng_SSE->Engine_sse::f4_volt_ptr;
				ArrayLib::ArrayNIJK<f4vector>& f4_curr = *m_Eng_SSE->Engine_sse::f4_curr_ptr;

				E_energy.v += f4_volt(0, pos[0], pos[1], pos[2]).v * f4_volt(0, pos[0], pos[1], pos[2]).v;
				E_energy.v += f4_volt(1, pos[0], pos[1], pos[2]).v * f4_volt(1, pos[0], pos[1], pos[2]).v;
				E_energy.v += f4_volt(2, pos[0], pos[1], pos[2]).v * f4_volt(2, pos[0], pos[1], pos[2]).v;

				H_energy.v += f4_curr(0, pos[0], pos[1], pos[2]).v * f4_curr(0, pos[0], pos[1], pos[2]).v;
				H_energy.v += f4_curr(1, pos[0], pos[1], pos[2]).v * f4_curr(1, pos[0], pos[1], pos[2]).v;
				H_energy.v += f4_curr(2, pos[0], pos[1], pos[2]).v * f4_curr(2, pos[0], pos[1], pos[2]).v;
			}
		}
	}

	double E_sum = E_energy.f[0]+E_energy.f[1]+E_energy.f[2]+E_energy.f[3];
	double H_sum = H_energy.f[0]+H_energy.f[1]+H_energy.f[2]+H_energy.f[3];

	// MPI: the z planes below m_EnergyStart are the lower rank's, the vectors span them
	for (pos[0]=m_EnergyStart[0]; pos[0]<m_Op->GetNumberOfLines(0)-1; ++pos[0])
		for (pos[1]=m_EnergyStart[1]; pos[1]<m_Op->GetNumberOfLines(1)-1; ++pos[1])
			for (pos[2]=0; pos[2]<m_EnergyStart[2]; ++pos[2])
				for (int n=0; n<3; ++n)
				{
					const double v = m_Eng_SSE->GetVolt(n, pos);
					const double i = m_Eng_SSE->GetCurr(n, pos);
					E_sum -= v*v;
					H_sum -= i*i;
				}

	return EPS0*E_sum + MUE0*H_sum;
}
