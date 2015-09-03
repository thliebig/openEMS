/*
*	Copyright (C) 2015 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "engine_ext_steadystate.h"
#include "operator_ext_steadystate.h"
#include "FDTD/engine_sse.h"
#include "FDTD/engine_interface_fdtd.h"

Engine_Ext_SteadyState::Engine_Ext_SteadyState(Operator_Ext_SteadyState* op_ext): Engine_Extension(op_ext)
{
	m_Op_SS = op_ext;
	m_Priority = ENG_EXT_PRIO_STEADYSTATE;

	for (size_t n=0;n<m_Op_SS->m_E_probe_dir.size();++n)
	{
		double* rec = new double[m_Op_SS->m_TS_period*2];
		m_E_records.push_back(rec);
	}
	m_last_max_diff = 1;
	last_total_energy = 0;
	m_Eng_Interface = NULL;
}

Engine_Ext_SteadyState::~Engine_Ext_SteadyState()
{
	for (size_t n=0;n<m_E_records.size();++n)
	{
		delete[] m_E_records.at(n);
		m_E_records.at(n) = NULL;
	}
	m_E_records.clear();
	delete m_Eng_Interface;
	m_Eng_Interface = NULL;
}

void Engine_Ext_SteadyState::Apply2Voltages()
{
	unsigned int p = m_Op_SS->m_TS_period;
	unsigned int TS = m_Eng->GetNumberOfTimesteps();
	unsigned int rel_pos = m_Eng->GetNumberOfTimesteps()%(2*p);
	for (size_t n=0;n<m_E_records.size();++n)
		m_E_records.at(n)[rel_pos] = m_Eng->GetVolt(m_Op_SS->m_E_probe_dir.at(n), m_Op_SS->m_E_probe_pos[0].at(n), m_Op_SS->m_E_probe_pos[1].at(n), m_Op_SS->m_E_probe_pos[2].at(n));
	if ((TS%(m_Op_SS->m_TS_period)==0) && (TS>=2*p))
	{
		bool no_valid = true;
		m_last_max_diff = 0;
		double curr_total_energy = m_Eng_Interface->CalcFastEnergy();
		if (last_total_energy>0)
		{
			m_last_max_diff = abs(curr_total_energy-last_total_energy)/last_total_energy;
			no_valid = false;
		}
		//cerr << curr_total_energy << "/" << last_total_energy << "=" << abs(curr_total_energy-last_total_energy)/last_total_energy << endl;
		last_total_energy = curr_total_energy;
		unsigned int old_pos = 0;
		unsigned int new_pos = p;
		if (rel_pos<=p)
		{
			new_pos = 0;
			old_pos = p;
		}
		//cerr << TS << "/" << rel_pos << ": one period complete, new_pos" << new_pos << " old pos: " << old_pos << endl;
		double *curr_pow = new double[m_E_records.size()];
		double *diff_pow = new double[m_E_records.size()];
		double max_pow = 0;
		for (size_t n=0;n<m_E_records.size();++n)
		{
			double *buf = m_E_records.at(n);
			curr_pow[n] = 0;
			diff_pow[n] = 0;
			for (unsigned int nt=0;nt<p;++nt)
			{
				curr_pow[n] += buf[nt+new_pos]*buf[nt+new_pos];
				diff_pow[n] += (buf[nt+old_pos]-buf[nt+new_pos])*(buf[nt+old_pos]-buf[nt+new_pos]);
			}
			max_pow = max(max_pow, curr_pow[n]);
		}
		for (size_t n=0;n<m_E_records.size();++n)
		{
			//cerr << "curr_pow: " << curr_pow[n] <<  " diff_pow: " << diff_pow[n] << " diff: " << diff_pow[n]/curr_pow[n] << endl;
			if (curr_pow[n]>max_pow*1e-2)
			{
				m_last_max_diff = max(m_last_max_diff, diff_pow[n]/curr_pow[n]);
				//cerr << m_last_max_diff << endl;
				no_valid = false;
			}
		}
		if ((no_valid) || (m_last_max_diff>1))
			m_last_max_diff = 1;
		delete[] curr_pow; curr_pow = NULL;
		//cerr << m_last_max_diff << endl;
	}
}

void Engine_Ext_SteadyState::Apply2Current()
{

}
