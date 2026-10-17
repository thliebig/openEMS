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

#include <cmath>

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
		CheckPeriod(rel_pos, [this](size_t n, unsigned int t) {return m_E_records[n][t];});
}

void Engine_Ext_SteadyState::CheckPeriod(unsigned int rel_pos, const std::function<double(size_t, unsigned int)>& record)
{
	const unsigned int p = m_Op_SS->m_TS_period;
	const size_t count = m_Op_SS->m_E_probe_dir.size();
	bool no_valid = true;
	m_last_max_diff = 0;
	double curr_total_energy = m_Eng_Interface->CalcFastEnergy();
	if (last_total_energy>0)
	{
		m_last_max_diff = std::fabs(curr_total_energy-last_total_energy)/last_total_energy;
		no_valid = false;
	}
	last_total_energy = curr_total_energy;
	unsigned int old_pos = 0;
	unsigned int new_pos = p;
	if (rel_pos<=p)
	{
		new_pos = 0;
		old_pos = p;
	}
	std::vector<double> curr_pow(count, 0), diff_pow(count, 0);
	double max_pow = 0;
	for (size_t n=0;n<count;++n)
	{
		for (unsigned int nt=0;nt<p;++nt)
		{
			const double b_new = record(n, nt+new_pos);
			const double b_old = record(n, nt+old_pos);
			curr_pow[n] += b_new*b_new;
			diff_pow[n] += (b_old-b_new)*(b_old-b_new);
		}
		max_pow = std::max(max_pow, curr_pow[n]);
	}
	for (size_t n=0;n<count;++n)
	{
		if (curr_pow[n]>max_pow*1e-2)
		{
			m_last_max_diff = std::max(m_last_max_diff, diff_pow[n]/curr_pow[n]);
			no_valid = false;
		}
	}
	if ((no_valid) || (m_last_max_diff>1))
		m_last_max_diff = 1;
}

void Engine_Ext_SteadyState::Apply2Current()
{

}
