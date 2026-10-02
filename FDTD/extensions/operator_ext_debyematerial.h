/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef OPERATOR_EXT_DEBYEMATERIAL_H
#define OPERATOR_EXT_DEBYEMATERIAL_H

#include "operator_ext_dispersive.h"
#include "FDTD/operator.h"

//! Operator extension for an n-term Debye dispersive material
/*!
  eps_r(w) = eps_inf + sum_k d_eps_k / (1 + j w tau_k)

  Each pole is an R-C branch (R_k = tau_k/C_k, C_k = eps0 d_eps_k A/l) across the
  cell capacitance. Unlike the Drude/Lorentz branch it carries no inductance, so
  its loop equation R_k I_k + V_Ck = V holds no time derivative and the branch
  current is not a state: the only state is the capacitor voltage V_Ck, which
  obeys

      tau_k dV_Ck/dt + V_Ck = V

  Integrating that with the trapezoidal rule gives

      V_Ck^(n+1) = c1_k V_Ck^n + c2_k (V^(n+1) + V^n)
      c1_k = (2 tau_k - dT)/(2 tau_k + dT),  c2_k = dT/(2 tau_k + dT)

  whose root c1_k is the same damped root the Drude/Lorentz branch gets, for any
  dT, with no limit on d_eps_k. The branch current the voltage update needs is
  just the increment, so the voltage decrement contributed by pole k is
  r_k (V_Ck^(n+1) - V_Ck^n), with r_k = C_k * vi / dT the branch to cell
  capacitance ratio and vi the cell's own voltage-update coefficient, which
  carries the cell losses and the cylindrical corrections.

  The state kept is that decrement itself, W_k = r_k V_Ck, rather than the
  capacitor voltage: it absorbs r_k, so only two coefficients per pole survive,

      W_k^(n+1) = c1_k W_k^n + d_k (V^(n+1) + V^n),   d_k = r_k c2_k
                                                         = C_k vi/(2 tau_k + dT)

  and the cell's total decrement is just sum_k (W_k^(n+1) - W_k^n). V^(n+1)
  appears on the right, so the implicit term is folded in per cell: with
  S = sum_k d_k, the single stored coefficient is 1/(1+S). See
  Engine_Ext_DebyeMaterial for the resulting update.
  */
class Operator_Ext_DebyeMaterial : public Operator_Ext_Dispersive
{
	friend class Engine_Ext_DebyeMaterial;
public:
	Operator_Ext_DebyeMaterial(Operator* op);
	virtual ~Operator_Ext_DebyeMaterial();

	virtual Operator_Extension* Clone(Operator* op);

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const {UNUSED(closedAlpha); UNUSED(R0_included); return true;}
	virtual bool IsCylindricalMultiGridSave(bool child) const {UNUSED(child); return true;}

	virtual std::string GetExtensionName() const
	{
		return std::string("Debye Dispersive Material Extension");
	}

	//! Number of Debye poles. Not the base's m_Order, which counts ADE
	//! application slots and is 1 here -- see below.
	virtual int GetDispersionOrder() {return m_PoleCount;}

	virtual void ShowStat(std::ostream &ostr) const;

protected:
	//! Copy constructor
	Operator_Ext_DebyeMaterial(Operator* op, Operator_Ext_DebyeMaterial* op_ext);

	void DeleteArrays();

	//! Number of Debye poles.
	/*!
	  The base's m_Order is deliberately 1 and does not hold this: it counts the
	  ADE corrections applied to a cell voltage, and the implicit solve couples
	  all poles into a single joint correction, carried in volt_ADE[0]. The cell
	  list is shared by every pole, so it lives in m_LM_pos[0]/m_LM_Count[0] and
	  is not replicated per pole. Operator_Ext_ConductingSheet uses m_Order the
	  same way, as a slot count rather than a material order.
	  */
	int m_PoleCount;

	//! How much of a pole's state survives one timestep: (2 tau - dT)/(2 tau + dT).
	//! Depends on tau alone, is always within (-1,1), and is the same damped root
	//! the Drude/Lorentz branch gets. Array setup: coeff[N_poles][direction][index]
	FDTD_FLOAT ***v_relax_ADE;

	//! How hard the cell voltage drives a pole: C_k vi/(2 tau_k + dT). Carries the
	//! pole strength d_eps_k, and is applied to V^n and to V^(n+1) in turn.
	//! Array setup: coeff[N_poles][direction][index]
	FDTD_FLOAT ***v_drive_ADE;

	//! Implicit solve factor 1/(1 + sum_k v_drive_ADE_k), which turns the raw
	//! updated voltage into V^(n+1). Array setup: coeff[direction][index]
	FDTD_FLOAT **v_solve_ADE;
};

#endif // OPERATOR_EXT_DEBYEMATERIAL_H
