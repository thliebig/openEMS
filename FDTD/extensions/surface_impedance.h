/*
*	Copyright (C) 2026 Roel Budé (roel.bude@proton.me)
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

#ifndef SURFACE_IMPEDANCE_H
#define SURFACE_IMPEDANCE_H

#include <complex>
#include <vector>

//! Surface impedance of a (rough) conductor and its passive rational fit, used by Operator_Ext_LossyMetal
namespace SurfaceImpedance
{

/*!
  Surface impedance of one face of a conductor slab of thickness \a t (m, may be infinite),
  referenced to the mean surface plane. The slab is driven symmetrically, i.e. the field
  vanishes in its centre plane: for a smooth face Zs = gamma/sigma * coth(gamma*t/2).
  With \a rq>0 the gradient roughness model is used: the conductivity rises across the
  surface as sigma(x) = sigma * Phi(x/rq), Phi the Gaussian CDF (G. Gold, K. Helmreich,
  "A Physical Surface Roughness Model and Its Applications", IEEE TMTT 2017), and the
  impedance is found by integrating Z' = -j w mu0 + sigma(x) Z^2 out of the metal.
  */
std::complex<double> FaceImpedance(double f, double sigma, double t, double rq);

//! Non-negative least squares (Lawson-Hanson): minimize |A x - b| with x>=0, A is m x n column major.
bool NNLS(const std::vector<double>& A, const std::vector<double>& b, int m, int n, std::vector<double>& x);

//! Rational fit of a face admittance: Y(w) = G + sum_k a_k/(1+j w/w_k), all G, a_k >= 0
struct Fit
{
	double G;
	std::vector<double> a;
	//! the impedance is referenced this far (m) outside the mean surface plane, see FitFace
	double d_ref;
	//! max. relative error of Re(Zs) and of Im(Zs) (relative to |Zs|) in the fitted band
	double err_re, err_im;
};

//! Corner frequencies w_k (rad/s) of the fit branches for a max. frequency \a f_max
std::vector<double> CornerFrequencies(double f_max, int num);

/*!
  Fit the admittance of FaceImpedance() over f_max/1000 ... f_max to G + sum_k a_k/(1+j w/w_k).
  Positive coefficients make the fit passive and causal; each term is a series R-L branch.
  A rough surface can have Im(Zs) < 0 in its mean plane, which no R-L network can represent.
  The fit is then done for the reference plane d_ref outside the mean plane, Zs + j w mu0 d_ref,
  with d_ref chosen from 0 ... 4 rq for the smallest error.
  */
Fit FitFace(double sigma, double t, double rq, double f_max, const std::vector<double>& w_k);

}

#endif // SURFACE_IMPEDANCE_H
