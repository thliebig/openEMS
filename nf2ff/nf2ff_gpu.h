/*
*	Copyright (C) 2026 Hugodonotexit
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

#ifndef NF2FF_GPU_H
#define NF2FF_GPU_H

#include <complex>
#include <string>

// GPU backend of the nf2ff surface integral, written in HIP and built with
// -DNF2FF_HIP=ON. It runs on AMD and on NVIDIA GPUs.
// Everything here is plain C++, only nf2ff_gpu.hip includes HIP headers.
namespace nf2ff_gpu
{

//! True if a GPU is usable. Probed once, device 0.
bool Available();

//! Name of the device in use, empty if there is none.
std::string DeviceName();

//! Why the last call failed, or why no device is usable.
std::string LastError();

/*! Surface integrals of one rectangular Cartesian plane at one frequency.

	The aperture phase separates per axis, exp(jk u.r') = exp(jk u_p p) exp(jk u_q q),
	so one factor per (angle, grid line) is enough and the sum over the
	plane is a complex matrix product followed by a row reduction.

	ax_n is the axis normal to the plane, ax_p/ax_q the in-plane axes.
	crd_p/crd_q are their grid lines and const_off the position of the plane,
	all relative to the phase center. G is (n_p, 4*n_q) row-major: block c
	holds the area weighted currents [J_b, J_c, M_b, M_c] at (p,q), with
	b=(ax_n+1)%3 and c=(ax_n+2)%3.

	out receives 4*n_th*n_ph values [Nt | Np | Lt | Lp], angle index tn*n_ph+pn.
	Returns false on failure, see LastError().
*/
bool PlaneSums(float k,
               unsigned int n_th, const float* theta,
               unsigned int n_ph, const float* phi,
               int ax_n, int ax_p, int ax_q,
               unsigned int n_p, const float* crd_p,
               unsigned int n_q, const float* crd_q,
               float const_off,
               const std::complex<float>* G,
               std::complex<float>* out);

} // namespace nf2ff_gpu

#endif // NF2FF_GPU_H
