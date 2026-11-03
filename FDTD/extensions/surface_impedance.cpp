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

#include "surface_impedance.h"

#include <cmath>
#include <limits>
#include <algorithm>

namespace SurfaceImpedance
{

typedef std::complex<double> cplx;

static const double s_pi = 3.14159265358979323846;
static const double s_mu0 = 4e-7*s_pi;

// tanh(z) without overflow for large Re(z)
static cplx SafeTanh(cplx z)
{
	if (std::real(z)>20)
		return cplx(1,0);
	return std::tanh(z);
}

cplx FaceImpedance(double f, double sigma, double t, double rq)
{
	double w = 2*s_pi*f;
	cplx gamma = std::sqrt(cplx(0, w*s_mu0*sigma));
	double half = t/2;
	if (rq<=0)
	{
		if (std::isinf(t))
			return gamma/sigma;
		return gamma/sigma/SafeTanh(gamma*half);
	}

	// gradient model, integrate the admittance Y = 1/Z = -H/E from inside the metal
	// to x0 outside the surface (x>0 is metal): Y' = j w mu0 Y^2 - sigma(x)
	double delta = std::sqrt(2/(w*s_mu0*sigma));
	double x0 = -6*rq;
	double xd = 6*rq; // the conductivity is uniform beyond this depth
	cplx Y;
	if (xd<half)
		Y = sigma/gamma*(std::isinf(t) ? cplx(1,0) : SafeTanh(gamma*(half-xd)));
	else
	{
		xd = half; // start in the centre plane, H=0
		Y = 0;
	}
	double h = std::min(rq, delta)/40;
	int n = (int)std::ceil((xd-x0)/h);
	h = (xd-x0)/n;
	cplx jwmu(0, w*s_mu0);
	double sq2 = std::sqrt(2.0);
	double x = xd;
	for (int i=0; i<n; ++i)
	{
		double s1 = sigma*0.5*(1+std::erf(x/rq/sq2));
		double s2 = sigma*0.5*(1+std::erf((x-h/2)/rq/sq2));
		double s3 = sigma*0.5*(1+std::erf((x-h)/rq/sq2));
		cplx k1 = jwmu*Y*Y - s1;
		cplx Y2 = Y - h/2*k1;
		cplx k2 = jwmu*Y2*Y2 - s2;
		Y2 = Y - h/2*k2;
		cplx k3 = jwmu*Y2*Y2 - s2;
		Y2 = Y - h*k3;
		cplx k4 = jwmu*Y2*Y2 - s3;
		Y -= h/6*(k1 + 2.0*k2 + 2.0*k3 + k4);
		x -= h;
	}
	// remove the air layer between x0 and the mean plane
	return 1.0/Y - jwmu*std::fabs(x0);
}

// least squares solution of A(:,cols) z = b by Householder QR, A is m x n column major
static bool LeastSquares(const std::vector<double>& A, const std::vector<double>& b, int m, const std::vector<int>& cols, std::vector<double>& z)
{
	int p = cols.size();
	std::vector<double> R(m*p);
	for (int j=0; j<p; ++j)
		for (int i=0; i<m; ++i)
			R[i+j*m] = A[i+cols[j]*m];
	std::vector<double> y(b);
	for (int k=0; k<p; ++k)
	{
		double norm = 0;
		for (int i=k; i<m; ++i)
			norm += R[i+k*m]*R[i+k*m];
		norm = std::sqrt(norm);
		if (norm==0)
			return false;
		double alpha = (R[k+k*m]>0) ? -norm : norm;
		std::vector<double> v(m, 0.0);
		for (int i=k; i<m; ++i)
			v[i] = R[i+k*m];
		v[k] -= alpha;
		double vnorm = 0;
		for (int i=k; i<m; ++i)
			vnorm += v[i]*v[i];
		if (vnorm==0)
			continue;
		for (int j=k; j<p; ++j)
		{
			double s = 0;
			for (int i=k; i<m; ++i)
				s += v[i]*R[i+j*m];
			s *= 2/vnorm;
			for (int i=k; i<m; ++i)
				R[i+j*m] -= s*v[i];
		}
		double s = 0;
		for (int i=k; i<m; ++i)
			s += v[i]*y[i];
		s *= 2/vnorm;
		for (int i=k; i<m; ++i)
			y[i] -= s*v[i];
	}
	z.assign(p, 0.0);
	for (int k=p-1; k>=0; --k)
	{
		double s = y[k];
		for (int j=k+1; j<p; ++j)
			s -= R[k+j*m]*z[j];
		if (R[k+k*m]==0)
			return false;
		z[k] = s/R[k+k*m];
	}
	return true;
}

bool NNLS(const std::vector<double>& A, const std::vector<double>& b, int m, int n, std::vector<double>& x)
{
	x.assign(n, 0.0);
	std::vector<bool> passive(n, false);
	std::vector<double> w(n);
	auto gradient = [&]()
	{
		std::vector<double> r(b);
		for (int j=0; j<n; ++j)
			if (x[j]!=0)
				for (int i=0; i<m; ++i)
					r[i] -= A[i+j*m]*x[j];
		for (int j=0; j<n; ++j)
		{
			w[j] = 0;
			for (int i=0; i<m; ++i)
				w[j] += A[i+j*m]*r[i];
		}
	};
	gradient();
	double wmax0 = 0;
	for (int j=0; j<n; ++j)
		wmax0 = std::max(wmax0, std::fabs(w[j]));
	double tol = 1e-12*wmax0;
	for (int iter=0; iter<3*n+10; ++iter)
	{
		int t = -1;
		double wmax = tol;
		for (int j=0; j<n; ++j)
			if (!passive[j] && w[j]>wmax)
			{
				wmax = w[j];
				t = j;
			}
		if (t<0)
			return true;
		passive[t] = true;
		for (int inner=0; inner<3*n+10; ++inner)
		{
			std::vector<int> cols;
			for (int j=0; j<n; ++j)
				if (passive[j])
					cols.push_back(j);
			std::vector<double> zp;
			if (!LeastSquares(A, b, m, cols, zp))
				return false;
			std::vector<double> z(n, 0.0);
			bool all_pos = true;
			for (size_t c=0; c<cols.size(); ++c)
			{
				z[cols[c]] = zp[c];
				if (zp[c]<=0)
					all_pos = false;
			}
			if (all_pos)
			{
				x = z;
				break;
			}
			double alpha = 1;
			for (int j=0; j<n; ++j)
				if (passive[j] && z[j]<=0)
					alpha = std::min(alpha, x[j]/(x[j]-z[j]));
			for (int j=0; j<n; ++j)
			{
				x[j] += alpha*(z[j]-x[j]);
				if (passive[j] && x[j]<=0)
				{
					passive[j] = false;
					x[j] = 0;
				}
			}
		}
		gradient();
	}
	return true;
}

std::vector<double> CornerFrequencies(double f_max, int num)
{
	// half a decade beyond the fitted band f_max/1000 ... f_max on either side
	std::vector<double> w_k(num);
	double lo = std::log10(f_max/3000), hi = std::log10(3*f_max);
	for (int k=0; k<num; ++k)
		w_k[k] = 2*s_pi*std::pow(10.0, lo + (hi-lo)*k/std::max(num-1,1));
	return w_k;
}

Fit FitFace(double sigma, double t, double rq, double f_max, const std::vector<double>& w_k)
{
	const int numF = 100;
	int N = w_k.size();
	int numCol = N+1;
	std::vector<double> f(numF);
	std::vector<cplx> Zm(numF);
	for (int i=0; i<numF; ++i)
	{
		f[i] = f_max/1000*std::pow(1000.0, (double)i/(numF-1));
		Zm[i] = FaceImpedance(f[i], sigma, t, rq);
	}

	Fit best;
	best.err_re = std::numeric_limits<double>::max();
	int numShift = (rq>0) ? 17 : 1;
	for (int s=0; s<numShift; ++s)
	{
		double d = 0.25*s*rq;
		// weighted rows: (fit - Y)/|Y|, real and imaginary parts, columns normalized
		std::vector<double> A(2*numF*numCol), b(2*numF);
		std::vector<cplx> Y(numF);
		for (int i=0; i<numF; ++i)
		{
			double w = 2*s_pi*f[i];
			Y[i] = 1.0/(Zm[i] + cplx(0, w*s_mu0*d));
			double aY = std::abs(Y[i]);
			for (int c=0; c<numCol; ++c)
			{
				cplx col = (c==0) ? cplx(1,0) : 1.0/cplx(1, w/w_k[c-1]);
				col /= aY;
				A[i+c*2*numF] = std::real(col);
				A[i+numF+c*2*numF] = std::imag(col);
			}
			b[i] = std::real(Y[i])/aY;
			b[i+numF] = std::imag(Y[i])/aY;
		}
		std::vector<double> scale(numCol);
		for (int c=0; c<numCol; ++c)
		{
			double s2 = 0;
			for (int i=0; i<2*numF; ++i)
				s2 += A[i+c*2*numF]*A[i+c*2*numF];
			scale[c] = std::sqrt(s2);
			for (int i=0; i<2*numF; ++i)
				A[i+c*2*numF] /= scale[c];
		}
		std::vector<double> x;
		if (!NNLS(A, b, 2*numF, numCol, x))
			continue;
		for (int c=0; c<numCol; ++c)
			x[c] /= scale[c];

		double err_re = 0, err_im = 0;
		for (int i=0; i<numF; ++i)
		{
			double w = 2*s_pi*f[i];
			cplx Yf = x[0];
			for (int k=0; k<N; ++k)
				Yf += x[k+1]/cplx(1, w/w_k[k]);
			cplx Zf = 1.0/Yf, Z = 1.0/Y[i];
			err_re = std::max(err_re, std::fabs(std::real(Zf)/std::real(Z)-1));
			err_im = std::max(err_im, std::fabs(std::imag(Zf)-std::imag(Z))/std::abs(Z));
		}
		if (err_re<best.err_re)
		{
			best.G = x[0];
			best.a.assign(x.begin()+1, x.end());
			best.d_ref = d;
			best.err_re = err_re;
			best.err_im = err_im;
		}
	}
	return best;
}

}
