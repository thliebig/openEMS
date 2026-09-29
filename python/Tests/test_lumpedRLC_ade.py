# -*- coding: utf-8 -*-
"""
Regression tests for the series lumped RLC/RL/RC extension.

Source under test:
    FDTD/extensions/operator_ext_lumpedRLC.cpp   (trapezoidal state-space coeffs)
    FDTD/extensions/engine_ext_lumpedRLC.cpp     (q/J update + implicit coupling)

The element is integrated with the trapezoidal rule in a first-order
state-space form (states: current J and charge q):

    J[n] = A*Vd[n] + B,      B = aV*Vd[n-1] + aQ*q[n-1] + aJ*J[n-1]
    q[n] = q[n-1] + (dT/2)*(J[n] + J[n-1])
    Vd[n] = Vraw - (dT/2Cd)*(J[n] + J[n-1])          (node coupling)

The old implementation used a second-order ADE recursion
    J[n] = ib0*(Vd[n]-Vd[n-2]) - (b1*ib0)*J[n-1] - (b2*ib0)*J[n-2].
That recursion is mathematically exact in FP64, but its two poles are
separated by only ~1.5e-5, so the root sensitivity 1/(p1-p2) ~ 7e4 turns the
~1e-7 FP32 rounding of `b1*ib0` / `b2*ib0` into a ~3e-4 pole displacement.
One pole landed at 1.00034 (outside the unit circle), producing the
~1.0004/step growth seen in long single-precision runs.  The state-space form
keeps the same accuracy but is well conditioned and stable in FP32.

Run:  python -m unittest python.Tests.test_lumpedRLC_ade
 or:  python python/Tests/test_lumpedRLC_ade.py
"""

import unittest
import numpy as np


# ---------------------------------------------------------------------------
# Coefficient expressions from operator_ext_lumpedRLC.cpp
# ---------------------------------------------------------------------------
def state_space_coefficients(dL, dC, dR, dT):
    """Trapezoidal state-space series coefficients (A, aV, aQ, aJ)."""
    if dL > 0.0 and dC > 0.0:                       # series RLC
        m = 1.0 + dT*dR/(2.0*dL) + dT*dT/(4.0*dL*dC)
        A = dT/(2.0*dL*m)
        return A, A, -2.0*A/dC, 2.0/m - 1.0
    if dL > 0.0:                                    # series RL
        den = dL/dT + dR/2.0
        A = 0.5/den
        return A, A, 0.0, (dL/dT - dR/2.0)/den
    if dC > 0.0:
        if dR > 0.0:                                # series RC
            K = 2.0*dR*dC + dT
            return 2.0*dC/K, -dT/(dR*K), -(2.0*dR*dC - dT)/(dR*dC*K), 0.0
        A = 2.0*dC/dT                               # series C only
        return A, -A, 0.0, -1.0
    return 1.0/dR, 0.0, 0.0, 0.0                    # series R only


def state_space_admittance(f, dL, dC, dR, dT):
    """J(z)/Vd(z) of the state-space update, evaluated at z = exp(j w dT)."""
    A, aV, aQ, aJ = state_space_coefficients(dL, dC, dR, dT)
    z = np.exp(1j * 2.0 * np.pi * np.asarray(f, dtype=float) * dT)
    num = (A + aV * z ** -1) * (1.0 - z ** -1)
    den = (1.0 - aJ * z ** -1) * (1.0 - z ** -1) \
        - aQ * (dT / 2.0) * z ** -1 * (1.0 + z ** -1)
    return num / den


def physical_admittance(f, dL, dC, dR):
    w = 2.0 * np.pi * np.asarray(f, dtype=float)
    Z = dR + 1j * w * dL
    if dC > 0.0:
        Z += 1.0 / (1j * w * dC)
    return 1.0 / Z


# ---------------------------------------------------------------------------
# Closed-loop transition matrix: 1D Yee grid + one series element.
# Mirrors engine_sse and Engine_Ext_LumpedRLC exactly.
# State = [V(1..N-2), I(0..N-2), q, J].
# ---------------------------------------------------------------------------
def closed_loop_matrix(dL, dC, dR, dT, Cd, dtype=np.float64, cour=0.5, N=11, k=5):
    T = lambda x: dtype(x)
    A, aV, aQ, aJ = state_space_coefficients(dL, dC, dR, dT)
    VI = T(dT) / T(Cd)
    IV = T(cour * cour) / VI
    nv, ni = N - 2, N - 1
    iV = lambda i: i - 1
    iI = lambda i: nv + i
    h2 = nv + ni
    iQ, iJ = h2, h2 + 1
    M = np.zeros((h2 + 2, h2 + 2), dtype=dtype)

    # base voltage update
    for i in range(1, N - 1):
        r = iV(i)
        M[r, r] = T(1.0)
        M[r, iI(i - 1)] += VI
        M[r, iI(i)] -= VI

    raw = M[iV(k)].copy()
    vcd = T(dT) / (2.0 * T(Cd))
    vvd = T(1.0) / (T(1.0) + vcd * T(A))

    # explicit part B = aV*Vd[n-1] + aQ*q + aJ*J, with Vd[n-1] = state V[k]
    B = np.zeros(h2 + 2, dtype=dtype)
    B[iV(k)] += T(aV)
    B[iQ] += T(aQ)
    B[iJ] += T(aJ)

    # Vd[n] = vvd*(Vraw - vcd*(B + J[n-1]))
    row = vvd * raw - vvd * vcd * B
    row[iJ] += -vvd * vcd * T(1.0)
    M[iV(k)] = row
    # J[n] = A*Vd[n] + B
    M[iJ] = T(A) * row + B
    # q[n] = q[n-1] + (dT/2)*(J[n] + J[n-1])
    d2 = T(dT / 2.0)
    M[iQ] = d2 * M[iJ].copy()
    M[iQ, iJ] += d2
    M[iQ, iQ] += T(1.0)

    # base current update
    for j in range(0, N - 1):
        r = iI(j)
        M[r, r] += T(1.0)
        if 1 <= j <= N - 2:
            M[r] += IV * M[iV(j)]
        if 1 <= j + 1 <= N - 2:
            M[r] -= IV * M[iV(j + 1)]
    return M


class TestSeriesRlcAde(unittest.TestCase):
    # Board decoupling capacitor from the stability report.
    R = 20e-3
    L = 0.5e-9
    C = 220e-9
    dT = 78.65e-15

    def test_transfer_function_matches_physical(self):
        f = np.array([1e6, 5e6, 15.17e6, 30e6, 100e6])
        Y_phys = physical_admittance(f, self.L, self.C, self.R)
        Y_code = state_space_admittance(f, self.L, self.C, self.R, self.dT)
        rel = np.abs(Y_code - Y_phys) / np.maximum(np.abs(Y_phys), 1e-30)
        self.assertLess(np.max(rel), 1e-5)

    def test_poles_carry_physical_damping(self):
        A, aV, aQ, aJ = state_space_coefficients(self.L, self.C, self.R, self.dT)
        # free evolution matrix on [J, q]
        d2 = self.dT / 2.0
        M = np.array([[aJ, aQ],
                      [d2 * (aJ + 1.0), 1.0 + d2 * aQ]])
        poles = np.linalg.eigvals(M)
        rad = -self.R / (2.0 * self.L)
        wd = np.sqrt(1.0 / (self.L * self.C) - rad * rad)
        phys = np.exp((rad + 1j * wd) * self.dT)
        self.assertAlmostEqual(np.abs(poles).max(), abs(phys), places=9)
        self.assertLessEqual(np.abs(poles).max(), 1.0)

    def test_fp32_poles_inside(self):
        """The FP32 state-space update stays inside the unit circle."""
        T = np.float32
        A, aV, aQ, aJ = state_space_coefficients(self.L, self.C, self.R, self.dT)
        d2 = T(self.dT) / T(2.0)
        M = np.array([[T(aJ), T(aQ)],
                      [d2 * (T(aJ) + T(1.0)), T(1.0) + d2 * T(aQ)]], dtype=T)
        poles = np.linalg.eigvals(M.astype(np.float64))
        self.assertLessEqual(
            np.abs(poles).max(), 1.0,
            msg='FP32 state-space pole outside unit circle: |p|={:.9f}'
                .format(np.abs(poles).max()))

    def test_closed_loop_is_stable_fp32(self):
        for Cd in [1e-15, 3.5e-15, 1e-14, 1e-13, 1e-12, 1e-11]:
            M = closed_loop_matrix(self.L, self.C, self.R, self.dT, Cd,
                                   dtype=np.float32)
            rho = np.max(np.abs(np.linalg.eigvals(M.astype(np.float64))))
            self.assertLessEqual(rho, 1.0 + 1e-6,
                                 msg='Cd={:.1e} rho={:.3e}'.format(Cd, rho))

    def test_all_topologies_stable_and_accurate(self):
        dT = 78.65e-15
        cases = [('RLC', 0.5e-9, 220e-9, 20e-3),
                 ('RL', 10e-9, 0.0, 10.0),
                 ('RC', 0.0, 1e-12, 50.0),
                 ('C', 0.0, 1e-12, 0.0),
                 ('L', 10e-9, 0.0, 0.0),
                 ('R', 0.0, 0.0, 50.0)]
        f = np.array([1e6, 1e8, 1e9])
        for name, L, C, R in cases:
            Y_phys = physical_admittance(f, L, C, R)
            Y_code = state_space_admittance(f, L, C, R, dT)
            rel = np.abs(Y_code - Y_phys) / np.maximum(np.abs(Y_phys), 1e-30)
            self.assertLess(np.max(rel), 1e-5, msg=name)


if __name__ == '__main__':
    unittest.main(verbosity=2)
