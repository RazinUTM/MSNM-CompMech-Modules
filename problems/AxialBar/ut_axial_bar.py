# Title: Code for the axial bar, uniform or tapered
# Attempt: to give the exact displacement, axial force & strain, for the
#          comparison (P2 of the ladder)
# Created on: 24/11/2026 by Claude (not yet double checked by Razin)
#   on the author's word of 24/11/2026 (the kernels that core/ has not are
#   written, modular, & kept in core/): the lines of the exact solution of
#   MAIN_AxialBar (README, open item 22 (c)), which now calls this file;
#   written from the equations below, (P2.4) of the catalogue for a = 0
# Last updated: 24/11/2026 by Claude
# Double checked on: -
# Ported to Python on: 24/11/2026 (line by line from UT_AxialBar.m)

import numpy as np


def ut_axial_bar(x,L,E,A0,a,q,P,ubar):
    # UDF for the displacement (u), the axial force (N) & the strain (du = du/dx)
    # at x of a bar of length L & modulus E, of area A(x) = A0*(1 - a*x/L) with
    # 0 <= a < 1, under the uniform load q along it & the load P at x = L, with
    # u = ubar held at x = 0 (q, P & N positive along +x & in tension)
    #   N  = P + q*(L - x),  du = N/(E*A),  u = ubar + integral from 0 to x of du
    #   a = 0:  u = ubar + (q*(L*x - x^2/2) + P*x)/(E*A0)
    #   a > 0:  u = ubar + L/(a*E*A0)*(q*x - (P + q*L - q*L/a)*log(1 - a*x/L))
    # (log1p for the log: the second form loses about eps/a of its figures)
    # u, N & du have the shape of x

    if a < 0 or a >= 1:
        raise ValueError('ut_axial_bar: the taper a must be 0 <= a < 1')
    x = np.asarray(x,dtype=float)

    # axial force & strain
    N = P + q*(L - x)
    du = N/(E*A0*(1 - a*x/L))

    # displacement
    if a == 0:
        u = ubar + (q*(L*x - x**2/2) + P*x)/(E*A0)
    else:
        u = ubar + L/(a*E*A0)*(q*x - (P + q*L - q*L/a)*np.log1p(-a*x/L))

    return u,N,du
