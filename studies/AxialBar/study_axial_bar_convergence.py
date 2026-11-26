# Title: Code for an Axial Bar by 1D Finite Elements: the errors against the mesh
# Attempt: to solve a uniform & a tapered bar by 2- & 3-noded elements on
#          five meshes, & to measure the rates of convergence against the
#          closed form (P2 of the ladder)
# Created on: 24/11/2026 by Claude (not yet double checked by Razin)
#   on the author's answers of 23/11/2026 on the choices of P2 (2 changed):
#   these are the lines of problems/AxialBar/MAIN_AxialBar as they stood on
#   24/11/2026, on the kernels KM_Line1D, SV_Partition, SF_Line1D &
#   UT_AxialBar, with its two loops, over the bars (Bar, one row [a p]
#   each) & over the meshes (Ne), its data, its guards & its rule for the
#   rates that are not given. They were that main file until this day,
#   with these numbers; the main file is now one bar on one mesh, & as it
#   stands it is the entry [a = 0.5, p = 1, Ne = 8] of this study.
#   The names are those of the main file & not the Loop1, Error1, Rate1 of
#   the other studies, so that an entry reads against the main file:
#     eN, eU, eE, eR   one row per bar, one column per mesh: the errors of
#                      the main file (at the nodes, L2 norm, energy norm; the
#                      reaction with the loads)
#     UL               u at x = L on every mesh, & the exact one (last column)
#     RateN, RateU, RateE   the observed order of eN, eU & eE between one
#                      mesh & the next (RateN is new with this file)
#   What it shows: with a constant EA the elements are exact at the nodes &
#   the 3-noded ones everywhere; on the tapered bar the L2 norm falls at the
#   order p+1, the energy norm at p & the error at the nodes at 2p.
#   With P = q*L*(1/a - 1) the strain of the tapered bar is uniform & the
#   elements are exact: such a bar has no rate to show.
#   The load is uniform: that the nodes of a uniform bar are exact "for any
#   load" (the catalogue) is shown for this load only.
#   differences from the lines of 24/11/2026: the shape functions of the
#   element by SF_FEMLG1D (the lines wrote [s 0] to call SF_FEMLG), the rate
#   at the nodes, & what eR checks, in its comment
# Last updated: 25/11/2026 by Claude (not yet double checked by Razin)
#   on the author's word of 25/11/2026 (his answers on the twelve choices of
#   P4: choice 12, the rates of P2 on the kernel):
#   (1) the rates are PL_Rate (core/), with the rule of before (no rate where
#       the error of the coarser mesh is below 1e-10); the same numbers in
#       every bit
# Double checked on: -
# Ported to Python on: 24/11/2026 (line by line from STUDY_AxialBarConvergence.m)
#   the nodes are counted from 0 (Elem); the gauss points & weights
#   are 1D arrays, so a product down the rows of a matrix is written with
#   [:,None]

# %% Project initialisation
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))
import warnings
import numpy as np
import matplotlib.pyplot as plt
from core.quadrature.gq_points import gq_points
from core.basis.sf_femlg1d import sf_femlg1d
from core.assembly.km_line1d import km_line1d
from core.basis.sf_line1d import sf_line1d
from core.solver.solve_prescribed import solve_prescribed
from problems.AxialBar.ut_axial_bar import ut_axial_bar
from problems.AxialBar.pl_axial_bar import pl_axial_bar
from core.post.pl_rate import pl_rate


# %% DATA INPUT

# the bar: x along the bar from the held end (N & m); q & P are positive along
# +x (away from the support) & N = E*A*du/dx is positive in tension
L = 2       # length (m)
E = 200e9   # Young's modulus (N/m^2)
A0 = 1e-3   # area at x = 0 (m^2)
q = 50e3    # distributed axial load, uniform (N/m)
P = 50e3    # load at the free end x = L (N)
ubar = 0    # displacement held at x = 0 (m): 0 as in (P2.4); with another
            # value every u moves by it through the term -K(f,c)*U(c), which
            # test_axial_bar holds on the main file, with ubar = 2e-3

# the bars solved, one row each [taper a , order p of the elements]:
# A(x) = A0*(1 - a*x/L) with 0 <= a < 1, a = 0 the uniform bar of the
# catalogue; p = 1 the linear element (2 nodes), p = 2 the quadratic (3 nodes).
# The closed form of the tapered bar loses about eps/a of its figures: a
# between 0 & about 1e-4 is not for the rates of the quadratic elements
Bar = np.array([[0,   1],
                [0,   2],
                [0.5, 1],
                [0.5, 2]])

# number of elements, one solution per entry (each twice the one before; two
# entries at the least: a rate is taken between one mesh & the next)
Ne = np.array([2, 4, 8, 16, 32])

ng = 2  # gauss points of an element for [k] & {r}: exact for these bars (degree 3); fewer than p make [K] singular
nge = 8 # gauss points of an element for the error norms (a higher rule: the exact u is not a polynomial)

# the solution drawn at the end [row of Bar , entry of Ne] (counted from 1;
# both must exist: change it with Bar & Ne)
Draw = [3, 1]

# the ranges of the input (the closed form & the elements hold only there)
if np.any((Bar[:,0] < 0) | (Bar[:,0] >= 1)):
    raise ValueError('study_axial_bar_convergence: the taper a must be 0 <= a < 1 (row %d of Bar)' % (np.flatnonzero((Bar[:,0] < 0) | (Bar[:,0] >= 1))[0]+1))
if ng < np.max(Bar[:,1]):
    raise ValueError('study_axial_bar_convergence: ng = %d is fewer than the order p = %d: [K] would be singular' % (ng,np.max(Bar[:,1])))
if np.any((Bar[:,0] > 0) & (Bar[:,0] < 1e-4)):
    warnings.warn('study_axial_bar_convergence: a taper below 1e-4 leaves the closed form about eps/a of its figures: not for the rates')


# %% SOLUTION OF EVERY BAR ON EVERY MESH
eN = np.zeros((Bar.shape[0],Ne.size)) # error at the nodes (largest, of the largest u)
eU = np.zeros((Bar.shape[0],Ne.size)) # error of the displacement over the bar (L2 norm)
eE = np.zeros((Bar.shape[0],Ne.size)) # error in the energy norm
eR = np.zeros((Bar.shape[0],Ne.size)) # reaction & loads, of the loads
UL = np.zeros((Bar.shape[0],Ne.size+1)) # u at x = L on every mesh, & exact
for b in range(Bar.shape[0]):
    # loop over the bars
    a = Bar[b,0]
    p = int(Bar[b,1])

    # AREA & AXIAL STIFFNESS: A(x) = A0*(1 - a*x/L), k = E*A for KM_Line1D;
    # the exact u, N & du/dx are those of UT_AxialBar
    Aex = lambda x: A0*(1 - a*x/L)
    EAx = lambda x: E*Aex(x)
    qx = lambda x: q + 0*x

    # SHAPE FUNCTIONS OF THE ELEMENT in s from -1 to 1, the xi of the catalogue
    # (the same for every element): its nodes from end to end, its gauss
    # points & weights for [k] & {r}, & those for the error norms
    sn = np.linspace(-1,1,p+1)
    G,W = gq_points(ng)
    Ns,dNs = sf_femlg1d(G,sn)[0:2]
    Ge,We = gq_points(nge)
    Nse,dNse = sf_femlg1d(Ge,sn)[0:2]

    for m in range(Ne.size):
        # loop over the meshes

        # NODES & ELEMENTS: equal elements, the nodes of each from end to end
        nN = p*Ne[m]+1
        xi = np.linspace(0,L,nN)
        Elem = p*np.arange(Ne[m])[:,None] + np.arange(p+1)

        # STIFFNESS MATRIX & FORCE VECTOR: the integrals of E*A*dNi/dx*dNj/dx &
        # Ni*q over the elements (KM_Line1D), & the load at x = L
        K,F = km_line1d(xi,Elem,Ns,dNs,W,EAx,qx)
        F[nN-1] = F[nN-1] + P

        # SUPPORT, SOLUTION & REACTION: the node at x = 0 held at ubar, by
        # partition, K(f,f)*U(f) = F(f) - K(f,c)*U(c) & R = K(c,:)*U - F(c)
        U,R = solve_prescribed(K,F,[0],[ubar])

        # COMPARE WITH THE EXACT SOLUTION
        # u & du/dx of the elements at the gauss points of the error norms
        # [xQ], with their weights [WxQ] (one column per element)
        xQ,WxQ,Uq,dUq = sf_line1d(xi,Elem,U,Nse,dNse,We)
        Uex,_,dUex = ut_axial_bar(xQ,L,E,A0,a,q,P,ubar)
        Uxi = ut_axial_bar(xi,L,E,A0,a,q,P,ubar)[0]
        # largest error at the nodes, of the largest displacement
        eN[b,m] = np.max(np.abs(U-Uxi))/np.max(np.abs(Uxi))
        # error of the displacement over the bar (L2 norm)
        eU[b,m] = np.sqrt(np.sum(WxQ*(Uq-Uex)**2)/np.sum(WxQ*Uex**2))
        # error in the energy norm, the integral of E*A*(du/dx - duex/dx)^2
        eE[b,m] = np.sqrt(np.sum(WxQ*E*Aex(xQ)*(dUq-dUex)**2)/np.sum(WxQ*E*Aex(xQ)*dUex**2))
        # equilibrium: the reaction & the applied loads add up to zero. This
        # is a check of the assembly & of the partition, not of the accuracy:
        # it holds for any U that solves the free equations, on any mesh
        eR[b,m] = np.abs(R[0] + P + q*L)/(np.abs(P) + np.abs(q*L))
        UL[b,m] = U[nN-1]

        # kept for the figure: u at the nodes, & the axial force E*A*du/dx
        # along every element from end to end (one column per element)
        if b == Draw[0]-1 and m == Draw[1]-1:
            sD = np.linspace(-1,1,11)
            NsD,dNsD = sf_femlg1d(sD,sn)[0:2]
            xD,_,_,dUD = sf_line1d(xi,Elem,U,NsD,dNsD,None)
            ND = E*Aex(xD)*dUD
            xDn = xi; UDn = U
            xDe = np.linspace(0,L,201)
            UDe,NDe,_ = ut_axial_bar(xDe,L,E,A0,a,q,P,ubar)
    UL[b,-1] = ut_axial_bar(L,L,E,A0,a,q,P,ubar)[0]


# %% RATES OF CONVERGENCE

# observed order between one mesh & the next, log(e1/e2)/log(l1/l2) with the
# element length l = L/Ne, & no rate where the error of the coarser mesh of a
# pair is already at the level of rounding (below 1e-10): the elements are
# exact there (the rounding grows with the mesh) (PL_Rate, the meshes along
# the 2nd dimension)
RateN = pl_rate(eN,L/Ne,1e-10)
RateU = pl_rate(eU,L/Ne,1e-10)
RateE = pl_rate(eE,L/Ne,1e-10)


# %% DISPLAY RESULTS
print('Reaction & loads: largest |R + P + q*L| = %.1e of the loads' % np.max(eR))
for b in range(Bar.shape[0]):
    print('Bar %d: A(L) = %g*A0, elements of order %d, u(L) exact = %.8e m' % (b+1,1-Bar[b,0],Bar[b,1],UL[b,-1]))
    print('   Ne   u(L) (m)         at the nodes (rate)   L2 norm (rate)     energy norm (rate)')
    for m in range(Ne.size):
        if m == 0:
            print('  %3d   %.8e   %.3e             %.3e          %.3e' % (Ne[m],UL[b,m],eN[b,m],eU[b,m],eE[b,m]))
        else:
            print('  %3d   %.8e   %.3e (%.2f)      %.3e (%.2f)   %.3e (%.2f)' % (Ne[m],UL[b,m],eN[b,m],RateN[b,m-1],eU[b,m],RateU[b,m-1],eE[b,m],RateE[b,m-1]))
print('=> Rates between the last two meshes, one bar per row [at the nodes , L2 norm , energy norm]')
for b in range(Bar.shape[0]):
    print('   %.3f   %.3f   %.3f' % (RateN[b,-1],RateU[b,-1],RateE[b,-1]))

# the errors against the element length, every bar that has a rate
l = L/Ne
plt.figure(); Leg = []
for b in range(Bar.shape[0]):
    if not np.isnan(RateN[b,-1]):
        plt.loglog(l,eN[b,:],'^:')
        Leg.append('a = %g, p = %d: at the nodes (%.2f)' % (Bar[b,0],Bar[b,1],RateN[b,-1]))
    if not np.isnan(RateU[b,-1]):
        plt.loglog(l,eU[b,:],'o-')
        Leg.append('a = %g, p = %d: L2 norm (%.2f)' % (Bar[b,0],Bar[b,1],RateU[b,-1]))
    if not np.isnan(RateE[b,-1]):
        plt.loglog(l,eE[b,:],'s--')
        Leg.append('a = %g, p = %d: energy norm (%.2f)' % (Bar[b,0],Bar[b,1],RateE[b,-1]))
plt.xlabel('element length l (m)'); plt.ylabel('relative error'); plt.grid(True)
plt.legend(Leg,loc='lower right')
plt.title('Axial bar: error against the element length (rate between the last two meshes)')

# the displacement & the axial force of the solution Draw: the force is not
# continuous from one element to the next
pl_axial_bar(xDn, UDn, xD, ND, xDe, UDe, NDe,
             "Axial bar: a = %g, %d elements of order %d" % (Bar[Draw[0]-1,0],Ne[Draw[1]-1],Bar[Draw[0]-1,1]))

plt.show()
