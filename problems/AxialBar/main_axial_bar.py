# Title: Code for an Axial Bar by 1D Finite Elements (linear or quadratic)
# Attempt: to solve a tapered bar under a distributed load & an end load by
#          2-noded (or 3-noded) elements on one mesh, & to compare with the
#          closed form (P2 of the ladder)
# Created on: 23/11/2026 by Claude (not yet double checked by Razin)
#   on the author's word of 23/11/2026 (the next problem of his ladder): P2
#   of "A Graded Catalogue of Computational Mechanics Problems for Civil
#   Engineering" (v0.1, 22/11/2026). The functions of an element are those
#   of SF_FEMLG on its nodes in s (one row of nodes, so they are the 1D
#   Lagrange functions) & its gauss points are those of GQ_Points; the
#   element matrices are assembled as in MAIN_SpringSystem & the node at
#   x = 0 is held by partition. Nothing is printed or drawn inside the loops.
#   With P = q*L*(1/a - 1) the strain of the tapered bar is uniform & the
#   elements are exact: such a bar has no rate to show.
#   The load is uniform: that the nodes of a uniform bar are exact "for any
#   load" (the catalogue) is shown for this load only
# Last updated: 23/11/2026 by Claude
# Last updated: 24/11/2026 by Claude (not yet double checked by Razin)
#   on the author's word of 24/11/2026 (the codes written his way, by his
#   skill "my way of writing codes"), after an AI review against it:
#   (1) no rate where the error of the coarser mesh of a pair is already at
#       the level of rounding, as in MAIN_SteadyConduction1D (the rule of the
#       first mesh printed rates of rounding once the meshes were fine)
#   (2) the file stops where a taper is not 0 <= a < 1 or ng is fewer than
#       p, & warns where a taper is below 1e-4
#   (3) a comment: the sign convention
# Last updated: 24/11/2026 by Claude (not yet double checked by Razin)
#   on the author's word of 24/11/2026 (use the kernels of core/; where there
#   is none, write one, modular, & keep it in core/):
#   (4) the exact solution is UT_AxialBar, the assembly KM_Line1D, the
#       partition SV_Partition & u, du/dx at the gauss points of the error
#       norms & of the figure SF_Line1D (core/), in place of their lines here;
#       no recorded number moved
# Last updated: 24/11/2026 by Claude (not yet double checked by Razin)
#   on the author's answers of 23/11/2026 on the choices of P2 (2 changed):
#   (5) one bar on one mesh (a = 0.5, p = 1, Ne = 8). The loops over the bars
#       & the meshes, the rates & the figure of the errors are now
#       studies/AxialBar/STUDY_AxialBarConvergence, which keeps the names of
#       this file: this file as it stands is its entry [a = 0.5, p = 1, Ne = 8]
#   (6) the shape functions by SF_FEMLG1D (the file wrote the points & the
#       nodes as [s 0] to call SF_FEMLG), so sn is a vector
#   (7) the comment at eR says what it checks
#   (8) the guards & the warning are written for the one bar (a, p); the
#       rule of the rates went with the rates to the study
#   no number moved: eN, eU, eE & eR are those of the file of before for this
#   bar & mesh in every figure (Python & GNU Octave)
# Double checked on: -
# Ported to Python on: 23/11/2026 (line by line from MAIN_AxialBar.m)
#   the nodes are counted from 0 (Elem); the gauss points & weights
#   are 1D arrays, so a product down the rows of a matrix is written with
#   [:,None]

# %% Project initialisation
import os
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


# %% DATA INPUT

# the bar: x along the bar from the held end (N & m); q & P are positive along
# +x (away from the support) & N = E*A*du/dx is positive in tension
L = 2       # length (m)
E = 200e9   # Young's modulus (N/m^2)
A0 = 1e-3   # area at x = 0 (m^2)
a = 0.5     # taper: A(x) = A0*(1 - a*x/L) with 0 <= a < 1; a = 0 is the uniform bar of the catalogue
q = 50e3    # distributed axial load, uniform (N/m)
P = 50e3    # load at the free end x = L (N)
ubar = 0    # displacement held at x = 0 (m): 0 as in (P2.4); with another
            # value every u moves by it through the term -K(f,c)*U(c) of
            # solve_prescribed, which test_axial_bar holds with ubar = 2e-3

# the elements
p = 1       # order of the elements: p = 1 the linear element (2 nodes); p = 2 gives the 3-noded (quadratic) element
Ne = 8      # number of elements

ng = 2  # gauss points of an element for [k] & {r}: exact for this bar (degree 3); fewer than p make [K] singular
nge = 8 # gauss points of an element for the error norms (a higher rule: the exact u is not a polynomial)

# the ranges of the input (the closed form & the elements hold only there)
if a < 0 or a >= 1:
    raise ValueError('main_axial_bar: the taper a must be 0 <= a < 1 (a = %g)' % a)
if ng < p:
    raise ValueError('main_axial_bar: ng = %d is fewer than the order p = %d: [K] would be singular' % (ng,p))
if a > 0 and a < 1e-4:
    warnings.warn('main_axial_bar: a taper below 1e-4 leaves the closed form about eps/a of its figures')


# %% AREA & AXIAL STIFFNESS

# A(x) = A0*(1 - a*x/L), k = E*A for KM_Line1D; the exact u, N & du/dx are
# those of UT_AxialBar
Aex = lambda x: A0*(1 - a*x/L)
EAx = lambda x: E*Aex(x)
qx = lambda x: q + 0*x


# %% SHAPE FUNCTIONS OF THE ELEMENT

# in s from -1 to 1, the xi of the catalogue (the same for every element): its
# nodes from end to end, its gauss points & weights for [k] & {r}, & those for
# the error norms
sn = np.linspace(-1,1,p+1)
G,W = gq_points(ng)
Ns,dNs = sf_femlg1d(G,sn)[0:2]
Ge,We = gq_points(nge)
Nse,dNse = sf_femlg1d(Ge,sn)[0:2]


# %% NODES & ELEMENTS

# equal elements, the nodes of each from end to end
nN = p*Ne+1
xi = np.linspace(0,L,nN)
Elem = p*np.arange(Ne)[:,None] + np.arange(p+1)


# %% STIFFNESS MATRIX & FORCE VECTOR

# the integrals of E*A*dNi/dx*dNj/dx & Ni*q over the elements (KM_Line1D), &
# the load at x = L
K,F = km_line1d(xi,Elem,Ns,dNs,W,EAx,qx)
F[nN-1] = F[nN-1] + P


# %% SUPPORT, SOLUTION & REACTION

# the node at x = 0 held at ubar, by partition, K(f,f)*U(f) = F(f) - K(f,c)*U(c)
# & R = K(c,:)*U - F(c)
U,R = solve_prescribed(K,F,[0],[ubar])


# %% COMPARE WITH THE EXACT SOLUTION

# u & du/dx of the elements at the gauss points of the error norms [xQ], with
# their weights [WxQ] (one column per element)
xQ,WxQ,Uq,dUq = sf_line1d(xi,Elem,U,Nse,dNse,We)
Uex,_,dUex = ut_axial_bar(xQ,L,E,A0,a,q,P,ubar)
Uxi = ut_axial_bar(xi,L,E,A0,a,q,P,ubar)[0]
# largest error at the nodes, of the largest displacement
eN = np.max(np.abs(U-Uxi))/np.max(np.abs(Uxi))
# error of the displacement over the bar (L2 norm)
eU = np.sqrt(np.sum(WxQ*(Uq-Uex)**2)/np.sum(WxQ*Uex**2))
# error in the energy norm, the integral of E*A*(du/dx - duex/dx)^2
eE = np.sqrt(np.sum(WxQ*E*Aex(xQ)*(dUq-dUex)**2)/np.sum(WxQ*E*Aex(xQ)*dUex**2))
# equilibrium: the reaction & the applied loads add up to zero. This is a
# check of the assembly & of the partition, not of the accuracy: it holds for
# any U that solves the free equations, on any mesh
eR = np.abs(R[0] + P + q*L)/(np.abs(P) + np.abs(q*L))


# %% DISPLAY RESULTS
print('Axial bar: A(L) = %g*A0, %d elements of order %d (%d nodes)' % (1-a,Ne,p,nN))
print('  u(L) = %.8e m (exact %.8e m)' % (U[nN-1],Uxi[nN-1]))
print('  R    = %.8g N; R + P + q*L = %.1e of the loads' % (R[0],eR))
print('=> Error: at the nodes %.3e, L2 norm %.3e, energy norm %.3e' % (eN,eU,eE))

# the displacement & the axial force against the exact ones: u at the nodes,
# & the axial force E*A*du/dx along every element from end to end (one column
# per element), which is not continuous from one element to the next
sD = np.linspace(-1,1,11)
NsD,dNsD = sf_femlg1d(sD,sn)[0:2]
xD,_,_,dUD = sf_line1d(xi,Elem,U,NsD,dNsD,None)
ND = E*Aex(xD)*dUD
xDe = np.linspace(0,L,201)
UDe,NDe,_ = ut_axial_bar(xDe,L,E,A0,a,q,P,ubar)

fig = pl_axial_bar(xi, U, xD, ND, xDe, UDe, NDe, f"Axial bar: a = {a:g}, {Ne} elements of order {p}")

# %% OPTIONAL EXPORT
output_dir = os.environ.get("P2_OUTPUT_DIR")
if output_dir:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(output / "p2_displacement_python.csv", np.column_stack((xi, U, Uxi)),
               delimiter=",", fmt="%.16g", comments="", header="x_m,fem_m,exact_m")
    elements = np.repeat(np.arange(1, Ne + 1), xD.shape[0])
    forces = np.column_stack((elements, xD.ravel(order="F"), ND.ravel(order="F"), P + q*(L - xD.ravel(order="F"))))
    np.savetxt(output / "p2_force_python.csv", forces, delimiter=",", fmt="%.16g",
               comments="", header="element,x_m,fem_N,exact_N")
    fig.savefig(output / "p2_results_python.png", dpi=160)

plt.show()
