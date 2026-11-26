# Title: Code for the matrix & the vector of a line of 1D elements
# Attempt: to assemble [K] & {F} of -d/dx(k*du/dx) = Q on 1D elements of
#          any order, as a kernel
# Created on: 24/11/2026 by Claude (not yet double checked by Razin)
#   on the author's word of 24/11/2026 (the kernels that core/ has not are
#   written, modular, & kept in core/): the loop over the elements that
#   MAIN_AxialBar (k = E*A, Q = q) & MAIN_SteadyConduction1D (k = k*A,
#   Q = Q*A) each wrote (README, open item 22 (b)), which now call this file;
#   the functions of the element are those of SF_FEMLG on one row of nodes,
#   its points those of GQ_Points
# Last updated: 24/11/2026 by Claude
# Double checked on: -
# Ported to Python on: 24/11/2026 (line by line from KM_Line1D.m)
#   the nodes are counted from 0 (Elem); the gauss points & weights are 1D
#   arrays, so a product down the rows of a matrix is written with [:,None]

import numpy as np
from .as_global import as_global


def km_line1d(xi,Elem,Ns,dNs,W,kx,Qx):
    # UDF for the matrix [K] & the vector {F} of -d/dx(k(x)*du/dx) = Q(x) on a
    # line of 1D elements, by gauss quadrature on each element
    #   K(i,j) = sum over the elements of the integral of k*dNi/dx*dNj/dx dx
    #   F(i)   = sum over the elements of the integral of Ni*Q dx
    # with x = N(s)*Xe, dx/ds = dN/ds*Xe & dN/dx = (dN/ds)/(dx/ds) on an element
    # xi: the nodes (x), nN of them; Elem: the nodes of each element from end to
    # end, one row per element; Ns, dNs: the functions of the element & their
    # derivatives in s (from -1 to 1) at its gauss points, one row per point, with
    # their weights W; kx, Qx: k(x) & Q(x), functions of an array of x.
    # Nothing of a face is added: a load at a node or a film coefficient is added
    # by the caller

    xi = np.asarray(xi,dtype=float).ravel()
    Elem = np.atleast_2d(Elem)   # one element given as one row
    W = np.asarray(W,dtype=float).ravel()
    nN = xi.size
    K = np.zeros((nN,nN))
    F = np.zeros(nN)
    for el in range(Elem.shape[0]):
        # loop over elements
        en = Elem[el,:]
        Xe = xi[en]        # nodes of the element
        J = dNs@Xe         # dx/ds at the gauss points (half the length)
        if np.any(J <= 0):
            raise ValueError('km_line1d: dx/ds <= 0 in element %d (its nodes must run from end to end, x increasing)' % (el+1))
        x = Ns@Xe          # the gauss points on the line
        B = dNs/J[:,None]  # dN/dx
        ke = B.T@((W*J*kx(x))[:,None]*B) # integral of k*dNi/dx*dNj/dx
        re = Ns.T@(W*J*Qx(x))            # integral of Ni*Q
        as_global(K,ke,en)
        F[en] = F[en] + re

    return K,F
