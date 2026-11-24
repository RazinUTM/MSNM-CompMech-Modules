# SPDX-License-Identifier: MIT
"""Partitioned solution of a linear system with prescribed DOFs."""

import numpy as np


def solve_prescribed(K, F, constrained, prescribed):
    """Return full displacements and reactions for one load case.

    K is a square stiffness matrix and F is a one-dimensional load vector.
    constrained contains unique zero-based DOF indices. prescribed gives
    their values in the same order. Reactions follow that order and include
    any applied loads at constrained DOFs.
    """
    K = np.asarray(K)
    F = np.asarray(F)
    c = np.asarray(constrained)
    values = np.asarray(prescribed)
    if F.ndim != 1 or K.shape != (F.size, F.size):
        raise ValueError("K must be square and F must match its size")
    if c.ndim != 1 or values.shape != c.shape:
        raise ValueError("Provide one prescribed value per constrained DOF")
    if c.size and (not np.issubdtype(c.dtype, np.integer)
                   or np.any(c < 0) or np.any(c >= F.size)
                   or np.unique(c).size != c.size):
        raise ValueError("Constrained DOFs must be unique, valid integer indices")
    c = c.astype(np.intp)
    f = np.setdiff1d(np.arange(F.size), c)
    U = np.zeros(F.size, dtype=np.result_type(K, F, values, float))
    U[c] = values
    if f.size:
        U[f] = np.linalg.solve(K[np.ix_(f, f)], F[f] - K[np.ix_(f, c)] @ U[c])
    return U, K[c, :] @ U - F[c]
