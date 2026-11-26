# SPDX-License-Identifier: MIT
"""Accumulate one element matrix into a preallocated global matrix."""

import numpy as np


def as_global(K, ke, dofs):
    """Add ke to K in place using zero-based global DOFs."""
    np.add.at(K, np.ix_(dofs, dofs), ke)
