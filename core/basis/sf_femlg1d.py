# Title: 1D SF Code considering Lagrange interpolation
# Attempt: to give the Lagrange SF of a line with their 1st & 2nd
#          derivatives, for the elements of a 1D problem
# Created on: 23/11/2026 by Claude (not yet double checked by Razin)
#   on the author's word of 23/11/2026 (his answers on P2 of the ladder): a
#   1D problem wrote the points & the nodes as [s 0] to call SF_FEMLG, which
#   is now done here, once. The numbers are those of SF_FEMLG: this file
#   calls it with a 2nd coordinate of zero, so there is one implementation
#   & the outputs are the same in every figure
# Last updated: 23/11/2026 by Claude
# Double checked on: -
# Ported to Python on: 23/11/2026 (line by line from SF_FEMLG1D.m)
#   points & nodes may be a number, a list or an array of any shape: they
#   are read as vectors (np.reshape), as points(:) & nodes(:) do

import numpy as np
from .sf_femlg import sf_femlg


def sf_femlg1d(points, nodes):
    # UDF for the 1D Lagrange functions (N) through the nodes, with their 1st
    # (dN) & 2nd (ddN) derivatives over the coordinate, at the points
    #   points = the coordinates of the points, a vector (row or column)
    #   nodes  = the coordinates of the nodes, a vector (row or column)
    # rows are points; columns are the nodes in ascending order (unique(nodes)):
    # column k is node k only when the nodes are given in that order. A node
    # given twice is the error of SF_FEMLG. The note of SF_FEMLG holds here too:
    # keep to 12 nodes or fewer

    points = np.reshape(np.asarray(points,dtype=float),-1)
    nodes = np.reshape(np.asarray(nodes,dtype=float),-1)

    # SF_FEMLG on one row of nodes [x 0] (the 3rd output, over y, is zero)
    N,dN,_,ddN = sf_femlg(np.c_[points, np.zeros(points.shape)],np.c_[nodes, np.zeros(nodes.shape)])[0:4]

    return N,dN,ddN
