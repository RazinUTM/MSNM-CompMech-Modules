# Title: 2D SF Code considering Lagrange interpolation
# Attempt: to test Lagrange SF & 1st derivatives (fitting)
# Created on: 05/07/2021 by Razin
# Last updated: 05/07/2021 by Razin
# Double checked on: 05/07/2021 by Razin (before the 12/11/2026 update)
# Ported to Python on: 12/11/2026 (line by line from SF_FEMLG.m)
# Last updated: 12/11/2026 by Claude (not yet double checked by Razin)
#   (1) 2nd derivatives dNNxx, dNNxy, dNNyy added; NN, dNNx, dNNy unchanged
# Last updated: 15/11/2026 by Claude (not yet double checked by Razin)
#   file & function named sf_femlg (snake_case in Python; was SF_FEMLG.py)
# Last updated: 18/11/2026 by Claude (not yet double checked by Razin)
#   (2) the 1D functions are worked out for all the points together (A holds
#       one row for each point) & every output is sized before it is filled;
#       the terms of each product & of each sum are those of before, in the
#       same order (3.5 s became 0.11 s for 11 x 11 nodes & 900 points; the six
#       outputs are those of SF_FEMLG.m in GNU Octave in every figure)
# Last updated: 22/11/2026 by Claude (not yet double checked by Razin)
#   (3) a note, no line of code changed. This is ONE Lagrange element through
#       all the nodes of a direction. Through equally spaced nodes its
#       functions grow at the edges (the largest |N| at the gauss points of a
#       plate is 28, 201, 1660 & 16000 for 11, 13, 15 & 17 nodes a side; along
#       an edge 11, 18, 99 & 557 for 12, 13, 16 & 19 nodes), so figures are lost
#       beyond about 12 nodes in a direction: in the plate buckling 13 nodes
#       still give k to 4 figures,
#       & 14 x 14 nodes give NaN (README, Findings 43e). No number of nodes
#       is an exact limit, so there is no warning: keep to 12 or fewer
# Last updated: 23/11/2026 by Claude (not yet double checked by Razin)
#   (4) an error where the nodes are not a full grid, each node once: the
#       columns are the nodes of the grid unique(x) x unique(y), so on other
#       nodes the columns were not the nodes & nothing said so; a note under
#       the function line says in which order the columns come

import numpy as np

# MATLAB operators
eye = lambda A: np.eye(*A.shape,dtype=bool) #logical(eye(size(A)))


def sf_femlg(points, nodes):
    # rows are points; columns are the nodes of the grid unique(x) x unique(y),
    # each sorted, x running first: column k is node k only when the nodes are
    # in that order (sortrows(nodes,[2 1]))

    if len(np.unique(nodes[:,0]))*len(np.unique(nodes[:,1])) != nodes.shape[0] or \
            len(np.unique(nodes,axis=0)) != nodes.shape[0]:
        raise ValueError('SF_FEMLG: the nodes must be a full grid: %d different x & %d different y make %d nodes, these are %d'
                         % (len(np.unique(nodes[:,0])),len(np.unique(nodes[:,1])),len(np.unique(nodes[:,0]))*len(np.unique(nodes[:,1])),nodes.shape[0]))

    # calculate SFunc for 1 element
    # in 1D (x-direction)
    Xp = points[:,0][:,None]
    Xm = np.unique(nodes[:,0])
    Xj = Xm[:,None]
    A = Xp-Xm    #numerator (one row for each point of interest)
    B = Xj-Xm   #denominator
    B[eye(B)] = 1
    Nx = np.zeros((len(Xp),len(Xm)))
    dNx = np.zeros((len(Xp),len(Xm)))
    ddNx = np.zeros((len(Xp),len(Xm)))
    for j in range(len(Xm)): #loops over the nodes of the direction
        dA = A.copy()
        dA[:,j] = 1
        Nx[:,j] = np.prod(dA/B[j,:],axis=1)

        #derivative over-x
        for k in range(len(Xm)): #1st term taken out of the product
            if k!=j:
                dA = A.copy()
                dA[:,[j,k]] = 1
                dNx[:,j] = dNx[:,j] + np.prod(dA/B[j,:],axis=1)

        #2nd derivative over-x
        for k in range(len(Xm)): #2nd term taken out of the product
            if k!=j:
                ddNxjk = np.zeros(len(Xp))
                for l in range(len(Xm)):
                    if l!=j and l!=k:
                        dA = A.copy()
                        dA[:,[j,k,l]] = 1
                        ddNxjk = ddNxjk + np.prod(dA/B[j,:],axis=1)
                ddNx[:,j] = ddNx[:,j] + ddNxjk

    # in 1D (y-direction)
    Yp = points[:,1][:,None]
    Ym = np.unique(nodes[:,1])
    Yj = Ym[:,None]
    A = Yp-Ym    #numerator (one row for each point of interest)
    B = Yj-Ym   #denominator
    B[eye(B)] = 1
    Ny = np.zeros((len(Yp),len(Ym)))
    dNy = np.zeros((len(Yp),len(Ym)))
    ddNy = np.zeros((len(Yp),len(Ym)))
    for j in range(len(Ym)): #loops over the nodes of the direction
        dA = A.copy()
        dA[:,j] = 1
        Ny[:,j] = np.prod(dA/B[j,:],axis=1)

        #derivative over-y
        for k in range(len(Ym)): #1st term taken out of the product
            if k!=j:
                dA = A.copy()
                dA[:,[j,k]] = 1
                dNy[:,j] = dNy[:,j] + np.prod(dA/B[j,:],axis=1)

        #2nd derivative over-y
        for k in range(len(Ym)): #2nd term taken out of the product
            if k!=j:
                ddNyjk = np.zeros(len(Yp))
                for l in range(len(Ym)):
                    if l!=j and l!=k:
                        dA = A.copy()
                        dA[:,[j,k,l]] = 1
                        ddNyjk = ddNyjk + np.prod(dA/B[j,:],axis=1)
                ddNy[:,j] = ddNy[:,j] + ddNyjk

    # in 2D (x-y directions)
    k = -1
    NN = np.zeros((points.shape[0],len(Xm)*len(Ym)))
    dNNx = np.zeros((points.shape[0],len(Xm)*len(Ym)))
    dNNy = np.zeros((points.shape[0],len(Xm)*len(Ym)))
    dNNxx = np.zeros((points.shape[0],len(Xm)*len(Ym)))
    dNNxy = np.zeros((points.shape[0],len(Xm)*len(Ym)))
    dNNyy = np.zeros((points.shape[0],len(Xm)*len(Ym)))
    for j in range(len(Ym)):
        for i in range(len(Xm)):
            k = k + 1
            NN[:,k] = Nx[:,i]*Ny[:,j]
            dNNx[:,k] = dNx[:,i]*Ny[:,j]
            dNNy[:,k] = Nx[:,i]*dNy[:,j]
            dNNxx[:,k] = ddNx[:,i]*Ny[:,j]
            dNNxy[:,k] = dNx[:,i]*dNy[:,j]
            dNNyy[:,k] = Nx[:,i]*ddNy[:,j]

    return NN,dNNx,dNNy,dNNxx,dNNxy,dNNyy
