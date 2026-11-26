# Ported to Python on: 11/11/2026 (line by line from GQ_Points.m)
# Last updated: 11/11/2026 by Claude (not yet double checked by Razin)
#   points & weights to full precision from numpy (the table had 5
#   significant digits, n = 1 to 9); same call [G,W] = GQ_Points(n) and
#   same order as the table (positive points first)
# Last updated: 15/11/2026 by Claude (not yet double checked by Razin)
#   file & function named gq_points (snake_case in Python; was GQ_Points.py)
# Last updated: 23/11/2026 by Claude (not yet double checked by Razin)
#   (1) import numpy as np is taken out: np is not used here
#   (2) a note, no line of code changed. n below 1 or not whole gives numpy's
#       error, not the message of GQ_Points.m; so does a whole n given as a
#       float (3.0), which GQ_Points takes

from numpy.polynomial.legendre import leggauss


def gq_points(n):

    G,W = leggauss(n)
    G = G[::-1]
    W = W[::-1]

    return G,W
