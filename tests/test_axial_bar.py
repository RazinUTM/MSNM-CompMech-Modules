# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Ahmad Razin Zainal Abidin and Universiti Teknologi Malaysia
"""Tests of the axial bar by 1D finite elements (P2 of the ladder).

The same checks as tests/Validation/TEST_AxialBar.m: the element matrices by the lines of
the main file against the closed forms of the catalogue (P2.3), and the main file, run in a
process of its own from its own folder (the fixture run_driver), against the exact solution
of the tapered bar on its one mesh; then each term alone, with lines of its input changed in
memory (run_driver_changed), and its guards. The main file as it stands has ubar = 0: the term
-K(f,c)*U(c) of its partition is zero there; the partition is sv_partition (core/), held with
values away from zero by tests/Kernels/test_line1d.py, as are km_line1d and sf_line1d, and
ut_axial_bar is held by tests/Kernels/test_exact.py. That the main file hands ubar to it is
held here, with ubar = 2e-3 on the uniform bar.

Since 24 November 2026 the main file solves one bar on one mesh (the author's answers on
P2). What this file held of the four bars and the five meshes (the uniform bar, the rates,
the values by hand on 2 elements) is in test_axial_bar_study.py, each with its limit.
"""

import numpy as np

from core.assembly.km_line1d import km_line1d
from core.quadrature.gq_points import gq_points
from core.basis.sf_femlg1d import sf_femlg1d


def test_element_matrices_against_the_closed_forms():
    # one element of length l with constant EA & q, by the lines of the main file:
    # sf_femlg1d gives the 1D Lagrange functions on the nodes of the element in s
    l, EA, q = 0.4, 3.0, 7.0
    G, W = gq_points(2)
    kex = [EA / l * np.array([[1, -1], [-1, 1]]), EA / (3 * l) * np.array([[7, -8, 1], [-8, 16, -8], [1, -8, 7]])]
    rex = [q * l / 2 * np.array([1, 1]), q * l / 6 * np.array([1, 4, 1])]
    for p in (1, 2):
        sn = np.linspace(-1, 1, p + 1)
        Ns, dNs = sf_femlg1d(G, sn)[0:2]
        Xe = 1.3 + np.linspace(0, l, p + 1)
        ke, re = km_line1d(Xe, np.arange(p + 1), Ns, dNs, W,
                           lambda x: EA + 0*x, lambda x: q + 0*x)
        # observed: 3.6e-16 & 1.6e-16 (linear), 5.3e-16 & 4.8e-16 (quadratic; 3.6e-16 &
        # 3.3e-16 in GNU Octave)
        # (in MATLAB on 23/11/2026, through SF_FEMLG: 3.6e-16 & 1.6e-16, 3.6e-16 & 3.3e-16)
        assert np.max(np.abs(ke - kex[p - 1])) < 1e-12 * np.max(np.abs(kex[p - 1]))
        assert np.max(np.abs(re - rex[p - 1])) < 1e-12 * np.max(np.abs(rex[p - 1]))


def test_axial_bar_against_the_exact_solution(run_driver):
    g, _ = run_driver("problems/AxialBar", "main_axial_bar.py",
                      ("a", "p", "Ne", "nN", "eN", "eU", "eE", "eR", "F", "U", "R", "Uxi", "WxQ", "xD", "ND"))
    eN, eU, eE, eR = g["eN"], g["eU"], g["eE"], g["eR"]
    F, U, Uxi = np.array(g["F"]), np.array(g["U"]), np.array(g["Uxi"])

    # the main file as it stands: the tapered bar (a = 0.5) on 8 linear elements, 9 nodes
    assert g["a"] == 0.5 and g["p"] == 1 and g["Ne"] == 8 and g["nN"] == 9 and U.shape == (9,)
    # the closed form by hand: the exact u at the nodes is 0 at the held end & at the free
    # end (a = 0.5) u(L) = 2e-3 - 1e-3*log(2) m. observed: 0 & 0
    assert Uxi[0] == 0 and abs(Uxi[-1] - (2e-3 - 1e-3 * np.log(2))) < 1e-15
    # the held end stays where it is held (ubar = 0)
    assert U[0] == 0
    # the reaction & the loads add up to zero: R = -(P + q*L) = -150000 N. observed:
    # 1.6e-15 of the loads (1.4e-15 in GNU Octave)
    assert eR < 1e-10
    assert abs(g["R"][0] + 150000) < 1e-10 * 150000
    # the weights of the error norms add up to the length of the bar. observed: 1.3e-15
    # (0 in GNU Octave)
    assert abs(np.sum(np.array(g["WxQ"])) - 2) < 1e-12
    # tapered bar: not exact at the nodes. observed: 3.7e-4
    assert eN > 1e-9
    # the energy norm by another route: the error of a Galerkin solution has the energy
    # a(e,e) = a(u,u) - F.U (ubar = 0), with a(u,u) = 100 + 50*log(2) J by hand & F.U of
    # the main file (134.633 J). observed: 4.6e-12 of it (3.0e-12 in GNU Octave)
    assert abs(eE - np.sqrt(1 - (F @ U) / (100 + 50 * np.log(2)))) < 1e-8 * eE
    # recorded on 23/11/2026 (Python & GNU Octave agree to 3e-11 of them): the L2 & the
    # energy norms of the error. They are those the file of before gave for this bar &
    # mesh, in every figure; test_axial_bar_study holds the same two as an entry of the study
    # (in MATLAB on 23/11/2026, the file of before: the two in every figure printed)
    # (in MATLAB on 24/11/2026, this file: the two in every figure printed; the reaction
    # 1.9e-16, the energy norm by F.U 1.3e-12, each term alone 1.6e-15 at most, ubar = 2e-3
    # 1.4e-16, the uniform strain 6.2e-16 at most)
    rec = np.array([7.503743828e-04, 1.344212226e-02])
    assert abs(eU - rec[0]) < 1e-7 * rec[0] and abs(eE - rec[1]) < 1e-7 * rec[1]
    # the figure: the axial force along every element, 11 points on each of the 8
    assert np.array(g["xD"]).shape == (11, 8) and np.array(g["ND"]).shape == (11, 8)


def test_axial_bar_each_term_alone(run_driver_changed):
    # in the data of the main file q = P = 50e3 & q*L/2 = P, so a load put for another would not
    # show: the distributed load alone, the end load alone & the two, on a uniform bar of 3 m
    # (a = 0), whose nodes the elements give exactly: u(L) against (P2.4),
    # u(L) = q*L^2/(2*E*A0) + P*L/(E*A0). observed: 3.2e-16 & 3.2e-16 at most
    # (4.8e-16 & 4.8e-16 in GNU Octave)
    for q, P in ((0, 70e3), (30e3, 0), (30e3, 70e3)):
        g, err = run_driver_changed("problems/AxialBar", "main_axial_bar.py",
                                    [("L = 2       #", "L = 3       #"), ("a = 0.5     #", "a = 0       #"),
                                     ("q = 50e3    #", "q = %r    #" % q), ("P = 50e3    #", "P = %r    #" % P)], ("U", "eN", "eR"))
        assert err is None, err
        uni = q * 3**2 / (2 * 2e8) + P * 3 / 2e8
        assert abs(g["U"][-1] - uni) < 1e-13 * uni and g["eN"] < 1e-10 and g["eR"] < 1e-10
    # & the support away from zero: with ubar = 2e-3 m every u of the uniform bar moves by it,
    # through the term -K(f,c)*U(c) of the partition, which is zero in the data of the main
    # file: u(L) = 2e-3 + 1e-3 m. observed: 1.2e-15
    # (1.4e-15 in GNU Octave)
    g, err = run_driver_changed("problems/AxialBar", "main_axial_bar.py",
                                [("a = 0.5     #", "a = 0       #"), ("ubar = 0    #", "ubar = 2e-3 #")], ("U", "eN", "eR"))
    assert err is None, err
    assert g["U"][0] == 2e-3 and abs(g["U"][-1] - 3e-3) < 1e-13 * 3e-3 and g["eN"] < 1e-10 and g["eR"] < 1e-10


def test_axial_bar_uniform_strain(run_driver_changed):
    # P = q*L*(1/a - 1) = 100e3 N: N/A is constant, the strain uniform, & the elements are exact
    # whatever their EA (the header of ut_axial_bar). observed: 1.5e-15, 1.7e-15
    # & 1.7e-15 (1.2e-15, 1.5e-15 & 1.7e-15 in GNU Octave)
    g, err = run_driver_changed("problems/AxialBar", "main_axial_bar.py", [("P = 50e3    #", "P = 100e3    #")], ("eN", "eU", "eE"))
    assert err is None, err
    assert g["eN"] < 1e-10 and g["eU"] < 1e-10 and g["eE"] < 1e-10


def test_axial_bar_guards(run_driver_changed):
    # a taper of 1, & fewer gauss points than the order, stop the main file with a message
    g1, m1 = run_driver_changed("problems/AxialBar", "main_axial_bar.py", [("a = 0.5     #", "a = 1       #")])
    g2, m2 = run_driver_changed("problems/AxialBar", "main_axial_bar.py", [("p = 1       #", "p = 2       #"), ("ng = 2  #", "ng = 1  #")])
    assert g1 is None and "main_axial_bar: the taper a must be 0 <= a < 1" in m1
    assert g2 is None and "[K] would be singular" in m2
