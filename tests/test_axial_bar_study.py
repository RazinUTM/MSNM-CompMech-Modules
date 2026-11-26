# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Ahmad Razin Zainal Abidin and Universiti Teknologi Malaysia
"""Test of the study of the axial bar (studies/AxialBar/; P2 of the ladder).

The same checks as tests/Validation/TEST_AxialBarStudy.m. The study is run as it stands, in a
process of its own from its own folder (the fixture run_driver), and so is the main file of
the problem. Each limit is stated with the value observed in Python and in GNU Octave 8.4.0
on 24 November 2026. What is held: the four bars on the five meshes against the exact
solution (the uniform bar exact at the nodes, the quadratic elements exact everywhere; the
rates; the values worked by hand; the four recorded errors), which test_axial_bar.py held
until the main file became one bar on one mesh, each with the limit it had there; and, new
with the study, the rate of the error at the nodes, and the main file as it stands as the
entry [a = 0.5, p = 1, Ne = 8] of the study; then each term alone, with lines of the input of
the study changed in memory (run_driver_changed), and its guards.
"""

import numpy as np

from core.quadrature.gq_points import gq_points

FOLDER = "studies/AxialBar"


def test_axial_bar_convergence(run_driver):
    g, _ = run_driver(FOLDER, "study_axial_bar_convergence.py",
                      ("eN", "eU", "eE", "eR", "UL", "RateN", "RateU", "RateE", "ND", "WxQ", "Bar", "Ne", "Draw"))
    eN, eU, eE, eR, UL = (np.array(g[k]) for k in ("eN", "eU", "eE", "eR", "UL"))
    RateN, RateU, RateE, ND = (np.array(g[k], dtype=float) for k in ("RateN", "RateU", "RateE", "ND"))

    # the study as it stands: rows of Bar = uniform linear, uniform quadratic, tapered
    # linear, tapered quadratic; five meshes; the solution drawn is the tapered bar on 2
    # linear elements
    assert g["Bar"] == [[0, 1], [0, 2], [0.5, 1], [0.5, 2]] and g["Ne"] == [2, 4, 8, 16, 32]
    assert g["Draw"] == [3, 1]
    assert eN.shape == (4, 5) and UL.shape == (4, 6) and RateN.shape == (4, 4)
    # the closed forms by hand: uniform, u(L) = q*L^2/(2*E*A) + P*L/(E*A) = 5e-4 + 5e-4 m;
    # tapered (a = 0.5), u(L) = 2e-3 - 1e-3*log(2) m
    assert np.max(np.abs(UL[0:2, -1] - 1e-3)) < 1e-15 and np.max(np.abs(UL[2:4, -1] - (2e-3 - 1e-3 * np.log(2)))) < 1e-15
    # u(L) of the elements: at the closed form on every mesh of the uniform bar; by hand for
    # the tapered bar on 2 linear elements, 125e3/175e6 = 1/1400 m at the middle node &
    # 1/1400 + 75e3/125e6 = 23/17500 m at the end. observed: 9.6e-17 & 2.2e-19 (1.5e-16 & 0
    # in GNU Octave)
    assert np.max(np.abs(UL[0:2, :-1] - 1e-3)) < 1e-13 and abs(UL[2, 0] - 23 / 17500) < 1e-13
    # the reaction & the loads add up to zero on every mesh. observed: 9.2e-14 (1.4e-13 in
    # GNU Octave)
    assert np.all(eR < 1e-10)
    # the weights of the error norms (of the last bar & mesh) add up to the length of the
    # bar. observed: 0 (8.9e-16 in GNU Octave)
    assert abs(np.sum(np.array(g["WxQ"])) - 2) < 1e-12
    # uniform bar: the elements are exact at the nodes, & the quadratic ones everywhere, so
    # that bar has no rate. observed: 9.6e-14, 1.0e-13 & 1.0e-13 (1.5e-13 in GNU Octave)
    # (in MATLAB on 23/11/2026, with these lines in the main file: 1.5e-13; the reaction 1.3e-13)
    assert np.all(eN[0:2, :] < 1e-10)
    assert np.all(eU[1, :] < 1e-10) and np.all(eE[1, :] < 1e-10)
    assert np.all(np.isnan(RateU[1, :])) and np.all(np.isnan(RateE[1, :]))
    # & the nodes of the uniform bar have no rate, with either element
    assert np.all(np.isnan(RateN[0:2, :]))
    # uniform bar, linear elements: the error between the nodes is the same curve on every
    # element, q/(2EA)*s*(l-s), so the rates are 2 & 1 to rounding. observed: 3.2e-11 &
    # 1.4e-15 from them (5.9e-11 in GNU Octave)
    assert np.all(np.abs(RateU[0, :] - 2) < 1e-8) and np.all(np.abs(RateE[0, :] - 1) < 1e-8)
    # & the errors themselves by hand: the L2 norm is l^2/(2*sqrt(204)) & the energy norm
    # l/sqrt(52), with l = L/Ne. observed: 2.6e-11 & 1.2e-15 of them (5.4e-11 & 4.0e-16 in
    # GNU Octave)
    Ne = np.array(g["Ne"])
    l = 2 / Ne
    assert np.all(np.abs(eU[0, :] - l**2 / (2 * np.sqrt(204))) < 1e-8 * eU[0, :])
    assert np.all(np.abs(eE[0, :] - l / np.sqrt(52)) < 1e-8 * eE[0, :])
    # tapered bar: no longer exact at the nodes. observed: 2.3e-5 at the least (linear),
    # 2.7e-9 (quadratic)
    assert np.all(eN[2:4, :] > 1e-9)
    # by hand on 2 linear elements: the largest error at the nodes is at x = L, 23/17500
    # against 2e-3 - 1e-3*log(2). observed: 1.7e-16 from it (0 in GNU Octave)
    assert abs(eN[2, 0] - (23 / 17500 - (2e-3 - 1e-3 * np.log(2))) / (2e-3 - 1e-3 * np.log(2))) < 1e-10
    # & the energy norm there by hand: a(e,e) = a(u,u) - F.U, with a(u,u) = 100 + 50*log(2) J
    # & F.U = 50e3/1400 + 75e3*23/17500 = 940/7 J. observed: 2.4e-14 of it
    assert abs(eE[2, 0] - np.sqrt(1 - (940 / 7) / (100 + 50 * np.log(2)))) < 1e-8 * eE[2, 0]
    # tapered bar: the rates between the last two meshes are p+1 in the L2 norm & p in the
    # energy norm. observed: 1.9989 & 0.9995, 2.9978 & 1.9986
    assert abs(RateU[2, -1] - 2) < 0.01 and abs(RateE[2, -1] - 1) < 0.01
    assert abs(RateU[3, -1] - 3) < 0.01 and abs(RateE[3, -1] - 2) < 0.01
    # & every error falls from one mesh to the next
    assert np.all(np.diff(eU[[0, 2, 3], :], axis=1) < 0) and np.all(np.diff(eE[[0, 2, 3], :], axis=1) < 0)
    # recorded on 23/11/2026 (Python & GNU Octave agree to 3e-11 of them): the tapered bar
    # on 8 elements, linear & quadratic
    # (in MATLAB on 23/11/2026, with these lines in the main file: the four in every figure
    # printed, & the rates)
    # (in MATLAB on 24/11/2026, the study: the four & every rate in every figure printed; the
    # main file against its entry 0; each term alone 1.7e-13 & 1.0e-15 at most)
    rec = np.array([[7.503743828e-04, 1.240235526e-05], [1.344212226e-02, 3.418021231e-04]])
    assert np.all(np.abs(eU[2:4, 2] - rec[0, :]) < 1e-7 * rec[0, :]) and np.all(np.abs(eE[2:4, 2] - rec[1, :]) < 1e-7 * rec[1, :])
    # the solution drawn (tapered, 2 linear elements): the axial force of the elements is
    # not continuous at the node between them; by hand it jumps from 750000/7 to 90000 N
    # there (the exact force is 100000 N)
    assert abs(ND[-1, 0] - 750000 / 7) < 1e-4 and abs(ND[0, 1] - 90000) < 1e-4

    # new with the study: the rate of the error at the nodes.
    # the rates that the study prints & draws are those of its errors. observed: 0
    dN = np.log(Ne[1:] / Ne[:-1])
    assert np.max(np.abs(RateN[2:4, :] - np.log(eN[2:4, :-1] / eN[2:4, 1:]) / dN)) < 1e-12
    assert np.max(np.abs(RateU[[0, 2, 3], :] - np.log(eU[[0, 2, 3], :-1] / eU[[0, 2, 3], 1:]) / dN)) < 1e-12
    assert np.max(np.abs(RateE[[0, 2, 3], :] - np.log(eE[[0, 2, 3], :-1] / eE[[0, 2, 3], 1:]) / dN)) < 1e-12
    # tapered bar: the error at the nodes falls on every mesh, at the order 2p: 2 with the
    # linear elements & towards 4 with the quadratic ones, whose rate still rises on these
    # meshes. observed: 1.9474, 1.9857, 1.9963, 1.9991 (linear); 3.3885, 3.6698, 3.8291,
    # 3.9132 (quadratic; 3.9131 in GNU Octave, where the error is 2.7e-9 of u)
    assert np.all(np.diff(eN[2:4, :], axis=1) < 0)
    assert abs(RateN[2, -1] - 2) < 0.01
    assert np.all(np.diff(RateN[3, :]) > 0) and abs(RateN[3, -1] - 4) < 0.15
    # so the nodes are better than the bar between them: the error at the nodes over the
    # L2 norm falls with the mesh for the quadratic elements. observed: 0.107 to 0.0137
    assert np.all(np.diff(eN[3, :] / eU[3, :]) < 0)


def test_main_file_is_an_entry_of_the_study(run_driver):
    # the main file as it stands (one bar on one mesh) is the entry of the study with its
    # taper a, its order p & its number of elements Ne: [a = 0.5, p = 1, Ne = 8], row 3 of
    # Bar & the third mesh
    s, _ = run_driver(FOLDER, "study_axial_bar_convergence.py", ("eN", "eU", "eE", "eR", "UL", "Bar", "Ne"))
    g, _ = run_driver("problems/AxialBar", "main_axial_bar.py", ("a", "p", "Ne", "eN", "eU", "eE", "eR", "U", "Uxi"))
    b = [row == [g["a"], g["p"]] for row in s["Bar"]].index(True)
    m = s["Ne"].index(g["Ne"])
    assert (b, m) == (2, 2)
    # the three errors & u(L), of the elements & exact. observed: 0, in every figure, in
    # Python & in GNU Octave (the same lines on the same kernels)
    for k in ("eN", "eU", "eE"):
        assert abs(g[k] - s[k][b][m]) < 1e-12 * g[k]
    assert abs(g["U"][-1] - s["UL"][b][m]) < 1e-12 * abs(g["U"][-1])
    assert abs(g["Uxi"][-1] - s["UL"][b][-1]) < 1e-12 * abs(g["Uxi"][-1])
    # the reaction with the loads is the same rounding in both. observed: no difference
    # (each is 1.6e-15; 1.4e-15 in GNU Octave)
    assert abs(g["eR"] - s["eR"][b][m]) < 1e-10


def test_axial_bar_study_each_term_alone(run_driver_changed):
    # in the data of the study q = P = 50e3 & q*L/2 = P, so a load put for another would not
    # show: the distributed load alone, the end load alone & the two, on a bar of 3 m. The
    # uniform bar against (P2.4), u(L) = q*L^2/(2*E*A0) + P*L/(E*A0); the tapered bar against an
    # independent integral of N/(E*A) by 40 gauss points. observed: 1.9e-14 & 7.9e-16 at most
    # (1.7e-13 & 4.5e-16 in GNU Octave; 2.4e-14 there until the assembly was KM_Line1D, which
    # takes E*A as one function: [K] moves by an ulp & u(L) of 32 quadratic elements by
    # 1.3e-13, rounding, so the limit is 1e-12 since 24/11/2026)
    t40, w40 = gq_points(40)
    s40 = 3 / 2 * (t40 + 1)
    for q, P in ((0, 70e3), (30e3, 0), (30e3, 70e3)):
        g, err = run_driver_changed(FOLDER, "study_axial_bar_convergence.py",
                                    [("L = 2       #", "L = 3       #"), ("q = 50e3    #", "q = %r    #" % q),
                                     ("P = 50e3    #", "P = %r    #" % P)], ("UL", "eR"))
        assert err is None, err
        UL, eR = np.array(g["UL"]), np.array(g["eR"])
        uni = q * 3**2 / (2 * 2e8) + P * 3 / 2e8
        tap = 3 / 2 * np.sum(w40 * (P + q * (3 - s40)) / (2e8 * (1 - 0.5 * s40 / 3)))
        assert np.all(np.abs(UL[0:2, :] - uni) < 1e-12 * uni) and np.all(np.abs(UL[2:4, -1] - tap) < 1e-13 * tap) and np.all(eR < 1e-10)


def test_axial_bar_study_uniform_strain(run_driver_changed):
    # P = q*L*(1/a - 1) = 100e3 N: N/A is constant, the strain uniform, & the elements are exact
    # whatever their EA (the header). observed: 9.9e-15 & 2.0e-14
    g, err = run_driver_changed(FOLDER, "study_axial_bar_convergence.py", [("P = 50e3    #", "P = 100e3    #")], ("eU", "eE", "RateU"))
    assert err is None, err
    assert np.all(np.array(g["eU"])[2:4, :] < 1e-10) and np.all(np.array(g["eE"])[2:4, :] < 1e-10)
    assert np.all(np.isnan(np.array(g["RateU"], dtype=float)[2:4, :]))


def test_axial_bar_study_guards(run_driver_changed):
    # a taper of 1 & fewer gauss points than the order stop the study, with a message
    g1, m1 = run_driver_changed(FOLDER, "study_axial_bar_convergence.py", [("                [0.5, 1],", "                [1,   1],")])
    g2, m2 = run_driver_changed(FOLDER, "study_axial_bar_convergence.py", [("ng = 2  #", "ng = 1  #")])
    assert g1 is None and "taper a must be 0 <= a < 1 (row 3 of Bar)" in m1
    assert g2 is None and "[K] would be singular" in m2
