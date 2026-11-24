# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Ahmad Razin Zainal Abidin and Universiti Teknologi Malaysia
"""Tests of the system of springs by the direct stiffness method (P1 of the ladder).

The same checks as TEST_SpringSystem.m, against the closed form for
zero and nonzero prescribed support displacements.
"""

from pathlib import Path
import runpy

import numpy as np


def test_spring_system_against_the_closed_form(monkeypatch):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    monkeypatch.setattr(plt, "show", lambda: None)
    monkeypatch.delenv("P1_OUTPUT_DIR", raising=False)
    driver = Path(__file__).resolve().parents[1] / "problems/SpringSystem/main_spring_system.py"
    try:
        g = runpy.run_path(str(driver))
        assert len(g["fig"].axes) == 4
    finally:
        plt.close("all")
    U, Uex, R = np.array(g["U"]), np.array(g["Uex"]), np.array(g["R"])
    eU, eC, eP, eR = (np.array(g[k]) for k in ("eU", "eC", "eP", "eR"))

    # before the support every row of [K] adds up to zero. observed: 0
    assert g["Krow"] < 1e-15
    # the closed form by hand for u1 = 0: u2 = 100/2000 = 0.05, u3 = 0.05 + 100/500 = 0.25
    assert np.max(np.abs(Uex[:, 0] - [0, 0.05, 0.25])) < 1e-15
    # by partition, for u1 = 0 and u1 = 0.01 m. observed: 0 and 2.7e-17
    assert np.all(eU < 1e-14)
    # the reaction is -F, and with the loads it adds up to zero. observed: 0
    assert np.max(np.abs(R - g["Rex"])) < 1e-12 and np.all(eR < 1e-14)
    # without the term -K(f,c)U(c) the answer is right only for u1 = 0: the check catches the
    # mistake, which leaves u off by u1. observed: 0 and 3.8e-2
    assert eC[0] < 1e-14 and abs(eC[1] - g["ubar"][1] / np.max(np.abs(Uex[:, 1]))) < 1e-12
    # the penalty is off by about F/(Alpha max(diag(K))) and conditions [K] worse.
    # observed: 1.6e-7 and 1.5e-7; 6.85 and 6.55e6
    assert np.all((eP > 1e-8) & (eP < 1e-6))
    assert abs(g["condKff"] - 6.85) < 0.01 and g["condKp"] > 1e5 * g["condKff"]
    np.testing.assert_allclose(U[:, 1], [0.01, 0.06, 0.26], atol=1e-14, rtol=0)
    np.testing.assert_allclose(g["extension"], [[0.05, 0.05], [0.20, 0.20]], atol=1e-14, rtol=0)
    np.testing.assert_allclose(g["spring_force"], 100, atol=1e-12, rtol=0)
    np.testing.assert_allclose(g["equilibrium_residual"], 0, atol=1e-12, rtol=0)
