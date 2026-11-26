# SPDX-License-Identifier: MIT
"""Plot the axial-bar solution and its closed-form reference."""

import matplotlib.pyplot as plt


def pl_axial_bar(x_nodes, displacement, x_element, force, x_exact,
                 exact_displacement, exact_force, title):
    fig, (ax_u, ax_n) = plt.subplots(2, 1, figsize=(8, 6), layout="constrained")
    ax_u.plot(x_exact, exact_displacement, "k-", label="Closed form")
    ax_u.plot(x_nodes, displacement, "ro", label="FEM at nodes")
    ax_u.set(xlabel="x (m)", ylabel="Displacement (m)", title=title)
    ax_u.legend()
    ax_n.plot(x_exact, exact_force, "k-", label="Closed form")
    segments = ax_n.plot(x_element, force, "r-")
    segments[0].set_label("FEM by element")
    ax_n.set(xlabel="x (m)", ylabel="Axial force (N)")
    ax_n.legend()
    return fig
