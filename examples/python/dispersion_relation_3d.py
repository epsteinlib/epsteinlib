"""Quantum spin wave dispersion relation for spins on a 3D square lattice.

This script calculates and plots the dispersion relation of quantum spin waves
on a 3D square lattice (Λ=ℤ³) with power-law long-range interactions. The
dispersion relation ω(k) is computed as a function of k₁ for k₂=k₃=0, using
the Epstein zeta function.

Key features:
1. Models ferromagnetic interactions between spins on a lattice.
2. Demonstrates different scaling behaviors:
   - Typical scaling ω(k) ~ k² for ν ≥ d+2 (green line, ν = 5)
   - Anomalous scaling for d < ν < d+2:
     * Linear dispersion ω(k) ~ k for ν = d+1 (orange line, ν = 4)
     * Square root behavior ω(k) ~ k^(1/2) for ν = 3.5 (blue line)

The dispersion relation is calculated using the Epstein zeta function.
"""

# SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
# SPDX-License-Identifier: AGPL-3.0-only

from typing import List

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from epsteinlib import epstein_zeta


def epstein_zeta_dispersion(nu: float, k: NDArray[np.float64]) -> float:
    """Calculate the Epstein zeta dispersion."""
    return np.real(epstein_zeta(nu, np.eye(len(k)), np.zeros(len(k)), k))


def epstein_zeta_dispersion_tab3d(
    nu: float, qrange: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Generate 3D dispersion table for a given nu and qrange."""
    return np.array(
        [
            (
                q,
                epstein_zeta_dispersion(nu, np.array([0, 0, 0]))
                - epstein_zeta_dispersion(nu, np.array([q, 0, 0])),
            )
            for q in qrange
        ]
    )


def plot_dispersion_relation(
    nurange: NDArray[np.float64], data: List[NDArray[np.float64]]
) -> None:
    """Plot the dispersion relation."""
    frame_label_size = 19
    frame_tick_size = frame_label_size
    k_ticks = np.arange(-0.4, 0.41, 0.2)
    y_ticks = np.arange(0, 41, 10)
    plot_colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]  # Blue, Orange, Green
    line_plot_style = {"color": "black", "alpha": 0.7, "linewidth": 2}

    fig, ax = plt.subplots(figsize=(12, 8))
    fig.suptitle(
        "Quantum Spin Wave Dispersion Relation on 3D Square Lattice",
        fontsize=frame_label_size + 2,
    )
    ax.set_xlabel(r"$k_1$", fontsize=frame_label_size)
    ax.set_ylabel(
        r"$\hbar \,\omega(\mathbf{k})/(JS)$", fontsize=frame_label_size
    )
    ax.set_xlim([-0.51, 0.51])
    ax.set_ylim([-1, 40])
    ax.set_xticks(k_ticks)
    ax.set_yticks(y_ticks)
    ax.tick_params(axis="both", which="major", labelsize=frame_tick_size)

    nu1range = np.linspace(0, 0.28, 100)
    nu2range = np.linspace(0, 0.37, 100)
    nu3range = np.linspace(-0.38, 0, 100)
    ax.plot(nu1range, [59 * nu for nu in nu1range], **line_plot_style)
    ax.plot(nu2range, [53 * np.sqrt(nu) for nu in nu2range], **line_plot_style)
    ax.plot(nu3range, [220 * nu**2 for nu in nu3range], **line_plot_style)

    for i, nu in enumerate(nurange):
        ax.plot(
            *zip(*data[i]),
            label=rf"$\nu = {nu}$",
            color=plot_colors[i],
            linewidth=2,
        )

    ax.text(-0.41, 35, r"$\sim k_1^2$", fontsize=frame_label_size)
    ax.text(0.4, 35, r"$\sim\sqrt{k_1}$", fontsize=frame_label_size)
    ax.text(0.32, 18, r"$\sim k_1 $", fontsize=frame_label_size)

    ax.legend(loc="upper center", fontsize=frame_label_size)
    plt.show()


def main() -> None:
    """Main function to compute data and plot the dispersion relation."""
    nurange = np.array([3.5, 4, 5])
    qrange = np.arange(-0.5, 0.51, 0.01)
    data = [epstein_zeta_dispersion_tab3d(nu, qrange) for nu in nurange]
    plot_dispersion_relation(nurange, data)


if __name__ == "__main__":
    main()
