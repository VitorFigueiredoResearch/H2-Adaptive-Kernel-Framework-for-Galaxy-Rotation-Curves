"""
generate_acceleration_space_plot.py
====================================
Generates the acceleration-space diagnostic figure for the H2 paper.

For each of the 80 validated SPARC galaxies, loads the SPARC rotation-curve
decomposition file (*_rotmod.dat) and computes:

    g_bar(r) = [V_gas(r)^2 + V_disk(r)^2 + V_bul(r)^2] / r
    g_tot(r) = V_obs(r)^2 / r

in units of m/s^2.  All data points with V_obs > 0, r > 0, and
V_obs/errV > 2 are retained.

The McGaugh et al. (2016) RAR is overplotted as:
    g_tot = g_bar / (1 - exp(-sqrt(g_bar / g_dagger)))
with g_dagger = 1.2e-10 m/s^2.

Points are colour-coded by the H2 dynamical regime classification
(baryon-dominated / balanced / DM-dominated) from fleet_summary_80galaxy.csv.

Output
------
acceleration_space_80galaxy.pdf  (in the same directory as this script, or
                                   the path given by --output)

Usage
-----
    python generate_acceleration_space_plot.py [--output PATH]

Run from the H2 repository root or any working directory.  Paths are
resolved relative to this script's location.
"""

import os
import sys
import argparse
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")           # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ---------------------------------------------------------------------------
# Physical constants / conversion factors
# ---------------------------------------------------------------------------
KPC_TO_M   = 3.08567758e19     # 1 kpc in metres
KMS_TO_MS  = 1.0e3             # 1 km/s in m/s
G_DAGGER   = 1.2e-10           # McGaugh+16 RAR characteristic acceleration [m/s^2]

# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

FLEET_CSV = os.path.join(
    REPO_ROOT,
    "H2_PUBLICATION_RELEASE",
    "fleet_expansion_80galaxy",
    "fleet_summary_80galaxy.csv",
)
SPARC_DIR = os.path.join(REPO_ROOT, "data", "sparc")

# ---------------------------------------------------------------------------
# Regime colours and labels
# ---------------------------------------------------------------------------
REGIME_STYLE = {
    "baryon-dom": dict(color="#E64B35", label="Baryon-dominated"),
    "balanced":   dict(color="#4DBBD5", label="Balanced"),
    "DM-dom":     dict(color="#00A087", label="DM-dominated"),
}
POINT_ALPHA  = 0.45
POINT_SIZE   = 4

# ---------------------------------------------------------------------------
# Helper: compute g [m/s^2] from V [km/s] and r [kpc]
# ---------------------------------------------------------------------------
def accel(V_kms, r_kpc):
    """Return centripetal acceleration V^2/r in m/s^2."""
    V_ms = V_kms * KMS_TO_MS
    r_m  = r_kpc * KPC_TO_M
    return V_ms**2 / r_m


# ---------------------------------------------------------------------------
# Helper: load SPARC rotmod file
# ---------------------------------------------------------------------------
def load_rotmod(galaxy, sparc_dir):
    """
    Return arrays (r_kpc, Vobs, errV, Vgas, Vdisk, Vbul) for a galaxy.
    Lines starting with '#' are treated as comments.
    Signed component velocities are used as unsigned magnitudes.
    """
    path = os.path.join(sparc_dir, galaxy + "_rotmod.dat")
    rows = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            vals = line.split()
            if len(vals) < 6:
                continue
            try:
                r, vobs, errv, vgas, vdisk, vbul = [float(v) for v in vals[:6]]
                rows.append((r, vobs, errv, abs(vgas), abs(vdisk), abs(vbul)))
            except ValueError:
                continue
    if not rows:
        return None
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4], arr[:, 5]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(output_path):

    # -- Load regime classification ----------------------------------------
    regime_map = {}
    with open(FLEET_CSV, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            regime_map[row["galaxy"]] = row["regime"]

    # -- Collect all acceleration pairs ------------------------------------
    data_by_regime = {r: {"gbar": [], "gtot": []} for r in REGIME_STYLE}

    n_loaded  = 0
    n_points  = 0
    n_skipped = 0

    for galaxy, regime in regime_map.items():
        result = load_rotmod(galaxy, SPARC_DIR)
        if result is None:
            n_skipped += 1
            continue
        r_kpc, Vobs, errV, Vgas, Vdisk, Vbul = result

        # Quality cuts: positive and finite velocities, S/N > 2
        mask = (
            (r_kpc  > 0.0) &
            (Vobs   > 0.0) &
            (errV   > 0.0) &
            (Vobs / errV > 2.0) &
            np.isfinite(Vobs) &
            np.isfinite(Vgas) &
            np.isfinite(Vdisk) &
            np.isfinite(Vbul)
        )
        r_    = r_kpc[mask]
        Vobs_ = Vobs[mask]
        Vgas_ = Vgas[mask]
        Vd_   = Vdisk[mask]
        Vb_   = Vbul[mask]

        if r_.size == 0:
            n_skipped += 1
            continue

        # Baryonic velocity magnitude (quadrature sum)
        Vbar = np.sqrt(Vgas_**2 + Vd_**2 + Vb_**2)
        mask2 = Vbar > 0.0
        r_    = r_[mask2]
        Vobs_ = Vobs_[mask2]
        Vbar  = Vbar[mask2]

        if r_.size == 0:
            n_skipped += 1
            continue

        gbar = accel(Vbar,  r_)
        gtot = accel(Vobs_, r_)

        # Only keep physically sensible pairs (g > 0)
        mask3 = (gbar > 0) & (gtot > 0)
        gbar = gbar[mask3]
        gtot = gtot[mask3]

        if regime not in data_by_regime:
            regime = "balanced"  # fallback

        data_by_regime[regime]["gbar"].extend(gbar.tolist())
        data_by_regime[regime]["gtot"].extend(gtot.tolist())
        n_loaded += 1
        n_points += gbar.size

    print(f"Loaded: {n_loaded} galaxies, {n_points} data points "
          f"({n_skipped} skipped).", flush=True)

    # -- McGaugh+16 RAR theory curve ---------------------------------------
    g_min = 1e-13
    g_max = 5e-9
    gbar_rar = np.logspace(np.log10(g_min), np.log10(g_max), 500)
    gtot_rar = gbar_rar / (1.0 - np.exp(-np.sqrt(gbar_rar / G_DAGGER)))

    # -- Figure setup -------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.set_xscale("log")
    ax.set_yscale("log")

    # 1:1 line (g_tot = g_bar, Newtonian / baryon-only limit)
    g_11 = np.array([g_min, g_max])
    ax.plot(g_11, g_11, color="0.6", lw=0.8, ls=":", zorder=1,
            label=r"$g_\mathrm{tot} = g_\mathrm{bar}$ (Newtonian)")

    # RAR curve
    ax.plot(gbar_rar, gtot_rar, color="k", lw=1.6, zorder=2,
            label=r"McGaugh et al. (2016) RAR")

    # Data points, one regime at a time (plot in order: DM-dom last → on top)
    plot_order = ["baryon-dom", "balanced", "DM-dom"]
    handles_data = []
    for regime in plot_order:
        style = REGIME_STYLE[regime]
        gbar_arr = np.array(data_by_regime[regime]["gbar"])
        gtot_arr = np.array(data_by_regime[regime]["gtot"])
        if gbar_arr.size == 0:
            continue
        sc = ax.scatter(
            gbar_arr, gtot_arr,
            s=POINT_SIZE,
            color=style["color"],
            alpha=POINT_ALPHA,
            linewidths=0,
            zorder=3 + plot_order.index(regime),
            label=style["label"],
            rasterized=True,   # reduce PDF file size
        )
        handles_data.append(sc)

    # H2 null-result annotation
    ax.text(
        0.03, 0.97,
        r"$|\Delta\sigma| < 10^{-12}\,\mathrm{km\,s}^{-1}$" + "\nfor all 74 measured galaxies",
        transform=ax.transAxes,
        fontsize=7.5,
        va="top", ha="left",
        color="0.3",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.7),
    )

    # Axes labels
    ax.set_xlabel(
        r"$g_\mathrm{bar} \equiv (V_\mathrm{gas}^2+V_\mathrm{disk}^2+V_\mathrm{bul}^2)/r$"
        r"  $[\mathrm{m\,s}^{-2}]$",
        fontsize=10,
    )
    ax.set_ylabel(
        r"$g_\mathrm{tot} \equiv V_\mathrm{obs}^2/r$"
        r"  $[\mathrm{m\,s}^{-2}]$",
        fontsize=10,
    )
    ax.set_title(
        r"Acceleration-space coverage of the 80-galaxy H2 validation sample",
        fontsize=9.5,
        pad=8,
    )

    ax.set_xlim(g_min, g_max)
    ax.set_ylim(g_min, g_max)
    ax.tick_params(which="both", direction="in", top=True, right=True)

    # Legend
    legend = ax.legend(
        fontsize=7.5,
        loc="lower right",
        framealpha=0.9,
        edgecolor="0.7",
        handletextpad=0.4,
        labelspacing=0.3,
    )

    # g_dagger marker
    ax.axvline(G_DAGGER, color="0.55", lw=0.7, ls="--", zorder=1)
    ax.text(
        G_DAGGER * 1.12, g_min * 4,
        r"$g^\dagger$",
        fontsize=8, color="0.45", va="bottom",
    )

    fig.tight_layout()

    # -- Save ---------------------------------------------------------------
    fig.savefig(output_path, dpi=300, bbox_inches="tight", format="pdf")
    print(f"Saved: {output_path}", flush=True)
    plt.close(fig)


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate acceleration-space figure for H2 paper."
    )
    parser.add_argument(
        "--output",
        default=os.path.join(SCRIPT_DIR, "acceleration_space_80galaxy.pdf"),
        help="Output PDF path (default: same directory as this script).",
    )
    args = parser.parse_args()
    main(args.output)
