from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from core.galaxy_io import load_galaxy_rc, get_Rd_star_kpc
from core.chi import compute_chi_from_rc


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--galaxy", default="NGC3198", help="Galaxy name, e.g. NGC3198")
    ap.add_argument("--no-smooth", action="store_true", help="Disable smoothing before differentiation")
    ap.add_argument("--sigma-idx", type=float, default=1.0, help="Gaussian smoothing sigma in index-space (fixed global)")
    ap.add_argument("--radius", type=int, default=4, help="Gaussian kernel half-width in samples")
    args = ap.parse_args()

    gal = args.galaxy
    smooth = (not args.no_smooth)

    g = load_galaxy_rc(gal)
    Rd = get_Rd_star_kpc(gal)

    res = compute_chi_from_rc(
        g.R_kpc,
        g.V_baryon_kms,
        Rd_star_kpc=Rd,
        smooth=smooth,
        sigma_idx=args.sigma_idx,
        radius=args.radius,
        eps=1e-30,
    )

    # write outputs
    out_dir = repo_root() / "data" / "derived" / "phase3"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Assemble CSV
    df = pd.DataFrame({
        "galaxy": gal,
        "R_kpc": res["raw"].r_kpc,
        "g_bar_raw": res["raw"].g_bar,
        "dgdr_raw": res["raw"].dgdr,
        "chi_raw": res["raw"].chi,
    })

    if "smooth" in res:
        df["g_bar_smooth"] = res["smooth"].g_bar
        df["dgdr_smooth"] = res["smooth"].dgdr
        df["chi_smooth"] = res["smooth"].chi

    csv_path = out_dir / f"chi_{gal}.csv"
    df.to_csv(csv_path, index=False)

    # Plot χ(r)
    plt.figure()
    plt.plot(res["raw"].r_kpc, res["raw"].chi, label="chi (raw)")
    if "smooth" in res:
        plt.plot(res["smooth"].r_kpc, res["smooth"].chi, label=f"chi (smooth, sigma_idx={args.sigma_idx})")
    plt.axvline(0.3 * float(np.max(res["raw"].r_kpc)), linestyle="--", label="~r_frac=0.3 (visual)")
    plt.xlabel("r (kpc)")
    plt.ylabel("chi (dimensionless)")
    plt.title(f"H2 Phase-3: chi(r) using Rd_star | {gal} (Rd_star={Rd:.3g} kpc)")
    plt.legend()

    png_path = out_dir / f"chi_{gal}.png"
    plt.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close()

    print(f"[Phase3] Galaxy: {gal}")
    print(f"[Phase3] Rd_star = {Rd:.6g} kpc")
    print(f"[Phase3] Points = {len(res['raw'].r_kpc)} | Rmax = {float(np.max(res['raw'].r_kpc)):.6g} kpc")
    print(f"[Phase3] Wrote CSV: {csv_path}")
    print(f"[Phase3] Wrote PNG: {png_path}")


if __name__ == "__main__":
    main()
