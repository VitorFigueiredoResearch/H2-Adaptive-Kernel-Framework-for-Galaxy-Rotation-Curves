# diagnostics/test3_inner_scatter.py
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def load_vobs_from_rotmod(rotmod_path: Path):
    """
    SPARC rotmod.dat: col0=R[kpc], col1=Vobs[km/s]
    (ignores comment lines starting with #)
    """
    Robs, Vobs = np.loadtxt(rotmod_path, comments="#", usecols=(0, 1), unpack=True)
    return Robs.astype(float), Vobs.astype(float)


def rms_scatter_log10(V_model, V_obs, eps=1e-12):
    # guard against zeros/negatives
    V_model = np.clip(V_model, eps, None)
    V_obs = np.clip(V_obs, eps, None)
    resid = np.log10(V_model) - np.log10(V_obs)
    return float(np.sqrt(np.mean(resid**2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--galaxy", default="NGC3198")
    ap.add_argument("--inner_frac", type=float, default=0.5,
                    help="Inner region defined as R < inner_frac * Rmax (Rmax from observed grid).")
    args = ap.parse_args()

    g = args.galaxy

    # Paths (match your repo)
    rotmod = Path("data/sparc") / f"{g}_rotmod.dat"
    h1_csv = Path("data/h1_frozen/per_galaxy") / f"rc_decomp_{g}_best.csv"
    h2_csv = Path("data/derived/phase4/h2_outputs") / f"rc_decomp_{g}_H2_adaptive.csv"

    if not rotmod.exists():
        raise FileNotFoundError(f"Missing rotmod file: {rotmod}")
    if not h1_csv.exists():
        raise FileNotFoundError(f"Missing H1 file: {h1_csv}")
    if not h2_csv.exists():
        raise FileNotFoundError(f"Missing H2 file: {h2_csv}")

    # Load observed
    R_obs, V_obs = load_vobs_from_rotmod(rotmod)

    # Load models
    h1 = pd.read_csv(h1_csv)
    h2 = pd.read_csv(h2_csv)

    # Column guards (your files use these)
    if "R_kpc" not in h1.columns:
        raise RuntimeError(f"{h1_csv} missing column R_kpc")
    if "R_kpc" not in h2.columns:
        raise RuntimeError(f"{h2_csv} missing column R_kpc")

    # H1 total column name may vary; try common options
    def get_vtotal(df, label):
        for c in ["V_total", "V_total_H1", "Vtot", "V_model"]:
            if c in df.columns:
                return df[c].to_numpy(float)
        raise RuntimeError(f"Could not find a total velocity column in {label}. Columns: {list(df.columns)}")

    V_h1_grid = get_vtotal(h1, "H1")
    V_h2_grid = get_vtotal(h2, "H2")

    R_h1 = h1["R_kpc"].to_numpy(float)
    R_h2 = h2["R_kpc"].to_numpy(float)

    # Interpolate onto observed grid (so comparisons are apples-to-apples)
    V_h1 = np.interp(R_obs, R_h1, V_h1_grid)
    V_h2 = np.interp(R_obs, R_h2, V_h2_grid)

    # Inner mask
    R_max = float(np.max(R_obs))
    R_cut = args.inner_frac * R_max
    mask = R_obs < R_cut

    N = int(np.sum(mask))
    if N < 3:
        raise RuntimeError(f"Too few inner points (N={N}). Try increasing --inner_frac.")

    sigma_h1 = rms_scatter_log10(V_h1[mask], V_obs[mask])
    sigma_h2 = rms_scatter_log10(V_h2[mask], V_obs[mask])
    delta = sigma_h2 - sigma_h1

    print("=" * 60)
    print("TEST-3: INNER SCATTER COMPARISON")
    print("=" * 60)
    print(f"Galaxy: {g}")
    print(f"rotmod: {rotmod}")
    print(f"H1:     {h1_csv}")
    print(f"H2:     {h2_csv}")
    print(f"Inner region: R < {R_cut:.2f} kpc (inner_frac={args.inner_frac}) | N={N}")
    print()
    print(f"sigma_inner(H1) = {sigma_h1:.4f} dex")
    print(f"sigma_inner(H2) = {sigma_h2:.4f} dex")
    print(f"Delta sigma (H2 - H1) = {delta:.4f} dex")
    print()
    if delta < -0.01:
        print("Result: ✅ H2 reduces inner scatter (meaningful improvement).")
    elif delta > 0.01:
        print("Result: ⚠️ H2 increases inner scatter (mechanism works, direction may be wrong).")
    else:
        print("Result: ➖ H2 is neutral within ±0.01 dex (mechanism works, leverage modest).")
    print("=" * 60)


if __name__ == "__main__":
    main()