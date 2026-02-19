from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from core.galaxy_io import load_galaxy_rc, get_Rd_star_kpc, get_h1_params_for_galaxy

# If you have core.chi.py already (you do), this import enables "compute if missing"
try:
    from core.chi import compute_chi_from_rc
except Exception:
    compute_chi_from_rc = None  # type: ignore


REPO_ROOT = Path(__file__).resolve().parents[1]
PHASE3_DIR = REPO_ROOT / "data" / "derived" / "phase3"


def leff_linear(L0_kpc: float, chi: np.ndarray, lmin_frac: float = 0.05) -> np.ndarray:
    """
    Diary-canonical first pass:
        L_eff = L0 / (1 + chi)

    Optional global clamp (frozen, not tuned) to avoid pathological collapse.
    """
    chi = np.asarray(chi, dtype=np.float64)
    Le = float(L0_kpc) / (1.0 + chi)
    Le = np.maximum(Le, lmin_frac * float(L0_kpc))
    return Le.astype(np.float64)


def _default_chi_csv_path(galaxy: str) -> Path:
    return PHASE3_DIR / f"chi_{galaxy}.csv"


def _default_out_csv_path(galaxy: str) -> Path:
    return PHASE3_DIR / f"leff_{galaxy}.csv"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--galaxy", default="NGC3198")
    ap.add_argument("--sigma-idx", type=float, default=1.0, help="Gaussian smoothing sigma in grid-index units (used only if computing chi)")
    ap.add_argument("--radius", type=int, default=4, help="Gradient stencil radius (used only if computing chi)")
    ap.add_argument("--lmin-frac", type=float, default=0.05, help="Clamp: minimum L_eff as fraction of L0")
    ap.add_argument("--prefer-existing", action="store_true", help="Use existing chi_<GAL>.csv if present (default behavior)")
    ap.add_argument("--force-recompute", action="store_true", help="Recompute chi even if chi CSV exists")
    args = ap.parse_args()

    gal = args.galaxy

    PHASE3_DIR.mkdir(parents=True, exist_ok=True)

    # Load metadata (R, etc.)
    rc = load_galaxy_rc(gal)
    Rd = get_Rd_star_kpc(gal)
    p = get_h1_params_for_galaxy(gal)
    L0 = float(p["L"])

    chi_csv = _default_chi_csv_path(gal)

    use_existing = (chi_csv.exists() and not args.force_recompute)
    if use_existing:
        df_chi = pd.read_csv(chi_csv)
        if "R_kpc" not in df_chi.columns or ("chi_raw" not in df_chi.columns and "chi_smooth" not in df_chi.columns and "chi" not in df_chi.columns):
            raise ValueError(
                f"{chi_csv} exists but does not have expected columns. "
                f"Need R_kpc and chi_raw/chi_smooth (or chi). Found: {list(df_chi.columns)}"
            )
        print(f"[Phase3-Leff] Loaded existing chi CSV: {chi_csv}")
    else:
        if compute_chi_from_rc is None:
            raise RuntimeError(
                "compute_chi_from_rc not available (core.chi import failed). "
                "Either fix core/chi.py importability, or generate chi_<GAL>.csv first."
            )

        # Compute chi (raw + smooth)
        res = compute_chi_from_rc(
            rc.R_kpc,
            rc.V_baryon_kms,
            Rd_star_kpc=Rd,
            smooth=True,
            sigma_idx=args.sigma_idx,
            radius=args.radius,
            eps=1e-30,
        )

        # We store both raw and smooth (preferred for diary 2.7)
        df_chi = pd.DataFrame({
            "galaxy": gal,
            "R_kpc": res["raw"].r_kpc,
            "chi_raw": res["raw"].chi,
            "chi_smooth": res["smooth"].chi,
        })

        df_chi.to_csv(chi_csv, index=False)
        print(f"[Phase3-Leff] Computed chi and wrote: {chi_csv}")

    # Decide which chi column to use
    if "chi_smooth" in df_chi.columns:
        chi_used = pd.to_numeric(df_chi["chi_smooth"], errors="coerce").to_numpy(dtype=np.float64)
        chi_used_label = f"chi_smooth (sigma_idx={args.sigma_idx})"
    elif "chi" in df_chi.columns:
        chi_used = pd.to_numeric(df_chi["chi"], errors="coerce").to_numpy(dtype=np.float64)
        chi_used_label = "chi"
    else:
        chi_used = pd.to_numeric(df_chi["chi_raw"], errors="coerce").to_numpy(dtype=np.float64)
        chi_used_label = "chi_raw"

    R = pd.to_numeric(df_chi["R_kpc"], errors="coerce").to_numpy(dtype=np.float64)

    # L_eff(r)
    Le = leff_linear(L0, chi_used, lmin_frac=args.lmin_frac)

    out_df = pd.DataFrame({
        "galaxy": gal,
        "R_kpc": R,
        "chi_used": chi_used,
        "chi_used_label": chi_used_label,
        "Rd_star_kpc": float(Rd),
        "L0_kpc": float(L0),
        "Leff_kpc": Le,
        "lmin_frac": float(args.lmin_frac),
        "mu": float(p.get("mu", np.nan)),
        "kernel": str(p.get("kernel", "")),
    })

    out_csv = _default_out_csv_path(gal)
    out_df.to_csv(out_csv, index=False)

    # --- Plot 1: chi(r)
    plt.figure()
    plt.plot(R, chi_used, label=chi_used_label)
    plt.xlabel("r (kpc)")
    plt.ylabel("chi (dimensionless)")
    plt.title(f"H2 Phase-3: chi(r) | {gal} (Rd_star={Rd:.3g} kpc)")
    plt.legend()
    chi_png = PHASE3_DIR / f"chi_used_{gal}.png"
    plt.savefig(chi_png, dpi=200, bbox_inches="tight")
    plt.close()

    # --- Plot 2: Leff(r)
    plt.figure()
    plt.plot(R, Le, label="L_eff(r)")
    plt.axhline(L0, linestyle="--", label=f"L0={L0:g} kpc")
    plt.xlabel("r (kpc)")
    plt.ylabel("L_eff (kpc)")
    plt.title(f"H2 Phase-3: L_eff(r)=L0/(1+chi) | {gal}")
    plt.legend()
    leff_png = PHASE3_DIR / f"leff_{gal}.png"
    plt.savefig(leff_png, dpi=200, bbox_inches="tight")
    plt.close()

    print(f"[Phase3-Leff] Galaxy: {gal}")
    print(f"[Phase3-Leff] Rd_star = {Rd:.6g} kpc | L0 = {L0:.6g} kpc | mu = {p.get('mu')} | kernel = {p.get('kernel')}")
    print(f"[Phase3-Leff] chi used: {chi_used_label}")
    print(f"[Phase3-Leff] Wrote CSV: {out_csv}")
    print(f"[Phase3-Leff] Wrote PNG: {chi_png}")
    print(f"[Phase3-Leff] Wrote PNG: {leff_png}")


if __name__ == "__main__":
    main()
