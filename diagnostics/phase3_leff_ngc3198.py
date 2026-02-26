from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from core.galaxy_io import load_galaxy_rc, get_Rd_star_kpc, get_h1_params_for_galaxy
from core.chi import compute_chi_from_rc

# Optional smoothing (SciPy). If not available, fall back to a tiny manual Gaussian.
try:
    from scipy.ndimage import gaussian_filter1d  # type: ignore
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False


def _as_galaxy_rc(obj):
    """
    Accept either:
      - GalaxyRC
      - dict-like wrapper with a 'raw' entry
    Return a GalaxyRC-like object with attributes r_kpc, v_baryon_kms, v_total_kms (or similar).
    """
    # common wrapper pattern used earlier: {"raw": GalaxyRC, "meta": ...}
    if isinstance(obj, dict) and "raw" in obj:
        return obj["raw"]
    return obj


def _get_rc_arrays(rc):
    """
    Robustly extract arrays from the GalaxyRC object (attribute names can vary).
    We try a few common names.
    """
    # radii
    for name in ("r_kpc", "R_kpc", "r"):
        if hasattr(rc, name):
            R = np.asarray(getattr(rc, name), dtype=np.float64)
            break
    else:
        raise AttributeError("GalaxyRC missing radius attribute (expected r_kpc/R_kpc/r).")

    # baryon velocity
    for name in ("v_baryon_kms", "V_baryon", "V_baryon_kms", "vbar_kms", "Vb"):
        if hasattr(rc, name):
            Vb = np.asarray(getattr(rc, name), dtype=np.float64)
            break
    else:
        raise AttributeError("GalaxyRC missing baryon velocity attribute (expected v_baryon_kms / V_baryon...).")

    # total velocity
    for name in ("v_total_kms", "V_total", "V_total_kms", "vtot_kms", "Vt"):
        if hasattr(rc, name):
            Vt = np.asarray(getattr(rc, name), dtype=np.float64)
            break
    else:
        raise AttributeError("GalaxyRC missing total velocity attribute (expected v_total_kms / V_total...).")

    if len(R) != len(Vb) or len(R) != len(Vt):
        raise ValueError(f"Length mismatch: len(R)={len(R)}, len(Vb)={len(Vb)}, len(Vt)={len(Vt)}")

    return R, Vb, Vt


def sigmoid_mask(r_frac: np.ndarray, r0: float = 0.70, k: float = 80.0) -> np.ndarray:
    """
    m(r): ~0 inner, ~1 outer. Smooth logistic hand-off.
    """
    x = k * (r_frac - r0)
    x = np.clip(x, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(-x))


def _smooth_1d(x: np.ndarray, sigma_idx: float) -> np.ndarray:
    if sigma_idx <= 0:
        return x.copy()

    if _HAS_SCIPY:
        return gaussian_filter1d(x, sigma=sigma_idx, mode="nearest")

    # fallback: small discrete Gaussian kernel (sigma in index units)
    # kernel half-width ~ 3*sigma
    hw = int(max(1, np.ceil(3.0 * sigma_idx)))
    grid = np.arange(-hw, hw + 1, dtype=np.float64)
    ker = np.exp(-0.5 * (grid / sigma_idx) ** 2)
    ker /= np.sum(ker)
    # pad edges
    xp = np.pad(x, (hw, hw), mode="edge")
    y = np.convolve(xp, ker, mode="valid")
    return y.astype(np.float64)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--galaxy", required=True, help="e.g. NGC3198")
    ap.add_argument("--alpha", type=float, default=1.0, help="L_eff = L0/(1+alpha*chi) (default=1.0)")
    ap.add_argument("--sigma_idx", type=float, default=1.0, help="chi smoothing in index units (default=1.0)")
    ap.add_argument("--taper", action="store_true", help="apply sigmoid taper to enforce L_eff -> L0 at large r")
    ap.add_argument("--taper_r0", type=float, default=0.70, help="sigmoid inflection in r_frac (default=0.70)")
    ap.add_argument("--taper_k", type=float, default=80.0, help="sigmoid steepness (default=80)")
    args = ap.parse_args()

    gal = args.galaxy.strip()

    # --- Load RC (GalaxyRC or wrapper) ---
    rc_any = load_galaxy_rc(gal)
    rc = _as_galaxy_rc(rc_any)
    R_kpc, Vb_kms, Vt_kms = _get_rc_arrays(rc)

    # --- H1 frozen params for this galaxy ---
    h1 = get_h1_params_for_galaxy(gal)  # expected keys: L, mu, mafe, kernel
    if "L" not in h1:
        raise KeyError(f"get_h1_params_for_galaxy('{gal}') did not return key 'L'. Got keys: {list(h1.keys())}")

    L0 = float(h1["L"])
    mu = float(h1.get("mu", np.nan))
    kernel = str(h1.get("kernel", "unknown"))
    Rd = float(get_Rd_star_kpc(gal))

    # --- Compute chi (core.chi returns dict[str, Chi1DResult]) ---
    # Signature: compute_chi_from_rc(r_kpc, V_baryon_kms, *, Rd_star_kpc=..., smooth=..., sigma_idx=...)
    chi_out = compute_chi_from_rc(
        R_kpc,
        Vb_kms,
        Rd_star_kpc=Rd,   # keyword-only
        smooth=False      # keep local smoothing step below (Phase-3 plotting)
    )
    chi_raw = np.asarray(chi_out["raw"].chi, dtype=np.float64)

    # smooth chi (NOT gbar) — cheap, stable, and sufficient for the Phase-3 diagnostic plots
    chi_smooth = _smooth_1d(chi_raw, sigma_idx=float(args.sigma_idx))
    chi_used = chi_raw  # Use stronger raw signal for deformation
    chi_amp = chi_used / np.max(chi_used)  # Normalize for deformation
    # --- L_adapt and L_eff ---
    L_adapt = L0 / (1.0 + float(args.alpha) * chi_amp)  # Use normalized chi

    # taper (optional)
    R_max = float(np.max(R_kpc))
    r_frac = (R_kpc / R_max) if R_max > 0 else np.zeros_like(R_kpc)

    if args.taper:
        m = sigmoid_mask(r_frac, r0=float(args.taper_r0), k=float(args.taper_k))
        L_eff = (1.0 - m) * L_adapt + m * L0
    else:
        m = np.zeros_like(r_frac)
        L_eff = L_adapt

    # --- Write outputs ---
    H2_ROOT = Path(__file__).resolve().parents[1]
    out_dir = H2_ROOT / "data" / "derived" / "phase3"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_csv = out_dir / f"leff_{gal}.csv"
    df = pd.DataFrame({
        "R_kpc": R_kpc,
        "r_frac": r_frac,
        "chi_raw": chi_raw,
        "chi_smooth": chi_smooth,
        "chi_used": chi_used,
        "mask_outer": m,              # 0 if taper off
        "L_adapt_kpc": L_adapt,
        "L_eff_kpc": L_eff,           # <-- Phase-4 needs this exact name
    })

    # scalar metadata (broadcast to all rows)
    df["galaxy"] = gal
    df["L0_kpc"] = L0
    df["alpha"] = float(args.alpha)
    df["sigma_idx"] = float(args.sigma_idx)
    df["taper_on"] = int(bool(args.taper))
    df["taper_r0"] = float(args.taper_r0)
    df["taper_k"] = float(args.taper_k)
    df["mu"] = mu
    df["kernel"] = kernel
    df["Rd_star_kpc"] = Rd
    df.to_csv(out_csv, index=False)

    print(f"[Phase3-Leff] Galaxy: {gal}")
    print(f"[Phase3-Leff] Rd_star = {Rd:.3g} kpc | L0 = {L0:g} kpc | mu = {mu:g} | kernel = {kernel}")
    print(f"[Phase3-Leff] chi used: chi_smooth (sigma_idx={args.sigma_idx:g})")
    print(f"[Phase3-Leff] taper: {'ON' if args.taper else 'OFF'} (r0={args.taper_r0:g}, k={args.taper_k:g})")
    print(f"[Phase3-Leff] Wrote CSV: {out_csv}")

    # --- Plots ---
    # Plot A: chi
    plt.figure(figsize=(8, 6))
    plt.plot(R_kpc, chi_raw, label="chi_raw")
    plt.plot(R_kpc, chi_smooth, label=f"chi_smooth (sigma_idx={args.sigma_idx:g})")
    plt.xlabel("r (kpc)")
    plt.ylabel("chi (dimensionless)")
    plt.title(f"H2 Phase-3: chi(r) | {gal} (Rd_star={Rd:.3g} kpc)")
    plt.legend()
    plt.tight_layout()
    chi_png = out_dir / f"chi_{gal}.png"
    plt.savefig(chi_png, dpi=200)
    plt.close()

    # Plot B: L_eff
    plt.figure(figsize=(8, 6))
    plt.plot(R_kpc, L_adapt, label="L_adapt = L0/(1+alpha*chi)")
    if args.taper:
        plt.plot(R_kpc, L_eff, label="L_eff (tapered)")
        plt.plot(R_kpc, m * L0, alpha=0.25, label="mask_outer * L0 (visual)")
    else:
        plt.plot(R_kpc, L_eff, label="L_eff (no taper)")
    plt.axhline(L0, linestyle="--", label=f"L0={L0:g} kpc")
    plt.xlabel("r (kpc)")
    plt.ylabel("L_eff (kpc)")
    plt.title(f"H2 Phase-3: L_eff(r) | {gal}")
    plt.legend()
    plt.tight_layout()
    leff_png = out_dir / f"leff_{gal}.png"
    plt.savefig(leff_png, dpi=200)
    plt.close()

    print(f"[Phase3-Leff] Wrote PNG: {chi_png}")
    print(f"[Phase3-Leff] Wrote PNG: {leff_png}")


if __name__ == "__main__":
    main()
