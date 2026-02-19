from __future__ import annotations

import argparse
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Paths (repo-local, deterministic)
# ============================================================
H2_ROOT = Path(__file__).resolve().parents[1]

PHASE3_DIR = H2_ROOT / "data" / "derived" / "phase3"
PHASE4_DIR = H2_ROOT / "data" / "derived" / "phase4"

BASIS_DIR = PHASE4_DIR / "basis_h1"
OUT_DIR = PHASE4_DIR / "h2_outputs"

# Vendored H1 snapshot
DEFAULT_H1_ROOT = H2_ROOT / "vendor" / "h1_src"
DEFAULT_H1_RUNNER = DEFAULT_H1_ROOT / "run_sparc_lite.py"
DEFAULT_H1_RESULTS = DEFAULT_H1_ROOT / "results"


# ============================================================
# Config
# ============================================================
@dataclass(frozen=True)
class BasisConfig:
    # Fixed basis multipliers (NO tuning, deterministic)
    multipliers: tuple[float, ...] = (1.00, 0.85, 0.70, 0.55)

    # Outer gate
    rfrac_outer: float = 0.70
    tol_outer_kms: float = 2.0


# ============================================================
# Utilities
# ============================================================
def _require_exists(p: Path, label: str) -> None:
    if not p.exists():
        raise FileNotFoundError(f"{label} not found: {p}")


def _ensure_float64_series(df: pd.DataFrame, cols: tuple[str, ...], name: str) -> None:
    for c in cols:
        if c not in df.columns:
            raise ValueError(f"{name} missing column '{c}'. Found: {list(df.columns)}")
        df[c] = pd.to_numeric(df[c], errors="raise").astype(np.float64)


def _load_rc_decomp_csv(p: Path) -> pd.DataFrame:
    df = pd.read_csv(p)
    _ensure_float64_series(df, ("R_kpc", "V_baryon", "V_total"), p.name)
    df = df.sort_values("R_kpc", kind="mergesort").reset_index(drop=True)
    return df


def _load_frozen_reference(galaxy: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Canonical reference for Phase-4:
      - R_kpc grid
      - V_baryon (invariant)
      - V_total_H1 (baseline truth for Test-1)
    """
    p = H2_ROOT / "data" / "h1_frozen" / "per_galaxy" / f"rc_decomp_{galaxy}_best.csv"
    _require_exists(p, "Frozen per-galaxy H1 baseline rc_decomp")
    df = _load_rc_decomp_csv(p)

    R = df["R_kpc"].to_numpy(dtype=np.float64)
    Vb = df["V_baryon"].to_numpy(dtype=np.float64)
    Vt = df["V_total"].to_numpy(dtype=np.float64)
    return R, Vb, Vt


def _load_phase3_leff(galaxy: str) -> pd.DataFrame:
    """
    Phase-3 contract: requires leff_<GAL>.csv with column L_eff_kpc.
    """
    leff_path = PHASE3_DIR / f"leff_{galaxy}.csv"
    _require_exists(leff_path, "Phase-3 Leff CSV")

    df = pd.read_csv(leff_path)

    # Required minimal columns
    required = ("R_kpc", "L_eff_kpc", "L0_kpc", "mu", "kernel")
    for c in required:
        if c not in df.columns:
            raise ValueError(
                f"{leff_path.name} missing required column '{c}'. Found: {list(df.columns)}"
            )

    _ensure_float64_series(df, ("R_kpc", "L_eff_kpc", "L0_kpc", "mu"), leff_path.name)
    df = df.sort_values("R_kpc", kind="mergesort").reset_index(drop=True)
    return df


def _pick_basis_Ls(L0_kpc: float, cfg: BasisConfig) -> list[float]:
    # Descending list (nice for display); we will sort ascending for interpolation.
    basis = sorted({float(L0_kpc * m) for m in cfg.multipliers}, reverse=True)
    return basis


def _basis_filename(galaxy: str, L_kpc: float) -> str:
    # Stable filename: treat L as integer kpc if near integer, else keep 3 decimals.
    if abs(L_kpc - round(L_kpc)) < 1e-9:
        tag = f"{int(round(L_kpc))}kpc"
    else:
        tag = f"{L_kpc:.3f}kpc"
    return f"rc_decomp_{galaxy}_L{tag}.csv"


def _run_h1_single_galaxy(
    *,
    h1_root: Path,
    galaxy: str,
    L_kpc: float,
    mu: float,
    kernel: str,
) -> Path:
    """
    Deterministic H1 basis run:
      - no survey loop
      - no kernel grid search
      - no downloads
      - directly calls predict_rc_for_params(gal, L, mu, kernel)
    Produces: <h1_root>/results/rc_decomp_<GAL>_best.csv
    """
    _require_exists(h1_root, "H1 root")
    runner = h1_root / "run_sparc_lite.py"
    _require_exists(runner, "H1 runner")

    results_dir = h1_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    expected = results_dir / f"rc_decomp_{galaxy}_best.csv"
    if expected.exists():
        expected.unlink()

    snippet = f"""
import sys
from pathlib import Path
import importlib.util
import numpy as np
import pandas as pd

# --- paths ---
h1_root = Path(r"{str(h1_root)}").resolve()
h2_data = Path(r"{str(H2_ROOT / 'data')}").resolve()
table_path = (h2_data / "galaxies.csv").resolve()

# --- import vendored H1 runner by path ---
sys.path.insert(0, str(h1_root))
sys.path.insert(0, str(h1_root / "src"))

p = (h1_root / "run_sparc_lite.py").resolve()
spec = importlib.util.spec_from_file_location("run_sparc_lite", str(p))
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# --- HARD SAFETY: no downloads ---
if hasattr(mod, "download_and_extract_data"):
    mod.download_and_extract_data = lambda: None

# --- CRITICAL: disable single-galaxy debug filter in vendor runner ---
if hasattr(mod, "TARGET_GALAXY"):
    mod.TARGET_GALAXY = None

# (Optional) try to steer module-level data dir variables if present
if hasattr(mod, "DATA_DIR"):
    mod.DATA_DIR = str(h2_data)
if hasattr(mod, "data_dir"):
    mod.data_dir = str(h2_data)

def norm(s):
    if s is None:
        return ""
    s = str(s).strip().upper()
    for ch in [" ", "_", "-", "\\t"]:
        s = s.replace(ch, "")
    return s

# --- Load galaxies table DIRECTLY (no vendor stubs, no renaming) ---
df = pd.read_csv(table_path)
if "name" not in df.columns:
    raise RuntimeError(
        "H2 galaxies.csv must contain a 'name' column. Columns: " + str(list(df.columns))
    )

gals = df.to_dict(orient="records")

target = norm("{galaxy}")

gal = None
for g in gals:
    if norm(g.get("name")) == target:
        gal = g
        break

if gal is None:
    raise RuntimeError(
        "Galaxy not found in H2 galaxies.csv: {galaxy}. "
        "First 20 names: " + str([g.get("name") for g in gals[:20]])
    )

# --- deterministic prediction (single param set) ---
res = mod.predict_rc_for_params(gal, {float(L_kpc)}, {float(mu)}, "{kernel}", beta=1.15)
if res is None:
    raise RuntimeError("predict_rc_for_params returned None")

R_pred, V_pred, V_b, V_k = res

# --- write standard format into vendor results folder ---
out = (h1_root / "results" / "rc_decomp_{galaxy}_best.csv")
out.parent.mkdir(parents=True, exist_ok=True)

np.savetxt(
    str(out),
    np.c_[R_pred, V_b, V_k, V_pred],
    delimiter=",",
    header="R_kpc,V_baryon,V_kernel,V_total",
    comments=""
)

print("WROTE", str(out))
"""

    import os

    env = os.environ.copy()
    env["DEBUG_DX"] = "1.0"

    proc = subprocess.run(
        ["python", "-c", snippet],
        cwd=str(h1_root),
        capture_output=True,
        text=True,
        env=env,
    )


    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr)
        raise RuntimeError(f"H1 basis run failed for {galaxy} at L={L_kpc} kpc (mu={mu}, kernel={kernel})")

    if not expected.exists():
        print(proc.stdout)
        print(proc.stderr)
        raise FileNotFoundError(f"H1 did not produce expected file: {expected}")

    return expected



def _gate_outer_deltaV(
    R_kpc: np.ndarray,
    V_h1: np.ndarray,
    V_h2: np.ndarray,
    *,
    rfrac_outer: float,
    tol_kms: float,
) -> tuple[float, bool]:
    rmax = float(np.max(R_kpc))
    rfrac = R_kpc / rmax if rmax > 0 else np.zeros_like(R_kpc)
    mask = rfrac >= rfrac_outer
    if not np.any(mask):
        return float("nan"), False
    dV = float(np.max(np.abs(V_h2[mask] - V_h1[mask])))
    return dV, (dV <= tol_kms)


def _save_outputs_and_plots(
    *,
    galaxy: str,
    R_kpc: np.ndarray,
    V_baryon: np.ndarray,
    V_h1: np.ndarray,
    V_h2: np.ndarray,
    L_eff: np.ndarray,
    meta: dict,
) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # CSV output (H2 adaptive decomposition)
    out_csv = OUT_DIR / f"rc_decomp_{galaxy}_H2_adaptive.csv"
    df_out = pd.DataFrame(
        {
            "R_kpc": R_kpc,
            "V_baryon": V_baryon,
            "V_total_H1": V_h1,
            "V_total_H2": V_h2,
            "dV_H2_minus_H1": (V_h2 - V_h1),
            "L_eff_kpc": L_eff,
        }
    )
    for k, v in meta.items():
        df_out[k] = v
    df_out.to_csv(out_csv, index=False)

    # Plot 1: rotation curves
    plt.figure(figsize=(8, 6))
    plt.plot(R_kpc, V_h1, label="H1 frozen V_total")
    plt.plot(R_kpc, V_h2, label="H2 adaptive V_total")
    plt.plot(R_kpc, V_baryon, label="V_baryon (canonical)")
    plt.xlabel("R (kpc)")
    plt.ylabel("V (km/s)")
    plt.title(f"H2 Phase-4: Adaptive convolution — {galaxy}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_DIR / f"phase4_rc_{galaxy}.png", dpi=180)
    plt.close()

    # Plot 2: deltaV
    plt.figure(figsize=(8, 5))
    plt.plot(R_kpc, V_h2 - V_h1)
    plt.axhline(0.0)
    plt.xlabel("R (kpc)")
    plt.ylabel("ΔV (km/s) [H2 − H1]")
    plt.title(f"H2 Phase-4: ΔV(r) — {galaxy}")
    plt.tight_layout()
    plt.savefig(OUT_DIR / f"phase4_deltaV_{galaxy}.png", dpi=180)
    plt.close()

    # Plot 3: L_eff
    plt.figure(figsize=(8, 5))
    plt.plot(R_kpc, L_eff)
    plt.axhline(float(meta["L0_kpc"]), linestyle="--", label="L0")
    plt.xlabel("R (kpc)")
    plt.ylabel("L_eff (kpc)")
    plt.title(f"H2 Phase-4: L_eff(r) — {galaxy}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_DIR / f"phase4_leff_{galaxy}.png", dpi=180)
    plt.close()

    print(f"[Phase4] Wrote CSV: {out_csv}")
    print(f"[Phase4] Wrote PNGs into: {OUT_DIR}")


# ============================================================
# Main
# ============================================================
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--galaxy", required=True, help="e.g. NGC3198")
    ap.add_argument(
        "--h1_root",
        default=str(DEFAULT_H1_ROOT),
        help="path to vendored H1 snapshot (default: H2/vendor/h1_src)",
    )
    ap.add_argument(
        "--no_run_h1",
        action="store_true",
        help="do not execute H1; use cached basis CSVs if present",
    )
    ap.add_argument(
        "--strict_radii",
        action="store_true",
        help="hard-fail if basis R_kpc mismatches canonical (recommended ON for reviewer-proof runs)",
    )
    args = ap.parse_args()

    galaxy = args.galaxy.strip()
    cfg = BasisConfig()

    # ----------------------------
    # Canonical baseline truth
    # ----------------------------
    R_ref, Vb_ref, Vt_h1_ref = _load_frozen_reference(galaxy)

    # ----------------------------
    # Phase-3 adaptive length
    # ----------------------------
    leff_df = _load_phase3_leff(galaxy)

    # Enforce Phase-3 radii == canonical radii
    R_leff = leff_df["R_kpc"].to_numpy(dtype=np.float64)
    if (len(R_leff) != len(R_ref)) or (not np.allclose(R_leff, R_ref, atol=1e-10, rtol=0.0)):
        raise RuntimeError(
            "Phase-3 radii mismatch frozen H1 radii.\n"
            "Fix: regenerate Phase-3 leff_<GAL>.csv using the frozen per-galaxy rc_decomp_<GAL>_best.csv grid.\n"
            f"Phase-3 n={len(R_leff)} vs frozen n={len(R_ref)}"
        )

    L0 = float(leff_df["L0_kpc"].iloc[0])
    mu = float(leff_df["mu"].iloc[0])
    kernel = str(leff_df["kernel"].iloc[0])
    L_eff = leff_df["L_eff_kpc"].to_numpy(dtype=np.float64)

    # ----------------------------
    # Basis grid
    # ----------------------------
    basis_L_desc = _pick_basis_Ls(L0, cfg)
    basis_L_asc = np.array(sorted(basis_L_desc), dtype=np.float64)  # ascending for np.interp

    print(f"[Phase4] Galaxy={galaxy}")
    print(f"[Phase4] Canonical: L0={L0:g} kpc | mu={mu:g} | kernel={kernel}")
    print(f"[Phase4] Basis L grid (kpc): {list(basis_L_desc)}")
    print(f"[Phase4] Leff min/max (kpc): {float(np.min(L_eff)):.6g} / {float(np.max(L_eff)):.6g}")

    # ----------------------------
    # Ensure dirs
    # ----------------------------
    BASIS_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ----------------------------
    # Build or reuse basis runs
    # ----------------------------
    h1_root = Path(args.h1_root).resolve()
    _require_exists(h1_root, "H1 root")
    _require_exists(h1_root / "run_sparc_lite.py", "H1 runner")

    basis_files: dict[float, Path] = {}

    for L in basis_L_asc:
        out = BASIS_DIR / _basis_filename(galaxy, float(L))
        basis_files[float(L)] = out

        if args.no_run_h1 and out.exists():
            continue

        produced = _run_h1_single_galaxy(
            h1_root=h1_root,
            galaxy=galaxy,
            L_kpc=float(L),
            mu=mu,
            kernel=kernel,
        )
        shutil.copy2(produced, out)

    # ----------------------------
    # Load basis V_total(L, r) with strict grid enforcement
    # ----------------------------
    basis_Vt = np.zeros((len(basis_L_asc), len(R_ref)), dtype=np.float64)

    for i, L in enumerate(basis_L_asc):
        p = basis_files[float(L)]
        _require_exists(p, f"Basis file for L={L:g} kpc")

        df = _load_rc_decomp_csv(p)

        # NOTE: we IGNORE df["V_baryon"] entirely (canonical Vb comes from frozen baseline)
        R_now = df["R_kpc"].to_numpy(dtype=np.float64)
        Vt_now = df["V_total"].to_numpy(dtype=np.float64)

        # Strict: basis radii MUST match canonical radii
        if (len(R_now) != len(R_ref)) or (not np.allclose(R_now, R_ref, atol=1e-10, rtol=0.0)):
            msg = (
                f"Basis radii mismatch at L={L:g} kpc.\n"
                f"basis n={len(R_now)} vs canonical n={len(R_ref)}\n"
                "This invalidates interpolation across L.\n"
                "Fix H1 runner to output identical R_kpc for all L and match frozen per_galaxy grid."
            )
            if args.strict_radii:
                raise RuntimeError(msg)
            else:
                # Non-strict fallback (not reviewer-proof): resample V_total onto canonical grid
                # Kept only for debugging, NOT recommended.
                Vt_now = np.interp(R_ref, R_now, Vt_now, left=Vt_now[0], right=Vt_now[-1])

        basis_Vt[i, :] = Vt_now

    # ----------------------------
    # Adaptive interpolation across L at each radius
    # ----------------------------
    L_min = float(basis_L_asc[0])
    L_max = float(basis_L_asc[-1])
    L_eff_clip = np.clip(L_eff, L_min, L_max)

    Vt_h2 = np.zeros_like(R_ref, dtype=np.float64)
    for j in range(len(R_ref)):
        Vt_h2[j] = np.interp(L_eff_clip[j], basis_L_asc, basis_Vt[:, j])

    # ----------------------------
    # Gate: outer stability vs frozen baseline
    # ----------------------------
    max_outer, pass_outer = _gate_outer_deltaV(
        R_ref,
        Vt_h1_ref,
        Vt_h2,
        rfrac_outer=cfg.rfrac_outer,
        tol_kms=cfg.tol_outer_kms,
    )

    print(f"[Phase4] TEST-1 Outer Stability (r_frac >= {cfg.rfrac_outer:.2f}):")
    print(f"[Phase4]   max |ΔV| outer = {max_outer:.6g} km/s  (tol={cfg.tol_outer_kms:.3g})")
    print(f"[Phase4]   PASS={pass_outer}")

    # ----------------------------
    # Save outputs + plots
    # ----------------------------
    meta = {
        "galaxy": galaxy,
        "L0_kpc": L0,
        "mu": mu,
        "kernel": kernel,
        "basis_L_min_kpc": L_min,
        "basis_L_max_kpc": L_max,
        "outer_rfrac": cfg.rfrac_outer,
        "outer_tol_kms": cfg.tol_outer_kms,
        "outer_max_abs_dV_kms": max_outer,
        "outer_pass": int(bool(pass_outer)),
        "strict_radii": int(bool(args.strict_radii)),
    }

    _save_outputs_and_plots(
        galaxy=galaxy,
        R_kpc=R_ref,
        V_baryon=Vb_ref,
        V_h1=Vt_h1_ref,
        V_h2=Vt_h2,
        L_eff=L_eff,
        meta=meta,
    )


if __name__ == "__main__":
    main()
