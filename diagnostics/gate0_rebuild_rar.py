from __future__ import annotations

import re
from pathlib import Path
import numpy as np
import pandas as pd

# -----------------------
# CONFIG: edit if needed
# -----------------------
H2_ROOT = Path(__file__).resolve().parents[1]
H1_PER_GALAXY = H2_ROOT / "data" / "h1_frozen" / "per_galaxy"
H15_RAR_POINTS = H2_ROOT / "data" / "h1p5_diagnostics" / "rar_points.csv"
OUT_REBUILT = H2_ROOT / "data" / "derived" / "rar_points_rebuilt.csv"

def parse_galaxy_name_from_filename(p: Path) -> str:
    m = re.match(r"rc_decomp_(.+?)_best\.csv$", p.name)
    if not m:
        raise ValueError(f"Unrecognized per-galaxy filename: {p.name}")
    return m.group(1)

def load_rc_decomp(p: Path) -> pd.DataFrame:
    df = pd.read_csv(p)
    needed = {"R_kpc", "V_baryon", "V_total"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"{p.name}: missing columns {missing}; found {list(df.columns)}")
    # enforce float64
    for c in ["R_kpc", "V_baryon", "V_total"]:
        df[c] = df[c].astype(np.float64)
    return df

def compute_accel_natural(V_kms: np.ndarray, R_kpc: np.ndarray) -> np.ndarray:
    # Natural units: (km/s)^2 / kpc
    R_kpc = np.maximum(R_kpc, 1e-30)
    return (V_kms * V_kms) / R_kpc

def rebuild_rar_points() -> pd.DataFrame:
    rows = []
    files = sorted(H1_PER_GALAXY.glob("rc_decomp_*_best.csv"))
    if not files:
        raise RuntimeError(f"No per-galaxy files found in {H1_PER_GALAXY}")

    print(f"[Gate0] Found {len(files)} files in {H1_PER_GALAXY}")

    for p in files:
        gal = parse_galaxy_name_from_filename(p)
        df = load_rc_decomp(p)

        gbar = compute_accel_natural(df["V_baryon"].to_numpy(), df["R_kpc"].to_numpy())
        gobs = compute_accel_natural(df["V_total"].to_numpy(),  df["R_kpc"].to_numpy())

        out = pd.DataFrame({
            "galaxy": gal,
            "R_kpc": df["R_kpc"].to_numpy(),
            "g_bar": gbar,
            "g_obs": gobs,
        })
        
        # Calculate r_frac just for info, but we won't merge on it
        r = df["R_kpc"].to_numpy()
        rmax = float(np.max(r)) if len(r) else np.nan
        out["r_frac"] = r / rmax
        
        rows.append(out)

    rar = pd.concat(rows, ignore_index=True)

    # Sort by galaxy and radius to ensure stable order
    rar.sort_values(["galaxy", "R_kpc"], inplace=True, kind="mergesort")
    
    # ---------------------------------------------------------
    # CRITICAL FIX: Add 'idx' (row order) for stable merging
    # ---------------------------------------------------------
    rar["idx"] = rar.groupby("galaxy").cumcount()
    
    rar.reset_index(drop=True, inplace=True)
    return rar

def try_compare_to_h15(rebuilt: pd.DataFrame) -> None:
    if not H15_RAR_POINTS.exists():
        print(f"[Gate0] Reference file not found: {H15_RAR_POINTS}")
        return

    ref = pd.read_csv(H15_RAR_POINTS)
    print(f"[Gate0] Loaded H1.5 reference: {H15_RAR_POINTS} with {len(ref)} rows")
    
    # ---------------------------------------------------------
    # DIAGNOSTIC 1: Find the Intruder Galaxy
    # ---------------------------------------------------------
    ref_gals = set(ref["name"].unique())
    reb_gals = set(rebuilt["galaxy"].unique())
    
    extra_in_rebuilt = reb_gals - ref_gals
    missing_in_rebuilt = ref_gals - reb_gals
    
    if extra_in_rebuilt:
        print(f"\n[Gate0] ⚠️  FOUND EXTRA GALAXIES in Rebuilt: {extra_in_rebuilt}")
        print("         (This explains the 5280 vs 5250 mismatch. Exclude these to pass.)\n")
    if missing_in_rebuilt:
        print(f"[Gate0] ⚠️  MISSING GALAXIES in Rebuilt: {missing_in_rebuilt}")

    # Map columns
    cols = set(ref.columns)
    candidates = {
        "galaxy": ["galaxy", "name", "gal", "GALAXY"],
        "g_bar":  ["g_bar", "gbar", "g_baryon", "gbar_m_s2"],
        "g_obs":  ["g_obs", "gobs", "g_total", "gobs_m_s2"],
        "r_frac": ["r_frac", "rfrac", "r_fraction"],
    }
    mapping = {}
    for std, alts in candidates.items():
        for a in alts:
            if a in cols:
                mapping[std] = a
                break

    if set(mapping.keys()) != {"galaxy", "g_bar", "g_obs", "r_frac"}:
        print("[Gate0] Schema mismatch. Cannot compare pointwise.")
        return

    # Standardize Reference
    ref_std = pd.DataFrame({
        "galaxy": ref[mapping["galaxy"]].astype(str),
        "g_bar": pd.to_numeric(ref[mapping["g_bar"]], errors="coerce").astype(np.float64),
        "g_obs": pd.to_numeric(ref[mapping["g_obs"]], errors="coerce").astype(np.float64),
        "r_frac": pd.to_numeric(ref[mapping["r_frac"]], errors="coerce").astype(np.float64),
    }).dropna()

    # Sort and add index for stable merge
    ref_std.sort_values(["galaxy", "r_frac"], inplace=True, kind="mergesort")
    ref_std["idx"] = ref_std.groupby("galaxy").cumcount()

    # Merge on (galaxy, idx) NOT r_frac
    merged = rebuilt.merge(ref_std, on=["galaxy", "idx"], how="inner", suffixes=("_reb", "_ref"))
    
    print(f"[Gate0] Merge matched rows: {len(merged)}")
    
    if len(merged) == 0:
        print("[Gate0] No overlap. Check galaxy naming.")
        return

    # Compute Deltas
    d_gbar = np.max(np.abs(merged["g_bar_reb"] - merged["g_bar_ref"]))
    d_gobs = np.max(np.abs(merged["g_obs_reb"] - merged["g_obs_ref"]))

    print(f"[Gate0] max |Δg_bar| = {d_gbar:.3e}")
    print(f"[Gate0] max |Δg_obs| = {d_gobs:.3e}")

    # STRICT PASS CONDITION
    count_match = (len(rebuilt) == 5250) and (len(ref) == 5250)
    precision_pass = (d_gbar < 1e-10) and (d_gobs < 1e-10)

    if count_match and precision_pass:
        print("[Gate0] ✅ PASS: Perfect Identity with H1.5 (Count + Precision).")
    elif precision_pass:
        print("[Gate0] ⚠️  PASS (Precision) but FAIL (Count). Remove the extra galaxy!")
    else:
        print("[Gate0] ❌ FAIL: Precision mismatch.")

def main():
    OUT_REBUILT.parent.mkdir(parents=True, exist_ok=True)
    rebuilt = rebuild_rar_points()
    rebuilt.to_csv(OUT_REBUILT, index=False)
    print(f"[Gate0] Wrote rebuilt RAR points: {OUT_REBUILT}")
    print(f"[Gate0] Rebuilt rows = {len(rebuilt)} | galaxies = {rebuilt['galaxy'].nunique()}")
    try_compare_to_h15(rebuilt)

if __name__ == "__main__":
    main()