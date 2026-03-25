"""
h2_metric_reconstruction.py
============================
Reconstruct H2 inner-region scatter in the harmonized metric for all 80
fleet galaxies.

Harmonized metric: RMS of log10(V_model) - log10(V_obs) [dex]
Inner region:      R < 0.5 * R_max (R_max from observed SPARC rotmod grid)
Source:            Directly from diagnostics/test3_inner_scatter.py

Run from H2_ROOT:
    python comparative_analysis/metric_harmonization/h2_metric_reconstruction.py
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths (all relative to H2_ROOT, which must be the working directory)
# ---------------------------------------------------------------------------
H2_ROOT = Path(".").resolve()
FLEET_CSV = H2_ROOT / "H2_PUBLICATION_RELEASE/fleet_expansion_80galaxy/fleet_summary_80galaxy.csv"
SPARC_DIR = H2_ROOT / "data/sparc"
H2_OUTPUTS_DIR = H2_ROOT / "data/derived/phase4/h2_outputs"
H1_DIR = H2_ROOT / "data/h1_frozen/per_galaxy"
OUT_DIR = H2_ROOT / "comparative_analysis/metric_harmonization"
OUT_CSV = OUT_DIR / "h2_metric_summary.csv"


# ---------------------------------------------------------------------------
# Core metric functions (matching test3_inner_scatter.py exactly)
# ---------------------------------------------------------------------------

def rms_scatter_log10(V_model, V_obs, eps=1e-12):
    """RMS of log10(V_model) - log10(V_obs). Matches test3_inner_scatter.py."""
    V_model = np.clip(V_model, eps, None)
    V_obs = np.clip(V_obs, eps, None)
    resid = np.log10(V_model) - np.log10(V_obs)
    return float(np.sqrt(np.mean(resid ** 2)))


def load_rotmod(rotmod_path):
    """
    Parse SPARC rotmod file.
    Columns: R[kpc], Vobs[km/s], errV, Vgas, Vdisk, Vbul
    Skips comment lines (starting with # or !) and blank lines.
    Returns (R, Vobs) as numpy arrays.
    """
    data = []
    with open(rotmod_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('!'):
                continue
            try:
                vals = [float(x) for x in line.split()]
                if len(vals) >= 6:
                    data.append(vals[:6])
            except ValueError:
                continue
    if not data:
        raise RuntimeError(f"No data rows parsed from {rotmod_path}")
    data = np.array(data)
    return data[:, 0], data[:, 1]  # R_kpc, Vobs_kms


def get_vtotal(df, label, prefer_h2=False):
    """
    Extract total velocity column from a model CSV.

    For H2 files (prefer_h2=True): tries V_total_H2 first, then V_total,
    then V_total_H1, Vtot, V_model, V_tot.
    For H1 files: tries V_total first, then V_total_H1, Vtot, V_model, V_tot.
    """
    if prefer_h2:
        candidates = ["V_total_H2", "V_total", "V_total_H1", "Vtot", "V_model", "V_tot"]
    else:
        candidates = ["V_total", "V_total_H1", "Vtot", "V_model", "V_tot"]
    for c in candidates:
        if c in df.columns:
            return df[c].to_numpy(float), c
    raise RuntimeError(
        f"No total velocity column found in {label}. "
        f"Available columns: {list(df.columns)}"
    )


def compute_inner_sigma(R_obs, V_obs, R_model, V_model_grid, inner_frac=0.5):
    """
    Interpolate V_model onto observed R grid, compute inner-region log10 RMS.
    Returns (sigma, n_inner, R_cut).
    """
    V_model_interp = np.interp(R_obs, R_model, V_model_grid)
    R_max = float(np.max(R_obs))
    R_cut = inner_frac * R_max
    mask = R_obs < R_cut
    n_inner = int(np.sum(mask))
    if n_inner == 0:
        return np.nan, 0, R_cut
    sigma = rms_scatter_log10(V_model_interp[mask], V_obs[mask])
    return sigma, n_inner, R_cut


# ---------------------------------------------------------------------------
# Main reconstruction
# ---------------------------------------------------------------------------

def process_galaxy(galaxy, regime):
    """
    Process one galaxy. Returns a dict with all output columns.
    """
    result = {
        "Galaxy": galaxy,
        "Regime": regime,
        "inner_region_definition": "R < 0.5 * R_max",
        "n_inner_points": np.nan,
        "h2_sigma_baseline": np.nan,
        "h1_sigma_baseline": np.nan,
        "h2_delta_sigma_harmonized": np.nan,
        "reconstruction_status": "failed",
        "notes": "",
    }

    # --- Load rotmod ---
    rotmod_path = SPARC_DIR / f"{galaxy}_rotmod.dat"
    if not rotmod_path.exists():
        result["notes"] = f"rotmod not found: {rotmod_path.name}"
        return result
    try:
        R_obs, V_obs = load_rotmod(rotmod_path)
    except Exception as e:
        result["notes"] = f"rotmod load error: {e}"
        return result

    # --- Try H2 model file ---
    h2_paths = [
        H2_OUTPUTS_DIR / f"rc_decomp_{galaxy}_H2_adaptive.csv",
        H2_ROOT / "data/derived/phase4" / f"rc_decomp_{galaxy}_H2_adaptive.csv",
    ]
    h2_csv = None
    for p in h2_paths:
        if p.exists():
            h2_csv = p
            break

    # --- Load H1 model file ---
    h1_path = H1_DIR / f"rc_decomp_{galaxy}_best.csv"
    h1_available = h1_path.exists()

    # --- Compute H1 sigma ---
    sigma_h1 = np.nan
    if h1_available:
        try:
            h1_df = pd.read_csv(h1_path)
            if "R_kpc" not in h1_df.columns:
                result["notes"] += f" H1 missing R_kpc."
                h1_available = False
            else:
                R_h1 = h1_df["R_kpc"].to_numpy(float)
                V_h1_grid, col_h1 = get_vtotal(h1_df, f"H1/{galaxy}", prefer_h2=False)
                sigma_h1, n_inner, _ = compute_inner_sigma(R_obs, V_obs, R_h1, V_h1_grid)
                result["h1_sigma_baseline"] = sigma_h1
                result["n_inner_points"] = n_inner
        except Exception as e:
            result["notes"] += f" H1 error: {e}."
            h1_available = False

    # --- Compute H2 sigma ---
    if h2_csv is not None:
        try:
            h2_df = pd.read_csv(h2_csv)
            if "R_kpc" not in h2_df.columns:
                result["notes"] += f" H2 missing R_kpc."
                # Fall through to H1-only or failed
            else:
                R_h2 = h2_df["R_kpc"].to_numpy(float)
                V_h2_grid, col_h2 = get_vtotal(h2_df, f"H2/{galaxy}", prefer_h2=True)
                sigma_h2, n_inner_h2, _ = compute_inner_sigma(R_obs, V_obs, R_h2, V_h2_grid)
                result["h2_sigma_baseline"] = sigma_h2
                # Use n_inner from H2 (consistent R grid from same rotmod)
                result["n_inner_points"] = n_inner_h2

                if not np.isnan(sigma_h1) and not np.isnan(sigma_h2):
                    result["h2_delta_sigma_harmonized"] = sigma_h2 - sigma_h1
                    result["reconstruction_status"] = "exact"
                    result["notes"] += f" H2_col={col_h2} H1_col={col_h1}"
                elif np.isnan(sigma_h1):
                    result["h2_delta_sigma_harmonized"] = np.nan
                    result["reconstruction_status"] = "exact_no_h1"
                    result["h2_sigma_baseline"] = sigma_h2
                    result["notes"] += f" H2 found but H1 failed; delta_sigma=NaN."
        except Exception as e:
            result["notes"] += f" H2 error: {e}."
            # Fall through to H1-only or failed

    # Determine final status
    if result["reconstruction_status"] == "failed":
        if h1_available and not np.isnan(sigma_h1):
            result["reconstruction_status"] = "H1_only"
            result["h2_delta_sigma_harmonized"] = np.nan
            result["notes"] = (result["notes"] +
                               " H2 file not found; H1 sigma computed only.").strip()
        elif h2_csv is None and not h1_available:
            result["reconstruction_status"] = "failed"
            result["notes"] = (result["notes"] +
                               " Both H2 and H1 files missing or failed.").strip()

    result["notes"] = result["notes"].strip()
    return result


def main():
    print("=" * 65)
    print("H2 METRIC RECONSTRUCTION")
    print("Harmonized metric: RMS log10(V_model/V_obs) [dex]")
    print("Inner region: R < 0.5 * R_max")
    print("=" * 65)

    # Load fleet
    fleet = pd.read_csv(FLEET_CSV)
    # Ensure output directory exists
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    records = []
    status_counts = {"exact": 0, "exact_no_h1": 0, "H1_only": 0, "H2_file_missing": 0, "failed": 0}

    for _, row in fleet.iterrows():
        galaxy = str(row["galaxy"]).strip()
        regime = str(row["regime"]).strip()

        result = process_galaxy(galaxy, regime)

        # Reclassify H1_only as H2_file_missing for clarity
        if result["reconstruction_status"] == "H1_only":
            result["reconstruction_status"] = "H2_file_missing"

        status = result["reconstruction_status"]
        status_counts[status] = status_counts.get(status, 0) + 1

        records.append(result)

        marker = {
            "exact": "[EXACT]",
            "exact_no_h1": "[H2_ONLY]",
            "H2_file_missing": "[H1_ONLY]",
            "failed": "[FAILED]",
        }.get(status, f"[{status}]")

        sigma_h2 = result["h2_sigma_baseline"]
        sigma_h1 = result["h1_sigma_baseline"]
        delta = result["h2_delta_sigma_harmonized"]
        n = result["n_inner_points"]

        print(
            f"  {marker:12s} {galaxy:<15s} {regime:<12s} "
            f"n={str(n) if not (isinstance(n, float) and np.isnan(n)) else 'NaN':>4s} "
            f"sigma_H2={sigma_h2:.4f} sigma_H1={sigma_h1:.4f} delta={delta:.4f}"
            if (not np.isnan(sigma_h2) if isinstance(sigma_h2, float) else True)
            and (not np.isnan(sigma_h1) if isinstance(sigma_h1, float) else True)
            else f"  {marker:12s} {galaxy:<15s} {regime:<12s} "
                 f"  {result['notes']}"
        )

    df_out = pd.DataFrame(records)
    df_out.to_csv(OUT_CSV, index=False, float_format="%.6f")

    print()
    print("=" * 65)
    print("RECONSTRUCTION STATUS SUMMARY")
    print("=" * 65)
    for status, count in sorted(status_counts.items()):
        print(f"  {status:<20s}: {count:3d} galaxies")
    print(f"  {'TOTAL':<20s}: {len(records):3d} galaxies")
    print()
    print(f"Output written to: {OUT_CSV}")
    print("=" * 65)

    # Quick stats for exact reconstructions
    exact = df_out[df_out["reconstruction_status"] == "exact"]
    if len(exact) > 0:
        print()
        print(f"EXACT H2 DELTA-SIGMA STATS ({len(exact)} galaxies):")
        delta_vals = exact["h2_delta_sigma_harmonized"].dropna()
        print(f"  Median delta_sigma: {delta_vals.median():.4f} dex")
        print(f"  Mean delta_sigma:   {delta_vals.mean():.4f} dex")
        print(f"  Std delta_sigma:    {delta_vals.std():.4f} dex")
        print(f"  Min delta_sigma:    {delta_vals.min():.4f} dex")
        print(f"  Max delta_sigma:    {delta_vals.max():.4f} dex")
        print(f"  Fraction |delta|<0.01: "
              f"{(delta_vals.abs() < 0.01).sum()}/{len(delta_vals)}")

    # H1-only sigma stats
    h1_only = df_out[df_out["reconstruction_status"] == "H2_file_missing"]
    if len(h1_only) > 0:
        print()
        print(f"H1-ONLY SIGMA BASELINE STATS ({len(h1_only)} galaxies):")
        sigma_vals = h1_only["h1_sigma_baseline"].dropna()
        print(f"  Median sigma_H1: {sigma_vals.median():.4f} dex")
        print(f"  Mean sigma_H1:   {sigma_vals.mean():.4f} dex")


if __name__ == "__main__":
    main()
