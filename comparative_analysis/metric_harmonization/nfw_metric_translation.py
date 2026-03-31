"""
nfw_metric_translation.py
==========================
Load existing NFW perturbation analysis results and reformat them to
the harmonized naming convention for direct comparison with H2.

The NFW analysis already uses the same metric as H2:
  RMS of log10(V_model) - log10(V_obs) [dex], inner region R < 0.5 * R_max

No metric conversion is needed. This script documents that fact and
produces a harmonized-named output CSV.

Run from H2_ROOT:
    python comparative_analysis/metric_harmonization/nfw_metric_translation.py
"""

import numpy as np
import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
H2_ROOT = Path(".").resolve()
NFW_CSV = H2_ROOT / "comparative_analysis/nfw/nfw_perturbation_summary.csv"
OUT_DIR = H2_ROOT / "comparative_analysis/metric_harmonization"
OUT_CSV = OUT_DIR / "nfw_metric_harmonized_summary.csv"


def main():
    print("=" * 65)
    print("NFW METRIC TRANSLATION")
    print("=" * 65)

    if not NFW_CSV.exists():
        raise FileNotFoundError(f"NFW summary not found: {NFW_CSV}")

    nfw = pd.read_csv(NFW_CSV)
    print(f"Loaded NFW summary: {len(nfw)} galaxies")
    print(f"NFW columns: {list(nfw.columns)}")

    # --- Confirm metric alignment ---
    # NFW columns:
    #   Galaxy, Regime, V200_kms, C200, rs_kpc,
    #   sigma_baseline,          <- RMS log10 dex of NFW baseline fit
    #   n_inner_points,          <- count of obs points in inner region
    #   inner_region_defined,    <- True if n_inner >= 3
    #   delta_sigma_mass_p10, delta_sigma_mass_p20, ...  <- perturbation deltas
    #   max_abs_delta_sigma,     <- max |delta_sigma| across all perturbations
    #   V_bar_over_V_NFW_inner_median

    print()
    print("METRIC ALIGNMENT CONFIRMATION")
    print("-" * 40)
    print("NFW uses: RMS of log10(V_model) - log10(V_obs) [dex]")
    print("H2  uses: RMS of log10(V_model) - log10(V_obs) [dex]")
    print("Inner region (NFW): R < 0.5 * R_max  [confirmed from nfw_perturbation_diagnostic.py]")
    print("Inner region (H2):  R < 0.5 * R_max  [from diagnostics/test3_inner_scatter.py]")
    print("=> NO metric conversion required. Direct comparison is valid.")
    print()

    # --- Build harmonized output ---
    out_records = []
    for _, row in nfw.iterrows():
        galaxy = str(row["Galaxy"]).strip()
        regime = str(row["Regime"]).strip()
        n_inner = int(row["n_inner_points"])
        inner_defined = bool(row["inner_region_defined"])
        sigma_baseline = float(row["sigma_baseline"])
        max_abs_delta = float(row["max_abs_delta_sigma"])

        # Determine metric_status
        if not inner_defined or n_inner < 3:
            metric_status = "insufficient_inner_points"
        else:
            metric_status = "valid"

        # Notes: record perturbation scenario achieving max
        # Identify which scenario gives max_abs_delta_sigma
        delta_cols = [c for c in nfw.columns if c.startswith("delta_sigma_")]
        max_scenario = "unknown"
        if delta_cols:
            row_deltas = {c: abs(float(row[c])) for c in delta_cols
                          if not pd.isna(row[c])}
            if row_deltas:
                max_scenario = max(row_deltas, key=row_deltas.get)

        notes = (
            f"NFW fit: V200={row['V200_kms']:.1f} km/s, "
            f"C200={row['C200']:.2f}, rs={row['rs_kpc']:.2f} kpc; "
            f"max perturbation scenario: {max_scenario}"
        )

        out_records.append({
            "Galaxy": galaxy,
            "Regime": regime,
            "n_inner_points": n_inner,
            "nfw_sigma_baseline": sigma_baseline,
            "nfw_max_abs_delta_sigma_harmonized": max_abs_delta,
            "metric_status": metric_status,
            "notes": notes,
        })

    df_out = pd.DataFrame(out_records)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(OUT_CSV, index=False, float_format="%.6f")

    print(f"Galaxies with valid inner region (n>=3): "
          f"{(df_out['metric_status'] == 'valid').sum()}/{len(df_out)}")
    print(f"Galaxies with insufficient inner points: "
          f"{(df_out['metric_status'] == 'insufficient_inner_points').sum()}/{len(df_out)}")
    print()

    # --- Summary stats ---
    valid = df_out[df_out["metric_status"] == "valid"]
    print(f"NFW MAX ABS DELTA-SIGMA STATS ({len(valid)} valid galaxies):")
    vals = valid["nfw_max_abs_delta_sigma_harmonized"]
    print(f"  Median max|delta_sigma|: {vals.median():.4f} dex")
    print(f"  Mean max|delta_sigma|:   {vals.mean():.4f} dex")
    print(f"  Std  max|delta_sigma|:   {vals.std():.4f} dex")
    print(f"  Min  max|delta_sigma|:   {vals.min():.4f} dex")
    print(f"  Max  max|delta_sigma|:   {vals.max():.4f} dex")
    print(f"  Fraction <= 0.01 dex:    {(vals <= 0.01).sum()}/{len(vals)}")
    print(f"  Fraction > 0.05 dex:     {(vals > 0.05).sum()}/{len(vals)}")
    print()

    # Regime-stratified
    print("BY REGIME:")
    for regime in ["baryon-dom", "balanced", "DM-dom"]:
        sub = valid[valid["Regime"] == regime]
        if len(sub) == 0:
            continue
        v = sub["nfw_max_abs_delta_sigma_harmonized"]
        print(f"  {regime:<12s}: n={len(sub):2d}  "
              f"median={v.median():.4f}  mean={v.mean():.4f}  "
              f"max={v.max():.4f} dex")

    print()
    print(f"Output written to: {OUT_CSV}")
    print("=" * 65)


if __name__ == "__main__":
    main()
