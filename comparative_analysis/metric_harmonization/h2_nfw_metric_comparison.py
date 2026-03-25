"""
h2_nfw_metric_comparison.py
============================
Direct H2 vs NFW comparison in the harmonized metric.

Loads:
  comparative_analysis/metric_harmonization/h2_metric_summary.csv
  comparative_analysis/metric_harmonization/nfw_metric_harmonized_summary.csv

Outputs:
  comparative_analysis/metric_harmonization/h2_nfw_metric_comparison_summary.csv

Run from H2_ROOT:
    python comparative_analysis/metric_harmonization/h2_nfw_metric_comparison.py
"""

import numpy as np
import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
H2_ROOT = Path(".").resolve()
OUT_DIR = H2_ROOT / "comparative_analysis/metric_harmonization"
H2_CSV = OUT_DIR / "h2_metric_summary.csv"
NFW_CSV = OUT_DIR / "nfw_metric_harmonized_summary.csv"
OUT_COMPARE = OUT_DIR / "h2_nfw_metric_comparison_summary.csv"


def regime_stats(df, col, label):
    """Print regime-stratified statistics for a column."""
    vals = df[col].dropna()
    if len(vals) == 0:
        print(f"  {label}: no valid data")
        return
    print(f"  {label} (n={len(vals)}):")
    print(f"    Median: {vals.median():.4f} dex")
    print(f"    Mean:   {vals.mean():.4f} dex")
    print(f"    Std:    {vals.std():.4f} dex")
    print(f"    Fraction <= 0.01 dex: {(vals.abs() <= 0.01).sum()}/{len(vals)}")
    print(f"    Fraction >  0.05 dex: {(vals.abs() > 0.05).sum()}/{len(vals)}")


def main():
    print("=" * 65)
    print("H2 vs NFW METRIC COMPARISON")
    print("Harmonized metric: RMS log10(V_model/V_obs) [dex]")
    print("=" * 65)

    # --- Load ---
    h2 = pd.read_csv(H2_CSV)
    nfw = pd.read_csv(NFW_CSV)

    print(f"H2 summary: {len(h2)} galaxies")
    print(f"NFW summary: {len(nfw)} galaxies")

    # --- Merge on Galaxy ---
    merged = pd.merge(h2, nfw, on="Galaxy", how="inner", suffixes=("_h2", "_nfw"))
    print(f"Merged on Galaxy: {len(merged)} galaxies in common")

    # --- Validity filter ---
    # Valid: n_inner >= 3 in both, reconstruction_status != "failed"
    # n_inner_points_h2 from h2 summary, n_inner_points from nfw summary
    # Rename nfw n_inner_points
    merged = merged.rename(columns={
        "n_inner_points_h2": "n_inner_h2",
        "n_inner_points": "n_inner_nfw",
        "Regime_h2": "Regime",
    })
    if "Regime_nfw" in merged.columns:
        merged = merged.drop(columns=["Regime_nfw"])

    # Valid NFW: n_inner_nfw >= 3
    merged["nfw_valid"] = merged["n_inner_nfw"] >= 3

    # H2 exact: reconstruction_status == "exact"
    merged["h2_exact"] = merged["reconstruction_status"] == "exact"

    # H2 baseline valid: H1 sigma is available (H2_file_missing still has h1_sigma_baseline)
    merged["h2_baseline_valid"] = (
        merged["reconstruction_status"].isin(["exact", "H2_file_missing"])
        & merged["h1_sigma_baseline"].notna()
        & (merged["n_inner_h2"].fillna(0) >= 3)
    )

    print()
    print("COVERAGE:")
    print(f"  Galaxies with exact H2 delta_sigma: {merged['h2_exact'].sum()}")
    print(f"  Galaxies with valid NFW delta_sigma: {merged['nfw_valid'].sum()}")
    print(f"  Galaxies exact H2 AND valid NFW: "
          f"{(merged['h2_exact'] & merged['nfw_valid']).sum()}")
    print(f"  Galaxies H2 baseline valid AND valid NFW: "
          f"{(merged['h2_baseline_valid'] & merged['nfw_valid']).sum()}")

    # --- Absolute H2 delta_sigma ---
    merged["h2_abs_delta_sigma"] = merged["h2_delta_sigma_harmonized"].abs()

    # --- Build output ---
    cols_out = [
        "Galaxy", "Regime",
        "n_inner_h2", "n_inner_nfw",
        "h1_sigma_baseline",        # H1 baseline sigma (available for all)
        "h2_sigma_baseline",        # H2 sigma (only for exact)
        "h2_delta_sigma_harmonized",
        "h2_abs_delta_sigma",
        "nfw_sigma_baseline",
        "nfw_max_abs_delta_sigma_harmonized",
        "reconstruction_status",
        "nfw_valid",
        "h2_exact",
        "h2_baseline_valid",
    ]
    # Add notes columns
    for c in ["notes_h2", "notes_nfw"]:
        if c in merged.columns:
            cols_out.append(c)

    df_out = merged[[c for c in cols_out if c in merged.columns]].copy()
    df_out.to_csv(OUT_COMPARE, index=False, float_format="%.6f")
    print()
    print(f"Comparison summary written to: {OUT_COMPARE}")

    # --- Print first 10 rows ---
    print()
    print("FIRST 10 ROWS OF COMPARISON SUMMARY:")
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 200)
    print(df_out.head(10).to_string(index=False))

    # =======================================================================
    # STATS: Exact H2 vs NFW (9 galaxies)
    # =======================================================================
    print()
    print("=" * 65)
    print("COMPARISON STATS: EXACT H2 vs NFW")
    print("(9 galaxies where H2 adaptive files exist)")
    print("=" * 65)

    exact_valid = df_out[df_out["h2_exact"] & df_out["nfw_valid"]].copy()
    print(f"N galaxies: {len(exact_valid)}")

    if len(exact_valid) > 0:
        print()
        print("H2 |delta_sigma| [dex]:")
        regime_stats(exact_valid, "h2_abs_delta_sigma", "  ALL")
        for regime in ["baryon-dom", "balanced", "DM-dom"]:
            sub = exact_valid[exact_valid["Regime"] == regime]
            if len(sub) > 0:
                regime_stats(sub, "h2_abs_delta_sigma", f"  {regime}")

        print()
        print("NFW max|delta_sigma| [dex]:")
        regime_stats(exact_valid, "nfw_max_abs_delta_sigma_harmonized", "  ALL")
        for regime in ["baryon-dom", "balanced", "DM-dom"]:
            sub = exact_valid[exact_valid["Regime"] == regime]
            if len(sub) > 0:
                regime_stats(sub, "nfw_max_abs_delta_sigma_harmonized", f"  {regime}")

        print()
        print("GALAXY-BY-GALAXY (exact H2):")
        print(f"  {'Galaxy':<15s} {'Regime':<12s} "
              f"{'H2 |delta|':>10s} {'NFW max|delta|':>14s} {'Ratio NFW/H2':>12s}")
        print("  " + "-" * 60)
        for _, row in exact_valid.iterrows():
            h2d = row["h2_abs_delta_sigma"]
            nfwd = row["nfw_max_abs_delta_sigma_harmonized"]
            ratio = nfwd / h2d if (not np.isnan(h2d) and h2d > 0) else np.nan
            print(f"  {row['Galaxy']:<15s} {row['Regime']:<12s} "
                  f"{h2d:>10.4f} {nfwd:>14.4f} {ratio:>12.2f}")

    # =======================================================================
    # STATS: H1 sigma baseline comparison (all valid NFW galaxies)
    # =======================================================================
    print()
    print("=" * 65)
    print("BASELINE SIGMA COMPARISON: H1 vs NFW (all valid galaxies)")
    print("(H1 sigma = what NFW baseline was fitted against)")
    print("=" * 65)

    baseline_valid = df_out[df_out["h2_baseline_valid"] & df_out["nfw_valid"]].copy()
    print(f"N galaxies: {len(baseline_valid)}")

    if len(baseline_valid) > 0:
        sigma_diff = baseline_valid["h1_sigma_baseline"] - baseline_valid["nfw_sigma_baseline"]
        print(f"  Median H1_sigma - NFW_sigma: {sigma_diff.median():.4f} dex")
        print(f"  Mean:                        {sigma_diff.mean():.4f} dex")
        print(f"  Std:                         {sigma_diff.std():.4f} dex")
        print()
        print("  BY REGIME:")
        for regime in ["baryon-dom", "balanced", "DM-dom"]:
            sub = baseline_valid[baseline_valid["Regime"] == regime]
            if len(sub) > 0:
                d = sub["h1_sigma_baseline"] - sub["nfw_sigma_baseline"]
                print(f"    {regime:<12s}: n={len(sub):2d}  "
                      f"median={d.median():.4f}  mean={d.mean():.4f} dex")

    # =======================================================================
    # KEY FINDING
    # =======================================================================
    print()
    print("=" * 65)
    print("KEY FINDING")
    print("=" * 65)

    if len(exact_valid) > 0:
        h2_med = exact_valid["h2_abs_delta_sigma"].median()
        nfw_med = exact_valid["nfw_max_abs_delta_sigma_harmonized"].median()
        ratio = nfw_med / h2_med if h2_med > 0 else float("inf")
        print(f"For the {len(exact_valid)} galaxies with exact H2 reconstructions:")
        print(f"  H2 median |delta_sigma|:   {h2_med:.4f} dex")
        print(f"  NFW median max|delta_sigma|: {nfw_med:.4f} dex")
        print(f"  Ratio NFW/H2:              {ratio:.1f}x")
        if nfw_med > h2_med:
            print(f"  => NFW produces {ratio:.1f}x larger inner scatter perturbations than H2")
        else:
            print(f"  => H2 produces larger inner scatter than NFW")

    # Full NFW population stats
    all_nfw_valid = df_out[df_out["nfw_valid"]]
    nfw_full_med = all_nfw_valid["nfw_max_abs_delta_sigma_harmonized"].median()
    nfw_frac_gt_01 = (all_nfw_valid["nfw_max_abs_delta_sigma_harmonized"] > 0.01).sum()
    print()
    print(f"For all {len(all_nfw_valid)} NFW-valid galaxies:")
    print(f"  NFW median max|delta_sigma|: {nfw_full_med:.4f} dex")
    print(f"  NFW fraction > 0.01 dex: {nfw_frac_gt_01}/{len(all_nfw_valid)}")
    print(f"  (H2 claims delta_sigma ~ 0; NFW shows systematic perturbation sensitivity)")

    print("=" * 65)


if __name__ == "__main__":
    main()
