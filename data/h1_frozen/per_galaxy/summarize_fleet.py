import os
import glob
import json
import warnings

import numpy as np
import pandas as pd


# ------------------------------------------------------------
# Configuration (easy to tweak later)
# ------------------------------------------------------------

RESULTS_DIR = os.path.abspath(os.path.dirname(__file__))
SUMMARY_CSV = os.path.join(RESULTS_DIR, "sparc_lite_summary.csv")

OUT_CSV = os.path.join(RESULTS_DIR, "fleet_summary_compact.csv")
OUT_JSON = os.path.join(RESULTS_DIR, "fleet_summary_compact.json")

KERNEL_ACTIVE_THRESHOLD = 0.10     # max(V_kernel) / max(V_baryon)
INNER_FRACTION = 0.25              # inner region = first 25% of radii
OVERSHOOT_THRESHOLD = 0.10         # 10% excess


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def safe_read_csv(path):
    try:
        return pd.read_csv(path)
    except Exception as e:
        warnings.warn(f"Failed to read {path}: {e}")
        return None


# ------------------------------------------------------------
# Main summarization
# ------------------------------------------------------------

def summarize_fleet():
    # Load fleet fit summary
    if not os.path.exists(SUMMARY_CSV):
        raise FileNotFoundError(f"Missing {SUMMARY_CSV}")

    summary_df = pd.read_csv(SUMMARY_CSV)
    summary_df = summary_df.set_index("name")

    records = []

    # Find all per-galaxy decomposition files
    pattern = os.path.join(RESULTS_DIR, "rc_decomp_*_best.csv")
    files = sorted(glob.glob(pattern))

    if len(files) == 0:
        raise RuntimeError("No rc_decomp_*_best.csv files found.")

    for path in files:
        fname = os.path.basename(path)
        gal = fname.replace("rc_decomp_", "").replace("_best.csv", "")

        if gal not in summary_df.index:
            warnings.warn(f"{gal}: missing from sparc_lite_summary.csv — skipped")
            continue

        df = safe_read_csv(path)
        if df is None:
            continue

        # Required columns check
        required = {"R_kpc", "V_baryon", "V_kernel", "V_total"}
        if not required.issubset(df.columns):
            warnings.warn(f"{gal}: missing required columns — skipped")
            continue

        # Extract arrays
        R = df["R_kpc"].to_numpy()
        Vb = df["V_baryon"].to_numpy()
        Vk = df["V_kernel"].to_numpy()
        Vt = df["V_total"].to_numpy()

        # Basic sanity
        if len(R) < 5 or np.all(Vb == 0):
            warnings.warn(f"{gal}: insufficient or degenerate data — skipped")
            continue

        # Amplitudes
        max_v_baryon = float(np.nanmax(Vb))
        max_v_kernel = float(np.nanmax(np.abs(Vk)))
        max_v_total = float(np.nanmax(Vt))

        ratio_kb = max_v_kernel / max_v_baryon if max_v_baryon > 0 else 0.0

        # Inner-region behavior
        n_inner = max(1, int(len(R) * INNER_FRACTION))
        inner_excess = np.mean(
            (Vt[:n_inner] - Vb[:n_inner]) / np.maximum(Vb[:n_inner], 1e-6)
        )

        inner_overshoot = inner_excess > OVERSHOOT_THRESHOLD

        # Outer behavior proxy (trend at last bins)
        outer_decay = np.mean(Vk[-3:]) < np.mean(Vk[:3])

        record = {
            "name": gal,
            "mafe": float(summary_df.loc[gal, "mafe"]),
            "best_L": float(summary_df.loc[gal, "best_L"]),
            "best_mu": float(summary_df.loc[gal, "best_mu"]),
            "max_v_baryon": max_v_baryon,
            "max_v_kernel": max_v_kernel,
            "max_v_total": max_v_total,
            "kernel_to_baryon_ratio": ratio_kb,
            "kernel_active": ratio_kb > KERNEL_ACTIVE_THRESHOLD,
            "inner_overshoot": bool(inner_overshoot),
            "outer_decay": bool(outer_decay),
        }

        records.append(record)

    if len(records) == 0:
        raise RuntimeError("No valid galaxies summarized.")

    out_df = pd.DataFrame.from_records(records)
    out_df = out_df.sort_values("mafe").reset_index(drop=True)

    # Write outputs
    out_df.to_csv(OUT_CSV, index=False)

    with open(OUT_JSON, "w") as f:
        json.dump(out_df.to_dict(orient="records"), f, indent=2)

    print(f"[OK] Fleet summary written:")
    print(f"  - {OUT_CSV}")
    print(f"  - {OUT_JSON}")
    print(f"  Galaxies summarized: {len(out_df)}")


# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------

if __name__ == "__main__":
    summarize_fleet()
