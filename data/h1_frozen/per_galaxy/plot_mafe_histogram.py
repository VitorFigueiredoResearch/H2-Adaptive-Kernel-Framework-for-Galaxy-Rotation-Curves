import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# -------------------------------
# Configuration
# -------------------------------
RESULTS_DIR = Path(".")
SUMMARY_FILE = RESULTS_DIR / "sparc_lite_summary.csv"
OUTPUT_FIG = RESULTS_DIR / "fig_mafe_histogram.pdf"

# acceptance thresholds (as used in paper)
THRESH_10 = 0.10
THRESH_15 = 0.15

# -------------------------------
# Load data
# -------------------------------
if not SUMMARY_FILE.exists():
    raise FileNotFoundError(f"Missing file: {SUMMARY_FILE}")

df = pd.read_csv(SUMMARY_FILE)

if "mafe" not in df.columns:
    raise ValueError("Expected column 'mafe' not found in sparc_lite_summary.csv")

mafe = df["mafe"].dropna().values

n_total = len(mafe)
median_mafe = np.median(mafe)
q16, q84 = np.percentile(mafe, [16, 84])
frac_10 = np.mean(mafe < THRESH_10)
frac_15 = np.mean(mafe < THRESH_15)

# -------------------------------
# Plot
# -------------------------------
plt.figure(figsize=(7, 5))

bins = np.linspace(0.0, max(mafe.max(), 0.3), 25)
plt.hist(mafe, bins=bins, edgecolor="black", alpha=0.75)

# Threshold lines
plt.axvline(THRESH_10, linestyle="--", linewidth=1.5, label=r"MAFE = 0.10")
plt.axvline(THRESH_15, linestyle="--", linewidth=1.5, label=r"MAFE = 0.15")

# Median line
plt.axvline(median_mafe, linestyle="-", linewidth=2.0, label=f"Median = {median_mafe:.3f}")

plt.xlabel("Median Absolute Fractional Error (MAFE)")
plt.ylabel("Number of Galaxies")
plt.title("Fleet-wide H1 Rotation Curve Performance (N = {})".format(n_total))
plt.legend(frameon=False)

plt.tight_layout()
plt.savefig(OUTPUT_FIG)
plt.close()

# -------------------------------
# Console summary (for paper text)
# -------------------------------
print("[OK] MAFE histogram generated:", OUTPUT_FIG)
print(f"  Galaxies: {n_total}")
print(f"  Median MAFE: {median_mafe:.4f}")
print(f"  16–84 percentile: [{q16:.4f}, {q84:.4f}]")
print(f"  Fraction MAFE < 0.10: {frac_10:.2%}")
print(f"  Fraction MAFE < 0.15: {frac_15:.2%}")
