import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ===============================
# Configuration (FROZEN)
# ===============================
INPUT_CSV = "fleet_summary_compact.csv"
OUTPUT_PDF = "fig_parameter_distributions.pdf"

# Histogram settings
NBINS_L = 20
NBINS_MU = 20

# ===============================
# Load data
# ===============================
df = pd.read_csv(INPUT_CSV)

L_vals = df["best_L"].values
mu_vals = df["best_mu"].values

n_gal = len(df)

# ===============================
# Figure
# ===============================
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

# --- Panel (a): L distribution ---
axes[0].hist(
    L_vals,
    bins=NBINS_L,
    edgecolor="black",
    alpha=0.8
)
axes[0].set_xlabel(r"$L$ [kpc]")
axes[0].set_ylabel("Number of galaxies")
axes[0].set_title(r"(a) Distribution of $L$")

# --- Panel (b): mu distribution (log scale) ---
axes[1].hist(
    mu_vals,
    bins=NBINS_MU,
    edgecolor="black",
    alpha=0.8
)
axes[1].set_xscale("log")
axes[1].set_xlabel(r"$\mu$")
axes[1].set_ylabel("Number of galaxies")
axes[1].set_title(r"(b) Distribution of $\mu$")

plt.tight_layout()
plt.savefig(OUTPUT_PDF)
plt.close()

# ===============================
# Console summary (diagnostic only)
# ===============================
print(f"[OK] Parameter distribution figure generated: {OUTPUT_PDF}")
print(f"  Galaxies included: {n_gal}")
print(f"  L range:  min={L_vals.min():.2f}, max={L_vals.max():.2f}")
print(f"  mu range: min={mu_vals.min():.2f}, max={mu_vals.max():.2f}")
