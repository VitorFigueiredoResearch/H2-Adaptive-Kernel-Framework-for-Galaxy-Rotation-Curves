import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(BASE_DIR, "..", "figures")

os.makedirs(FIG_DIR, exist_ok=True)

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv("fleet_summary_compact.csv")

# Drop any pathological rows (safety)
df = df.replace([np.inf, -np.inf], np.nan).dropna(
    subset=["best_L", "best_mu", "max_v_baryon"]
)

# Proxy for baryonic mass:
# v_baryon^2 ~ GM/R → monotonic proxy is sufficient for scaling
mass_proxy = df["max_v_baryon"]**2

L = df["best_L"]
mu = df["best_mu"]

# -----------------------------
# Create figure
# -----------------------------
plt.figure(figsize=(15, 4))

# Panel 1: L vs mass proxy
plt.subplot(1, 3, 1)
plt.scatter(mass_proxy, L, s=12, alpha=0.7)
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"Baryonic mass proxy ($V_{\mathrm{baryon,max}}^2$)")
plt.ylabel(r"Kernel scale $L$")
plt.title(r"$L$ vs baryonic mass proxy")

# Panel 2: mu vs mass proxy
plt.subplot(1, 3, 2)
plt.scatter(mass_proxy, mu, s=12, alpha=0.7)
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"Baryonic mass proxy ($V_{\mathrm{baryon,max}}^2$)")
plt.ylabel(r"Kernel amplitude $\mu$")
plt.title(r"$\mu$ vs baryonic mass proxy")

# Panel 3: L vs mu
plt.subplot(1, 3, 3)
plt.scatter(L, mu, s=12, alpha=0.7)
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$L$")
plt.ylabel(r"$\mu$")
plt.title(r"$\mu$ vs $L$")

plt.tight_layout()

# -----------------------------
# Save
# -----------------------------
outpath = os.path.join(FIG_DIR, "fig_scaling_relations.pdf")
plt.savefig(outpath, bbox_inches="tight")
print(f"[OK] Scaling relations figure saved: {outpath}")
