import pandas as pd
import numpy as np

g = "NGC3198"

p3 = pd.read_csv(f"data/derived/phase3/leff_{g}.csv").sort_values("R_kpc")
p4 = pd.read_csv(f"data/derived/phase4/h2_outputs/rc_decomp_{g}_H2_adaptive.csv").sort_values("R_kpc")

# chi
if "chi_used" not in p3.columns:
    raise RuntimeError(f"chi_used not found in phase3 file. Columns={list(p3.columns)}")

# align
if len(p3) != len(p4) or not np.allclose(p3.R_kpc.values, p4.R_kpc.values, atol=1e-10, rtol=0):
    raise RuntimeError("R grids mismatch between phase3 and phase4")

R = p4.R_kpc.values
rmax = R.max()

# inner+mid region where adaptation acts (same mask you used)
mask = (R > 0) & ((R / rmax) < 0.3)

chi = p3["chi_used"].values[mask]

Vh1 = p4.V_total_H1.values[mask]
Vb  = p4.V_baryon.values[mask]
Rk  = R[mask]

# build H1 RAR residual proxy: log10(g_obs_H1) - log10(g_bar)
# constants cancel, but keep it explicit
gobs = (Vh1 * 1e3) ** 2 / (Rk * 3.085677581e19)
gbar = (Vb  * 1e3) ** 2 / (Rk * 3.085677581e19)

resid = np.log10(gobs) - np.log10(gbar)

# Pearson r
chi0 = chi - np.mean(chi)
res0 = resid - np.mean(resid)
r = float(np.sum(chi0 * res0) / np.sqrt(np.sum(chi0 * chi0) * np.sum(res0 * res0)))

print("N points:", len(chi))
print("Pearson r (chi vs H1 RAR residual proxy):", r)
print("chi range:", float(np.min(chi)), float(np.max(chi)))
print("resid range:", float(np.min(resid)), float(np.max(resid)))
