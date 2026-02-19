import pandas as pd
import numpy as np

g = "NGC3198"

p3 = pd.read_csv(f"data/derived/phase3/leff_{g}.csv")
p4 = pd.read_csv(f"data/derived/phase4/h2_outputs/rc_decomp_{g}_H2_adaptive.csv")

# Find chi column
chi_col = None
for c in ["chi_used", "chi", "chi_raw", "chi_smooth"]:
    if c in p3.columns:
        chi_col = c
        break

if chi_col is None:
    raise RuntimeError(
        f"No chi column found in phase3 leff file. Columns: {list(p3.columns)}"
    )

# Align on R
p3 = p3.sort_values("R_kpc")
p4 = p4.sort_values("R_kpc")

if len(p3) != len(p4) or not np.allclose(
    p3.R_kpc.values, p4.R_kpc.values, atol=1e-10, rtol=0
):
    raise RuntimeError("R grids mismatch between phase3 and phase4")

chi = p3[chi_col].values
dV = p4.V_total_H2.values - p4.V_total_H1.values

# Use inner+mid region where adaptation acts
R = p4.R_kpc.values
rmax = R.max()
mask = (R / rmax) < 0.7

chi = chi[mask]
dV = dV[mask]

# Pearson r
chi = chi - np.mean(chi)
dV = dV - np.mean(dV)

r = float(np.sum(chi * dV) / np.sqrt(np.sum(chi * chi) * np.sum(dV * dV)))

print("chi column:", chi_col)
print("Pearson r:", r)
print("N points:", len(chi))
