import argparse
import pandas as pd
import numpy as np
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Test-2: Chi-ΔV Correlation")
    parser.add_argument("--galaxy", required=True, help="Galaxy name")
    args = parser.parse_args()
    
    galaxy = args.galaxy
    
    p3 = pd.read_csv(f"data/derived/phase3/leff_{galaxy}.csv")
    p4 = pd.read_csv(f"data/derived/phase4/h2_outputs/rc_decomp_{galaxy}_H2_adaptive.csv")
    
    chi_col = 'chi_used' if 'chi_used' in p3.columns else 'chi_raw'
    
    R_p3 = p3['R_kpc'].values
    chi = p3[chi_col].values
    
    R_p4 = p4['R_kpc'].values
    dV = p4['V_total_H2'].values - p4['V_total_H1'].values
    
    if not np.allclose(R_p3, R_p4, atol=0.01):
        dV = np.interp(R_p3, R_p4, dV)
        R = R_p3
    else:
        R = R_p3
    
    mask = (R / R.max()) < 0.7
    chi_masked = chi[mask] - chi[mask].mean()
    dV_masked = dV[mask] - dV[mask].mean()
    
    r = np.sum(chi_masked * dV_masked) / np.sqrt(np.sum(chi_masked**2) * np.sum(dV_masked**2))
    
    print(f"Galaxy: {galaxy}")
    print(f"chi column: {chi_col}")
    print(f"Pearson r: {r:.4f}")
    print(f"N points: {mask.sum()}")
    print(f"max|ΔV| (masked): {np.max(np.abs(dV[mask])):.4f} km/s")
    
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())