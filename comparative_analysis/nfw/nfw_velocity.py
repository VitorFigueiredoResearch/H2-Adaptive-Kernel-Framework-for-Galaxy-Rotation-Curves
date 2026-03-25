"""
NFW circular velocity module.
Li et al. (2020) convention: V200 in km/s, C200 dimensionless, rs in kpc.
H0 = 73 km/s/Mpc -> r200[kpc] = V200 / 0.73
"""
import numpy as np

H0_KMS_MPC = 73.0  # km/s/Mpc


def r200_kpc(V200_kms):
    """r200 in kpc from V200 in km/s. Derived from M200 = (4/3)pi 200 rho_crit r200^3.

    At r200: V200^2 / r200 = G M200 / r200^2 = (4pi/3) G 200 rho_crit r200
    -> V200^2 = (4pi/3) G 200 rho_crit r200^2
    For H0=73 km/s/Mpc: r200[kpc] = V200[km/s] / 0.73
    """
    return V200_kms / (10.0 * H0_KMS_MPC / 1000.0)  # /1000 converts Mpc->kpc units


def nfw_vcirc(r_kpc, V200_kms, C200, rs_kpc):
    """
    NFW circular velocity at radius r_kpc.
    V^2_NFW(r) = (V200^2 * r200 / r) * f(r/rs) / f(C200)
    where f(x) = ln(1+x) - x/(1+x)
    Uses rs_kpc directly from catalog.
    Returns V_NFW in km/s.

    Parameters
    ----------
    r_kpc : array-like
        Radii in kpc
    V200_kms : float
        Circular velocity at r200 in km/s
    C200 : float
        NFW concentration parameter (r200/rs)
    rs_kpc : float
        NFW scale radius in kpc (from catalog directly)

    Returns
    -------
    V_NFW : ndarray
        Circular velocity in km/s at each radius
    """
    r_kpc = np.asarray(r_kpc, dtype=float)
    r200 = r200_kpc(V200_kms)
    x = r_kpc / rs_kpc
    c = C200

    def f(t):
        return np.log(1.0 + t) - t / (1.0 + t)

    fc = f(c)
    if fc <= 0:
        return np.zeros_like(r_kpc)

    fx = f(x)
    # V^2(r) = V200^2 * (r200/r) * f(x)/f(c)
    V2 = V200_kms**2 * (r200 / r_kpc) * fx / fc
    V2 = np.where(V2 > 0, V2, 0.0)
    return np.sqrt(V2)


def load_nfw_catalog(catalog_path):
    """
    Load the Li et al. (2020) NFW catalog in MRT fixed-width format.

    Returns
    -------
    catalog : dict
        Keys are galaxy names, values are dicts with V200, e_V200, C200, e_C200,
        rs, e_rs, log_rhos, e_log_rhos, log_M200, e_log_M200, chi2, and other params.
    """
    catalog = {}
    with open(catalog_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.rstrip('\n')
            if len(line) < 14:
                continue
            name = line[:14].strip()
            if not name:
                continue
            try:
                vals = line[14:].split()
                if len(vals) < 18:
                    continue
                chi_str = vals[18] if len(vals) > 18 else 'nan'
                chi_val = float('inf') if chi_str.lower() == 'inf' else float(chi_str)
                catalog[name] = {
                    'Ydisk':      float(vals[0]),
                    'e_Ydisk':    float(vals[1]),
                    'Ybul':       float(vals[2]),
                    'e_Ybul':     float(vals[3]),
                    'D':          float(vals[4]),
                    'e_D':        float(vals[5]),
                    'inc':        float(vals[6]),
                    'e_inc':      float(vals[7]),
                    'V200':       float(vals[8]),
                    'e_V200':     float(vals[9]),
                    'C200':       float(vals[10]),
                    'e_C200':     float(vals[11]),
                    'rs':         float(vals[12]),
                    'e_rs':       float(vals[13]),
                    'log_rhos':   float(vals[14]),
                    'e_log_rhos': float(vals[15]),
                    'log_M200':   float(vals[16]),
                    'e_log_M200': float(vals[17]),
                    'chi2':       chi_val,
                }
            except (ValueError, IndexError):
                continue
    return catalog


def validate_ngc3198(catalog_path=None):
    """Quick validation for NGC3198 from Li+2020 NFW_LCDM catalog."""
    import os
    if catalog_path is None:
        catalog_path = os.path.join(os.path.dirname(__file__),
                                    '../../data/nfw/Fits/ByModel/Table/parameter_NFW_LCDM.mrt')

    catalog = load_nfw_catalog(catalog_path)
    entry = catalog.get('NGC3198')
    if entry is None:
        print("NGC3198 not found in catalog")
        return

    V200 = entry['V200']
    C200 = entry['C200']
    rs   = entry['rs']
    r200 = r200_kpc(V200)

    print("NGC3198 validation:")
    print(f"  V200={V200:.2f} km/s, C200={C200:.2f}, rs={rs:.2f} kpc")
    print(f"  r200={r200:.2f} kpc (catalog rs*C200={(rs*C200):.2f} kpc)")
    print(f"  rs check: r200/C200 = {r200/C200:.2f} vs catalog rs={rs:.2f} kpc")

    r_test = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 30.0])
    v_test = nfw_vcirc(r_test, V200, C200, rs)
    print("  Radial profile:")
    for r, v in zip(r_test, v_test):
        print(f"    V_NFW({r:5.1f} kpc) = {v:.1f} km/s")
    print("  (Expected: rising from ~50 km/s at 1 kpc to ~130-145 km/s peak, then flat)")


def validate_d564_8(catalog_path=None):
    """Validation for D564-8 (pre-confirmed in data_inspection_report.txt)."""
    import os
    if catalog_path is None:
        catalog_path = os.path.join(os.path.dirname(__file__),
                                    '../../data/nfw/Fits/ByModel/Table/parameter_NFW_LCDM.mrt')

    catalog = load_nfw_catalog(catalog_path)
    entry = catalog.get('D564-8')
    if entry is None:
        print("D564-8 not found in catalog")
        return

    V200 = entry['V200']
    C200 = entry['C200']
    rs   = entry['rs']
    r200 = r200_kpc(V200)
    rs_calc = r200 / C200

    print("D564-8 verification:")
    print(f"  V200={V200:.2f} km/s -> r200={r200:.2f} kpc")
    print(f"  C200={C200:.2f} -> rs_calc={rs_calc:.2f} kpc")
    print(f"  Catalog rs={rs:.2f} kpc")
    print(f"  Residual: {abs(rs_calc - rs):.3f} kpc", "PASS" if abs(rs_calc - rs) < 0.05 else "FAIL")


if __name__ == '__main__':
    import os
    catalog_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '../../data/nfw/Fits/ByModel/Table/parameter_NFW_LCDM.mrt')
    validate_d564_8(catalog_path)
    print()
    validate_ngc3198(catalog_path)
