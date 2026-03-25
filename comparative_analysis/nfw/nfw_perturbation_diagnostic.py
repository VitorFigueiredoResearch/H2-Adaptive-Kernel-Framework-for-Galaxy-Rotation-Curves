"""
NFW Perturbation Diagnostic
=============================
Tests scatter sensitivity to ±10%/±20% perturbations in NFW best-fit parameters
(V200, C200) for Li et al. (2020) fits to H2 fleet galaxies.

Scatter metric: RMS of log10(V_tot) - log10(V_obs) over inner region (R < 0.5*R_max).
Output: nfw_perturbation_summary.csv and diagnostic_report.txt
Figures: figures/ngc3198_validation.png, figures/nfw_scatter_sensitivity.png
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Ensure we can import nfw_velocity from the same directory
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)
from nfw_velocity import r200_kpc, nfw_vcirc, load_nfw_catalog

# -----------------------------------------------------------------------
# Paths (relative to repo root; this script is run from repo root)
# -----------------------------------------------------------------------
REPO_ROOT = os.path.dirname(os.path.dirname(THIS_DIR))
NFW_CATALOG_PATH  = os.path.join(REPO_ROOT, 'data', 'nfw', 'Fits', 'ByModel', 'Table',
                                  'parameter_NFW_LCDM.mrt')
FLEET_CSV_PATH    = os.path.join(REPO_ROOT, 'H2_PUBLICATION_RELEASE',
                                  'fleet_expansion_80galaxy', 'fleet_summary_80galaxy.csv')
SPARC_DIR         = os.path.join(REPO_ROOT, 'data', 'sparc')
OUT_DIR           = THIS_DIR
FIGURES_DIR       = os.path.join(THIS_DIR, 'figures')

os.makedirs(FIGURES_DIR, exist_ok=True)


# -----------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------

def norm_name(s):
    """Normalize galaxy name for matching: uppercase, strip spaces/underscores."""
    return s.strip().upper().replace('_', '').replace(' ', '')


def load_sparc_rotmod(galaxy_name):
    """
    Load SPARC rotation curve data for a galaxy.
    Returns (R, Vobs, errV, Vgas, Vdisk, Vbul) as numpy arrays.
    Returns None if file not found or parsing fails.
    """
    fname = os.path.join(SPARC_DIR, f'{galaxy_name}_rotmod.dat')
    if not os.path.exists(fname):
        return None
    data = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('!'):
                continue
            try:
                vals = [float(x) for x in line.split()]
                if len(vals) >= 6:
                    data.append(vals[:6])
            except (ValueError, TypeError):
                continue
    if len(data) == 0:
        return None
    data = np.array(data)
    R, Vobs, errV, Vgas, Vdisk, Vbul = (data[:, i] for i in range(6))
    # Only keep rows with Vobs > 0 and R > 0
    mask = (Vobs > 0) & (R > 0)
    if mask.sum() < 3:
        return None
    return R[mask], Vobs[mask], errV[mask], Vgas[mask], Vdisk[mask], Vbul[mask]


def compute_vbar(Vgas, Vdisk, Vbul):
    """Baryonic circular velocity: sqrt(Vgas^2 + Vdisk^2 + Vbul^2)."""
    return np.sqrt(Vgas**2 + Vdisk**2 + Vbul**2)


def compute_vtot(Vbar, V_NFW):
    """Total circular velocity: sqrt(Vbar^2 + V_NFW^2)."""
    return np.sqrt(Vbar**2 + V_NFW**2)


def inner_scatter(Vtot, Vobs, inner_mask):
    """
    RMS of log10(Vtot) - log10(Vobs) over inner region points.
    Returns NaN if no valid points.
    """
    if inner_mask.sum() < 1:
        return np.nan
    Vt = Vtot[inner_mask]
    Vo = Vobs[inner_mask]
    valid = (Vt > 0) & (Vo > 0)
    if valid.sum() < 1:
        return np.nan
    residuals = np.log10(Vt[valid]) - np.log10(Vo[valid])
    return np.sqrt(np.mean(residuals**2))


def process_galaxy(galaxy_name, regime, nfw_entry, rotmod_data):
    """
    Compute baseline scatter and perturbation delta-sigmas for one galaxy.

    Returns a dict with all output columns, or None if analysis fails.
    """
    R, Vobs, errV, Vgas, Vdisk, Vbul = rotmod_data

    V200 = nfw_entry['V200']
    C200 = nfw_entry['C200']
    rs   = nfw_entry['rs']

    # Physical parameter checks
    if V200 <= 0 or C200 <= 0 or rs <= 0:
        return None

    Vbar = compute_vbar(Vgas, Vdisk, Vbul)

    # Baseline NFW velocity
    V_NFW_base = nfw_vcirc(R, V200, C200, rs)
    Vtot_base  = compute_vtot(Vbar, V_NFW_base)

    # Inner region definition: R < 0.5 * R_max
    R_max = R.max()
    R_cut = 0.5 * R_max
    inner = R < R_cut
    n_inner = int(inner.sum())
    inner_defined = n_inner >= 3

    sigma_base = inner_scatter(Vtot_base, Vobs, inner)

    # Median Vbar/V_NFW in inner region
    if inner.sum() > 0 and np.any(V_NFW_base[inner] > 0):
        ratio_inner = Vbar[inner] / np.where(V_NFW_base[inner] > 0, V_NFW_base[inner], np.nan)
        vbar_over_vnfw = float(np.nanmedian(ratio_inner))
    else:
        vbar_over_vnfw = np.nan

    # ---- Perturbations ----
    # Perturbation scales
    scales = {10: 0.10, 20: 0.20}

    delta_sigmas = {}
    for scale_pct, frac in scales.items():
        for sign, label_sign in [(+1, 'p'), (-1, 'm')]:
            # Mass perturbation: change V200 (V200 ~ M^1/3 r^-1/2 but we perturb V200 directly)
            V200_p = V200 * (1.0 + sign * frac)
            V_NFW_p = nfw_vcirc(R, V200_p, C200, rs)
            Vtot_p  = compute_vtot(Vbar, V_NFW_p)
            sig_p   = inner_scatter(Vtot_p, Vobs, inner)
            delta_sigmas[f'delta_sigma_mass_{label_sign}{scale_pct}'] = (sig_p - sigma_base
                                                                          if not np.isnan(sig_p) else np.nan)

            # Concentration perturbation: change C200, keep V200
            C200_p = C200 * (1.0 + sign * frac)
            V_NFW_p = nfw_vcirc(R, V200, C200_p, rs)
            Vtot_p  = compute_vtot(Vbar, V_NFW_p)
            sig_p   = inner_scatter(Vtot_p, Vobs, inner)
            delta_sigmas[f'delta_sigma_c_{label_sign}{scale_pct}'] = (sig_p - sigma_base
                                                                        if not np.isnan(sig_p) else np.nan)

            # Joint perturbation: change both V200 and C200
            V200_p = V200 * (1.0 + sign * frac)
            C200_p = C200 * (1.0 + sign * frac)
            V_NFW_p = nfw_vcirc(R, V200_p, C200_p, rs)
            Vtot_p  = compute_vtot(Vbar, V_NFW_p)
            sig_p   = inner_scatter(Vtot_p, Vobs, inner)
            delta_sigmas[f'delta_sigma_both_{label_sign}{scale_pct}'] = (sig_p - sigma_base
                                                                           if not np.isnan(sig_p) else np.nan)

    # max |delta_sigma| across all perturbations
    ds_vals = [v for v in delta_sigmas.values() if not np.isnan(v)]
    max_abs_ds = float(np.max(np.abs(ds_vals))) if ds_vals else np.nan

    return {
        'Galaxy': galaxy_name,
        'Regime': regime,
        'V200_kms': V200,
        'C200': C200,
        'rs_kpc': rs,
        'sigma_baseline': sigma_base,
        'n_inner_points': n_inner,
        'inner_region_defined': inner_defined,
        **delta_sigmas,
        'max_abs_delta_sigma': max_abs_ds,
        'V_bar_over_V_NFW_inner_median': vbar_over_vnfw,
    }


# -----------------------------------------------------------------------
# Validation plot for NGC3198
# -----------------------------------------------------------------------

def make_ngc3198_validation_plot(catalog, rotmod_dir):
    """
    Plot NGC3198 validation: Vobs, Vbar, baseline V_NFW, ±10%/±20% perturbed V_NFW curves.
    Mark inner region boundary.
    """
    galaxy = 'NGC3198'
    entry = catalog.get(galaxy)
    if entry is None:
        print(f"Warning: {galaxy} not in catalog, skipping validation plot.")
        return

    rotmod_data = load_sparc_rotmod(galaxy)
    if rotmod_data is None:
        print(f"Warning: rotmod data for {galaxy} not found, skipping validation plot.")
        return

    R, Vobs, errV, Vgas, Vdisk, Vbul = rotmod_data
    Vbar = compute_vbar(Vgas, Vdisk, Vbul)

    V200 = entry['V200']
    C200 = entry['C200']
    rs   = entry['rs']

    R_max = R.max()
    R_cut = 0.5 * R_max

    # Fine radial grid for smooth curves
    r_grid = np.linspace(0.01, R_max * 1.05, 500)

    V_NFW_base = nfw_vcirc(r_grid, V200, C200, rs)

    fig, ax = plt.subplots(figsize=(8, 5))

    # Perturbed curves
    for frac, ls in [(0.10, '--'), (0.20, ':')]:
        pct = int(frac * 100)
        V_nfw_pp = nfw_vcirc(r_grid, V200 * (1 + frac), C200, rs)
        V_nfw_pm = nfw_vcirc(r_grid, V200 * (1 - frac), C200, rs)
        ax.plot(r_grid, V_nfw_pp, ls=ls, color='steelblue', alpha=0.7,
                label=f'V_NFW mass +{pct}%')
        ax.plot(r_grid, V_nfw_pm, ls=ls, color='orange', alpha=0.7,
                label=f'V_NFW mass -{pct}%')

    # Baseline NFW
    ax.plot(r_grid, V_NFW_base, '-', color='navy', lw=2, label='V_NFW baseline')

    # Vbar at data points
    ax.plot(R, Vbar, 's', color='green', ms=4, label='V_bar (data points)')

    # Vobs with error bars
    ax.errorbar(R, Vobs, yerr=errV, fmt='o', color='black', ms=4, lw=1,
                label='V_obs (observed)')

    # Inner region boundary
    ax.axvline(R_cut, color='red', ls='-', lw=1.5, alpha=0.8,
               label=f'Inner region cut (R_cut={R_cut:.1f} kpc)')

    # Shade inner region
    ax.axvspan(0, R_cut, alpha=0.06, color='red', label='Inner region')

    ax.set_xlabel('R [kpc]', fontsize=12)
    ax.set_ylabel('V [km/s]', fontsize=12)
    ax.set_title(f'NGC3198 NFW validation\nV200={V200:.1f} km/s, C200={C200:.2f}, rs={rs:.2f} kpc',
                 fontsize=11)
    ax.legend(fontsize=7, loc='lower right', ncol=2)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)

    out_path = os.path.join(FIGURES_DIR, 'ngc3198_validation.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_path}")


# -----------------------------------------------------------------------
# Scatter sensitivity figure
# -----------------------------------------------------------------------

def make_scatter_sensitivity_figure(df_results):
    """
    Top panel: max|Dsigma| vs median(Vbar/V_NFW) inner, color by regime
    Bottom panel: histogram of max|Dsigma|
    """
    df = df_results[df_results['inner_region_defined']].copy()
    df = df.dropna(subset=['max_abs_delta_sigma', 'V_bar_over_V_NFW_inner_median'])

    regime_colors = {
        'baryon-dom': 'red',
        'balanced':   'blue',
        'DM-dom':     'green',
    }
    default_color = 'gray'

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 9))
    fig.suptitle('NFW Parameter Perturbation Sensitivity\nInner Region Scatter Response',
                 fontsize=13)

    # Top panel
    for regime, color in regime_colors.items():
        sub = df[df['Regime'] == regime]
        ax1.scatter(sub['V_bar_over_V_NFW_inner_median'], sub['max_abs_delta_sigma'],
                    color=color, alpha=0.7, s=50, label=regime, zorder=3)

    # Any unclassified
    known_regimes = set(regime_colors.keys())
    other = df[~df['Regime'].isin(known_regimes)]
    if len(other) > 0:
        ax1.scatter(other['V_bar_over_V_NFW_inner_median'], other['max_abs_delta_sigma'],
                    color=default_color, alpha=0.7, s=50, label='other', zorder=3)

    ax1.set_xlabel('Median V_bar / V_NFW (inner region)', fontsize=11)
    ax1.set_ylabel('max |Δσ| [dex]', fontsize=11)
    ax1.set_title('Scatter Sensitivity vs Baryon Dominance', fontsize=10)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Annotation
    n_analyzed = len(df)
    ax1.text(0.98, 0.97, f'N = {n_analyzed}', transform=ax1.transAxes,
             ha='right', va='top', fontsize=10,
             bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='gray', alpha=0.8))

    # Bottom panel: histogram
    ax2.hist(df['max_abs_delta_sigma'], bins=20, color='steelblue', edgecolor='white',
             alpha=0.85)
    ax2.set_xlabel('max |Δσ| [dex]', fontsize=11)
    ax2.set_ylabel('Count', fontsize=11)
    ax2.set_title('Distribution of Maximum Scatter Change\n(±10%/±20% NFW parameter perturbations)',
                  fontsize=10)
    ax2.grid(True, alpha=0.3, axis='y')

    # Threshold lines
    for threshold, label, color in [(0.01, '0.01 dex', 'green'), (0.05, '0.05 dex', 'orange')]:
        ax2.axvline(threshold, color=color, ls='--', lw=1.5, label=label)
    ax2.legend(fontsize=9)

    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'nfw_scatter_sensitivity.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_path}")


# -----------------------------------------------------------------------
# Main analysis
# -----------------------------------------------------------------------

def main():
    print("=" * 60)
    print("NFW Perturbation Diagnostic")
    print("=" * 60)

    # Load NFW catalog
    print(f"\nLoading NFW catalog: {NFW_CATALOG_PATH}")
    catalog = load_nfw_catalog(NFW_CATALOG_PATH)
    print(f"  NFW catalog entries: {len(catalog)}")

    # Load H2 fleet
    print(f"\nLoading H2 fleet: {FLEET_CSV_PATH}")
    fleet_df = pd.read_csv(FLEET_CSV_PATH)
    fleet_df.columns = [c.strip() for c in fleet_df.columns]
    # Build regime lookup
    regime_lookup = {}
    for _, row in fleet_df.iterrows():
        gname = str(row['galaxy']).strip()
        regime = str(row['regime']).strip() if 'regime' in row else 'unknown'
        regime_lookup[norm_name(gname)] = (gname, regime)
    print(f"  H2 fleet galaxies: {len(fleet_df)}")

    # Build NFW name lookup (normalized)
    nfw_norm = {norm_name(k): k for k in catalog.keys()}

    # Diagnostic tracking
    not_in_nfw = []
    no_rotmod = []
    low_inner_pts = []
    bad_params = []
    excluded = []
    results = []

    print("\nProcessing galaxies...")
    for norm_g, (gname, regime) in sorted(regime_lookup.items()):
        # Match to NFW catalog
        nfw_key = nfw_norm.get(norm_g)
        if nfw_key is None:
            not_in_nfw.append(gname)
            excluded.append(gname)
            continue

        nfw_entry = catalog[nfw_key]

        # Physical parameter check
        if nfw_entry['V200'] <= 0 or nfw_entry['C200'] <= 0 or nfw_entry['rs'] <= 0:
            bad_params.append(gname)
            excluded.append(gname)
            continue

        # Load SPARC rotmod
        rotmod_data = load_sparc_rotmod(gname)
        if rotmod_data is None:
            no_rotmod.append(gname)
            excluded.append(gname)
            continue

        # Process
        result = process_galaxy(gname, regime, nfw_entry, rotmod_data)
        if result is None:
            excluded.append(gname)
            continue

        if result['n_inner_points'] < 3:
            low_inner_pts.append(gname)
            # Still include in CSV but mark inner_region_defined=False

        results.append(result)

    print(f"  Successfully analyzed: {len(results)}")
    print(f"  Excluded: {len(excluded)}")

    # Save summary CSV
    if results:
        df_out = pd.DataFrame(results)
        # Ensure column order
        col_order = [
            'Galaxy', 'Regime', 'V200_kms', 'C200', 'rs_kpc', 'sigma_baseline',
            'n_inner_points', 'inner_region_defined',
            'delta_sigma_mass_p10', 'delta_sigma_mass_p20',
            'delta_sigma_mass_m10', 'delta_sigma_mass_m20',
            'delta_sigma_c_p10',    'delta_sigma_c_p20',
            'delta_sigma_c_m10',    'delta_sigma_c_m20',
            'delta_sigma_both_p10', 'delta_sigma_both_p20',
            'delta_sigma_both_m10', 'delta_sigma_both_m20',
            'max_abs_delta_sigma', 'V_bar_over_V_NFW_inner_median',
        ]
        # Add any missing columns as NaN
        for col in col_order:
            if col not in df_out.columns:
                df_out[col] = np.nan
        df_out = df_out[col_order]

        csv_path = os.path.join(OUT_DIR, 'nfw_perturbation_summary.csv')
        df_out.to_csv(csv_path, index=False, float_format='%.6f')
        print(f"\nSaved: {csv_path}")
        print("\nFirst 10 rows:")
        print(df_out.head(10).to_string())
    else:
        df_out = pd.DataFrame()
        print("No results to save.")

    # Write diagnostic report
    report_path = os.path.join(OUT_DIR, 'diagnostic_report.txt')
    with open(report_path, 'w') as rpt:
        rpt.write("NFW Perturbation Diagnostic Report\n")
        rpt.write("=" * 50 + "\n\n")

        rpt.write(f"H2 fleet galaxies:          {len(fleet_df)}\n")
        rpt.write(f"NFW catalog entries:         {len(catalog)}\n")
        rpt.write(f"Successfully analyzed:       {len(results)}\n")
        rpt.write(f"Total excluded:              {len(excluded)}\n\n")

        rpt.write("Exclusion reasons:\n")
        rpt.write(f"  Not in NFW catalog:        {len(not_in_nfw)}\n")
        rpt.write(f"  No SPARC rotmod file:      {len(no_rotmod)}\n")
        rpt.write(f"  Non-physical params:       {len(bad_params)}\n\n")

        if not_in_nfw:
            rpt.write("H2 galaxies NOT in NFW catalog:\n")
            for g in sorted(not_in_nfw):
                rpt.write(f"  {g}\n")
            rpt.write("\n")

        if no_rotmod:
            rpt.write("H2 galaxies with missing rotmod file:\n")
            for g in sorted(no_rotmod):
                rpt.write(f"  {g}\n")
            rpt.write("\n")

        if bad_params:
            rpt.write("Galaxies with non-physical NFW parameters (V200<=0 or C200<=0 or rs<=0):\n")
            for g in sorted(bad_params):
                rpt.write(f"  {g}\n")
            rpt.write("\n")

        if low_inner_pts:
            rpt.write("Galaxies with n_inner < 3 (flagged, excluded from sensitivity figures):\n")
            for g in sorted(low_inner_pts):
                rpt.write(f"  {g}\n")
            rpt.write("\n")

        if len(results) > 0:
            df_valid = df_out[df_out['inner_region_defined']].dropna(subset=['max_abs_delta_sigma'])
            n_valid = len(df_valid)
            rpt.write(f"Galaxies with n_inner >= 3 (used in sensitivity figures): {n_valid}\n\n")

            if n_valid > 0:
                rpt.write("Summary statistics for max|Delta_sigma| [dex]:\n")
                rpt.write(f"  Mean:   {df_valid['max_abs_delta_sigma'].mean():.4f}\n")
                rpt.write(f"  Median: {df_valid['max_abs_delta_sigma'].median():.4f}\n")
                rpt.write(f"  Std:    {df_valid['max_abs_delta_sigma'].std():.4f}\n")
                rpt.write(f"  Min:    {df_valid['max_abs_delta_sigma'].min():.4f}\n")
                rpt.write(f"  Max:    {df_valid['max_abs_delta_sigma'].max():.4f}\n\n")

                for regime in ['baryon-dom', 'balanced', 'DM-dom']:
                    sub = df_valid[df_valid['Regime'] == regime]
                    if len(sub) > 0:
                        rpt.write(f"  Regime '{regime}' (N={len(sub)}): "
                                  f"mean={sub['max_abs_delta_sigma'].mean():.4f}, "
                                  f"median={sub['max_abs_delta_sigma'].median():.4f}\n")
                rpt.write("\n")

    print(f"Saved: {report_path}")

    # Make validation plot for NGC3198
    print("\nGenerating NGC3198 validation plot...")
    make_ngc3198_validation_plot(catalog, SPARC_DIR)

    # Make scatter sensitivity figure
    if len(results) > 0:
        print("Generating scatter sensitivity figure...")
        make_scatter_sensitivity_figure(df_out)

    print("\nDone.")
    return df_out


if __name__ == '__main__':
    df = main()
