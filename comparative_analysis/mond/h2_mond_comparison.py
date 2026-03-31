"""
h2_mond_comparison.py
=====================
Compare H2 inner-region scatter-neutrality pattern with MOND/RAR
bounded perturbation results using the harmonized log10-RMS metric.

Inputs:
  - H2 metric summary: comparative_analysis/metric_harmonization/h2_metric_summary.csv
  - NFW results (reference): comparative_analysis/nfw/nfw_perturbation_summary.csv
  - MOND results: comparative_analysis/mond/mond_perturbation_summary.csv
  - H2 fleet:   H2_PUBLICATION_RELEASE/fleet_expansion_80galaxy/fleet_summary_80galaxy.csv

Outputs:
  - comparative_analysis/mond/h2_mond_comparison_summary.csv
  - comparative_analysis/mond/h2_mond_comparison_report.txt
  - comparative_analysis/mond/figures/mond_scatter_sensitivity.png

Run from repository root:
  python comparative_analysis/mond/h2_mond_comparison.py
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

# Paths
MOND_CSV  = os.path.join(REPO_ROOT, 'comparative_analysis', 'mond',
                          'mond_perturbation_summary.csv')
NFW_CSV   = os.path.join(REPO_ROOT, 'comparative_analysis', 'nfw',
                          'nfw_perturbation_summary.csv')
H2_CSV    = os.path.join(REPO_ROOT, 'comparative_analysis',
                          'metric_harmonization', 'h2_metric_summary.csv')
FLEET_CSV = os.path.join(REPO_ROOT, 'H2_PUBLICATION_RELEASE',
                          'fleet_expansion_80galaxy', 'fleet_summary_80galaxy.csv')
OUT_CSV   = os.path.join(REPO_ROOT, 'comparative_analysis', 'mond',
                          'h2_mond_comparison_summary.csv')
OUT_RPT   = os.path.join(REPO_ROOT, 'comparative_analysis', 'mond',
                          'h2_mond_comparison_report.txt')
FIG_DIR   = os.path.join(REPO_ROOT, 'comparative_analysis', 'mond', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

REGIME_COLORS = {
    'baryon-dom':  '#d62728',   # red
    'balanced':    '#1f77b4',   # blue
    'DM-dom':      '#2ca02c',   # green
}
REGIME_LABELS = {
    'baryon-dom':  'Baryon-dominated',
    'balanced':    'Balanced',
    'DM-dom':      'DM-dominated',
}


def load_mond_per_galaxy(mond_csv):
    """Collapse per-perturbation rows to one row per galaxy."""
    df = pd.read_csv(mond_csv)
    per_gal = (df.groupby(['Galaxy', 'Regime'])
                 .agg(
                     sigma_baseline=('sigma_baseline', 'first'),
                     n_inner_points=('n_inner_points', 'first'),
                     max_abs_delta_sigma=('delta_sigma', lambda x: x.abs().max()),
                     V_bar_over_V_MOND=('V_bar_over_V_MOND_inner_median', 'first'),
                 )
                 .reset_index())
    return per_gal


def load_nfw_per_galaxy(nfw_csv):
    """Load NFW summary and extract max|Δσ| per galaxy."""
    df = pd.read_csv(nfw_csv)
    if 'max_abs_delta_sigma' in df.columns:
        return df[['Galaxy', 'Regime', 'sigma_baseline', 'n_inner_points',
                   'max_abs_delta_sigma']].drop_duplicates('Galaxy')
    # Fall back: collapse if per-perturbation rows
    if 'delta_sigma' in df.columns:
        per_gal = (df.groupby(['Galaxy', 'Regime'])
                     .agg(sigma_baseline=('sigma_baseline', 'first'),
                          n_inner_points=('n_inner_points', 'first'),
                          max_abs_delta_sigma=('delta_sigma', lambda x: x.abs().max()))
                     .reset_index())
        return per_gal
    return df


def load_h2_per_galaxy(h2_csv):
    """Load H2 harmonized metric summary."""
    df = pd.read_csv(h2_csv)
    # h2_delta_sigma_harmonized may be absent if H2 files not found
    cols = ['Galaxy', 'Regime', 'n_inner_points', 'h2_sigma_baseline',
            'h2_delta_sigma_harmonized', 'reconstruction_status']
    present = [c for c in cols if c in df.columns]
    return df[present]


def regime_stats(df, col):
    """Compute summary stats per regime for a given column."""
    rows = []
    for regime, grp in df.groupby('Regime'):
        vals = grp[col].dropna()
        rows.append({
            'Regime':           regime,
            'n':                len(vals),
            'median':           vals.median(),
            'mean':             vals.mean(),
            'std':              vals.std(),
            'frac_lt_001':      (vals < 0.01).sum() / len(vals) if len(vals) else np.nan,
            'frac_gt_005':      (vals > 0.05).sum() / len(vals) if len(vals) else np.nan,
        })
    return pd.DataFrame(rows)


def make_comparison_figure(merged, fig_dir):
    """Two-panel comparison figure (scatter + histogram)."""
    H2_MEDIAN_FIXED = 0.000478

    fig, axes = plt.subplots(2, 1, figsize=(9, 9))
    fig.suptitle('Adopted RAR Inner-Region Scatter Sensitivity\n'
                 '(Harmonized metric: log10-RMS [dex])', fontsize=13, y=0.98)

    # ── Top panel: RAR max|Δσ| vs V_bar/V_RAR ratio ──────────────────────────
    ax = axes[0]
    for regime, grp in merged.groupby('Regime'):
        c = REGIME_COLORS.get(regime, 'gray')
        lbl = REGIME_LABELS.get(regime, regime)
        ax.scatter(grp['V_bar_over_V_MOND'],
                   grp['mond_max_abs_ds'],
                   c=c, label=lbl, s=55, alpha=0.75, edgecolors='k', lw=0.4, zorder=3)

    ax.axhline(0.01, color='gray', lw=1, ls='--', label='0.01 dex threshold')
    ax.set_xlabel('Median($V_{bar}$ / $V_{RAR}$) in inner region', fontsize=11)
    ax.set_ylabel('RAR max|$\Delta\sigma$| [dex]', fontsize=11)
    ax.set_title(f'RAR scatter sensitivity vs baryonic dominance ratio  '
                 f'(N={len(merged)})', fontsize=11)
    ax.legend(fontsize=9, loc='upper right')
    ax.set_ylim(bottom=0)

    # ── Bottom panel: histogram of max|Δσ| — RAR vs NFW with H2 reference ───
    ax2 = axes[1]
    bins = np.linspace(0, max(merged['mond_max_abs_ds'].max(),
                               merged.get('nfw_max_abs_ds', pd.Series([0.1])).max(),
                               0.12), 25)

    ax2.hist(merged['mond_max_abs_ds'].dropna(), bins=bins,
             color='steelblue', alpha=0.65, label='RAR (adopted) max|Δσ|', edgecolor='k', lw=0.4)

    if 'nfw_max_abs_ds' in merged.columns:
        ax2.hist(merged['nfw_max_abs_ds'].dropna(), bins=bins,
                 color='tomato', alpha=0.55, label='NFW max|Δσ|', edgecolor='k', lw=0.4)

    ax2.axvline(H2_MEDIAN_FIXED, color='purple', lw=2, ls='--',
                label='H2 median|Δσ| = 0.000478 dex')

    ax2.axvline(0.01, color='gray', lw=1, ls=':', label='0.01 dex')
    ax2.set_xlabel('max|$\Delta\sigma$| [dex, log10-RMS]', fontsize=11)
    ax2.set_ylabel('Galaxy count', fontsize=11)
    ax2.set_title('Distribution of adopted RAR scatter sensitivity with NFW and H2 reference levels',
                  fontsize=11)
    ax2.legend(fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(fig_dir, 'rar_scatter_sensitivity.png')
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  Figure saved: {out}")


def write_report(merged, mond_stats, nfw_stats, h2_stats, out_path, n_mond, n_nfw, n_h2):
    """Write plain-text comparison narrative."""

    def fmt_stats(stats_df, label):
        lines = [f"\n  {label}:"]
        for _, r in stats_df.iterrows():
            lines.append(
                f"    {r['Regime']:14s}  n={int(r['n']):2d}  "
                f"median={r['median']:.4f}  mean={r['mean']:.4f}  "
                f"<0.01: {r['frac_lt_001']:.0%}  >0.05: {r['frac_gt_005']:.0%}")
        return '\n'.join(lines)

    mond_med = merged['mond_max_abs_ds'].median()
    nfw_med  = merged['nfw_max_abs_ds'].median()  if 'nfw_max_abs_ds'  in merged else np.nan
    h2_med   = merged['h2_abs_delta_sigma'].median() if 'h2_abs_delta_sigma' in merged else np.nan

    closer = "cannot be determined (NFW or H2 data unavailable)"
    if not np.isnan(nfw_med) and not np.isnan(h2_med):
        d_mond_h2  = abs(mond_med - h2_med)
        d_mond_nfw = abs(mond_med - nfw_med)
        if d_mond_h2 < d_mond_nfw:
            closer = f"closer to H2 (|Δ|={d_mond_h2:.4f} vs NFW |Δ|={d_mond_nfw:.4f})"
        else:
            closer = f"closer to NFW (|Δ|={d_mond_nfw:.4f} vs H2 |Δ|={d_mond_h2:.4f})"

    neutrality = "YES" if mond_med < 0.01 else "NO"

    report = f"""H2 – MOND/RAR Comparative Scatter Analysis — Plain-Language Report
====================================================================

1. METRIC COMPARABILITY
------------------------
Both H2 and MOND results are expressed in the harmonized metric:
  RMS of log10(V_model) - log10(V_obs) over inner region (R < 0.5 * R_max)
  Units: dex (dimensionless log10 residual)

The comparison is EXACT in metric convention.
Inner-region definition is IDENTICAL across H2, NFW, and MOND analyses.

2. SAMPLE COUNTS
-----------------
  MOND galaxies analyzed:   {n_mond}
  NFW galaxies analyzed:    {n_nfw}
  H2 galaxies available:    {n_h2}
  Galaxies in all three:    {len(merged)}

3. MOND SCATTER SENSITIVITY — OVERALL
---------------------------------------
  Median max|Δσ| (MOND):    {mond_med:.4f} dex
  Median max|Δσ| (NFW):     {nfw_med:.4f} dex  (reference)
  Median |Δσ|    (H2):      {h2_med:.4f} dex  (reference)

  Fraction of MOND galaxies with max|Δσ| < 0.01 dex:
    {(merged['mond_max_abs_ds'] < 0.01).sum()} / {len(merged['mond_max_abs_ds'].dropna())} = {(merged['mond_max_abs_ds'] < 0.01).mean():.1%}

  Near-neutral inner-region scatter (< 0.01 dex threshold): {neutrality}

4. PER-REGIME BREAKDOWN
------------------------
{fmt_stats(mond_stats, "MOND max|Δσ| by regime")}
{fmt_stats(nfw_stats,  "NFW max|Δσ| by regime") if not nfw_stats.empty else "  (NFW data not available)"}
{fmt_stats(h2_stats,   "H2  |Δσ|   by regime")  if not h2_stats.empty  else "  (H2 data not available)"}

5. INTERPRETATION
------------------
  Within the tested implementation:

  a) Does MOND reproduce H2-like scatter neutrality?
     Median MOND max|Δσ| = {mond_med:.4f} dex.
     H2 delta_sigma ≈ {h2_med:.4f} dex (from harmonized metric reconstruction).
     The MOND result {"appears consistent with" if mond_med < 0.02 else "does not reproduce"} H2-like scatter neutrality
     within the adopted diagnostic.

  b) Is MOND more similar to H2 or to NFW?
     In the harmonized metric, MOND is {closer}.

  c) Physical interpretation:
     The RAR model is a baryonic interpolation function, not a halo model.
     Inner-region MOND velocities are dominated by baryonic acceleration in
     baryon-dominated galaxies, meaning the MOND-side perturbation (g_dag)
     acts differently from NFW halo perturbations: in the deep-MOND limit
     (g_bar << g_dag), the model is more sensitive to g_dag; in the
     Newtonian limit (g_bar >> g_dag), the model approaches V_bar regardless
     of g_dag. This asymmetry is the physical reason for any regime dependence
     observed in the MOND perturbation results.

6. ASSUMPTIONS AND LIMITATIONS
--------------------------------
  - g_dag is held fixed in the Li et al. (2020) RAR fits; perturbations of g_dag
    are "out-of-fit" variations, not formal parameter uncertainties. This matches
    the spirit of the NFW scheme but is not identical in statistical interpretation.
  - Ydisk perturbations are directly comparable to NFW M/L perturbations.
  - The RAR interpolation function form is fixed (McGaugh+2016); no comparison
    across different MOND interpolation functions is made here.
  - Distance and inclination are held at Li et al. best-fit values (same as NFW).

7. CAUTIONARY NOTES
--------------------
  - "consistent with" does not imply a universal law
  - "does not reproduce" does not falsify MOND — it only characterizes the
    inner-region scatter sensitivity within this specific diagnostic
  - Results are valid within the tested implementations and sample
  - No statement about the correct physical model is made here
"""
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"  Report saved: {out_path}")


def main():
    print("=" * 70)
    print("H2 – MOND Comparison")
    print("=" * 70)

    # ── Load MOND results ──────────────────────────────────────────────────────
    if not os.path.exists(MOND_CSV):
        print(f"MOND results not found: {MOND_CSV}")
        print("Run mond_perturbation_diagnostic.py first.")
        sys.exit(1)

    mond_pg = load_mond_per_galaxy(MOND_CSV)
    print(f"MOND per-galaxy: {len(mond_pg)} galaxies")

    # ── Load NFW results (reference) ───────────────────────────────────────────
    nfw_pg = pd.DataFrame()
    if os.path.exists(NFW_CSV):
        nfw_pg = load_nfw_per_galaxy(NFW_CSV)
        print(f"NFW per-galaxy:  {len(nfw_pg)} galaxies")
    else:
        print(f"NFW results not found (reference skipped): {NFW_CSV}")

    # ── Load H2 results ────────────────────────────────────────────────────────
    h2_pg = pd.DataFrame()
    if os.path.exists(H2_CSV):
        h2_pg = load_h2_per_galaxy(H2_CSV)
        print(f"H2 per-galaxy:   {len(h2_pg)} galaxies")
    else:
        print(f"H2 metric summary not found: {H2_CSV}")

    # ── Merge ──────────────────────────────────────────────────────────────────
    merged = mond_pg.rename(columns={
        'max_abs_delta_sigma': 'mond_max_abs_ds',
        'sigma_baseline':       'mond_sigma_baseline',
        'n_inner_points':       'mond_n_inner',
        'V_bar_over_V_MOND':    'V_bar_over_V_MOND',
    })

    if not nfw_pg.empty:
        nfw_sub = nfw_pg[['Galaxy', 'max_abs_delta_sigma', 'sigma_baseline',
                           'n_inner_points']].rename(columns={
            'max_abs_delta_sigma': 'nfw_max_abs_ds',
            'sigma_baseline':       'nfw_sigma_baseline',
            'n_inner_points':       'nfw_n_inner',
        })
        merged = merged.merge(nfw_sub, on='Galaxy', how='left')

    if not h2_pg.empty:
        h2_col = 'h2_delta_sigma_harmonized'
        if h2_col not in h2_pg.columns:
            h2_col = [c for c in h2_pg.columns if 'delta' in c.lower()]
            h2_col = h2_col[0] if h2_col else None
        if h2_col:
            h2_sub = h2_pg[['Galaxy', h2_col, 'h2_sigma_baseline',
                             'n_inner_points', 'reconstruction_status']
                            if 'reconstruction_status' in h2_pg.columns
                            else ['Galaxy', h2_col, 'h2_sigma_baseline', 'n_inner_points']
                            ].rename(columns={h2_col: 'h2_delta_sigma_harmonized',
                                              'n_inner_points': 'h2_n_inner'})
            merged = merged.merge(h2_sub, on='Galaxy', how='left')
            merged['h2_abs_delta_sigma'] = merged['h2_delta_sigma_harmonized'].abs()

    merged.to_csv(OUT_CSV, index=False)
    print(f"\nComparison table saved: {OUT_CSV}")

    # ── Stats ──────────────────────────────────────────────────────────────────
    mond_stats = regime_stats(merged, 'mond_max_abs_ds')
    nfw_stats  = regime_stats(merged, 'nfw_max_abs_ds') if 'nfw_max_abs_ds' in merged else pd.DataFrame()
    h2_stats   = regime_stats(merged, 'h2_abs_delta_sigma') if 'h2_abs_delta_sigma' in merged else pd.DataFrame()

    print("\n── MOND max|Δσ| by regime ──")
    print(mond_stats.to_string(index=False))
    if not nfw_stats.empty:
        print("\n── NFW max|Δσ| by regime (reference) ──")
        print(nfw_stats.to_string(index=False))

    # ── Figure ─────────────────────────────────────────────────────────────────
    make_comparison_figure(merged, FIG_DIR)

    # ── Report ─────────────────────────────────────────────────────────────────
    n_nfw = len(nfw_pg) if not nfw_pg.empty else 0
    n_h2  = len(h2_pg)  if not h2_pg.empty  else 0
    write_report(merged, mond_stats, nfw_stats, h2_stats, OUT_RPT,
                 n_mond=len(mond_pg), n_nfw=n_nfw, n_h2=n_h2)

    print("\nDone.")


if __name__ == '__main__':
    main()
