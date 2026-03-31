"""
generate_three_way_comparison.py
=================================
Generates the three-way inner-region scatter sensitivity comparison figure
for the MNRAS paper:

  "Inner-Region Scatter Response to Bounded Perturbations:
   A Comparative Analysis of 74 SPARC Galaxies"

Figure design
-------------
Upper panel : Scatter sensitivity (max |Δσ|, dex) as a function of the
              median inner-region baryonic-dominance ratio V_bar/V_NFW,
              for NFW and RAR on the common 74-galaxy set, and for H2 on
              the full reconstructed 74-galaxy set (Tier A + Tier B + Tier C).
              NFW points (circles) show the strong anti-correlation trend.
              RAR points (diamonds) occupy an intermediate, regime-flat level.
              H2 points (squares) are shown for ALL 74 galaxies.

Lower panel : Box plots of the max |Δσ| distributions:
              NFW (N=74), RAR (N=74), H2 (N=74).
              Sample sizes are labelled on the x-axis.

Data sources
------------
  comparative_analysis/nfw/h2_nfw_comparison_summary.csv
      - NFW max_abs_delta_sigma, V_bar_over_V_NFW_inner_median, h2_regime
      - Used for 74-galaxy NFW data and for V_bar/V_NFW x-positions
  comparative_analysis/mond/h2_mond_comparison_summary.csv
      - RAR (MOND) mond_max_abs_ds for 74 galaxies
  h2_74galaxy_combined_metrics.csv   ← UPDATED SOURCE (was h2_full74_explicit_summary.csv)
      - H2 abs_delta_sigma for ALL 74 galaxies (Tier A + B + C)
      - Combined authoritative dataset produced by Phase 4 of the
        H2 74-galaxy reconstruction feasibility test (2026-03-29)
      - Method: compute_h2_full74.py V_total_H2 explicit load

The same six galaxies excluded in the paper (F567-2, NGC4389, NGC6789,
UGC00634, UGC05999, UGC09992) are already absent from these files.

Spearman results (from comparative_validation/rar_spearman_result.txt)
------------------------------------------------------------------
  NFW: ρ = −0.899,  p = 1.8×10⁻²⁷,  N = 74  (strong anti-correlation)
  RAR: ρ = +0.049,  p = 0.681,        N = 74  (null — not significant)

Usage
-----
    python generate_three_way_comparison.py [--output PATH]

Output
------
    three_way_comparison.png  (saved to the paper folder by default, or
                                to the path given by --output)

Script version: 3.0  (2026-03-29)
Change from v2.0:
  - H2 now uses ALL 74 galaxies (Tier A + Tier B + Tier C)
  - H2 source changed from h2_full74_explicit_summary.csv (Tier A+B only)
    to h2_74galaxy_combined_metrics.csv (full 74-galaxy authoritative table)
  - Removed N=30 hard assertion; now asserts N=74
  - H2 median updated from 0.000860 to 0.000478 dex (combined 74-galaxy result)
  - All N-labels updated from N=30 to N=74
  - Figure title updated to reflect N=74 H2 sample
  - H2 x-positions for scatter panel: merged from NFW CSV (V_bar/V_NFW)
    for all 74 galaxies; galaxies without NFW x-position are excluded
    from scatter panel only (reported in load summary)
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Paths (relative to repo root)
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..')
)
NFW_CSV = os.path.join(
    REPO_ROOT,
    'comparative_analysis', 'nfw', 'h2_nfw_comparison_summary.csv'
)
MOND_CSV = os.path.join(
    REPO_ROOT,
    'comparative_analysis', 'mond', 'h2_mond_comparison_summary.csv'
)
# UPDATED: full 74-galaxy combined metrics (replaces h2_full74_explicit_summary.csv)
H2_COMBINED_CSV = os.path.join(
    REPO_ROOT,
    'h2_74galaxy_combined_metrics.csv'
)

# Default output: paper folder
PAPER_DIR = os.path.abspath(
    os.path.join(REPO_ROOT, '..', '..', 'paper')
)
DEFAULT_OUTPUT = os.path.join(PAPER_DIR, 'three_way_comparison.png')

# ---------------------------------------------------------------------------
# Colour / style constants (MNRAS-compatible, colour-blind friendly)
# ---------------------------------------------------------------------------
COL_NFW  = '#D55E00'   # reddish-orange
COL_RAR  = '#0072B2'   # deep blue
COL_H2   = '#009E73'   # green

ALPHA_SCATTER = 0.75
MARKER_SIZE   = 22      # scatter marker area (s parameter)


def load_data(nfw_path, mond_path, h2_combined_path):
    """
    Load and prepare validated CSV outputs.

    Returns
    -------
    nfw_df   : DataFrame, 74 rows, NFW data with V_bar/V_NFW x-positions
    rar_df   : DataFrame, 74 rows, RAR data
    h2_df    : DataFrame, up to 74 rows, H2 all tiers merged with NFW x-positions
    merged74 : DataFrame, 74-row merge of NFW+RAR for joint plotting
    """
    nfw_df   = pd.read_csv(nfw_path)
    mond_df  = pd.read_csv(mond_path)
    h2_raw   = pd.read_csv(h2_combined_path)

    # --- H2: all 74 galaxies (Tier A + B + C) ---
    # h2_74galaxy_combined_metrics.csv has column 'abs_delta_sigma'
    h2_all = h2_raw.copy()
    assert len(h2_all) == 74, (
        f"Expected 74 H2 rows in combined metrics, got {len(h2_all)}"
    )
    assert h2_all['abs_delta_sigma'].isna().sum() == 0, (
        "Unexpected NaN in abs_delta_sigma for combined H2 data"
    )

    # --- Merge H2 with NFW to get V_bar/V_NFW x-positions ---
    # NFW CSV key column: 'h2_galaxy'; H2 combined key: 'Galaxy'
    h2_with_x = pd.merge(
        h2_all[['Galaxy', 'Regime', 'Tier', 'abs_delta_sigma']],
        nfw_df[['h2_galaxy', 'V_bar_over_V_NFW_inner_median', 'h2_regime']],
        left_on='Galaxy', right_on='h2_galaxy',
        how='inner'
    )
    n_h2_scatter = len(h2_with_x)
    n_h2_dropped = len(h2_all) - n_h2_scatter

    # --- Merge NFW and RAR for the 74-galaxy joint dataset ---
    merged74 = pd.merge(
        nfw_df[['h2_galaxy', 'h2_regime', 'max_abs_delta_sigma',
                'V_bar_over_V_NFW_inner_median']],
        mond_df[['Galaxy', 'mond_max_abs_ds']],
        left_on='h2_galaxy', right_on='Galaxy',
        how='inner'
    )

    print(f"[INFO] NFW dataset           : {len(nfw_df)} galaxies")
    print(f"[INFO] RAR dataset           : {len(mond_df)} galaxies")
    print(f"[INFO] H2 combined (all tiers): {len(h2_all)} galaxies")
    print(f"[INFO] H2 merged with NFW x  : {n_h2_scatter} galaxies (scatter panel)")
    if n_h2_dropped > 0:
        print(f"[INFO] H2 dropped (no NFW x) : {n_h2_dropped} (not in scatter panel)")
    print(f"[INFO] 74-galaxy merged set  : {len(merged74)} galaxies")
    print(f"")
    print(f"[INFO] NFW  median max|Δσ|   : {merged74['max_abs_delta_sigma'].median():.4f} dex  (N=74)")
    print(f"[INFO] RAR  median max|Δσ|   : {merged74['mond_max_abs_ds'].median():.4f} dex  (N=74)")
    print(f"[INFO] H2   median |Δσ|      : {h2_all['abs_delta_sigma'].median():.6f} dex  (N=74 all tiers)")
    print(f"[INFO] H2   median |Δσ|      : {h2_with_x['abs_delta_sigma'].median():.6f} dex  (N={n_h2_scatter} scatter panel)")

    return merged74, h2_with_x, h2_all


def make_figure(merged74, h2_scatter_df, h2_all_df, output_path):
    """
    Create and save the two-panel comparison figure.

    Upper panel : scatter sensitivity vs V_bar/V_NFW.
                  NFW: 74 galaxies. RAR: 74 galaxies. H2: all 74 galaxies
                  (scatter panel uses galaxies with NFW x-position available).
    Lower panel : box plots — NFW (N=74), RAR (N=74), H2 (N=74).
    """
    fig, axes = plt.subplots(
        2, 1,
        figsize=(3.46, 5.5),     # MNRAS single-column: 3.46 in wide
        gridspec_kw={'height_ratios': [1.5, 1.0],
                     'hspace': 0.42}
    )
    ax_scatter, ax_box = axes

    x_nfw  = merged74['V_bar_over_V_NFW_inner_median'].values
    y_nfw  = merged74['max_abs_delta_sigma'].values
    y_rar  = merged74['mond_max_abs_ds'].values

    x_h2   = h2_scatter_df['V_bar_over_V_NFW_inner_median'].values
    y_h2_scatter = h2_scatter_df['abs_delta_sigma'].values
    y_h2_all     = h2_all_df['abs_delta_sigma'].values    # for box plot

    # -------------------------------------------------------------------
    # Upper panel: scatter plot
    # -------------------------------------------------------------------

    # NFW — 74 galaxies (circles)
    ax_scatter.scatter(
        x_nfw, y_nfw,
        s=MARKER_SIZE, marker='o',
        color=COL_NFW, alpha=ALPHA_SCATTER,
        edgecolors='none', zorder=4,
        label='NFW'
    )

    # RAR — 74 galaxies (diamonds)
    ax_scatter.scatter(
        x_nfw, y_rar,
        s=MARKER_SIZE, marker='D',
        color=COL_RAR, alpha=ALPHA_SCATTER,
        edgecolors='none', zorder=4,
        label='RAR'
    )

    # H2 — all 74 galaxies (squares); x-position from NFW merge
    ax_scatter.scatter(
        x_h2, y_h2_scatter,
        s=MARKER_SIZE * 1.1, marker='s',
        color=COL_H2, alpha=0.90,
        edgecolors='none', zorder=5,
        label=r'H2 ($N{=}74$)'
    )

    # Median reference lines
    ax_scatter.axhline(
        np.median(y_nfw), color=COL_NFW, lw=1.4, ls='--', alpha=0.80, zorder=2
    )
    ax_scatter.axhline(
        np.median(y_rar), color=COL_RAR, lw=1.4, ls='--', alpha=0.80, zorder=2
    )
    ax_scatter.axhline(
        np.median(y_h2_all), color=COL_H2, lw=1.2, ls=':', alpha=0.90, zorder=2
    )

    # Linear trend for NFW (visual guide)
    p = np.polyfit(x_nfw, y_nfw, 1)
    x_fit = np.linspace(x_nfw.min(), x_nfw.max(), 100)
    ax_scatter.plot(
        x_fit, np.polyval(p, x_fit),
        color=COL_NFW, lw=1.0, ls='-', alpha=0.45, zorder=1
    )

    ax_scatter.set_xlabel(
        r'$V_{\rm bar}/V_{\rm NFW}$ (inner median)', fontsize=7
    )
    ax_scatter.set_ylabel(r'$\max\,|\Delta\sigma|\;(\rm dex)$', fontsize=7)
    ax_scatter.tick_params(labelsize=6.5)
    ax_scatter.set_ylim(-0.005, 0.130)
    ax_scatter.set_xlim(0.20, 3.40)

    # Spearman annotations (validated; from comparative_validation/rar_spearman_result.txt)
    ax_scatter.annotate(
        r'NFW: $\rho = -0.899$',
        xy=(3.30, 0.118), fontsize=5.8,
        color=COL_NFW, ha='right', va='top'
    )
    ax_scatter.annotate(
        r'RAR: $\rho = +0.049$ (n.s.)',
        xy=(3.30, 0.104), fontsize=5.8,
        color=COL_RAR, ha='right', va='top'
    )

    # Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_NFW,
               markersize=5, label='NFW ($N=74$)'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor=COL_RAR,
               markersize=5, label='RAR ($N=74$)'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor=COL_H2,
               markersize=5, label='H2 ($N=74$)'),
        Line2D([0], [0], color='grey', lw=1.2, ls='--', label='median'),
    ]
    ax_scatter.legend(
        handles=legend_elements, fontsize=5.8, loc='upper left',
        framealpha=0.88, edgecolor='none'
    )

    # Panel label
    ax_scatter.text(
        0.97, 0.97, '(a)', transform=ax_scatter.transAxes,
        fontsize=7, va='top', ha='right', fontweight='bold'
    )

    # -------------------------------------------------------------------
    # Lower panel: box plots
    # NOTE: All three implementations now use N=74.
    # -------------------------------------------------------------------
    data_box    = [y_h2_all, y_rar, y_nfw]
    colours_box = [COL_H2, COL_RAR, COL_NFW]
    tick_lbls   = [
        f'H2\n($N=74$)',
        f'RAR\n($N=74$)',
        f'NFW\n($N=74$)',
    ]

    bp = ax_box.boxplot(
        data_box,
        patch_artist=True,
        medianprops=dict(color='white', lw=1.8),
        whiskerprops=dict(lw=1.0),
        capprops=dict(lw=1.0),
        flierprops=dict(marker='+', markersize=4, alpha=0.6),
        widths=0.55,
        zorder=3
    )

    # Apply colours
    for patch, col in zip(bp['boxes'], colours_box):
        patch.set_facecolor(col)
        patch.set_alpha(0.75)
    for whisker in bp['whiskers']:
        whisker.set_color('grey')
    for cap in bp['caps']:
        cap.set_color('grey')

    # Set tick labels via ax.set_xticks (avoids deprecated boxplot 'labels' kwarg)
    ax_box.set_xticks([1, 2, 3])
    ax_box.set_xticklabels(tick_lbls, fontsize=6)

    # Annotate medians
    medians_box = [np.median(d) for d in data_box]
    offsets = [0.004, 0.003, 0.003]
    for i, (med, col, off) in enumerate(zip(medians_box, colours_box, offsets), start=1):
        ax_box.text(
            i, med + off, f'{med:.4f}',
            ha='center', va='bottom', fontsize=5.5, color='black'
        )

    ax_box.set_ylabel(r'$\max\,|\Delta\sigma|\;(\rm dex)$', fontsize=7)
    ax_box.tick_params(axis='y', labelsize=6.5)
    ax_box.set_ylim(-0.005, 0.135)

    # Panel label
    ax_box.text(
        0.97, 0.97, '(b)', transform=ax_box.transAxes,
        fontsize=7, va='top', ha='right', fontweight='bold'
    )

    # -------------------------------------------------------------------
    # Figure title and save
    # -------------------------------------------------------------------
    fig.suptitle(
        'Inner-region scatter sensitivity: NFW/RAR/H2 ($N=74$)',
        fontsize=6.8, y=0.99
    )

    fig.savefig(output_path, dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"[INFO] Figure saved: {output_path}")
    plt.close(fig)

    # --- Summary printout ---
    print(f"\n[SUMMARY]")
    print(f"  NFW  median max|Δσ| = {np.median(y_nfw):.4f} dex  (N=74)")
    print(f"  RAR  median max|Δσ| = {np.median(y_rar):.4f} dex  (N=74)")
    print(f"  H2   median |Δσ|    = {np.median(y_h2_all):.6f} dex  (N=74 all tiers)")
    print(f"  NFW  Spearman rho = -0.899  (validated)")
    print(f"  RAR  Spearman rho = +0.049  (null; validated)")


def main():
    parser = argparse.ArgumentParser(
        description='Generate three-way scatter sensitivity comparison figure (v3.0).'
    )
    parser.add_argument(
        '--output', default=DEFAULT_OUTPUT,
        help=f'Output PNG path (default: {DEFAULT_OUTPUT})'
    )
    args = parser.parse_args()

    merged74, h2_scatter_df, h2_all_df = load_data(NFW_CSV, MOND_CSV, H2_COMBINED_CSV)
    make_figure(merged74, h2_scatter_df, h2_all_df, args.output)


if __name__ == '__main__':
    main()
