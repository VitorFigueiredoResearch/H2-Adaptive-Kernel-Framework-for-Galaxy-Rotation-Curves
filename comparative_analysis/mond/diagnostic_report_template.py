"""
diagnostic_report_template.py
==============================
Generates comparative_analysis/mond/diagnostic_report.txt
after the main perturbation diagnostic has run.

Run from repository root AFTER mond_perturbation_diagnostic.py:
  python comparative_analysis/mond/diagnostic_report_template.py
"""

import os, sys
import numpy as np
import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
MOND_CSV  = os.path.join(REPO_ROOT, 'comparative_analysis', 'mond',
                          'mond_perturbation_summary.csv')
FLEET_CSV = os.path.join(REPO_ROOT, 'H2_PUBLICATION_RELEASE',
                          'fleet_expansion_80galaxy', 'fleet_summary_80galaxy.csv')
RAR_TABLE = os.path.join(REPO_ROOT, 'data', 'nfw', 'Fits', 'ByModel',
                          'Table', 'parameter_RAR.mrt')
OUT_RPT   = os.path.join(REPO_ROOT, 'comparative_analysis', 'mond',
                          'diagnostic_report.txt')

sys.path.insert(0, os.path.join(REPO_ROOT, 'comparative_analysis', 'mond'))
from mond_velocity import _parse_rar_table


def main():
    fleet   = pd.read_csv(FLEET_CSV)
    rar_map = _parse_rar_table(RAR_TABLE)

    all_gals    = set(fleet['galaxy'].tolist())
    rar_matched = set(rar_map.keys()) & all_gals
    no_match    = all_gals - rar_matched

    if os.path.exists(MOND_CSV):
        df = pd.read_csv(MOND_CSV)
        succeeded = set(df['Galaxy'].unique())
        inner_ok  = df[df['n_inner_points'] >= 3]['Galaxy'].unique()
        low_inner = df[df['n_inner_points'] <  3]['Galaxy'].unique()
    else:
        succeeded = set()
        inner_ok  = []
        low_inner = []

    excluded = all_gals - succeeded

    lines = [
        "MOND Comparative Analysis — Diagnostic Report",
        "=" * 60,
        f"Generated: 2026-03-24",
        "",
        "SAMPLE SUMMARY",
        "-" * 40,
        f"  H2 fleet galaxies:          {len(all_gals)}",
        f"  RAR catalog matches:        {len(rar_matched)}",
        f"  Galaxies successfully run:  {len(succeeded)}",
        f"  Galaxies excluded total:    {len(all_gals) - len(succeeded)}",
        f"  Galaxies with n_inner >= 3: {len(inner_ok)}",
        f"  Galaxies with n_inner <  3: {len(low_inner)}",
        "",
        "NO-MATCH GALAXIES (in H2 fleet but not in RAR catalog)",
        "-" * 40,
    ]
    if no_match:
        for g in sorted(no_match):
            lines.append(f"  {g}")
    else:
        lines.append("  None — 100% match rate")

    lines += [
        "",
        "LOW INNER-REGION POINT COUNT (n_inner < 3, excluded from scatter metric)",
        "-" * 40,
    ]
    if len(low_inner):
        for g in sorted(low_inner):
            n = df[df['Galaxy'] == g]['n_inner_points'].iloc[0]
            lines.append(f"  {g}  (n_inner={n})")
    else:
        lines.append("  None")

    lines += [
        "",
        "METRIC AND CONVENTION NOTES",
        "-" * 40,
        "  Scatter metric: RMS of log10(V_model) - log10(V_obs) [dex]",
        "  Inner region:   R < 0.5 * R_max",
        "  These are IDENTICAL to H2 and NFW analysis conventions.",
        "",
        "  g_dag perturbations: ±10%, ±20% around 1.2e-10 m/s^2",
        "  Ydisk perturbations: ±10%, ±20% around Li et al. best-fit value",
        "  Combined (both simultaneous): ±10%, ±20%",
        "",
        "  No interpolation-function switching was performed (Li et al. 2020",
        "  uses only the McGaugh+2016 form; no 'simple' vs 'standard' variant).",
        "",
        "ASSUMPTIONS REQUIRED",
        "-" * 40,
        "  1. g_dag perturbations are 'out-of-fit' variations (g_dag is fixed",
        "     in the Li et al. RAR fits, not a free parameter). This is",
        "     intentional and mirrors the NFW scheme (bounded but not refit).",
        "  2. Ydisk/Ybul held at Li et al. best-fit for g_dag perturbations.",
        "  3. g_dag held at default for Ydisk perturbations.",
        "  4. Distance and inclination are not perturbed (same as NFW choice).",
        "",
        "NARRATIVE/METRIC MISMATCH",
        "-" * 40,
        "  H2 public wording ('delta_sigma_kms') uses km/s language in some",
        "  narrative contexts, but the operational H2 diagnostic (test3_inner_scatter.py)",
        "  uses log10-RMS in dex. The harmonized comparison uses the OPERATIONAL",
        "  metric (dex) for mathematical consistency.",
    ]

    text = '\n'.join(lines)
    with open(OUT_RPT, 'w', encoding='utf-8') as f:
        f.write(text)
    print(f"Diagnostic report saved: {OUT_RPT}")


if __name__ == '__main__':
    main()
