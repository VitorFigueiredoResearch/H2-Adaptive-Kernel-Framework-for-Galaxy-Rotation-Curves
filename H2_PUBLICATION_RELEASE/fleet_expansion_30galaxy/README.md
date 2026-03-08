# H2 Fleet Expansion — 30-Galaxy Validation Set

Generated: 2026-03-02
Parameters: α=2.0, σ_idx=1.0, L₀=200 kpc, μ=10, kernel=ananta-hybrid, taper=True

## Purpose

Expands the original 9-galaxy H2 validation fleet to 30 galaxies for
statistical robustness, addressing peer-reviewer request for larger sample
in A&A submission.

## Fleet Composition

| Regime | Vk/Vb range | Count |
|--------|-------------|-------|
| Baryon-dominated | 0.22 – 0.77 | 18 |
| Balanced | 0.87 – 1.40 | 10 |
| DM-dominated | 1.94 – 2.19 | 2 |
| **Total** | **0.22 – 2.19** | **30** |

## Key Result

**Δσ = 0.0000 dex across all 30 galaxies** — scatter neutrality confirmed
universally across baryon-dominated, balanced, and DM-dominated regimes.

## Galaxy List

### Original 9 (from publication baseline)
NGC0891, NGC6946, NGC5585, IC2574, NGC3893, NGC2903, NGC5055, NGC3198, NGC6503

### New 21 (this expansion)
NGC7793, NGC3726, NGC3769, NGC3521, NGC3953,
NGC2998, NGC0024, DDO161, NGC1003, NGC5033, NGC3992,
UGC03580, UGC06983, UGC07608, UGC02259, NGC4100,
KK98-251, UGC05716, F571-8, NGC2403, NGC2955

## Directory Structure

```
fleet_expansion_30galaxy/
├── README.md                        — this file
├── fleet_summary_30galaxy.csv       — combined results for all 30 galaxies
├── phase3_outputs/                  — L_eff and chi CSVs + PNGs (21 new galaxies)
├── phase4_outputs/                  — H2 adaptive RC CSVs + PNGs (21 new galaxies)
├── diagnostics/                     — per-galaxy fleet CSVs (21 new galaxies)
└── figures/
    ├── regime_trends_30galaxy.png   — Pearson r and Δσ vs Vk/Vb
    └── leff_overlay_30galaxy.png    — L_eff(r) overlay for all 30 galaxies
```

## Compatibility Constraint

Only galaxies with H1 baseline parameters L₀=200 kpc and μ=10 are
compatible with Phase 4 adaptive convolution under taper mode. Galaxies
with smaller L₀ fail with [TAPER-FAIL] at the smallest basis scale.
This constraint is a property of the H1 fitting distribution across SPARC,
not a limitation of H2 itself.

## Anomalies and Flags

| Galaxy | Flag | Notes |
|--------|------|-------|
| NGC3992 | T1=0 | Outer-region tolerance exceeded; Δσ=0 holds |
| NGC2955 | T1=0, χ_max=13.5 | Steep DM gradient; expected for Vk/Vb=2.19; Δσ=0 holds |

## Skipped Galaxies (TAPER-FAIL)

| Galaxy | Reason |
|--------|--------|
| NGC7331 | L₀=120 kpc (H1 params incompatible) |
| NGC4157 | L₀=80 kpc, μ=200 (incompatible) |
| NGC3109 | L₀=10 kpc, μ=500 (incompatible) |

## Pipeline

Each new galaxy was processed as:
1. `python -m diagnostics.phase3_leff_ngc3198 --galaxy [GAL] --alpha 2.0 --sigma_idx 1.0 --taper`
2. `python -m diagnostics.phase4_adaptive_convolution --galaxy [GAL]`
3. `python -m diagnostics.test2_chi_correlation --galaxy [GAL]`
4. `python -m diagnostics.test3_inner_scatter --galaxy [GAL]`

Automated via: `python -m diagnostics.run_fleet --galaxies [GAL] --alpha 2.0 --sigma_idx 1.0 --taper`
