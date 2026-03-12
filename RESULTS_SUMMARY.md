# H2 A Null-Response Diagnostic of Adaptive Nonlocal Kernel Modulation in 80 SPARC Galaxies  Summary

> **Status:** 80-galaxy validation complete · Manuscript Submitted to The Astrophysical Journal 
> **Date:** 2026-03-09  
> **Sample:** 80 SPARC galaxies (canonical-baseline-compatible subset)

---

## Overview

This validation demonstrates that bounded local adaptive kernel modulation produces regime-dependent local response while maintaining scatter neutrality in the inner regions across 80 SPARC galaxies spanning the full range from extreme baryon domination to dark-matter domination.

**Key result:** Inner-region scatter remains invariant under adaptive kernel modulation despite strong regime-dependent correlations between local baryonic structure and velocity residuals.

---

## Sample Characteristics

| Parameter | Value |
|-----------|-------|
| **Total galaxies** | 80 |
| **Baryon-dominated** | 39 (Vobs/Vb < 0.8) |
| **Balanced** | 30 (0.8 ≤ Vobs/Vb < 1.5) |
| **DM-dominated** | 11 (Vobs/Vb ≥ 1.5) |
| **Vobs/Vb range** | 0.027 to 2.187 |
| **Canonical baseline** | L₀ = 200 kpc, μ = 10, β = 1.15 |

All 80 galaxies share identical frozen H1 baseline parameters, ensuring homogeneous adaptive conditions across the validation sample.

---

## Primary Results

### Test-2: χ–ΔV Correlation (Regime-Dependent Response)

**Pearson correlation range:** r = −0.371 to +0.981

**Regime-dependent medians (with 95% bootstrap CIs):**
- **Baryon-dominated** (n=39): median r = +0.73 [0.65, 0.80]
- **Balanced** (n=30): median r = +0.77 [0.68, 0.85]
- **DM-dominated** (n=11): median r = +0.83 [0.71, 0.91]

**Peak correlations:**
- Maximum positive: UGC07577 (r = +0.981)
- Minimum (anti-correlation): NGC2955 (r = −0.371)

**Interpretation:** The adaptive kernel responds systematically to local baryonic structure across all dynamical regimes, with strongest coupling in DM-dominated systems.

---

### Test-3: Inner-Region Scatter (Universal Neutrality)

**Δσ measurement results:**
- **Measured systems:** 74/80 (92.5%)
- **Δσ ≈ 0 km/s:** 74/74 (100% of measurable systems)
- **Undefined (inner_coverage):** 6/80 (7.5%)

**Systems with undefined Δσ due to insufficient inner-region data:**
F567-2, UGC09992, UGC05999, NGC6789, NGC4389, UGC00634

**Stop rule:** |Δσ| > 10⁻¹² (numerical precision threshold)  
**Violations:** 0 across all 80 galaxies

**Interpretation:** Despite regime-dependent adaptive coupling (Test-2), the inner-region scatter metric remains invariant, consistent with quadrature suppression in baryon-dominated inner regions.

---

### Geometric Floor Identification

**Empirical operating floor:** Leff,min = L₀/3 ≈ 66.67 kpc

**Leff,min range:** 66.67 to 96.95 kpc

**Systems reaching floor directly:** ~60% of sample  
**Systems approaching asymptotically:** ~40% of sample

**Interpretation:** The discretized eight-scale basis imposes a universal structural boundary on adaptive scale modulation, independent of galaxy dynamical regime.

---

## High-MAFE Stress Tests

Several systems with poor H1 baseline fits were deliberately retained to test universality:

| Galaxy | MAFE | Regime | Δσ Result |
|--------|------|--------|-----------|
| UGC00191 | 0.489 | Balanced | 0.0000 km/s |
| UGC08550 | 0.898 | Balanced | 0.0000 km/s |
| DDO154 | 1.154 | DM-dom | 0.0000 km/s |
| UGC05005 | 0.951 | DM-dom | 0.0000 km/s |

**Full-fleet MAFE summary:**  
Median = 0.21, Range = 0.07 to 1.67

**Result:** Scatter neutrality persists even when H1 baseline quality is poor, demonstrating that the null result is not an artifact of cherry-picked well-fit systems.

---

## Key Findings

1. **Regime-dependent adaptive response confirmed:** Pearson correlations between local baryonic structure (χ) and velocity residuals (ΔV) vary systematically across regimes, with strongest coupling in DM-dominated systems (median r = +0.83).

2. **Universal scatter neutrality:** Inner-region scatter remains unchanged (Δσ ≈ 0 km/s) for all 74 measurable systems, despite strong local adaptive response.

3. **Geometric floor identified:** All 80 systems reach or approach Leff,min = L₀/3 = 66.67 kpc, establishing a structural boundary imposed by the discretized kernel basis.

4. **Independence from baseline quality:** High-MAFE stress tests (including one system with MAFE > 1.0) maintain scatter neutrality, ruling out selection bias.

5. **Quadrature suppression mechanism:** The combination of active local response (variable Pearson r) with invariant global scatter is consistent with suppression by quadrature velocity coupling when the nonlocal contribution remains subdominant in inner regions.

---

## Statistical Robustness

**Sample size:** 80 galaxies spanning factor of ~81 in Vobs/Vb ratio

**Regime coverage:**
- Extreme baryon domination (Vobs/Vb = 0.027)
- Extreme DM domination (Vobs/Vb = 2.187)
- Continuous coverage across balanced regime

**Incremental validation trajectory:**
- 9 galaxies → 30 galaxies → 50 galaxies → 80 galaxies
- Staged expansion with consistent results at each phase
- Zero stop-rule violations across all phases

---

## Limitations

1. **Canonical-baseline-compatible subset:** The 80-galaxy sample represents systems with homogeneous H1 baseline parameters (L₀=200, μ=10). This subset selection was necessary to preserve controlled adaptive conditions and avoid per-subset recalibration.

2. **Inner-coverage undefined cases:** Six galaxies lack sufficient inner-region radial data points for scatter computation. This is a data geometry limitation, not a model failure (all six show good Pearson r values).

3. **Phenomenological interpretation of L₀=200 kpc:** The fixed scale L₀ = 200 kpc should be interpreted as a phenomenological coupling scale of the present implementation rather than a localized interaction range tied to optical disk size.

---

## Validation Products

Complete validation outputs available in `H2_PUBLICATION_RELEASE/`:

### Fleet Expansion Data
- `fleet_expansion_30galaxy/` — Initial 30-galaxy expansion
- `fleet_expansion_50galaxy/` — Intermediate 50-galaxy state  
- `fleet_expansion_80galaxy/` — Final 80-galaxy validation
  - `fleet_summary_80galaxy.csv` — Complete per-galaxy results
  - `figures/regime_trends_80galaxy.png` — Pearson r and Δσ vs regime
  - `figures/leff_overlay_80galaxy.png` — All 80 Leff(r) profiles
  - `README.md` — Detailed validation documentation

### Per-Galaxy Diagnostics
- `data/derived/phase3/` — Leff(r) profiles for all galaxies
- `data/derived/phase4/` — H2 rotation curves and comparisons
- `data/derived/fleet/` — Summary statistics and trend analysis

---

## Publication Status

**Manuscript:** Submitted to Astronomy & Astrophysics  
**Title:** "H2 Adaptive Kernel Diagnostics: Scatter Neutrality and Geometric Scaling Limits in 80 SPARC Galaxies"  
**DOI:** Pending acceptance  
**Preprint:** Available upon request

**Repository DOI:** [10.5281/zenodo.18793233](https://doi.org/10.5281/zenodo.18793233)

---

## Citation

If you use this validation data or methodology, please cite:

> Figueiredo, V. M. F. (2026). H2 Adaptive Kernel Framework for Galaxy Rotation Curves (v2.0.0).  
> Zenodo. https://doi.org/10.5281/zenodo.18793233

And the SPARC database:

> Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves. AJ, 152, 157.

---

## Contact

**Vítor M. F. Figueiredo**  
Email: vitor.figueiredo.research@protonmail.com  
ORCID: [0009-0004-7358-4622](https://orcid.org/0009-0004-7358-4622)

---

*Last updated: 2026-03-09*


