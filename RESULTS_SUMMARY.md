# Inner-Region Scatter Response to Bounded Perturbations: A Comparative Analysis of 74 SPARC Galaxies — Results Summary

> **Status:** Manuscript submitted to Monthly Notices of the Royal Astronomical Society (MNRAS)
> **Date:** 2026-03-28
> **Sample:** 74 SPARC galaxies (common comparison set)

---

## Overview

This paper presents a three-way inner-region scatter comparison across 74 SPARC galaxies between the NFW profile (Li et al. 2020), the Radial Acceleration Relation (RAR; Li et al. 2020), and the H2 adaptive kernel framework (Figueiredo 2026). The central question: does inner-region rotation-curve scatter respond to bounded perturbations applied to each modelling framework, and does the pattern of response differ between frameworks?

**Key result:** NFW perturbations produce a strong, statistically significant negative Spearman correlation between baryonic dominance and inner-region scatter (ρ = −0.899, p = 1.8×10⁻²⁷). RAR perturbations produce a null response (ρ = +0.049, p = 0.681). H2 perturbations produce near-zero log₁₀-dex scatter across all 30 explicitly archived galaxies (median |Δσ| = 0.000860 dex).

---

## Sample

| Parameter | Value |
|-----------|-------|
| **Common comparison set** | 74 galaxies |
| **Source catalogue** | SPARC (Lelli et al. 2016) |
| **NFW/RAR fits** | Li et al. (2020) |
| **H2 explicit archived** | 30 galaxies (Tier A+B) |
| **H2 Tier C (no archive)** | 44 galaxies — not plotted |
| **Excluded from 80-galaxy fleet** | 6 galaxies (inner-coverage undefined) |

The 6 excluded galaxies (F567-2, UGC09992, UGC05999, NGC6789, NGC4389, UGC00634) lack sufficient inner-region radial data for scatter computation and are excluded from all three analyses.

---

## Primary Results

### NFW — Perturbation Response

**Metric:** max|Δσ| (maximum absolute inner-region scatter across perturbation grid)
**Median (N=74):** 0.0548 dex
**Spearman correlation (scatter vs. V_bar/V_NFW):** ρ = −0.899, p = 1.8×10⁻²⁷

NFW inner-region scatter is strongly anti-correlated with baryonic dominance ratio: baryon-dominated galaxies show the largest scatter response; DM-dominated galaxies show the smallest. This indicates NFW inner-region fits are most sensitive to perturbation in the regime where baryons dominate the inner potential.

---

### RAR — Perturbation Response

**Metric:** max|Δσ| (maximum absolute inner-region scatter across perturbation grid)
**Median (N=74):** 0.0321 dex
**Spearman correlation (scatter vs. V_bar/V_NFW):** ρ = +0.049, p = 0.681 (not significant)

RAR inner-region scatter shows no statistically significant correlation with baryonic dominance. The perturbation response is regime-independent, consistent with the RAR's self-consistent embedding of baryonic structure via the acceleration relation.

---

### H2 — Perturbation Response

**Metric:** |Δσ| (log₁₀-dex inner-region scatter)
**Median (N=30 explicit Tier A+B):** 0.000860 dex
**Tier C (N=44):** No archived log₁₀-dex values; not plotted

H2 inner-region scatter remains near the numerical floor for all 30 explicitly archived galaxies. This null response persists across baryon-dominated, balanced, and DM-dominated systems, consistent with the quadrature suppression mechanism described in the H2 paper.

The 44 Tier C galaxies have fleet-summary Δσ values of 0 km/s in velocity units but no archived log₁₀-dex output. Their inclusion in the plotted distribution would require a full H1-basis-pipeline rerun. They are documented but excluded from the comparative figure (see `comparative_analysis/referee_resolution/tier_c_audit_report.txt`).

---

## Three-Way Comparison Summary

| Model | N | Median max\|Δσ\| (dex) | Spearman ρ | p-value |
|-------|---|------------------------|------------|---------|
| NFW | 74 | 0.0548 | −0.899 | 1.8×10⁻²⁷ |
| RAR | 74 | 0.0321 | +0.049 | 0.681 |
| H2 | 30 | 0.000860 | — | — |

The H2 median is approximately 37× below the RAR median and 64× below the NFW median. The H2 distribution does not span the regime-dominance axis in the same way as NFW/RAR (different perturbation scheme and different archived galaxy count), so a direct Spearman comparison to regime is not the primary H2 result.

---

## Key Findings

1. **NFW scatter is regime-sensitive:** Strong negative correlation (ρ = −0.899) between baryonic dominance and inner-region scatter implies NFW fits are most perturbed in baryon-dominated systems.

2. **RAR scatter is regime-independent:** Null Spearman result (ρ = +0.049, n.s.) consistent with RAR's self-consistent acceleration-based structure.

3. **H2 scatter remains near-zero:** Median 0.000860 dex across 30 explicit archived galaxies, consistent with the bounded kernel suppression mechanism reported in the H2 single-framework paper.

4. **Inner-region definition:** R < 0.5×R_max; metric is log₁₀-RMS residual in dex units. All three frameworks use the same spatial definition for comparability.

---

## Reproducibility

The main comparative figure (`three_way_comparison.png`) is fully reproducible from this repository:

```bash
python scripts/Figures/generate_three_way_comparison.py --output three_way_comparison.png
```

All required input CSVs are present in `comparative_analysis/`. The script loads:
- `comparative_analysis/nfw/h2_nfw_comparison_summary.csv` (74 rows)
- `comparative_analysis/mond/h2_mond_comparison_summary.csv` (74 rows)
- `comparative_analysis/referee_resolution/h2_full74_explicit_summary.csv` (filtered to 30 Tier A+B rows)

See `comparative_analysis/referee_resolution/rar_spearman_result.txt` for the Spearman computation details.

---

## Validation Products

| Path | Contents |
|------|----------|
| `comparative_analysis/referee_resolution/` | H2 tier breakdown, Spearman results, Tier C audit |
| `comparative_analysis/final_validation/` | Go/no-go memo, metric harmonization check, MOND defensibility check |
| `comparative_analysis/metric_harmonization/` | Metric choice rationale, harmonized summary CSVs |
| `comparative_analysis/nfw/` | NFW perturbation summary, diagnostic reports |
| `comparative_analysis/mond/` | RAR perturbation summary, diagnostic reports |
| `scripts/Figures/generate_three_way_comparison.py` | Main figure script (v2.0) |
| `H2_PUBLICATION_RELEASE/` | Archived H2 Tier A+B per-galaxy phase4 outputs |

---

## Publication Status

**Manuscript:** Submitted to Monthly Notices of the Royal Astronomical Society (MNRAS)
**Title:** "Inner-Region Scatter Response to Bounded Perturbations: A Comparative Analysis of 74 SPARC Galaxies"
**Author:** Figueiredo, V. M. F. (2026)
**DOI:** Pending acceptance
**Preprint:** Available upon request

---

## Citation

If you use this repository, data, or methodology, please cite:

> Figueiredo, V. M. F. (2026). H2-SPARC-Comparative-Galaxy-Dynamics (v2.0.0). Zenodo. https://doi.org/10.5281/zenodo.18793233

And the SPARC database and fit catalogues:

> Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves. AJ, 152, 157.

> Li, P., Lelli, F., McGaugh, S., & Schombert, J. (2020). Fitting the Radial Acceleration Relation to Individual SPARC Galaxies. ApJS, 247, 31.

---

## Contact

**Vítor M. F. Figueiredo**
Email: vitor.figueiredo.research@protonmail.com
ORCID: [0009-0004-7358-4622](https://orcid.org/0009-0004-7358-4622)

---

*Last updated: 2026-03-28*
