# SPARC-Comparative-Galaxy-Dynamics
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18793232.svg)](https://doi.org/10.5281/zenodo.18793232)

Comparative galaxy-dynamics repository built around the SPARC database, with harmonized inner-region scatter diagnostics across three bounded model implementations:

- **H2 adaptive nonlocal kernel**
- **NFW dark-matter halo perturbations**
- **MOND/RAR bounded sensitivity tests**

This repository evolved from the original H2 adaptive-kernel framework into a broader comparative study of **inner-region scatter sensitivity** in galaxy rotation-curve modelling.

---

## Repository Status

**Manuscript:** *"Inner-Region Scatter Response to Bounded Perturbations: A Comparative Analysis of 74 SPARC Galaxies"* — submitted to Monthly Notices of the Royal Astronomical Society (MNRAS).

**Current scientific status:** comparative analysis complete; manuscript under review.

Validated comparison classes:

- **H2 adaptive kernel**
- **NFW halo perturbation analysis**
- **Adopted RAR (Radial Acceleration Relation) implementation**

Common validated comparison set:

- **74 SPARC galaxies** (from 80-galaxy H2 fleet; 6 excluded for insufficient inner-region coverage)
- same inner-region definition (R < 0.5 × R\_max)
- same harmonized log₁₀-RMS scatter metric (dex)

This repository supports a comparative manuscript studying how different bounded model constructions respond to perturbations in the inner regions of galaxy rotation curves.

---

## Scientific Scope

This repository is **not** a general theory-of-gravity claim and is **not** presented as a replacement for ΛCDM or MOND.

It is a **comparative structural diagnostics repository** asking a narrower question:

> How sensitive is the adopted inner-region scatter metric to bounded perturbations across different galaxy rotation-curve model classes?

The current comparative result is:

- **H2** shows near-neutral inner-region scatter response
- **MOND/RAR** shows intermediate scatter sensitivity
- **NFW** shows larger scatter sensitivity

All three are evaluated with:
- the same SPARC galaxy basis
- the same inner-region definition
- the same harmonized log10-residual RMS metric

---

## Main Comparative Result

For the current common 74-galaxy comparison set:

- **H2:** median \(|\Delta\sigma| \approx 0.0009\) dex (30 explicit archived galaxies; Tier A+B)
- **MOND/RAR:** median max\(|\Delta\sigma| \approx 0.0321\) dex (N=74)
- **NFW:** median max\(|\Delta\sigma| \approx 0.0548\) dex (N=74)

Interpretation:

- the H2 bounded kernel implementation remains near-neutral in the adopted inner-region scatter diagnostic
- MOND/RAR bounded sensitivity tests do not reproduce full H2-like neutrality
- bounded NFW perturbations produce the largest scatter response among the three tested implementations
- in the NFW case, scatter sensitivity tracks baryonic dominance strongly

This repository documents those comparisons in a reproducible form.

---

## Repository Components

### 1. H2 validation program
The original H2 program tested a bounded adaptive nonlocal kernel on an 80-galaxy SPARC fleet.

Key H2 products include:
- 80-galaxy fleet summary
- regime diagnostics
- \(L_{\mathrm{eff}}\) profiles
- scatter-neutrality diagnostic outputs

### 2. NFW comparative analysis
Located in:

```text
comparative_analysis/nfw/
````

Key outputs:

* `nfw_perturbation_summary.csv`
* `diagnostic_report.txt`
* `h2_nfw_comparison_summary.csv`
* `h2_nfw_comparison_report.txt`
* `figures/ngc3198_validation.png`
* `figures/nfw_scatter_sensitivity.png`

This module tests whether bounded perturbations of NFW halo parameters reproduce the same inner-region scatter neutrality pattern seen in H2.

### 3. MOND/RAR comparative analysis

Located in:

```text
comparative_analysis/mond/
```

Key outputs:

* `data_inspection_report.txt`
* `metric_and_region_report.txt`
* `perturbation_scheme_report.txt`
* `mond_velocity.py`
* `mond_perturbation_diagnostic.py`
* `mond_perturbation_summary.csv`
* `diagnostic_report.txt`
* `h2_mond_comparison.py`
* `h2_mond_comparison_summary.csv`
* `h2_mond_comparison_report.txt`
* `figures/validation_NGC3198.png`
* `figures/mond_scatter_sensitivity.png`

This module tests whether bounded RAR/MOND-side perturbations reproduce the H2 null-response pattern in the same inner-region metric.

---

## Data Sources

This repository uses or compares against the following public data sources:

### SPARC database

Lelli, McGaugh & Schombert (2016)

* galaxy rotation curves
* baryonic decomposition
* structural galaxy information

### Li et al. dark-halo and RAR fit catalogues

Li, Lelli, McGaugh & Schombert (2020)

Used here for:

* NFW halo parameter tables
* RAR/MOND-equivalent fit tables

---

## Comparative Analysis Design

All comparative analyses are built around the same operational diagnostic choices wherever possible:

* **same galaxy sample logic**
* **same inner-region definition**
* **same scatter metric**
* **same bounded perturbation philosophy**
* **same SPARC observational basis**

### Harmonized scatter metric

The comparison uses:

* RMS of `log10(V_model) - log10(V_obs)`
* reported in **dex**

### Inner region

The inner region is defined operationally as:

```text
R < 0.5 × R_max
```

for direct comparability across H2, NFW, and MOND/RAR analyses.

---

## Installation

### Requirements

* Python 3.8+
* numpy
* scipy
* matplotlib
* pandas

### Setup

```bash
git clone https://github.com/VitorFigueiredoResearch/SPARC-Comparative-Galaxy-Dynamics
cd SPARC-Comparative-Galaxy-Dynamics
pip install -r requirements.txt
```

Windows users may need UTF-8 enabled for command-line execution.

---

## Quick Start

### Reproduce original H2 diagnostic example

```bash
python -m diagnostics.phase3_leff_ngc3198 --galaxy NGC3198 --alpha 2.0 --sigma_idx 1.0 --taper
python -m diagnostics.phase4_adaptive_convolution --galaxy NGC3198
python plot_rc_comparison.py --galaxy NGC3198
python -m diagnostics.test2_chi_correlation --galaxy NGC3198
python -m diagnostics.test3_inner_scatter --galaxy NGC3198
```

### Run NFW comparative analysis

```bash
python comparative_analysis/nfw/nfw_velocity.py
python comparative_analysis/nfw/nfw_perturbation_diagnostic.py
python comparative_analysis/nfw/h2_nfw_comparison.py
```

### Run MOND comparative analysis

```bash
python comparative_analysis/mond/mond_velocity.py
python comparative_analysis/mond/mond_perturbation_diagnostic.py
python comparative_analysis/mond/h2_mond_comparison.py
```

---

## Repository Structure

```text
SPARC-Comparative-Galaxy-Dynamics/
├── comparative_analysis/
│   ├── nfw/
│   │   ├── README.md
│   │   ├── nfw_velocity.py
│   │   ├── nfw_perturbation_diagnostic.py
│   │   ├── nfw_perturbation_summary.csv
│   │   ├── diagnostic_report.txt
│   │   ├── h2_nfw_comparison.py
│   │   ├── h2_nfw_comparison_summary.csv
│   │   ├── h2_nfw_comparison_report.txt
│   │   └── figures/
│   │       ├── ngc3198_validation.png
│   │       └── nfw_scatter_sensitivity.png
│   └── mond/
│       ├── data_inspection_report.txt
│       ├── metric_and_region_report.txt
│       ├── perturbation_scheme_report.txt
│       ├── mond_velocity.py
│       ├── mond_perturbation_diagnostic.py
│       ├── mond_perturbation_summary.csv
│       ├── diagnostic_report.txt
│       ├── h2_mond_comparison.py
│       ├── h2_mond_comparison_summary.csv
│       ├── h2_mond_comparison_report.txt
│       └── figures/
│           ├── validation_NGC3198.png
│           └── mond_scatter_sensitivity.png
├── data/
│   ├── sparc/
│   ├── h1_frozen/
│   ├── nfw/
│   │   └── Fits/
│   │       ├── ByGalaxy/
│   │       └── ByModel/
│   └── mond/
├── diagnostics/
├── H2_PUBLICATION_RELEASE/
├── paper/ or papers/
├── requirements.txt
├── USER_GUIDE.md
├── QUICKSTART.md
├── CITATION.cff
└── LICENSE.md
```

---

## Reproducibility

The comparative outputs are designed to be:

* script-generated
* path-traceable
* based on public source catalogues where possible
* transparent about exclusions and assumptions

Important reproducibility notes:

* the common comparison set uses **74 galaxies**
* **6 galaxies** are excluded consistently where inner-region coverage is insufficient for stable scatter evaluation
* metric harmonization and exclusion consistency are documented in the comparison-specific reports

**H2 archive coverage note:** The H2 explicit scatter metric is available for **30 of the 74 comparison galaxies** (Tier A: 9 pilot galaxies; Tier B: 21 expansion galaxies). The remaining 44 galaxies (Tier C) are valid H2 fleet members but have no archived log₁₀-dex phase4 output. Their fleet summaries record Δσ = 0 km/s in velocity units, but this cannot be converted to log₁₀-dex without re-running the H1 basis pipeline. The 44 Tier C galaxies are not treated as near-zero values — they are simply not plotted in H2 scatter comparisons. Full documentation in `comparative_analysis/comparative_validation/tier_c_audit_report.txt`.

---

## What This Repository Does Not Claim

This repository does **not** claim:

* a proof against ΛCDM
* a proof against MOND
* a universal law of galaxy dynamics
* a final theory of gravity

It provides:

* a reproducible comparative testbed
* a harmonized inner-region scatter diagnostic
* side-by-side bounded perturbation comparisons across three implementations

---

## Citation and Use

If you use this repository, please cite:

1. the SPARC database
2. the Li et al. fit catalogues used for the NFW and RAR comparisons
3. the repository archive / DOI once updated for the comparative version

You may also use the metadata in `CITATION.cff`.

---

## License

MIT License. See `LICENSE.md`.

---

## Contact

**Vítor M. F. Figueiredo**
Independent Researcher
ORCID: 0009-0004-7358-4622

```


