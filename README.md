# H2 Adaptive Kernel Framework for Galaxy Rotation Curves

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18793233.svg)](https://doi.org/10.5281/zenodo.18793233)

State-driven adaptive nonlocal gravity framework with regime-dependent kernel deformation.

**Status:** Companion code for paper submission · 9-galaxy validation complete

---

## Overview

H2 extends the H1 frozen nonlocal kernel framework by introducing adaptive kernel length
scales driven by local baryonic density gradients (the χ field). The kernel scale L_eff(r)
varies radially in response to the local baryon-to-dark-matter balance, producing
regime-dependent velocity corrections while preserving the H1 baseline in well-fit regions.

This repository contains:

- Complete implementation of the H2 adaptive kernel framework
- Full diagnostic pipeline (Test-1 outer stability, Test-2 χ correlation, Test-3 scatter)
- 9-galaxy validation results spanning three dynamical regimes
- Reproducible workflow documentation

---

## Key Findings

Across 9 SPARC galaxies spanning baryon-dominated to DM-dominated regimes
(Vk/Vb = 0.71 to 1.86):

| Finding | Value |
|---------|-------|
| **Universal inner-scatter neutrality** | Δσ = 0.0000 dex (all 9 galaxies) |
| **Regime-dependent χ–ΔV correlation** | r ranges from −0.37 to +0.95 |
| **Geometric L_eff floor** | min L_eff = L₀/3 ≈ 66.67 kpc (regime-independent) |
| **Baryon-dominated regime** | max\|ΔV\| ≤ 0.43 km/s — H2 ≈ H1 as expected |
| **DM-dominated regime** | max\|ΔV\| ≈ 2.5 km/s — active adaptive correction |

See `data/derived/fleet/REGIME_ANALYSIS.md` and `RESULTS_SUMMARY.md` for full analysis.

---

## Installation

### Prerequisites

- Python 3.8 or higher
- Required packages: `numpy`, `scipy`, `matplotlib`, `pandas`

### Setup

```bash
git clone https://github.com/VitorFigueiredoResearch/H2-Adaptive-Kernel-Framework-for-Galaxy-Rotation-Curves
cd H2-Adaptive-Framework
pip install -r requirements.txt
```

**Windows users:** See `USER_GUIDE.md §2.4` for the UTF-8 encoding requirement.

---

## Quick Start

Reproduce the NGC3198 validation results:

```bash
# Step 1 — Chi field and L_eff profile
python -m diagnostics.phase3_leff_ngc3198 --galaxy NGC3198 --alpha 2.0 --sigma_idx 1.0 --taper

# Step 2 — Adaptive convolution (Windows: set PYTHONUTF8=1 first)
python -m diagnostics.phase4_adaptive_convolution --galaxy NGC3198

# Step 3 — Rotation curve comparison plot
python plot_rc_comparison.py --galaxy NGC3198

# Step 4 — Test-2: chi-ΔV correlation
python -m diagnostics.test2_chi_correlation --galaxy NGC3198

# Step 5 — Test-3: inner scatter
python -m diagnostics.test3_inner_scatter --galaxy NGC3198
```

For the full 9-galaxy fleet:

```bash
set PYTHONUTF8=1
python -m diagnostics.run_fleet \
  --galaxies NGC0891,NGC6946,NGC5585,IC2574,NGC3893,NGC2903,NGC5055,NGC3198,NGC6503 \
  --alpha 2.0 --sigma_idx 1.0 --taper
```

See `USER_GUIDE.md` for complete workflow documentation and troubleshooting.
See `QUICKSTART.md` for a condensed 1-page pipeline reference.

---

## Repository Structure

```
H2/
├── core/                        # Core physics library
│   ├── chi.py                   # χ field computation
│   ├── leff.py                  # L_eff adaptive profile
│   ├── galaxy_io.py             # SPARC data loading
│   └── ...
├── kernels/                     # Kernel implementations
├── diagnostics/                 # Diagnostic pipeline
│   ├── phase3_leff_ngc3198.py   # Step 1: χ and L_eff
│   ├── phase4_adaptive_convolution.py  # Step 2: H2 rotation curve
│   ├── test2_chi_correlation.py        # Test-2: χ–ΔV Pearson r
│   ├── test3_inner_scatter.py          # Test-3: inner Δσ
│   └── run_fleet.py             # Multi-galaxy runner
├── tests/                       # Unit tests
├── data/
│   ├── sparc/                   # SPARC rotation curves (175 galaxies, input)
│   ├── h1_frozen/per_galaxy/    # H1 frozen baselines (input)
│   ├── galaxies.csv             # Galaxy structural parameters
│   └── derived/
│       ├── phase3/              # L_eff profiles (generated)
│       ├── phase4/              # H2 adaptive results (generated)
│       └── fleet/               # Summary tables and trend plots
├── paper_materials/             # LaTeX-ready tables
├── plot_rc_comparison.py        # RC comparison plot generator
├── generate_basis_minimal.py    # H1 basis file generator
├── run_diagnostics_utf8.bat     # Windows UTF-8 wrapper
├── requirements.txt
├── USER_GUIDE.md                # Full workflow documentation
├── QUICKSTART.md                # Condensed pipeline reference
├── RESULTS_SUMMARY.md           # 9-galaxy validation results
└── LICENSE.md
```

---

## Data Sources and Citations

### SPARC Database (mandatory citation)

This work uses rotation curve data from the Spitzer Photometry and Accurate Rotation
Curves (SPARC) database:

> **Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016)**
> *SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves*
> The Astronomical Journal, 152, 157.
> DOI: [10.3847/0004-6256/152/6/157](https://doi.org/10.3847/0004-6256/152/6/157)
> Data: [astroweb.cwru.edu/SPARC](http://astroweb.cwru.edu/SPARC/)

**If you use this code or the included data, you must cite the SPARC paper above.**

### H1 Nonlocal Baseline

The H1 frozen kernel parameters and baselines were developed as part of the H1 nonlocal
gravity framework. [Add reference if/when published.]

---

## Software Dependencies

This code uses the following scientific Python libraries. Please cite them in published work:

- **NumPy:** Harris, C. R., et al. (2020). Array programming with NumPy. *Nature*, 585, 357–362.
  DOI: [10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)

- **SciPy:** Virtanen, P., et al. (2020). SciPy 1.0: Fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272.
  DOI: [10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2)

- **Matplotlib:** Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, 9, 90–95.
  DOI: [10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55)

- **pandas:** McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 51–56.

---

## License

MIT License. See `LICENSE.md` for details.

---

## Citation

If you use this framework in published work, please cite:

> Figueiredo, V. M. F. (2026). H2 Adaptive Kernel Framework for Galaxy Rotation Curves (v1.0). 
> Zenodo. https://doi.org/10.5281/zenodo.18793233

And the SPARC database (mandatory — see above).

---

## Contact

Questions or issues? Open a GitHub issue or contact:
[vitor.figueiredo.research@protonmail.com]


