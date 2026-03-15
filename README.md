# H2 Adaptive Kernel Framework for Galaxy Rotation Curves

![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18911959.svg)](https://doi.org/10.5281/zenodo.18911959)
![License](https://img.shields.io/github/license/VitorFigueiredoResearch/H2-Adaptive-Kernel-Framework-for-Galaxy-Rotation-Curves)


Adaptive nonlocal kernel framework for galaxy rotation curves, with bounded local scale modulation and an 80-galaxy SPARC validation sample.

**Status:** Final paper-support repository · 80-galaxy validation complete

---
## Paper Status
**Currently in peer review** at *The Astrophysical Journal*  
- **Manuscript ID:** AAS74932  
- **Submitted:** March 11, 2026  
- **Status:** Peer review (as of March 13, 2026)
  
## Overview

H2 is a controlled extension of the frozen H1 nonlocal kernel framework. It introduces a locally adaptive effective kernel scale, $L_{\mathrm{eff}}(r)$, modulated by the baryonic gradient field $\chi(r)$, while preserving the canonical H1 baseline parameters:

- $L_0 = 200$ kpc
- $\mu = 10$
- $\beta = 1.15$

The purpose of H2 is not to propose a replacement cosmology, but to test whether bounded local kernel adaptation can reduce the inner-region scatter metric relative to the frozen H1 baseline.

This repository contains:

- The full H2 implementation
- The diagnostic pipeline used in the manuscript
- The 80-galaxy validation sample
- Reproducible figures, tables, and summary outputs
- Manuscript-support materials for the H2 paper

---

## Main Result
The final validation sample contains 80 SPARC galaxies...

[H2 Acceleration Space Results](scripts/Figures/acceleration_space_80galaxy.pdf)
*Figure 1: Acceleration space distribution for the 80-galaxy validation fleet, showing the transition from baryon-dominated to dark-matter-dominated regimes.*
The final validation sample contains **80 SPARC galaxies** from a **canonical-baseline-compatible subset** of the SPARC database:

- **39 baryon-dominated** systems
- **30 balanced** systems
- **11 dark-matter-dominated** systems

**Regime proxy range:** $V_{\mathrm{obs}}/V_b = 0.027$ to $2.187$

**Primary findings:**

- $\Delta\sigma$ is **consistent with zero** for all galaxies with sufficient inner-region coverage (74/74 measured systems)
- **6 galaxies** show undefined $\Delta\sigma$ due to insufficient inner-region sampling (data geometry limitation, not model failure)
- Pearson correlation $r$ spans **−0.371 to +0.981** (regime-dependent adaptive response)
- All validated systems reach or approach an **empirical operating floor** near $L_0/3 \approx 66.67$ kpc

This repository documents that bounded local adaptation produces regime-dependent local response while leaving the adopted inner-region scatter metric unchanged within the validated sample.

---

## Validation Products

Key outputs from the 80-galaxy validation:

- `H2_PUBLICATION_RELEASE/fleet_expansion_80galaxy/fleet_summary_80galaxy.csv`
- `H2_PUBLICATION_RELEASE/fleet_expansion_80galaxy/figures/regime_trends_80galaxy.png`
- `H2_PUBLICATION_RELEASE/fleet_expansion_80galaxy/figures/leff_overlay_80galaxy.png`
- `H2_PUBLICATION_RELEASE/fleet_expansion_80galaxy/README.md`

These files summarize the final 80-galaxy validation state used in the manuscript.

---

## Installation

### Requirements

- Python 3.8 or higher
- `numpy`
- `scipy`
- `matplotlib`
- `pandas`

### Setup
```bash
git clone https://github.com/VitorFigueiredoResearch/H2-Adaptive-Kernel-Framework-for-Galaxy-Rotation-Curves
cd H2-Adaptive-Kernel-Framework-for-Galaxy-Rotation-Curves
pip install -r requirements.txt
```

**Windows users:** UTF-8 encoding must be enabled for command-line execution. See `USER_GUIDE.md` for details or use the provided `run_diagnostics_utf8.bat` wrapper script.

---

## Quick Start

### Single-Galaxy Diagnostic Example

Reproduce NGC3198 validation results:
```bash
# Step 1: Chi field and L_eff profile
python -m diagnostics.phase3_leff_ngc3198 --galaxy NGC3198 --alpha 2.0 --sigma_idx 1.0 --taper

# Step 2: Adaptive convolution
python -m diagnostics.phase4_adaptive_convolution --galaxy NGC3198

# Step 3: Rotation curve comparison plot
python plot_rc_comparison.py --galaxy NGC3198

# Step 4: Test-2 (χ–ΔV correlation)
python -m diagnostics.test2_chi_correlation --galaxy NGC3198

# Step 5: Test-3 (inner scatter)
python -m diagnostics.test3_inner_scatter --galaxy NGC3198
```

### Fleet Execution

Process the original 9-galaxy baseline validation:
```bash
python -m diagnostics.run_fleet \
  --galaxies NGC0891,NGC6946,NGC5585,IC2574,NGC3893,NGC2903,NGC5055,NGC3198,NGC6503 \
  --alpha 2.0 --sigma_idx 1.0 --taper
```

For the complete workflow documentation, see:

- `USER_GUIDE.md` — Full pipeline documentation
- `QUICKSTART.md` — Condensed reference guide

---

## Repository Structure
```text
H2/
├── core/                         # Core physics and profile logic
│   ├── chi.py                    # χ field computation
│   ├── leff.py                   # L_eff adaptive profile
│   ├── galaxy_io.py              # SPARC data loading
│   └── ...
├── kernels/                      # Kernel implementation
├── diagnostics/                  # Phase 3, Phase 4, Test-2, Test-3, fleet tools
│   ├── phase3_leff_ngc3198.py
│   ├── phase4_adaptive_convolution.py
│   ├── test2_chi_correlation.py
│   ├── test3_inner_scatter.py
│   └── run_fleet.py
├── tests/                        # Unit tests
├── data/
│   ├── sparc/                    # SPARC input rotation curves
│   ├── h1_frozen/per_galaxy/     # Frozen H1 baselines
│   ├── galaxies.csv              # Galaxy structural parameters
│   └── derived/                  # Generated outputs
│       ├── phase3/               # L_eff profiles
│       ├── phase4/               # H2 adaptive results
│       └── fleet/                # Fleet summary outputs
├── H2_PUBLICATION_RELEASE/
│   ├── fleet_expansion_30galaxy/ # 30-galaxy validation
│   ├── fleet_expansion_50galaxy/ # 50-galaxy expansion
│   └── fleet_expansion_80galaxy/ # 80-galaxy final validation
├── papers/H2/h2/                 # Manuscript LaTeX materials
├── plot_rc_comparison.py
├── generate_basis_minimal.py
├── run_diagnostics_utf8.bat      # Windows UTF-8 wrapper
├── requirements.txt
├── USER_GUIDE.md
├── QUICKSTART.md
├── RESULTS_SUMMARY.md
├── CITATION.cff
└── LICENSE.md
```

---

## Scientific Scope

This repository supports a diagnostic paper about the structural limits of bounded adaptive nonlocal kernels in galaxy dynamics.

**It does NOT claim:**

- A replacement for ΛCDM
- A MOND alternative
- Elimination of dark matter
- A general theory of gravity

**It DOES provide:**

- A reproducible adaptive-kernel implementation
- A controlled 80-galaxy validation sample from the SPARC database
- A documented constraint on the adopted inner-region scatter metric
- Evidence of regime-dependent adaptive response with invariant scatter

---

## Data Sources and Citations

### SPARC Database (Mandatory Citation)

This work uses rotation-curve data from the Spitzer Photometry and Accurate Rotation Curves (SPARC) database:

> **Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016)**  
> *SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves*  
> The Astronomical Journal, 152, 157  
> DOI: [10.3847/0004-6256/152/6/157](https://doi.org/10.3847/0004-6256/152/6/157)  
> Data: [http://astroweb.cwru.edu/SPARC/](http://astroweb.cwru.edu/SPARC/)

**If you use this code or the included SPARC data, you MUST cite the SPARC paper above.**

---

## Software Dependencies

This code relies on the standard scientific Python ecosystem:

- **NumPy:** Harris, C. R., et al. (2020). Array programming with NumPy. *Nature*, 585, 357–362. DOI: [10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)

- **SciPy:** Virtanen, P., et al. (2020). SciPy 1.0: Fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272. DOI: [10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2)

- **Matplotlib:** Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, 9, 90–95. DOI: [10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55)

- **pandas:** McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 51–56.

See the manuscript acknowledgments and bibliography for complete software citations.

---

## License

MIT License. See `LICENSE.md` for details.

---

## Citation

If you use this repository in published work, please cite:

> **Figueiredo, V. M. F. (2026)**  
> *H2 Adaptive Kernel Framework for Galaxy Rotation Curves* (v1.0)  
> Zenodo  
> https://doi.org/10.5281/zenodo.18793233

And the SPARC database (mandatory — see above).

You may also use the repository citation metadata in `CITATION.cff`.

---

## Contact

For questions related to the repository or manuscript:

**Vítor M. F. Figueiredo**  
Email: vitor.figueiredo.research@protonmail.com  
ORCID: [0009-0004-7358-4622](https://orcid.org/0009-0004-7358-4622)

---

## Acknowledgments

This research made use of:
- The SPARC galaxy rotation curve database
- The Python scientific computing ecosystem (NumPy, SciPy, Matplotlib, pandas)
- Large language models (Anthropic Claude, OpenAI GPT, Google Gemini) for code debugging, documentation structuring, and manuscript preparation assistance

All scientific framework design, data analysis, diagnostic methodology, and research conclusions are the original work of the author.

---




