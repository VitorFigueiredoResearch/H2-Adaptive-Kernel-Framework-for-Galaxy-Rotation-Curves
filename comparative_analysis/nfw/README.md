# NFW Bounded Perturbation Analysis

## Purpose
Tests whether the scatter-neutrality pattern found in H2 (bounded kernel adaptation → Δσ ≈ 0 in inner regions) extends to standard NFW halo fits. We apply ±10%/±20% perturbations to Li et al. (2020) NFW best-fit parameters and measure the resulting change in inner-region scatter.

## Data Sources
- NFW catalog: Li, Lelli, McGaugh & Schombert (2020), ApJS 247, 31
  - File: `data/nfw/Fits/ByModel/Table/parameter_NFW_LCDM.mrt`
  - Parameters: V200 [km/s], C200 [dimensionless], rs [kpc]
  - Cosmology: H0=73 km/s/Mpc, 200×ρ_crit overdensity
- SPARC rotation curves: Lelli, McGaugh & Schombert (2016), AJ 152, 157
  - Files: `data/sparc/{GALAXY}_rotmod.dat`
- H2 fleet: 80 galaxies from fleet_summary_80galaxy.csv

## NFW Parameter Convention
- r200 [kpc] = V200 [km/s] / 0.73  (verified against catalog rs values)
- rs from catalog directly; C200 = r200/rs
- NFW: V²(r) = V200² × (r200/r) × f(r/rs) / f(C200), f(x)=ln(1+x)-x/(1+x)

## Inner Region Definition
- R < 0.5 × R_max (identical to H2 diagnostic, inner_frac=0.5)
- Scatter σ = RMS[log10(V_model) - log10(V_obs)] over inner region
- Galaxies with n_inner < 3 are excluded from analysis

## How to Run
```bash
cd [repo root]
python comparative_analysis/nfw/nfw_velocity.py          # validates NFW module
python comparative_analysis/nfw/nfw_perturbation_diagnostic.py  # main analysis
python comparative_analysis/nfw/h2_nfw_comparison.py     # comparison with H2
```

## Output Files
- `nfw_perturbation_summary.csv` — per-galaxy Δσ for all perturbations
- `diagnostic_report.txt` — coverage, exclusions, quality flags
- `h2_nfw_comparison_summary.csv` — cross-comparison table
- `h2_nfw_comparison_report.txt` — narrative results
- `figures/ngc3198_validation.png` — single-galaxy validation plot
- `figures/nfw_scatter_sensitivity.png` — main diagnostic figure

## Results
(See `h2_nfw_comparison_report.txt` after running)
