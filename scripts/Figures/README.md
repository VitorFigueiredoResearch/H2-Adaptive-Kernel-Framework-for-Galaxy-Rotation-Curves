# scripts/Figures — Figure Generation Scripts

This directory contains standalone Python scripts that generate publication
figures for the comparative paper:

*"Inner-Region Scatter Response to Bounded Perturbations: A Comparative Analysis of 74 SPARC Galaxies"*
(Figueiredo, V. M. F., 2026, MNRAS submitted)

---

## `generate_acceleration_space_plot.py`

### Purpose

Generates `acceleration_space_80galaxy.pdf`: an acceleration-space diagnostic
figure showing all radial data points from the 80-galaxy H2 validation sample
in `g_tot` vs `g_bar` space, overlaid on the McGaugh et al. (2016) radial
acceleration relation (RAR).

### What the figure shows

| Element | Description |
|---|---|
| Data points | Each point = one radial measurement from a SPARC rotmod file |
| Red points | Baryon-dominated galaxies (Vobs/Vb < 0.8, 39 galaxies) |
| Blue points | Balanced galaxies (0.8 ≤ Vobs/Vb < 1.5, 30 galaxies) |
| Green points | DM-dominated galaxies (Vobs/Vb ≥ 1.5, 11 galaxies) |
| Black curve | McGaugh et al. (2016) RAR: g†=1.2×10⁻¹⁰ m/s² |
| Dotted line | 1:1 Newtonian limit (g_tot = g_bar) |
| Dashed line | Characteristic acceleration g† |
| Annotation | H2 null result: \|Δσ\| < 10⁻¹² km/s for 74 measurable systems |

### Acceleration definitions

```
g_bar(r) = [V_gas(r)² + V_disk(r)² + V_bul(r)²] / r   [m/s²]
g_tot(r) = V_obs(r)² / r                                [m/s²]
```

Velocities from SPARC `*_rotmod.dat` files (columns: Rad, Vobs, errV, Vgas,
Vdisk, Vbul). Quality cut: Vobs/errV > 2 and all quantities positive/finite.

### McGaugh+16 RAR formula

```
g_tot = g_bar / (1 - exp(-sqrt(g_bar / g†)))
g† = 1.2 × 10⁻¹⁰ m/s²
```

### Data inputs

| File | Location |
|---|---|
| `fleet_summary_80galaxy.csv` | `H2_PUBLICATION_RELEASE/fleet_expansion_80galaxy/` |
| `*_rotmod.dat` (80 files) | `data/sparc/` |

### Usage

Run from any working directory (paths auto-resolved from script location):

```bash
# Default output: scripts/figures/acceleration_space_80galaxy.pdf
PYTHONUTF8=1 python scripts/figures/generate_acceleration_space_plot.py

# Custom output path
PYTHONUTF8=1 python scripts/figures/generate_acceleration_space_plot.py \
    --output path/to/output.pdf
```

### Output

`acceleration_space_80galaxy.pdf` — 6.5 × 5.5 inch, 300 dpi, vector PDF
suitable for journal submission.

### Dependencies

- Python ≥ 3.9
- numpy, matplotlib (standard scientific stack; see `requirements.txt`)

### Citation

If using this figure or script, cite:

> Figueiredo, V. M. F. (2026). *Inner-Region Scatter Response to Bounded
> Perturbations: A Comparative Analysis of 74 SPARC Galaxies.* MNRAS (submitted).

> McGaugh, S. S., Lelli, F., & Schombert, J. M. (2016). *Radial Acceleration
> Relation in Rotationally Supported Galaxies.* PRL, 117, 201101.

> Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). *SPARC: Mass Models
> for 175 Disk Galaxies.* AJ, 152, 157.

---

## `generate_three_way_comparison.py` (v2.0)

### Purpose

Generates `three_way_comparison.png`: the main comparative figure for the MNRAS paper showing inner-region scatter sensitivity distributions for all three model classes side by side.

### What the figure shows

| Panel | Content |
|---|---|
| (a) Scatter plot | max\|Δσ\| (dex) vs V\_bar/V\_NFW for NFW (74), RAR (74), and H2 (30 explicit) |
| (b) Box plots | Distribution summaries for NFW (N=74), RAR (N=74), H2 (N=30 explicit) |

Key annotations:
- NFW Spearman ρ = −0.899 (p < 10⁻²⁵, strong anti-correlation with baryonic dominance)
- RAR Spearman ρ = +0.049 (p = 0.68, null — not significant)
- Median reference lines per implementation
- H2 plotted as 30 explicit Tier A+B galaxies only

**H2 archive coverage:** The 44 Tier C galaxies (no archived phase4 log₁₀-dex output) are intentionally not plotted. See `comparative_analysis/referee_resolution/tier_c_audit_report.txt`.

### Data inputs

| File | Location |
|---|---|
| `h2_nfw_comparison_summary.csv` | `comparative_analysis/nfw/` |
| `h2_mond_comparison_summary.csv` | `comparative_analysis/mond/` |
| `h2_full74_explicit_summary.csv` | `comparative_analysis/referee_resolution/` |

### Usage

From the repository root:

```bash
python scripts/Figures/generate_three_way_comparison.py
# Default output: papers/H2_MNRAS/three_way_comparison.png

python scripts/Figures/generate_three_way_comparison.py --output path/to/output.png
```

### Output

`three_way_comparison.png` — 3.46 × 5.5 inch, 300 dpi, MNRAS single-column width.

### Validated medians (v2.0)

| Implementation | Median max\|Δσ\| | N |
|---|---|---|
| NFW | 0.0548 dex | 74 |
| RAR (adopted) | 0.0321 dex | 74 |
| H2 (explicit) | 0.000860 dex | 30 |

### Dependencies

- Python ≥ 3.9
- numpy, matplotlib, pandas (standard scientific stack; see `requirements.txt`)

### Citation

If using this figure or script, cite:

> Figueiredo, V. M. F. (2026). *Inner-Region Scatter Response to Bounded
> Perturbations: A Comparative Analysis of 74 SPARC Galaxies.* MNRAS (submitted).

> Li, P., Lelli, F., McGaugh, S. S., & Schombert, J. M. (2020). *A Comprehensive
> Catalog of Dark Matter Halo Models for SPARC Galaxies.* ApJS, 247, 31.

> Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). *SPARC: Mass Models
> for 175 Disk Galaxies.* AJ, 152, 157.

> Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). *SPARC: Mass Models
> for 175 Disk Galaxies.* AJ, 152, 157.
