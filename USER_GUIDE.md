# H2 Adaptive Kernel Framework — User Guide

> **Version:** 9-galaxy validated build · **Date:** 2026-02-24
> **Python:** 3.12 · **Platform tested:** Windows 10/11 (Anaconda), Linux

## Setup

Throughout this guide, `<H2_ROOT>` refers to your H2 installation directory.

| Platform | Example path |
|----------|-------------|
| Windows  | `C:\Users\YourName\Documents\H2-Adaptive-Framework` |
| Linux / macOS | `~/H2-Adaptive-Framework` |

All commands must be run from `<H2_ROOT>`. Relative paths like `data/sparc/...` will break if you run from a subdirectory.

## Citation Requirements

If you use this framework or data in published work, you **must** cite:

1. **SPARC Database:**
   Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016), AJ, 152, 157.
   DOI: [10.3847/0004-6256/152/6/157](https://doi.org/10.3847/0004-6256/152/6/157)

2. **H2 Framework:** 

3. **Python scientific stack:** See `README.md` for complete software citations.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Prerequisites and Setup](#2-prerequisites-and-setup)
3. [Running the Pipeline](#3-running-the-pipeline)
4. [Understanding the Tests](#4-understanding-the-tests)
5. [Interpreting Results](#5-interpreting-results)
6. [Advanced: Generating New Basis Files](#6-advanced-generating-new-basis-files)
7. [Understanding Output Files](#7-understanding-output-files)
8. [Repository Structure](#8-repository-structure)
9. [Expected Results by Galaxy Type](#9-expected-results-by-galaxy-type)
10. [Troubleshooting](#10-troubleshooting)

---

## 1. Overview

H2 is a diagnostic framework for testing **adaptive nonlocal gravity kernels** against observed galaxy rotation curves from the SPARC database. It extends a frozen H1 baseline by computing a spatially varying effective length scale L_eff(r) driven by a local baryonic structure field χ.

### What the pipeline does

```
SPARC observed data  ──┐
H1 frozen baseline   ──┤──► Phase-3 (χ + L_eff) ──► Phase-4 (adaptive V_H2)
H1 per-galaxy params ──┘                                    │
                                                            ▼
                                            Test-2 (χ–ΔV correlation)
                                            Test-3 (inner scatter Δσ)
                                            RC comparison plot
```

### What it does NOT do

- H2 does **not** fit free parameters — L0, mu, and kernel type are **frozen** from H1 per galaxy.
- H2 does **not** replace H1. It builds on top of the frozen H1 solution.
- H2 does **not** modify the SPARC data.

---

## 2. Prerequisites and Setup

### 2.1 Required files per galaxy

Before running any galaxy `X`, confirm both input files exist:

| File | Path | Description |
|------|------|-------------|
| SPARC data | `data/sparc/X_rotmod.dat` | Observed rotation curve |
| H1 frozen baseline | `data/h1_frozen/per_galaxy/rc_decomp_X_best.csv` | Frozen H1 solution |

175 galaxies currently have both files. Check with:

```bash
# From H2 root:
dir data\sparc\NGC0891_rotmod.dat
dir data\h1_frozen\per_galaxy\rc_decomp_NGC0891_best.csv
```

### 2.2 Python dependencies

```
Python 3.12
pandas  2.2.2
numpy   1.26.4
scipy   1.13.1
matplotlib 3.9.2
```

All dependencies are available in the project's Anaconda environment.

### 2.3 Working directory

**All commands must be run from the H2 root directory:**

```
C:\Users\username\Documents\GitHub\H2
```

If you run from a subdirectory, relative paths like `data/sparc/...` will fail.

### 2.4 Windows UTF-8 requirement

Several diagnostic scripts print Unicode characters (Δ, σ, etc.) which crash on Windows with the default cp1252 encoding. **Always use one of these options:**

**Option A — Environment variable (recommended for one-off runs):**
```batch
set PYTHONUTF8=1
python -m diagnostics.phase4_adaptive_convolution --galaxy NGC3198
```

**Option B — Wrapper script (recommended for regular use):**
```batch
run_diagnostics_utf8 phase4_adaptive_convolution --galaxy NGC3198
run_diagnostics_utf8 test2_chi_correlation --galaxy NGC3198
run_diagnostics_utf8 test3_inner_scatter --galaxy NGC3198
```

> `run_diagnostics_utf8.bat` is included in the H2 root. It sets `PYTHONUTF8=1` automatically.

Phase-3 and `plot_rc_comparison.py` do **not** require this workaround.

---

## 3. Running the Pipeline

The full pipeline for any galaxy `X` consists of 5 commands in order:

### Step 1 — Phase-3: χ field and L_eff profile

```bash
python -m diagnostics.phase3_leff_ngc3198 --galaxy X --alpha 2.0 --sigma_idx 1.0 --taper
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `--galaxy` | *(required)* | Galaxy name, e.g. `NGC0891` |
| `--alpha` | `1.0` | Sensitivity: L_eff = L0 / (1 + α·χ). **Use 2.0 for validated runs.** |
| `--sigma_idx` | `1.0` | Gaussian smoothing width for χ, in radial index units |
| `--taper` | off | Enforce L_eff → L0 at large r (sigmoid taper). **Always use for validated runs.** |
| `--taper_r0` | `0.70` | Sigmoid inflection point in fractional radius (0–1). canonical (do not change). |
| `--taper_k` | `80` | Sigmoid steepness. canonical (do not change). |

**Expected output:**
```
[Phase3-Leff] Galaxy: NGC0891
[Phase3-Leff] Rd_star = 2.55 kpc | L0 = 200 kpc | mu = 10 | kernel = ananta-hybrid
[Phase3-Leff] chi used: chi_smooth (sigma_idx=1)
[Phase3-Leff] taper: ON (r0=0.7, k=80)
[Phase3-Leff] Wrote CSV: data\derived\phase3\leff_NGC0891.csv
[Phase3-Leff] Wrote PNG: data\derived\phase3\chi_NGC0891.png
[Phase3-Leff] Wrote PNG: data\derived\phase3\leff_NGC0891.png
```

> **Note:** `ananta-hybrid` is the internal identifier for the power-law nonlocal kernel used in this framework. The string appears in log output only; it does not affect computed results.

---

### Step 2 — Phase-4: Adaptive convolution

```batch
set PYTHONUTF8=1
python -m diagnostics.phase4_adaptive_convolution --galaxy X
```

Or equivalently:
```batch
run_diagnostics_utf8 phase4_adaptive_convolution --galaxy X
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `--galaxy` | *(required)* | Galaxy name |
| `--no_run_h1` | off | Skip H1 re-runs; use cached basis CSVs if present (faster, but verify cache is current) |
| `--strict_radii` | off | Hard-fail if basis R grid mismatches canonical grid. Recommended for reviewer-proof runs. |

**Expected output:**
```
[Phase4] Galaxy=NGC0891
[Phase4] Canonical: L0=200 kpc | mu=10 | kernel=ananta-hybrid
[Phase4] Basis L grid (kpc): [200.0, 170.0, 140.0, 110.0, 80.0, 50.0, 30.0, 10.0]
[Phase4] Leff min/max (kpc): 66.6667 / 200
[Phase4] TEST-1 Outer Stability (r_frac >= 0.70):
[Phase4]   max |ΔV| outer = 0.0193 km/s  (tol=2)
[Phase4]   PASS=True
[Phase4] Wrote CSV: data\derived\phase4\h2_outputs\rc_decomp_NGC0891_H2_adaptive.csv
[Phase4] Wrote PNGs into: data\derived\phase4\h2_outputs
```

> **Note:** `ananta-hybrid` is the internal identifier for the power-law nonlocal kernel used in this framework. The string appears in log output only; it does not affect computed results.

**Test-1 must show `PASS=True`.** If it shows `PASS=False`, see [Troubleshooting §10.3](#issue-3-phase-4-test-1-outer-stability-fails).

---

### Step 3 — Rotation curve comparison plot

```bash
python plot_rc_comparison.py --galaxy X
```

**Expected output:**
```
✓ Loaded observed data: 18 points
✓ Loaded H1 frozen: 30 points
✓ Loaded H2 adaptive: 30 points
✓ Saved plot: data\derived\phase4\h2_outputs\rc_comparison_NGC0891.png

ΔV statistics (V_H2 - V_H1):
  max|ΔV| = 0.22 km/s
  mean ΔV = 0.07 km/s
  std ΔV  = 0.07 km/s
```

---

### Step 4 — Test-2: χ–ΔV correlation

```batch
run_diagnostics_utf8 test2_chi_correlation --galaxy X
```

**Expected output format:**
```
Galaxy: NGC0891
chi column: chi_used
Pearson r: -0.0440
N points: 21
max|ΔV| (masked): 0.2183 km/s
```

---

### Step 5 — Test-3: Inner scatter comparison

```batch
run_diagnostics_utf8 test3_inner_scatter --galaxy X
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `--galaxy` | `NGC3198` | Galaxy name |
| `--inner_frac` | `0.5` | Inner region = R < inner_frac × R_max |

**Expected output format:**
```
Galaxy: NGC0891
Inner region: R < 8.55 kpc (inner_frac=0.5) | N=8

sigma_inner(H1) = 0.1152 dex
sigma_inner(H2) = 0.1152 dex
Delta sigma (H2 - H1) = 0.0000 dex

Result: ➖ H2 is neutral within ±0.01 dex (mechanism works, leverage modest).
```

---

### Full pipeline as a single batch sequence

```batch
@rem Full H2 pipeline for galaxy X (replace X with galaxy name)
@rem Run from: <H2_ROOT>  (your H2 installation directory)

set GALAXY=NGC0891
set PYTHONUTF8=1

python -m diagnostics.phase3_leff_ngc3198 --galaxy %GALAXY% --alpha 2.0 --sigma_idx 1.0 --taper
python -m diagnostics.phase4_adaptive_convolution --galaxy %GALAXY%
python plot_rc_comparison.py --galaxy %GALAXY%
python -m diagnostics.test2_chi_correlation --galaxy %GALAXY%
python -m diagnostics.test3_inner_scatter --galaxy %GALAXY%
```

---

## 4. Understanding the Tests

### Test-1 (inside Phase-4): Outer Stability

**What it measures:** The maximum velocity difference |V_H2 − V_H1| in the outer region (r/r_max ≥ 0.70), where the taper is designed to enforce L_eff → L0.

**Why it matters:** If the adaptive kernel causes large deviations in the well-constrained outer disk, the L_eff profile is numerically unstable. The test ensures the taper is doing its job.

**Pass criterion:** max|ΔV|_outer < 2.0 km/s

**Interpretation:** Values well below 1 km/s (e.g., 0.019 km/s for NGC0891) confirm the taper is working correctly. Values near the 2 km/s limit warrant investigation of the Phase-3 L_eff profile.

---

### Test-2: χ–ΔV Correlation (Pearson r)

**What it measures:** The Pearson correlation coefficient between the χ field (from Phase-3) and the velocity residual ΔV = V_H2 − V_H1 at each radial point.

**Why it matters:** If the H2 mechanism is physically meaningful, the adaptive kernel correction should be driven by χ. A high positive r means the velocity change is coherently related to the local baryonic structure. A near-zero r means χ has no leverage — H2 is not adding information in that regime.

**Interpretation guide:**

| r value | Meaning |
|---------|---------|
| r > 0.8 | Strong χ–correction coupling. H2 is actively modifying the curve in response to structure. |
| 0.4 < r < 0.8 | Moderate coupling. H2 correction is present but partial. |
| -0.2 < r < 0.4 | Weak / no coupling. H2 is near-neutral for this galaxy. |
| r < -0.2 | Anti-correlation. Unexpected; investigate Phase-3 chi sign. |

**Output files used:** `data/derived/phase3/leff_X.csv` (χ column), `data/derived/phase4/h2_outputs/rc_decomp_X_H2_adaptive.csv` (ΔV column).

---

### Test-3: Inner Scatter Comparison (Δσ)

**What it measures:** The RMS log10 scatter of V_model/V_obs in the inner half of the radial grid for both H1 and H2, and their difference Δσ = σ_H2 − σ_H1 in dex.

**Why it matters:** A good adaptive mechanism should not degrade the rotation curve fit quality. Δσ = 0 means H2 neither improves nor worsens the inner scatter relative to H1. Positive Δσ would mean H2 is introducing scatter (a problem). Negative Δσ would mean H2 is actively improving the fit.

**Interpretation guide:**

| Δσ value | Meaning |
|----------|---------|
| Δσ = 0.0000 dex | H2 is neutral — does not change inner fit quality. Expected for all validated galaxies. |
| Δσ > +0.01 dex | H2 is degrading the inner fit. Investigate L_eff profile. |
| Δσ < -0.01 dex | H2 is improving the inner fit. Noteworthy positive result. |

**Note on N:** The number of inner points N depends on the galaxy's observed radial grid density and `--inner_frac`. It is not fixed across galaxies (NGC3198: N=31, IC2574: N=15, NGC0891: N=8).

---

## 5. Interpreting Results

### The three-number summary

For each galaxy, the key results are:

```
Test-2 r  = [Pearson correlation, −1 to +1]
Test-2 N  = [number of radial points]
max|ΔV|   = [peak H2−H1 difference, km/s]
Test-3 Δσ = [inner scatter change, dex]
```

### Regime-dependent interpretation

The galaxy's **Vk/Vb ratio** (kinematic-to-baryonic velocity) determines which regime it falls in and what to expect:

| Regime | Vk/Vb | Expected r | Expected max\|ΔV\| | Expected Δσ |
|--------|-------|-----------|-------------------|------------|
| DM-dominated | > 1.3 | 0.4–0.7 | 1–3 km/s | 0.0000 |
| Balanced | ~1.0 | 0.5–0.95 | 0.1–0.5 km/s | 0.0000 |
| Baryon-dominated | < 0.8 | −0.2 to +0.2 | 0.1–0.5 km/s | 0.0000 |

### Validated reference values (3-galaxy baseline)

| Galaxy | Vk/Vb | r | N | max\|ΔV\| | Δσ |
|--------|-------|---|---|----------|-----|
| NGC3198 | 1.42 | +0.572 | 21 | 2.14 km/s | 0.0000 |
| IC2574  | 1.00 | +0.945 | 21 | 0.15 km/s | 0.0000 |
| NGC0891 | 0.71 | −0.044 | 21 | 0.22 km/s | 0.0000 |

If a new galaxy produces Test-2 r well outside the expected range for its regime, or Δσ ≠ 0.0000, investigate before reporting.

---

## 6. Advanced: Generating New Basis Files

Phase-4 requires pre-computed H1 basis runs at 8 L values. For the 175 supported galaxies, these are generated automatically when Phase-4 first runs and cached in `data/derived/phase4/basis_h1/`.

If you need to test a new galaxy or regenerate a damaged cache:

```bash
python generate_basis_minimal.py --galaxy X --L_values 10 30 50 80 110 140 170 200
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `--galaxy` | *(required)* | Galaxy name |
| `--L_values` | *(required)* | Space-separated L values in kpc. Must include all 8 standard values. |
| `--mu` | `10.0` | Coupling parameter. Do not change — canonical (do not change). |
| `--kernel` | `ananta-hybrid` | Kernel functional form (power-law nonlocal). Fixed by reference model — do not change. |
| `--beta` | `1.15` | Kernel beta. Do not change — canonical (do not change). |

**Output:** 8 files in `data/derived/phase4/basis_h1/`:
```
rc_decomp_X_L10kpc.csv
rc_decomp_X_L30kpc.csv
...
rc_decomp_X_L200kpc.csv
```

> **Warning:** The basis files are derived from the vendored H1 runner (`vendor/h1_src/run_sparc_lite.py`). Do not modify that file.

---

## 7. Understanding Output Files

### Phase-3 outputs (`data/derived/phase3/`)

| File | Description |
|------|-------------|
| `leff_X.csv` | Per-radial-point table: R_kpc, chi_raw, chi_smooth (=chi_used), L_eff_kpc, L0_kpc, alpha, sigma_idx, taper_on, mask_outer, Rd_star_kpc, mu, kernel |
| `chi_X.png` | Plot of raw and smoothed χ vs radius |
| `leff_X.png` | Plot of L_eff(r) vs radius, showing taper effect |

**Key columns in `leff_X.csv`:**
- `chi_used`: The smoothed χ value used for L_eff computation
- `L_eff_kpc`: The adaptive kernel length scale at each radius
- `mask_outer`: Taper weight (0 = full adaptation, 1 = forced back to L0)

---

### Phase-4 outputs (`data/derived/phase4/h2_outputs/`)

| File | Description |
|------|-------------|
| `rc_decomp_X_H2_adaptive.csv` | Full rotation curve table: R_kpc, V_baryon, V_total_H1, V_total_H2, dV (=H2−H1), L_eff_used, outer_pass |
| `rc_comparison_X.png` | Three-curve plot: observed V_obs, H1 model, H2 model |
| `phase4_rc_X.png` | Phase-4 diagnostic: H1 vs H2 velocity curves |
| `phase4_deltaV_X.png` | ΔV(r) profile across the full radial range |
| `phase4_leff_X.png` | L_eff(r) as used by Phase-4 |

**Key columns in `rc_decomp_X_H2_adaptive.csv`:**
- `V_total_H1`: Frozen H1 model velocity
- `V_total_H2`: H2 adaptive model velocity
- `dV`: V_H2 − V_H1 (signed difference)
- `L_eff_used`: Effective kernel length scale applied at each radius

---

### Basis files (`data/derived/phase4/basis_h1/`)

| File | Description |
|------|----------------|
| `rc_decomp_X_L{N}kpc.csv` | H1 rotation curve computed at fixed L = N kpc. 8 files per galaxy. |

These are intermediate files used by Phase-4 to interpolate V(r, L_eff). They are generated automatically and cached.

---

### Summary outputs (`data/derived/fleet/`, `data/results_snapshot_3galaxy/`)

| File | Description |
|------|----------------|
| `data/derived/fleet/3galaxy_summary.csv` | Machine-readable results table for all 3 validated galaxies |
| `data/results_snapshot_3galaxy/` | Commit-ready snapshot: summary CSV + 3 RC comparison PNGs |

---

## 8. Repository Structure

```
H2/
├── data/
│   ├── sparc/                  # SPARC observed rotation curves (INPUT — do not modify)
│   │                           # 175 files: X_rotmod.dat
│   ├── h1_frozen/
│   │   └── per_galaxy/         # H1 frozen per-galaxy baselines (INPUT — do not modify)
│   │                           # 175 files: rc_decomp_X_best.csv
│   ├── galaxies.csv            # Galaxy structural parameters (Rd_star, Mstar, etc.)
│   └── derived/                # All generated outputs (safe to regenerate)
│       ├── phase3/             # Phase-3: chi profiles + L_eff CSVs and PNGs
│       ├── phase4/
│       │   ├── basis_h1/       # Cached H1 basis runs at 8 L values per galaxy
│       │   └── h2_outputs/     # Phase-4 H2 adaptive results + plots
│       └── fleet/              # Multi-galaxy summary tables
│
├── diagnostics/                # Diagnostic pipeline scripts (do not modify physics)
│   ├── phase3_leff_ngc3198.py  # Step 1: chi field and L_eff computation
│   ├── phase4_adaptive_convolution.py  # Step 2: adaptive V_H2 convolution
│   ├── test2_chi_correlation.py        # Step 4: chi–ΔV Pearson test
│   ├── test3_inner_scatter.py          # Step 5: inner scatter Δσ test
│   └── run_fleet.py            # Multi-galaxy fleet runner
│
├── core/                       # Core physics library (do not modify)
│   └── galaxy_io.py            # Galaxy data loading utilities
│
├── vendor/
│   └── h1_src/                 # Vendored H1 source (do not modify)
│       └── run_sparc_lite.py   # H1 rotation curve runner
│
├── kernels/                    # Kernel definitions (do not modify)
│
├── plot_rc_comparison.py       # Step 3: rotation curve comparison plot
├── generate_basis_minimal.py   # Advanced: regenerate H1 basis files
├── run_diagnostics_utf8.bat    # Windows UTF-8 wrapper for diagnostics
├── RESULTS_SUMMARY.md          # 3-galaxy validation results and findings
├── QUICKSTART.md               # Condensed pipeline guide
└── USER_GUIDE.md               # This file
```

---

## 9. Expected Results by Galaxy Type

Based on the 3-galaxy validation (NGC3198, IC2574, NGC0891). These ranges are indicative — use them as sanity checks, not hard thresholds.

### Balanced regime (Vk/Vb ≈ 1.0, e.g. IC2574)

The baryon and dark matter contributions are roughly equal. The χ field has moderate structure and L_eff varies meaningfully.

- **Test-2 r:** 0.50 to 0.95 (moderate to strong positive correlation)
- **max|ΔV|:** 0.1 to 0.5 km/s (small adaptive correction)
- **Test-3 Δσ:** 0.0000 dex (neutral)

### DM-dominated (Vk/Vb > 1.3, e.g. NGC3198)

Dark matter dominates; baryons set up a relatively smooth χ field. H2 produces larger absolute corrections but with moderate correlation.

- **Test-2 r:** 0.40 to 0.70 (moderate correlation)
- **max|ΔV|:** 1 to 3 km/s (larger magnitude corrections)
- **Test-3 Δσ:** 0.0000 dex (neutral)

### Baryon-dominated (Vk/Vb < 0.8, e.g. NGC0891)

Baryons already explain most of the rotation curve. The χ field may be large but L_eff has little room to improve the fit. H2 is nearly identical to H1.

- **Test-2 r:** −0.20 to +0.20 (weak or no correlation)
- **max|ΔV|:** 0.1 to 0.5 km/s (small magnitude)
- **Test-3 Δσ:** 0.0000 dex (neutral)

> **Important:** Δσ = 0.0000 dex is expected for all regimes. It is a property of the adaptive mechanism's design (it cannot increase inner scatter), not a coincidence.

---

## 10. Troubleshooting

### Issue 1: Wrong working directory

**Symptom:**
```
FileNotFoundError: data/sparc/NGC0891_rotmod.dat
```
or
```
ModuleNotFoundError: No module named 'diagnostics'
```

**Solution:** Always run from the H2 root:
```batch
cd C:\Users\username\Documents\GitHub\H2
python -m diagnostics.phase3_leff_ngc3198 --galaxy NGC0891
```

---

### Issue 2: UnicodeEncodeError (Windows)

**Symptom:**
```
UnicodeEncodeError: 'charmap' codec can't encode character '\u0394'
```

**Affects:** `phase4_adaptive_convolution.py`, `test2_chi_correlation.py`, `test3_inner_scatter.py`

**Solution:** Set UTF-8 mode before running:
```batch
set PYTHONUTF8=1
python -m diagnostics.phase4_adaptive_convolution --galaxy NGC0891
```
Or use the wrapper:
```batch
run_diagnostics_utf8 phase4_adaptive_convolution --galaxy NGC0891
```

---

### Issue 3: Phase-4 missing basis files

**Symptom:**
```
[Phase4] Missing basis file for L=10 kpc
```
or `FileNotFoundError` pointing to `data/derived/phase4/basis_h1/`.

**Cause:** Phase-4 generates basis files on first run. If the cache is missing or corrupted, they need to be regenerated.

**Solution:**
```bash
python generate_basis_minimal.py --galaxy X --L_values 10 30 50 80 110 140 170 200
```

Then re-run Phase-4. Phase-4 can also regenerate them itself if `--no_run_h1` is **not** set (which is the default).

---

### Issue 4: Test-2 returns identical results for different galaxies

**Symptom:** All galaxies show `r=0.572, N=21` regardless of which `--galaxy` is passed.

**Cause:** Old bug (fixed 2026-02-24) where the galaxy name was hardcoded instead of read from `args`.

**Verify the fix is present:**
```bash
python -c "import inspect, diagnostics.test2_chi_correlation as t; print(inspect.getsource(t.main))" | grep "args ="
```
Should show `args = parser.parse_args()`, not a hardcoded string.

If you have the old version, update from the repository.

---

### Issue 5: Phase-4 Test-1 Outer Stability FAIL

**Symptom:**
```
[Phase4]   max |ΔV| outer = 3.47 km/s  (tol=2)
[Phase4]   PASS=False
```

**Cause:** The taper is not sufficiently suppressing L_eff adaptation in the outer region. Most common causes:
1. Phase-3 was run **without** `--taper`
2. `--taper_r0` was changed from the default 0.70

**Solution:**
1. Re-run Phase-3 with `--taper` flag included
2. Confirm `--taper_r0` is at default (0.70) — this is a canonical (do not change) parameter
3. Re-run Phase-4

---

### Issue 6: Phase-3 "galaxy not found in galaxies.csv"

**Symptom:**
```
KeyError: Galaxy 'NGC9999' not in galaxies.csv
```

**Cause:** The galaxy name is not in `data/galaxies.csv`, which lists all 175 supported galaxies.

**Solution:** Check the exact name used in SPARC (names are case-sensitive and must match the `_rotmod.dat` filename prefix). List all supported galaxies:
```bash
python -c "import pandas as pd; print(pd.read_csv('data/galaxies.csv')['name'].tolist())"
```

---

### Issue 7: Phase-4 radii mismatch with `--strict_radii`

**Symptom:**
```
RuntimeError: R_kpc mismatch: basis n=29 vs canonical n=30
```

**Cause:** A cached basis file has a different radial grid than the canonical H1 frozen baseline. This can happen if a basis file was generated with a different version of the H1 runner.

**Solution:** Delete the affected basis files and re-run Phase-4 (without `--no_run_h1`):
```bash
del data\derived\phase4\basis_h1\rc_decomp_X_*.csv
python -m diagnostics.phase4_adaptive_convolution --galaxy X
```

---

## Known Limitations

1. **`TARGET_GALAXY` vendor patch:** `vendor/h1_src/run_sparc_lite.py` contains a hardcoded `TARGET_GALAXY = "NGC3198"` debug variable. Phase-4 neutralises this at runtime (`mod.TARGET_GALAXY = None`), but the dependency is fragile if the vendor file is updated.

2. **Windows UTF-8 encoding:** Print statements containing Δ/σ characters fail without `PYTHONUTF8=1`. This is a display-only issue — no computation is affected.

3. **3-galaxy validation scope:** The regime-dependent r trend is established on 3 galaxies (one per regime). A fleet run on ≥20 galaxies is needed to confirm the correlation slope and scatter statistically.
