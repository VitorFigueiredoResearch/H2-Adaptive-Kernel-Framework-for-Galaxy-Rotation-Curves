# H2 Adaptive Kernel Framework — Quick Start

> Full documentation: [USER_GUIDE.md](USER_GUIDE.md)

---

## Minimal Setup

**1. Open a terminal in the H2 root directory:**
```batch
cd C:\Users\username\Documents\GitHub\H2
```

**2. Verify both input files exist for your galaxy (replace `X`):**
```batch
dir data\sparc\X_rotmod.dat
dir data\h1_frozen\per_galaxy\rc_decomp_X_best.csv
```

**3. Set UTF-8 encoding (required on Windows for steps 2, 4, 5):**
```batch
set PYTHONUTF8=1
```
*Or just use `run_diagnostics_utf8 <module> <args>` — it sets this automatically.*

---

## 5-Step Pipeline

Run these in order, replacing `X` with the galaxy name (e.g. `NGC0891`):

```batch
set GALAXY=X
set PYTHONUTF8=1

:: Step 1 — Chi field and L_eff profile (no UTF-8 needed)
python -m diagnostics.phase3_leff_ngc3198 --galaxy %GALAXY% --alpha 2.0 --sigma_idx 1.0 --taper

:: Step 2 — Adaptive convolution (generates H2 rotation curve)
python -m diagnostics.phase4_adaptive_convolution --galaxy %GALAXY%

:: Step 3 — Rotation curve comparison plot (no UTF-8 needed)
python plot_rc_comparison.py --galaxy %GALAXY%

:: Step 4 — Test-2: chi correlation
python -m diagnostics.test2_chi_correlation --galaxy %GALAXY%

:: Step 5 — Test-3: inner scatter
python -m diagnostics.test3_inner_scatter --galaxy %GALAXY%
```

---

## What Each Step Produces

| Step | Output file(s) | What to check |
|------|---------------|---------------|
| 1 · Phase-3 | `data/derived/phase3/leff_X.csv` | Runs without error |
| 2 · Phase-4 | `data/derived/phase4/h2_outputs/rc_decomp_X_H2_adaptive.csv` | `PASS=True` in terminal |
| 3 · Plot | `data/derived/phase4/h2_outputs/rc_comparison_X.png` | File exists; curves visible |
| 4 · Test-2 | *(printed to terminal)* | Note r and N values |
| 5 · Test-3 | *(printed to terminal)* | `Delta sigma = 0.0000 dex` |

---

## Reading Your Results

After steps 4 and 5, you have three numbers:

```
Test-2 r    = Pearson correlation (chi vs ΔV)
max|ΔV|     = Peak H2−H1 velocity difference [km/s]
Test-3 Δσ  = Inner scatter change [dex]
```

Compare against the validated 3-galaxy reference:

| Galaxy | Vk/Vb | r | max\|ΔV\| | Δσ |
|--------|-------|---|----------|-----|
| NGC3198 (DM-dominated) | 1.42 | +0.572 | 2.14 km/s | 0.0000 |
| IC2574  (balanced)     | 1.00 | +0.945 | 0.15 km/s | 0.0000 |
| NGC0891 (baryon-dom.)  | 0.71 | −0.044 | 0.22 km/s | 0.0000 |

**Quick interpretation:**
- **r near +1 and your galaxy has Vk/Vb ≈ 1:** Good, expected for balanced regime
- **r near 0 and Vk/Vb < 0.8:** Good, expected for baryon-dominated regime
- **Δσ = 0.0000 dex for any galaxy:** Good — H2 is not hurting the fit
- **Δσ ≠ 0 or |r| anomalous for the regime:** Investigate (see USER_GUIDE §10)

---

## Common Issues

| Symptom | Fix |
|---------|-----|
| `ModuleNotFoundError: diagnostics` | Wrong directory — run from H2 root |
| `UnicodeEncodeError: charmap` | Add `set PYTHONUTF8=1` or use `run_diagnostics_utf8` |
| `FileNotFoundError: leff_X.csv` in Test-2/3 | Phase-3 not run yet — run Step 1 first |
| `FileNotFoundError: rc_decomp_X_H2_adaptive.csv` in Test-2/3 | Phase-4 not run yet — run Step 2 first |
| Phase-4 `PASS=False` | Re-run Phase-3 with `--taper` flag; then re-run Phase-4 |
| Phase-4 missing basis files | Run `python generate_basis_minimal.py --galaxy X --L_values 10 30 50 80 110 140 170 200` |

---

## Where to Get Help

- **Full pipeline details and argument reference:** [USER_GUIDE.md §3](USER_GUIDE.md#3-running-the-pipeline)
- **Test definitions and interpretation:** [USER_GUIDE.md §4–5](USER_GUIDE.md#4-understanding-the-tests)
- **All troubleshooting:** [USER_GUIDE.md §10](USER_GUIDE.md#10-troubleshooting)
- **Advanced: generating basis files for new galaxies:** [USER_GUIDE.md §6](USER_GUIDE.md#6-advanced-generating-new-basis-files)
- **Validated results and findings:** [RESULTS_SUMMARY.md](RESULTS_SUMMARY.md)

