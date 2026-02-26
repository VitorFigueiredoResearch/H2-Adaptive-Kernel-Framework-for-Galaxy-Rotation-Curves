# H2 Adaptive Framework — Pre-Upload Verification Checklist

> Complete this checklist before uploading to GitHub or submitting as supplementary material.
> Generated: 2026-02-25

---

## Code Quality
- [x] No hardcoded personal paths (`C:\Users\Vitor`) in any script — verified by grep
- [x] No debugging print statements with personal notes — verified by grep
- [x] No AI conversation fragments (Claude/GPT/Gemini) in source — verified by grep
- [x] Validated pipeline scripts have professional docstrings
- [x] `requirements.txt` created and imports verified (numpy, scipy, matplotlib, pandas)
- [ ] **TODO:** Optionally add professional headers to all scripts (see USER_GUIDE §4)

## Documentation
- [x] `README.md` complete: overview, key findings, install, quick start, citations, license
- [x] `USER_GUIDE.md` uses generic `<H2_ROOT>` paths (no personal paths)
- [x] `USER_GUIDE.md` has 12/12 completeness checks passing
- [x] SPARC citation present in README.md and USER_GUIDE.md
- [x] Software dependency citations (NumPy, SciPy, Matplotlib, pandas) in README.md
- [x] `LICENSE.md` present (MIT, 2026)
- [x] `QUICKSTART.md` present (condensed 1-page reference)
- [x] `RESULTS_SUMMARY.md` present (9-galaxy findings)
- [x] `data/derived/fleet/REGIME_ANALYSIS.md` present (regime interpretation)
- [ ] **TODO:** Fill in H2 paper citation in README.md once published
- [ ] **TODO:** Fill in H1 baseline reference in README.md if/when published

## Data
- [x] 9 validation galaxies present in `data/derived/phase3/` (leff CSVs + PNGs)
- [x] 9 validation galaxies present in `data/derived/phase4/h2_outputs/` (CSVs + PNGs)
- [x] 72 H1 basis files present in `data/derived/phase4/basis_h1/` (9 × 8 L values)
- [x] Summary tables included: `10galaxy_summary.csv`, `3galaxy_summary.csv`
- [x] Key fleet plots included: `regime_trends.png`, `leff_overlay_9galaxies.png`
- [x] SPARC input data included: `data/sparc/` (175 galaxies, text files)
- [x] H1 frozen baselines included: `data/h1_frozen/per_galaxy/` (175 galaxies)
- [x] Large intermediate files excluded: `vendor/`, `__pycache__/`, `figs/`, `results/`

## Reproducibility
- [x] `requirements.txt` tested — all imports successful
- [x] NGC3198 5-step pipeline documented in README.md Quick Start
- [x] File paths are relative (no hardcoded absolute paths in physics code)
- [x] Windows UTF-8 workaround documented (`run_diagnostics_utf8.bat` included)
- [x] Known failure case documented (NGC5907 — H1 parameter incompatibility)
- [x] `FILE_MANIFEST.txt` lists all 550 included files

## Legal
- [x] SPARC citation mandatory statement in README.md
- [x] SPARC citation in USER_GUIDE.md Citation Requirements section
- [x] MIT License file present (`LICENSE.md`)
- [x] No proprietary or personal data included
- [x] No passwords, API keys, or credentials in any file

## Size Check
- [x] Total folder size: **8.3 MB** (well under 100 MB target)
- [x] Total files: **550**
- [x] All files necessary — nothing redundant except intermediate fleet CSVs
         (`fleet_new7.csv`, `fleet_orig3.csv` — safe to delete before upload if preferred)

---

## Upload Instructions (Manual GitHub Upload)

Since git authentication is unavailable, upload the `H2_PUBLICATION_RELEASE/` folder
contents manually:

1. Go to your GitHub repository page
2. Click **Add file → Upload files**
3. Drag and drop the contents of `H2_PUBLICATION_RELEASE/` (not the folder itself)
4. Commit message suggestion: `"Add 9-galaxy validation results and publication materials"`
5. Verify the file tree matches `FILE_MANIFEST.txt`

Alternatively, zip the folder and attach as supplementary material directly to the journal submission system.

---

## Final Sign-Off

| Item | Status |
|------|--------|
| Code audit complete | ✅ |
| Documentation complete | ✅ |
| Citations present | ✅ |
| 9-galaxy data present | ✅ |
| Size acceptable (8.3 MB) | ✅ |
| Ready for upload | ✅ |
