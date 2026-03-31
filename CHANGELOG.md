# Changelog

## Version 2.2.2 (March 2026)

### Added
- `tools/` — reproducibility utilities and documentation for bounded workflow verification
- `tools/h2_diagnostic_tool.py` — standalone inner-region scatter metric verification utility for compatible model/data inputs
- `tools/TOOL_FEASIBILITY_REPORT.md` — feasibility audit defining the bounded scope of the public diagnostic utility
- `tools/TOOL_INTEGRATION_REPORT.md` — tooling architecture report documenting the layered reproducibility structure
- `tools/TOOL_SAFETY_HARDENING_REPORT.md` — release-safety audit, README hardening, and legacy-script warning status
- `tools/sparc_extractor/` — SPARC/H1 input reconstruction utility with validator and documentation for upstream reproducibility
- updated `comparative_analysis/mond/figures/rar_scatter_sensitivity.png` with manuscript-synchronised RAR-facing labels and fixed H2 reference level
- updated three-way comparison assets reflecting full-sample H2 coverage on the common 74-galaxy basis

### Changed
- manuscript synchronized to the final comparative MNRAS framing:
  *Inner-Region Scatter Response to Bounded Perturbations: A Comparative Analysis of 74 SPARC Galaxies*
- H2 full-sample median inner-region response updated to  
  `|Δσ|_H2 = 0.000478 dex` on the common 74-galaxy comparison sample
- comparative text, captions, tables, and figures updated to remove legacy 30-galaxy subset wording
- Figure 3 updated to show NFW, adopted RAR, and H2 on a symmetric `N = 74` footing
- repository documentation revised to distinguish:
  - upstream SPARC/H1 input reconstruction
  - downstream H2 metric verification
  - the non-portable internal H1+H2 middle pipeline

### Fixed
- removed stale figure annotations and manuscript references tied to the earlier 30-galaxy archived H2 subset
- corrected RAR figure wording to use adopted RAR language instead of legacy MOND-facing labels
- added warning shielding for legacy `test3_inner_scatter.py` so it is not mistaken for the authoritative manuscript path
- hardened public-facing tool documentation to clarify assumptions, limitations, and release-safe usage boundaries

### Notes
- The public tools are provided as reproducibility and diagnostic utilities only. They do not constitute a general external-data implementation of the full H1+H2 pipeline.
- Release of bundled SPARC/CDS source tables remains subject to verification of upstream data-usage and redistribution terms.

---

### Note on v1.1.3
The ApJ submission phase (Manuscript #AAS74932) preceded the current comparative reframing. The earlier H2_ApJ_SUBMISSION_FINAL.tex manuscript is preserved in `paper/` for archival reference.

---

## Version 1.1.3 (March 15, 2026)

### Added
- `scripts/figures/` folder with acceleration-space plot generation script
- Updated manuscript to ApJ submission version (`H2_ApJ_FINAL_READY.tex`)
- Added peer review status to README

### Changed
- Paper status: Now in peer review at ApJ (Manuscript #AAS74932)
- Removed outdated A&A manuscript version

### Fixed
- Complete reproducibility for all figures in paper

---

## Version 1.0.1 (February 2026)

### Initial Release
- Complete H2 validation pipeline
- 80-galaxy SPARC sample analysis
- All data and pre-computed outputs included
- Comprehensive documentation (README, QUICKSTART, USER_GUIDE)
