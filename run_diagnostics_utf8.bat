@echo off
:: run_diagnostics_utf8.bat
:: Wrapper to run H2 diagnostic modules on Windows without Unicode encoding errors.
:: Usage: run_diagnostics_utf8 <module_name> [args...]
:: Example: run_diagnostics_utf8 phase3_leff_ngc3198 --galaxy NGC3198 --alpha 2.0 --sigma_idx 1.0 --taper
:: Example: run_diagnostics_utf8 phase4_adaptive_convolution --galaxy NGC3198
:: Example: run_diagnostics_utf8 test2_chi_correlation --galaxy NGC3198
:: Example: run_diagnostics_utf8 test3_inner_scatter --galaxy NGC3198
set PYTHONUTF8=1
python -m diagnostics.%*
