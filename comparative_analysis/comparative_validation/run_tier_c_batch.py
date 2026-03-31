"""
run_tier_c_batch.py
-------------------
Batch processor for the 40 remaining Tier C galaxies.
Skips galaxies that already have a Phase4 output CSV (idempotent).
Runs Phase3 + Phase4 + Test3 for each galaxy.
Logs progress every 10 galaxies.

Usage:
    python run_tier_c_batch.py
"""

import subprocess
import sys
import time
import csv
import traceback
from pathlib import Path
from datetime import datetime

# ── All 40 Tier C galaxies assigned to Phase 3 (4 done in Phase 2 included) ─

ALL_44 = [
    # BARYON-DOMINATED (17)
    "UGC07261", "UGC11557", "UGC09037", "NGC4085", "UGC06628", "NGC2976",
    "UGC07577", "UGC06818", "NGC0100", "DDO064",
    "UGC07603", "UGC07151", "NGC4088", "UGC12732",
    "UGC07089", "UGC12632", "NGC4214",
    # BALANCED (19)
    "UGC00731", "NGC4559", "D564-8", "UGC07323", "UGC00191",
    "UGC04278", "UGC05918", "UGC11820", "UGC02916",
    "UGC11914", "F583-1", "UGC06917", "UGC08550",
    "UGC06930", "DDO170", "F568-3", "NGC4183",
    "NGC1090", "UGCA444",
    # DM-DOMINATED (8)
    "NGC0247", "UGC06399", "F563-1", "F565-V2",
    "DDO154", "F568-V1", "UGC08286", "UGC05005",
]

REPO       = Path(__file__).parent
PHASE4_DIR = REPO / "data" / "derived" / "phase4" / "h2_outputs"
PHASE3_DIR = REPO / "data" / "derived" / "phase3"
SPARC_DIR  = REPO / "data" / "sparc"
H1_DIR     = REPO / "data" / "h1_frozen" / "per_galaxy"

LOG_FILE   = REPO / "tier_c_batch_log.txt"
RESULT_CSV = REPO / "tier_c_batch_results.csv"

PYTHON = sys.executable


def already_done(galaxy):
    """True if Phase4 output CSV already exists."""
    return (PHASE4_DIR / f"rc_decomp_{galaxy}_H2_adaptive.csv").exists()


def run_cmd(cmd, timeout=600):
    t0 = time.time()
    proc = subprocess.run(
        cmd, capture_output=True, text=True, timeout=timeout, cwd=str(REPO)
    )
    elapsed = time.time() - t0
    output = proc.stdout + proc.stderr
    return proc.returncode, output, elapsed


def compute_delta_sigma(galaxy, inner_frac=0.5, eps=1e-12):
    import numpy as np
    import pandas as pd

    rotmod = SPARC_DIR / f"{galaxy}_rotmod.dat"
    h1_csv = H1_DIR / f"rc_decomp_{galaxy}_best.csv"
    h2_csv = PHASE4_DIR / f"rc_decomp_{galaxy}_H2_adaptive.csv"

    if not h2_csv.exists():
        return None, None, None, None
    try:
        Robs, Vobs = np.loadtxt(rotmod, comments="#", usecols=(0, 1), unpack=True)
        h1 = pd.read_csv(h1_csv)
        h2 = pd.read_csv(h2_csv)

        def get_vtot(df):
            for c in ["V_total", "V_total_H1", "Vtot"]:
                if c in df.columns:
                    return df[c].to_numpy(float)
            raise RuntimeError("No total velocity column")

        Vh1 = np.interp(Robs, h1["R_kpc"].to_numpy(), get_vtot(h1))
        Vh2 = np.interp(Robs, h2["R_kpc"].to_numpy(), get_vtot(h2))

        Rmax = float(np.max(Robs))
        mask = Robs < inner_frac * Rmax
        N = int(np.sum(mask))
        if N < 3:
            return None, None, None, N

        Vm   = np.clip(Vobs[mask], eps, None)
        Vh1m = np.clip(Vh1[mask], eps, None)
        Vh2m = np.clip(Vh2[mask], eps, None)

        sig_h1 = float(np.sqrt(np.mean((np.log10(Vh1m) - np.log10(Vm)) ** 2)))
        sig_h2 = float(np.sqrt(np.mean((np.log10(Vh2m) - np.log10(Vm)) ** 2)))
        return sig_h1, sig_h2, abs(sig_h2 - sig_h1), N
    except Exception:
        return None, None, None, None


def check_test1_pass(galaxy):
    import pandas as pd
    h2_csv = PHASE4_DIR / f"rc_decomp_{galaxy}_H2_adaptive.csv"
    if not h2_csv.exists():
        return None
    try:
        df = pd.read_csv(h2_csv)
        return bool(int(df["outer_pass"].iloc[0]))
    except Exception:
        return None


def log(fh, msg):
    ts = datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def main():
    results = []
    failed  = []
    skipped = []

    with open(LOG_FILE, "a") as logf:          # append so prior log is preserved
        log(logf, "")
        log(logf, "=" * 70)
        log(logf, "PHASE 3 BATCH RUN (RESUMED) — ALL 44 TIER C GALAXIES")
        log(logf, f"Started: {datetime.now().isoformat()}")
        log(logf, "=" * 70)

        total = len(ALL_44)

        for idx, galaxy in enumerate(ALL_44, start=1):

            # ── SKIP if already completed ─────────────────────────────
            if already_done(galaxy):
                s1, s2, delta, N = compute_delta_sigma(galaxy)
                t1 = check_test1_pass(galaxy)
                row = {
                    "galaxy": galaxy, "idx": idx,
                    "phase3_ok": True, "phase4_ok": True,
                    "test1_pass": t1,
                    "sigma_h1": round(s1, 6) if s1 is not None else None,
                    "sigma_h2": round(s2, 6) if s2 is not None else None,
                    "abs_delta_sigma": round(delta, 6) if delta is not None else None,
                    "n_inner": N,
                    "phase3_s": 0, "phase4_s": 0,
                    "error": "SKIPPED_ALREADY_DONE",
                }
                results.append(row)
                skipped.append(galaxy)
                log(logf, f"[{idx:02d}/{total}] {galaxy:<14} — SKIPPED (Phase4 CSV exists)"
                    + (f" | |Δσ|={delta:.6f} N={N}" if delta is not None else ""))
                continue

            log(logf, "")
            log(logf, f"[{idx:02d}/{total}] ── {galaxy} ─────────────────────────")

            row = {
                "galaxy": galaxy, "idx": idx,
                "phase3_ok": False, "phase4_ok": False,
                "test1_pass": None,
                "sigma_h1": None, "sigma_h2": None,
                "abs_delta_sigma": None, "n_inner": None,
                "phase3_s": None, "phase4_s": None,
                "error": "",
            }

            # ── PHASE 3 ───────────────────────────────────────────────
            try:
                rc, out, t = run_cmd([
                    PYTHON, "-m", "diagnostics.phase3_leff_ngc3198",
                    "--galaxy", galaxy,
                    "--alpha", "2.0",
                    "--sigma_idx", "1.0",
                    "--taper",
                ])
                row["phase3_s"] = round(t, 1)
                if rc != 0:
                    raise RuntimeError(f"exit {rc}: {out[-300:]}")
                if not (PHASE3_DIR / f"leff_{galaxy}.csv").exists():
                    raise RuntimeError("leff CSV missing after run")
                row["phase3_ok"] = True
                log(logf, f"  Phase3: ✅ {t:.1f}s")
            except Exception as e:
                row["error"] = f"PHASE3: {e}"
                log(logf, f"  Phase3: ❌ {e}")
                results.append(row)
                failed.append(galaxy)
                continue

            # ── PHASE 4 ───────────────────────────────────────────────
            try:
                rc, out, t = run_cmd([
                    PYTHON, "-m", "diagnostics.phase4_adaptive_convolution",
                    "--galaxy", galaxy,
                ])
                row["phase4_s"] = round(t, 1)
                if rc != 0:
                    raise RuntimeError(f"exit {rc}: {out[-400:]}")
                if not (PHASE4_DIR / f"rc_decomp_{galaxy}_H2_adaptive.csv").exists():
                    raise RuntimeError("H2 CSV missing after run")
                row["phase4_ok"] = True
                row["test1_pass"] = check_test1_pass(galaxy)
                log(logf, f"  Phase4: ✅ {t:.1f}s | Test-1: {'PASS' if row['test1_pass'] else 'FAIL'}")
            except Exception as e:
                row["error"] = f"PHASE4: {e}"
                log(logf, f"  Phase4: ❌ {e}")
                results.append(row)
                failed.append(galaxy)
                continue

            # ── TEST-3 METRIC ─────────────────────────────────────────
            s1, s2, delta, N = compute_delta_sigma(galaxy)
            row["sigma_h1"]        = round(s1, 6) if s1 is not None else None
            row["sigma_h2"]        = round(s2, 6) if s2 is not None else None
            row["abs_delta_sigma"] = round(delta, 6) if delta is not None else None
            row["n_inner"]         = N
            log(logf, f"  Test-3: σ_H1={s1:.6f}  σ_H2={s2:.6f}  |Δσ|={delta:.6f}  N={N}")

            results.append(row)

            # ── PROGRESS CHECKPOINT every 10 ─────────────────────────
            newly_done = [r for r in results if r["phase4_ok"] and r["error"] != "SKIPPED_ALREADY_DONE"]
            if (idx % 10 == 0):
                log(logf, "")
                log(logf, f"  ──── CHECKPOINT {idx}/{total} ──────────────────────────────")
                log(logf, f"  Newly processed: {len(newly_done)}")
                log(logf, f"  Skipped (cached): {len(skipped)}")
                log(logf, f"  Failed: {len(failed)}")
                if failed:
                    log(logf, f"  Failures: {', '.join(failed)}")
                log(logf, "")

        # ── WRITE RESULT CSV ──────────────────────────────────────────
        fields = ["idx","galaxy","phase3_ok","phase4_ok","test1_pass",
                  "sigma_h1","sigma_h2","abs_delta_sigma","n_inner",
                  "phase3_s","phase4_s","error"]
        with open(RESULT_CSV, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(results)

        # ── FINAL SUMMARY ─────────────────────────────────────────────
        all_done = [r for r in results if r["phase4_ok"]]
        log(logf, "")
        log(logf, "=" * 70)
        log(logf, "PHASE 3 BATCH COMPLETE")
        log(logf, f"Finished: {datetime.now().isoformat()}")
        log(logf, f"Total:    {total} galaxies")
        log(logf, f"Success:  {len(all_done)} ({len(skipped)} cached + {len(all_done)-len(skipped)} newly run)")
        log(logf, f"Failed:   {len(failed)}")
        if failed:
            log(logf, f"Failed galaxies: {', '.join(failed)}")
        log(logf, f"Results:  {RESULT_CSV}")
        log(logf, "=" * 70)

    print(f"\nBatch complete. Results: {RESULT_CSV}")


if __name__ == "__main__":
    main()
