import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _repo_root() -> Path:
    # diagnostics/.. => repo root
    return Path(__file__).resolve().parents[1]


def pearson_r(x, y) -> float:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    x0 = x - np.mean(x)
    y0 = y - np.mean(y)
    den = float(np.sqrt(np.sum(x0 * x0) * np.sum(y0 * y0)))
    return float(np.sum(x0 * y0) / den) if den > 0 else float("nan")


def read_rotmod_obs(galaxy: str, repo: Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Read SPARC rotmod observed RC:
      col0 = Rad [kpc]
      col1 = Vobs [km/s]
    """
    # Search recursively under H2/data/sparc for the rotmod file
    sparc_root = repo / "data" / "sparc"
    if not sparc_root.exists():
        raise RuntimeError(f"SPARC root not found: {sparc_root}")

    hits = list(sparc_root.rglob(f"{galaxy}_rotmod.dat"))
    if len(hits) == 0:
        # case-insensitive fallback
        hits = [p for p in sparc_root.rglob("*rotmod*.dat") if galaxy.upper() in p.name.upper()]

    if len(hits) == 0:
        raise RuntimeError(
            f"Could not find {galaxy}_rotmod.dat under {sparc_root} (recursive search)."
        )

    # If multiple, pick the shortest path (usually the canonical one)
    hits = sorted(hits, key=lambda p: len(str(p)))
    path = hits[0]


    R = []
    V = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            try:
                r = float(parts[0])
                v = float(parts[1])
            except Exception:
                continue
            if np.isfinite(r) and np.isfinite(v) and (r >= 0) and (v > 0):
                R.append(r)
                V.append(v)

    if len(R) < 5:
        raise RuntimeError(f"Too few valid obs points read from {path} (N={len(R)}).")

    R = np.asarray(R, dtype=float)
    V = np.asarray(V, dtype=float)

    # Sort just in case
    idx = np.argsort(R)
    return R[idx], V[idx]


def main():
    repo = _repo_root()
    # allow running as "python diagnostics/..." (not only -m)
    if str(repo) not in sys.path:
        sys.path.insert(0, str(repo))

    g = "NGC3198"

    p3_path = repo / "data" / "derived" / "phase3" / f"leff_{g}.csv"
    p4_path = repo / "data" / "derived" / "phase4" / "h2_outputs" / f"rc_decomp_{g}_H2_adaptive.csv"

    if not p3_path.exists():
        raise RuntimeError(f"Missing Phase-3 file: {p3_path}")
    if not p4_path.exists():
        raise RuntimeError(f"Missing Phase-4 file: {p4_path}")

    p3 = pd.read_csv(p3_path).sort_values("R_kpc").reset_index(drop=True)
    p4 = pd.read_csv(p4_path).sort_values("R_kpc").reset_index(drop=True)

    # Canonical grid from phase4
    R = p4["R_kpc"].to_numpy(float)
    rmax = float(np.max(R))

    # Choose chi column (prefer chi_used, else chi_raw)
    chi_col = None
    for c in ("chi_used", "chi_raw", "chi", "chi_smooth"):
        if c in p3.columns:
            chi_col = c
            break
    if chi_col is None:
        raise RuntimeError(
            f"No chi column found in {p3_path.name}. "
            f"Expected one of: chi_used, chi_raw, chi, chi_smooth"
        )

    R3 = p3["R_kpc"].to_numpy(float)
    chi3 = p3[chi_col].to_numpy(float)

    # Align chi to phase4 grid (should match; interpolate if not)
    if (len(R3) != len(R)) or (not np.allclose(R3, R, atol=1e-10, rtol=0.0)):
        chi = np.interp(R, R3, chi3, left=chi3[0], right=chi3[-1])
    else:
        chi = chi3

    # Observed RC from rotmod.dat
    R_obs, V_obs = read_rotmod_obs(g, repo)
    Vobs = np.interp(R, R_obs, V_obs, left=np.nan, right=np.nan)

    # H1 prediction from phase4 file
    if "V_total_H1" not in p4.columns:
        raise RuntimeError(f"Phase-4 file missing V_total_H1 column: {p4_path}")

    Vh1 = p4["V_total_H1"].to_numpy(float)

    # Inner region mask: r_frac < 0.30, and require finite Vobs
    mask = (R > 0) & ((R / rmax) < 0.30) & np.isfinite(Vobs) & (Vobs > 0)

    chi_m = chi[mask]
    dV_m = (Vh1 - Vobs)[mask]
    absdV = np.abs(dV_m)

    print("galaxy:", g)
    print("chi column:", chi_col)
    print("N points:", int(np.sum(mask)))
    print("Pearson r (chi vs signed dV=Vh1-Vobs):", pearson_r(chi_m, dV_m))
    print("Pearson r (chi vs |dV|):", pearson_r(chi_m, absdV))
    print("chi range:", float(np.min(chi_m)), float(np.max(chi_m)))
    print("|dV| range:", float(np.min(absdV)), float(np.max(absdV)))


if __name__ == "__main__":
    main()
