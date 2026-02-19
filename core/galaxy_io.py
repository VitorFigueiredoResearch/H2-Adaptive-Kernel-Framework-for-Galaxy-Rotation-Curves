from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


# -----------------------------
# Paths (relative to repo root)
# -----------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]

# H1 frozen per-galaxy RC decompositions
H1_PER_GALAXY_DIR = REPO_ROOT / "data" / "h1_frozen" / "per_galaxy"

# Frozen parameter table (your stated location)
H1_PARAMS_JSON = H1_PER_GALAXY_DIR / "all_galaxy_params.json"

# Galaxy structural table with Rd_star (kpc), etc.
# (You already have this from H1/H1.5)
GALAXIES_CSV = REPO_ROOT / "data" / "h1p5_diagnostics" / "galaxies.csv"


@dataclass(frozen=True)
class GalaxyRC:
    name: str
    R_kpc: np.ndarray
    V_baryon_kms: np.ndarray
    V_total_kms: np.ndarray


def _parse_name_from_filename(p: Path) -> str:
    """
    Expected per-galaxy filenames:
        rc_decomp_<NAME>_best.csv
    """
    m = re.match(r"rc_decomp_(.+?)_best\.csv$", p.name)
    if not m:
        raise ValueError(f"Unrecognized per-galaxy filename: {p.name}")
    return m.group(1)


def list_galaxies() -> List[str]:
    files = sorted(H1_PER_GALAXY_DIR.glob("rc_decomp_*_best.csv"))
    names = [_parse_name_from_filename(p) for p in files]
    return names


def load_galaxy_rc(name: str) -> GalaxyRC:
    """
    Load H1 frozen per-galaxy decomposition for a given galaxy.
    Required columns:
      R_kpc, V_baryon, V_total
    """
    p = H1_PER_GALAXY_DIR / f"rc_decomp_{name}_best.csv"
    if not p.exists():
        # help: try case-insensitive search
        matches = list(H1_PER_GALAXY_DIR.glob("rc_decomp_*_best.csv"))
        cand = [q for q in matches if _parse_name_from_filename(q).lower() == name.lower()]
        if cand:
            p = cand[0]
            name = _parse_name_from_filename(p)
        else:
            raise FileNotFoundError(f"Per-galaxy file not found for '{name}': {p}")

    df = pd.read_csv(p)
    needed = {"R_kpc", "V_baryon", "V_total"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"{p.name}: missing columns {missing}; found {list(df.columns)}")

    R = df["R_kpc"].astype(np.float64).to_numpy()
    Vb = df["V_baryon"].astype(np.float64).to_numpy()
    Vt = df["V_total"].astype(np.float64).to_numpy()

    if len(R) == 0:
        raise ValueError(f"{p.name}: file has zero rows")

    # enforce sorted radius (stable gradients)
    order = np.argsort(R)
    R = R[order]
    Vb = Vb[order]
    Vt = Vt[order]

    return GalaxyRC(name=name, R_kpc=R, V_baryon_kms=Vb, V_total_kms=Vt)


def _load_galaxies_table() -> pd.DataFrame:
    if not GALAXIES_CSV.exists():
        raise FileNotFoundError(
            f"Missing galaxies.csv at {GALAXIES_CSV}. "
            f"Place your H1/H1.5 galaxies.csv there."
        )
    df = pd.read_csv(GALAXIES_CSV)

    if "name" not in df.columns or "Rd_star" not in df.columns:
        raise ValueError(
            f"{GALAXIES_CSV.name} must contain columns ['name','Rd_star',...]. "
            f"Found: {list(df.columns)}"
        )

    # normalize name column to str
    df["name"] = df["name"].astype(str)
    df["Rd_star"] = pd.to_numeric(df["Rd_star"], errors="coerce")

    return df


def get_Rd_star_kpc(name: str) -> float:
    """
    Return stellar scale length Rd_star (kpc) from galaxies.csv.
    """
    df = _load_galaxies_table()

    # exact match first
    row = df.loc[df["name"] == name]
    if row.empty:
        # try case-insensitive
        row = df.loc[df["name"].str.lower() == name.lower()]
        if row.empty:
            raise KeyError(f"Galaxy '{name}' not found in {GALAXIES_CSV}")

    val = float(row.iloc[0]["Rd_star"])
    if not np.isfinite(val) or val <= 0:
        raise ValueError(f"Invalid Rd_star for '{name}': {val}")
    return val


def get_h1_params_for_galaxy(name: str) -> Dict[str, float]:
    """
    Load frozen H1 params from all_galaxy_params.json.
    Expected structure:
      { "CamB": {"L": 50.0, "mu": 10.0, "mafe": ..., "kernel": ...}, ... }
    """
    if not H1_PARAMS_JSON.exists():
        raise FileNotFoundError(f"Missing all_galaxy_params.json at {H1_PARAMS_JSON}")

    with open(H1_PARAMS_JSON, "r", encoding="utf-8") as f:
        data = json.load(f)

    if name in data:
        d = data[name]
    else:
        # try case-insensitive
        key = None
        for k in data.keys():
            if str(k).lower() == name.lower():
                key = k
                break
        if key is None:
            raise KeyError(f"Galaxy '{name}' not found in {H1_PARAMS_JSON}")
        d = data[key]
        name = key  # normalized

    # Must contain at least L and mu; keep others if present
    if "L" not in d or "mu" not in d:
        raise ValueError(f"{name}: missing 'L' or 'mu' in all_galaxy_params.json entry: {d}")

    out = dict(d)
    out["L"] = float(out["L"])
    out["mu"] = float(out["mu"])
    return out
