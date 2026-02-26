"""
H2 extraction of the H1 frozen kernel grid builder.

SOURCE (H1):
- kernels.py: U_plummer, U_exp_core, U_ananta_hybrid
- run_sparc_lite.py: build_U_grid, get_U_grid, U_CACHE

Design goal:
- identical numerical behavior to H1 for kernel grid construction
- no edits to H1 files required
"""

from __future__ import annotations
import numpy as np

# -----------------------------
# 1) Analytic kernel shapes (from src/kernels.py)
# -----------------------------
def U_plummer(r: np.ndarray, L: float) -> np.ndarray:
    """Legacy Plummer (1/r)"""
    L = float(L)
    r_safe = np.maximum(r, 1e-6)
    return 1.0 / np.sqrt(r_safe**2 + L**2)

def U_exp_core(r: np.ndarray, L: float) -> np.ndarray:
    """Legacy Exponential"""
    L = float(L)
    return np.exp(-r / L) / (r + 1e-6)

def U_ananta_hybrid(r: np.ndarray, L: float, beta: float = 1.0) -> np.ndarray:
    """
    Shape-only Ananta hybrid kernel. No amplitude normalization here.
    NOTE: beta is accepted for compatibility but ignored for amplitude.
    """
    L = float(L)
    eps = 0.01 * L
    r_safe = np.sqrt(r**2 + eps**2)
    U = 0.5 * np.log(1.0 + (r_safe / L) ** 2)
    return U.astype(np.float64)


# -----------------------------
# 2) Grid builder + taper + renorm + DC-guard (from run_sparc_lite.py)
# -----------------------------
U_CACHE: dict[tuple, np.ndarray] = {}

def build_U_grid(
    n: int,
    Lbox: float,
    L: float,
    kernel: str,
    beta: float = 1.0,
    *,
    logger_fix=print,
    logger_debug=lambda *a, **k: None,
) -> np.ndarray:
    """
    Build discrete kernel grid U(x,y,z) with H1 frozen rules:
      - analytic kernel shape
      - spherical taper/cut
      - renormalize to enforce ∫U d^3r = 1/L
      - DC guard: subtract mean; set U[0]=0
    Returns float32 array.
    """
    axis = np.linspace(-Lbox, Lbox, n, endpoint=False, dtype=np.float32)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    r = np.sqrt(x * x + y * y + z * z)

    # --- analytic kernel (no normalization here) ---
    if kernel == "plummer":
        U = U_plummer(r, L)
    elif kernel == "exp-core":
        U = U_exp_core(r, L)
    elif kernel == "ananta-hybrid":
        U = U_ananta_hybrid(r, L, beta=beta)
    else:
        raise ValueError(f"Unknown kernel='{kernel}'")

    logger_debug("[DBG-K] U.dtype/min/max/mean:",
                 U.dtype, float(U.min()), float(U.max()), float(np.mean(U)))

    # --- spherical taper/cut (H1 rule) ---
    R = r
    R_cut = min(3.0 * float(L), 0.45 * float(Lbox))
    taper = np.ones_like(U, dtype=np.float64)

    r0 = 0.85 * R_cut
    mask = (R > r0) & (R <= R_cut)
    if np.any(mask):
        frac = (R[mask] - r0) / (R_cut - r0)
        taper[mask] = 0.5 * (1.0 + np.cos(np.pi * frac))
    taper[R > R_cut] = 0.0

    U = U * taper

    nonzero_frac = np.count_nonzero(taper) / taper.size
    logger_fix(f"[TAPER] nonzero fraction = {nonzero_frac:.6f}")
    if nonzero_frac < 0.02:
        raise RuntimeError(f"[TAPER-FAIL] taper removed too much kernel: nonzero={nonzero_frac:.6f}")

    # --- renormalize: enforce ∫U d^3r = 1/L exactly (H1 rule) ---
    dx = float(axis[1] - axis[0])
    cell_vol = dx**3
    current_integral = float(np.sum(U) * cell_vol)
    desired_integral = 1.0 / max(1e-12, float(L))

    if (not np.isfinite(current_integral)) or abs(current_integral) < 1e-30:
        raise RuntimeError(f"Bad kernel integral {current_integral:.3e} for L={L} at dx={dx}")

    scale = desired_integral / current_integral
    U *= scale

    logger_fix(
        f"[FIX] kernel renormalized: integral {current_integral:.3e} "
        f"-> {float(np.sum(U) * cell_vol):.3e} (scale={scale:.6e})"
    )

    # --- DC guard & single-point zero (H1 rule) ---
    U -= float(np.mean(U))
    U.flat[0] = 0.0

    return U.astype(np.float32)


def get_U_grid(
    n: int,
    Lbox: float,
    L: float,
    kernel: str,
    beta: float = 1.0,
    *,
    logger_fix=print,
    logger_debug=lambda *a, **k: None,
) -> np.ndarray:
    """
    Cached wrapper for build_U_grid. Cache key matches H1 logic.
    """
    key = (kernel, float(L), int(n), round(float(Lbox), 2), float(beta))
    if key not in U_CACHE:
        U_CACHE[key] = build_U_grid(
            n, Lbox, L, kernel, beta=beta,
            logger_fix=logger_fix, logger_debug=logger_debug
        )
    return U_CACHE[key]
