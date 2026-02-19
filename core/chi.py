from __future__ import annotations

from dataclasses import dataclass
import numpy as np

Array = np.ndarray


@dataclass(frozen=True)
class Chi1DResult:
    r_kpc: Array
    g_bar: Array
    dgdr: Array
    chi: Array


def compute_accel_natural_from_rc(V_kms: Array, R_kpc: Array) -> Array:
    """
    Natural acceleration units used throughout your H1.5 diagnostics:
        g ~ (km/s)^2 / kpc

    This matches your rar_points magnitude (e.g., hundreds to thousands),
    and is fine because χ is made dimensionless via Rd_star.
    """
    V_kms = np.asarray(V_kms, dtype=np.float64)
    R_kpc = np.asarray(R_kpc, dtype=np.float64)
    R_kpc = np.maximum(R_kpc, 1e-30)
    return (V_kms * V_kms) / R_kpc


def gaussian_smooth_1d(
    y: Array,
    *,
    sigma_idx: float = 1.0,
    radius: int = 4,
) -> Array:
    """
    Fixed 1D Gaussian smoothing in index-space (NOT tuned per galaxy).
    This is the minimal "pre-diff noise control" analogue of Diary 2.7,
    adapted to 1D radial sampling.

    - sigma_idx = 1.0 means ~one sample spacing.
    - radius sets kernel half-width in samples (default 4 -> 9-tap kernel).
    """
    y = np.asarray(y, dtype=np.float64)

    if sigma_idx <= 0:
        return y.copy()

    radius = int(radius)
    if radius < 1:
        radius = 1

    x = np.arange(-radius, radius + 1, dtype=np.float64)
    k = np.exp(-0.5 * (x / sigma_idx) ** 2)
    k /= np.sum(k)

    # edge padding to avoid boundary artifacts
    ypad = np.pad(y, (radius, radius), mode="edge")
    ys = np.convolve(ypad, k, mode="same")
    ys = ys[radius:-radius]
    return ys


def chi_from_gbar_1d(
    r_kpc: Array,
    g_bar: Array,
    *,
    Rd_star_kpc: float,
    eps: float = 1e-30,
) -> Chi1DResult:
    """
    Dimensionless stiffness:
        χ(r) = |dg_bar/dr| / (g_bar / Rd_star)
             = |dg_bar/dr| * Rd_star / g_bar

    Units:
      dg/dr has units ( (km/s)^2/kpc ) / kpc = (km/s)^2/kpc^2
      multiply by Rd_star (kpc) -> (km/s)^2/kpc
      divide by g_bar -> dimensionless

    eps is a constant floor (spatially constant) to avoid division by zero.
    """
    r_kpc = np.asarray(r_kpc, dtype=np.float64)
    g_bar = np.asarray(g_bar, dtype=np.float64)

    if not np.all(np.isfinite(r_kpc)):
        raise ValueError("r_kpc contains non-finite values.")
    if not np.all(np.isfinite(g_bar)):
        raise ValueError("g_bar contains non-finite values.")
    if Rd_star_kpc <= 0 or not np.isfinite(Rd_star_kpc):
        raise ValueError(f"Invalid Rd_star_kpc: {Rd_star_kpc}")

    # derivative on non-uniform grid
    dgdr = np.gradient(g_bar, r_kpc)

    denom = (g_bar / float(Rd_star_kpc)) + float(eps)
    chi = np.abs(dgdr) / denom

    return Chi1DResult(r_kpc=r_kpc, g_bar=g_bar, dgdr=dgdr, chi=chi)


def compute_chi_from_rc(
    r_kpc: Array,
    V_baryon_kms: Array,
    *,
    Rd_star_kpc: float,
    smooth: bool = True,
    sigma_idx: float = 1.0,
    radius: int = 4,
    eps: float = 1e-30,
) -> dict[str, Chi1DResult]:
    """
    Convenience wrapper:
    - builds g_bar from (V_baryon, r)
    - optionally smooths g_bar before differentiation (Diary 2.7 analogue)
    - returns dict with keys: 'raw', 'smooth' (smooth optional)
    """
    g_bar_raw = compute_accel_natural_from_rc(V_baryon_kms, r_kpc)
    out: dict[str, Chi1DResult] = {}

    out["raw"] = chi_from_gbar_1d(r_kpc, g_bar_raw, Rd_star_kpc=Rd_star_kpc, eps=eps)

    if smooth:
        g_bar_s = gaussian_smooth_1d(g_bar_raw, sigma_idx=sigma_idx, radius=radius)
        out["smooth"] = chi_from_gbar_1d(r_kpc, g_bar_s, Rd_star_kpc=Rd_star_kpc, eps=eps)

    return out
