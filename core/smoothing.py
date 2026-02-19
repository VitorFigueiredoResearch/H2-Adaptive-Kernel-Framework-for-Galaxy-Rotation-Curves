"""
H2/core/smoothing.py

Periodic Gaussian smoothing implemented in Fourier space.

Key properties:
- DC mode preserved exactly (H(k=0)=1) -> constant fields remain constant.
- Periodic boundaries (FFT-consistent with H1 convolution).
- No edge artifacts.
- Uses real FFT (rfft2/irfft2) for efficiency and strictly-real output.

Diary discipline:
- sigma should be fixed globally (e.g., sigma_cells = 1.0) and never tuned per galaxy.
"""

from __future__ import annotations
import numpy as np

Array = np.ndarray


def gaussian_smooth_periodic(
    f: Array,
    dx: float,
    *,
    sigma_phys: float | None = None,
    sigma_cells: float | None = 1.0,
    check_finite: bool = True,
) -> Array:
    """
    Apply periodic Gaussian smoothing to a 2D real field f.

    Parameters
    ----------
    f : 2D array
        Input real field.
    dx : float
        Grid spacing (physical units).
    sigma_phys : float, optional
        Gaussian sigma in physical units. If provided, used directly.
    sigma_cells : float, optional
        Gaussian sigma in grid-cell units, i.e., sigma_phys = sigma_cells * dx.
        Default is 1.0 cell (Diary-safe). Ignored if sigma_phys is provided.
    check_finite : bool
        If True, raise on NaNs/Infs.

    Returns
    -------
    f_smooth : 2D float64 array
        Smoothed field.
    """
    f = np.asarray(f)
    if f.ndim != 2:
        raise ValueError("gaussian_smooth_periodic expects a 2D array")

    dx = float(dx)
    if dx <= 0:
        raise ValueError("dx must be > 0")

    if check_finite and (not np.isfinite(f).all()):
        raise ValueError("Input field contains NaN/Inf")

    # Work in float64 for numerical stability and reproducible null tests.
    f64 = f.astype(np.float64, copy=False)

    if sigma_phys is None:
        if sigma_cells is None:
            raise ValueError("Provide either sigma_phys or sigma_cells")
        sigma_phys = float(sigma_cells) * dx
    else:
        sigma_phys = float(sigma_phys)

    if sigma_phys < 0:
        raise ValueError("sigma must be >= 0")

    # sigma = 0 -> identity (exact)
    if sigma_phys == 0.0:
        return f64.copy()

    ny, nx = f64.shape

    # Fourier frequencies (cycles per unit length) -> convert to angular wavenumber.
    kx = 2.0 * np.pi * np.fft.rfftfreq(nx, d=dx)  # shape: (nx//2 + 1,)
    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=dx)   # shape: (ny,)

    # Build k^2 grid with broadcasting: (ny, 1) + (1, nkx)
    k2 = (ky[:, None] ** 2) + (kx[None, :] ** 2)

    # Gaussian transfer function: H(k) = exp(-0.5 * (sigma*k)^2)
    # Guarantees H(0)=1 exactly -> DC preserved -> constant fields preserved.
    H = np.exp(-0.5 * (sigma_phys ** 2) * k2)

    # Real FFT -> multiply -> inverse
    F = np.fft.rfft2(f64)
    F *= H
    out = np.fft.irfft2(F, s=(ny, nx))

    # Output is real float64 by construction. Numerical noise may produce -0.0; harmless.
    return out.astype(np.float64, copy=False)
