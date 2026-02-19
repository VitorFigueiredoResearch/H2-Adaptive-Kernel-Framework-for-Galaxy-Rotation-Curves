"""
H2/core/gradients.py

Gradient operators for grid-based diagnostics (χ construction).

Hard requirements:
- Null-field control: if f(x)=const, then grad(f)=0 everywhere (up to eps_mach).
- No hidden "feature creation" from regularization: eps must be spatially constant.
- Provide both periodic and non-periodic boundary handling.

Recommended default for FFT-grid quantities: bc="periodic"
"""

from __future__ import annotations
import numpy as np


Array = np.ndarray


def _ensure_float64(a: Array) -> Array:
    a = np.asarray(a)
    if a.dtype != np.float64:
        return a.astype(np.float64, copy=False)
    return a


def grad_scalar(
    f: Array,
    dx: float,
    *,
    bc: str = "periodic",
) -> tuple[Array, Array]:
    """
    Compute ∇f on a 2D grid: returns (df/dx, df/dy) with the same shape as f.

    bc:
      - "periodic": central differences using wrap (np.roll), FFT-consistent.
      - "neumann":  central differences interior, one-sided at edges (zero-flux).

    Notes:
      - For constant f, returns exact zeros (float arithmetic exactness).
    """
    f = _ensure_float64(f)
    dx = float(dx)
    if dx <= 0:
        raise ValueError("dx must be > 0")

    if bc == "periodic":
        dfdx = (np.roll(f, -1, axis=1) - np.roll(f, +1, axis=1)) / (2.0 * dx)
        dfdy = (np.roll(f, -1, axis=0) - np.roll(f, +1, axis=0)) / (2.0 * dx)
        return dfdx, dfdy

    if bc == "neumann":
        # interior: central
        dfdx = np.empty_like(f)
        dfdy = np.empty_like(f)

        dfdx[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dx)
        dfdy[1:-1, :] = (f[2:, :] - f[:-2, :]) / (2.0 * dx)

        # edges: one-sided (consistent with zero-flux notion)
        dfdx[:, 0] = (f[:, 1] - f[:, 0]) / dx
        dfdx[:, -1] = (f[:, -1] - f[:, -2]) / dx

        dfdy[0, :] = (f[1, :] - f[0, :]) / dx
        dfdy[-1, :] = (f[-1, :] - f[-2, :]) / dx

        return dfdx, dfdy

    raise ValueError(f"Unknown bc='{bc}' (use 'periodic' or 'neumann')")


def grad_mag(
    f: Array,
    dx: float,
    *,
    bc: str = "periodic",
) -> Array:
    """Return |∇f| for scalar field f on a 2D grid."""
    dfdx, dfdy = grad_scalar(f, dx, bc=bc)
    return np.sqrt(dfdx * dfdx + dfdy * dfdy)


def grad_log_scalar(
    f: Array,
    dx: float,
    *,
    eps: float = 1e-30,
    bc: str = "periodic",
) -> tuple[Array, Array]:
    """
    Compute ∇ log(f + eps), with eps spatially constant.

    Null-field safety:
      If f is constant, log(f+eps) is constant -> grad = 0 exactly.

    Important:
      eps must NOT depend on x. If you make eps(x), you will create artificial structure.
    """
    f = _ensure_float64(f)
    eps = float(eps)
    if eps < 0:
        raise ValueError("eps must be >= 0")
    g = np.log(f + eps)
    return grad_scalar(g, dx, bc=bc)


def grad_log_mag(
    f: Array,
    dx: float,
    *,
    eps: float = 1e-30,
    bc: str = "periodic",
) -> Array:
    """Return |∇ log(f + eps)|."""
    gx, gy = grad_log_scalar(f, dx, eps=eps, bc=bc)
    return np.sqrt(gx * gx + gy * gy)


def divergence(
    vx: Array,
    vy: Array,
    dx: float,
    *,
    bc: str = "periodic",
) -> Array:
    """
    Compute ∇·v for a 2D vector field v=(vx, vy).
    """
    vx = _ensure_float64(vx)
    vy = _ensure_float64(vy)
    dx = float(dx)

    if bc == "periodic":
        dvxdx = (np.roll(vx, -1, axis=1) - np.roll(vx, +1, axis=1)) / (2.0 * dx)
        dvydy = (np.roll(vy, -1, axis=0) - np.roll(vy, +1, axis=0)) / (2.0 * dx)
        return dvxdx + dvydy

    if bc == "neumann":
        dvxdx = np.empty_like(vx)
        dvydy = np.empty_like(vy)

        dvxdx[:, 1:-1] = (vx[:, 2:] - vx[:, :-2]) / (2.0 * dx)
        dvxdx[:, 0] = (vx[:, 1] - vx[:, 0]) / dx
        dvxdx[:, -1] = (vx[:, -1] - vx[:, -2]) / dx

        dvydy[1:-1, :] = (vy[2:, :] - vy[:-2, :]) / (2.0 * dx)
        dvydy[0, :] = (vy[1, :] - vy[0, :]) / dx
        dvydy[-1, :] = (vy[-1, :] - vy[-2, :]) / dx

        return dvxdx + dvydy

    raise ValueError(f"Unknown bc='{bc}'")


def laplacian(
    f: Array,
    dx: float,
    *,
    bc: str = "periodic",
) -> Array:
    """
    Compute ∇²f on a 2D grid.

    For constant f, laplacian = 0 exactly.
    """
    f = _ensure_float64(f)
    dx = float(dx)
    dx2 = dx * dx

    if bc == "periodic":
        return (
            (np.roll(f, -1, axis=0) + np.roll(f, +1, axis=0)
             + np.roll(f, -1, axis=1) + np.roll(f, +1, axis=1)
             - 4.0 * f) / dx2
        )

    if bc == "neumann":
        # Use central interior, copy edges for simple zero-flux handling
        out = np.empty_like(f)
        out[1:-1, 1:-1] = (
            (f[2:, 1:-1] + f[:-2, 1:-1] + f[1:-1, 2:] + f[1:-1, :-2] - 4.0 * f[1:-1, 1:-1]) / dx2
        )
        # edges: clamp (simple, stable)
        out[0, :] = out[1, :]
        out[-1, :] = out[-2, :]
        out[:, 0] = out[:, 1]
        out[:, -1] = out[:, -2]
        return out

    raise ValueError(f"Unknown bc='{bc}'")


def dx_from_box(n: int, Lbox: float) -> float:
    """
    Convenience: for H1/H2 grids built on axis linspace(-Lbox, Lbox, n, endpoint=False),
    the spacing is dx = 2*Lbox / n.
    """
    n = int(n)
    if n <= 0:
        raise ValueError("n must be > 0")
    Lbox = float(Lbox)
    if Lbox <= 0:
        raise ValueError("Lbox must be > 0")
    return (2.0 * Lbox) / n
