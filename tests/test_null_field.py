import numpy as np

from core.gradients import grad_scalar, grad_log_scalar
from core.smoothing import gaussian_smooth_periodic


def test_A_constant_field_grad_zero():
    n = 128
    dx = 1.0
    f = np.full((n, n), 3.14159, dtype=np.float64)
    gx, gy = grad_scalar(f, dx, bc="periodic")
    assert np.max(np.abs(gx)) < 1e-12
    assert np.max(np.abs(gy)) < 1e-12


def test_B_constant_field_smooth_then_grad_zero():
    n = 128
    dx = 1.0
    f = np.full((n, n), 2.0, dtype=np.float64)
    fs = gaussian_smooth_periodic(f, dx, sigma_cells=1.0)
    # smoothing must preserve constants
    assert np.max(np.abs(fs - f)) < 1e-12
    gx, gy = grad_scalar(fs, dx, bc="periodic")
    assert np.max(np.abs(gx)) < 1e-12
    assert np.max(np.abs(gy)) < 1e-12


def test_C_constant_field_smooth_then_gradlog_zero():
    n = 128
    dx = 1.0
    f = np.full((n, n), 2.0, dtype=np.float64)
    fs = gaussian_smooth_periodic(f, dx, sigma_cells=1.0)
    # use a constant eps (global) so we don't manufacture structure
    gx, gy = grad_log_scalar(fs, dx, eps=1e-30, bc="periodic")
    assert np.max(np.abs(gx)) < 1e-12
    assert np.max(np.abs(gy)) < 1e-12
