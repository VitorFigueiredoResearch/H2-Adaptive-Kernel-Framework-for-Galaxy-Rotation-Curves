import numpy as np
import pytest
from core.chi import compute_chi_field

def test_chi_honesty_constant_field():
    """
    If the galaxy mass is perfectly flat (constant), 
    the stiffness chi must be zero everywhere.
    """
    n = 128
    dx = 1.0
    rd = 3.0  # arbitrary scale length
    
    # Create a constant baryonic acceleration field
    g_bar = np.full((n, n), 1.2e-10, dtype=np.float64)
    
    # Compute chi
    chi = compute_chi_field(g_bar, dx, rd_kpc=rd)
    
    # The max value of chi should be effectively zero (machine precision)
    max_chi = np.max(np.abs(chi))
    
    print(f"\n[Test] Constant field max chi: {max_chi:.3e}")
    assert max_chi < 1e-12

def test_chi_scaling_invariance():
    """
    Test that chi is dimensionless. If we double the field strength 
    but keep the shape the same, chi (the relative gradient) should not change.
    """
    n = 128
    dx = 1.0
    rd = 2.5
    
    # Create a Gaussian 'blob'
    y, x = np.indices((n, n))
    r2 = (x - n//2)**2 + (y - n//2)**2
    field1 = np.exp(-r2 / (2 * 10**2))
    field2 = field1 * 10.0  # Same shape, 10x stronger
    
    chi1 = compute_chi_field(field1, dx, rd_kpc=rd)
    chi2 = compute_chi_field(field2, dx, rd_kpc=rd)
    
    # The stiffness depends on the SHAPE (gradient), not the total mass.
    np.testing.assert_allclose(chi1, chi2, atol=1e-12)
    print("[Test] Chi is correctly shape-dependent, not mass-dependent.")