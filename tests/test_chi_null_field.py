from core.chi import chi_null_field_check

def test_chi_null_field():
    passed, chi_max = chi_null_field_check(tol=1e-6)
    assert passed, f"chi_max={chi_max}"
