import numpy as np

def L_eff_linear(L0: float, chi: np.ndarray) -> np.ndarray:
    """
    Linear stiffness-controlled deformation.
    Diary 3.5 canonical form.

    L_eff = L0 / (1 + chi)

    Parameters
    ----------
    L0 : float
        Frozen H1 kernel scale (kpc)
    chi : ndarray
        Dimensionless stiffness field

    Returns
    -------
    ndarray
        Effective kernel scale at each radius
    """
    chi = np.asarray(chi, dtype=float)
    return L0 / (1.0 + chi)
