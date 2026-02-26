from __future__ import annotations
import numpy as np

Array = np.ndarray


def leff_linear(L0: float, chi: Array, alpha: float = 1.0, Lmin_frac: float = 0.05) -> Array:
    """
    Canonical H2 first-pass:
        L_eff = L0 / (1 + alpha * chi)

    - Lmin_frac prevents pathological collapse (frozen global clamp, not tuned).
    - Keep alpha fixed initially (e.g. 1.0) to avoid "tuning" accusations.
    """
    L0 = float(L0)
    alpha = float(alpha)
    chi = np.asarray(chi, dtype=np.float64)

    Le = L0 / (1.0 + alpha * chi)

    # global clamp: minimum fraction of L0
    Le_min = Lmin_frac * L0
    Le = np.maximum(Le, Le_min)

    return Le.astype(np.float64)
