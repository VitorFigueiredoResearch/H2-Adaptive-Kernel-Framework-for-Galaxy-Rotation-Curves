import numpy as np
from kernels.deformation import L_eff_linear

def effective_kernel_scale(L0, chi):
    """
    Adapter layer: H1 frozen kernel + H2 stiffness control.
    """
    return L_eff_linear(L0, chi)
