import numpy as np

def choose_box_and_grid(R_obs_max, L):
    """
    Choose FFT box size and grid resolution.
    This function is copied verbatim from H1 and must not be altered.
    """
    R_box = max(4.0 * L, 1.2 * R_obs_max)
    n = int(np.ceil(2 * R_box / (0.02 * L)))
    n = int(2 ** np.ceil(np.log2(n)))
    Lbox = R_box
    return n, Lbox
