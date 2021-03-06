"""
Prepare data for animation or other usages.
"""
import numpy as np


def sol_to_xyz(sol, n_steps, n_part):
    """
    6D to 3D for animation => array is flatten.
    :param sol: 6D from calc.traces()
    :return: 3D
    """
    ret = np.random.rand(n_steps*n_part, 3)
    for step in range(n_steps):
        for part in range(n_part):
            ret[step*n_part + part] = sol[part][step][:3]
    return ret


def field_to_surface(num, x, z, m):
    """
    For drawing in animate.py a 2D sketch of magnetic field without directions, only
    amplitudes.
    TypeError: Shape of x does not match that of z: found (200, 200) instead of (200, 400).
    """
    xx, zz = np.meshgrid(x, z)
    return xx, m, zz
