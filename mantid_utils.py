from typing import Tuple

import numpy as np

from typing_utils import Array

vec3float = Array["3", float]
vec2float = Array["2", float]
vec3int = Array["3", int]
vec2int = Array["2", int]

"""
Methods for manipulating data required for Mantid.
"""


def rotations_to_spherical_coordinates(
    zeroth_pixel_origin: vec2float, rotations: Tuple[Tuple[float, vec3int], ...]
) -> vec2float:

    """
    Corrects rotations for zeroth_pixel_origin and returns
    the gam and nu angles in spherical coordinates.
    """

    gam = rotations[0][0]

    try:
        nu = rotations[1][0]
    except IndexError:
        nu = 0

    xstart, ystart = zeroth_pixel_origin
    if np.sign(xstart) == np.sign(ystart):
        if abs(nu) > 0:
            nu *= -1
        else:
            gam *= -1

    elif abs(nu) > 0:
        gam -= 180

    return gam, nu


def spherical_coordinates_to_rotations(
    zeroth_pixel_origin: vec2float, gam: float, nu: float
) -> Tuple[Tuple[float, vec3int], ...]:

    """
    Corrects gam and nu for zeroth_pixel_origin and returns
    Euler angle rotations
    """
    xstart, ystart = zeroth_pixel_origin
    if np.sign(xstart) == np.sign(ystart):
        if abs(nu) > 0:
            nu *= -1
        else:
            gam *= -1

    elif abs(gam - 180) < 1e-1:
        gam -= 180

    elif abs(nu) > 0:
        gam += 180

    return ((gam, (0, 1, 0)), (nu, (1, 0, 0)))


def panel_idx_to_name(idx: int) -> str:
    return f"bank{idx}"


def panel_name_to_idx(name: str) -> int:
    return int("".join(filter(str.isdigit, name)))
