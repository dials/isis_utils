from math import radians
from typing import Tuple

import numpy as np

import coords_converter
from typing_utils import Array

vec3float = Array["3", float]
vec2float = Array["2", float]
vec3int = Array["3", int]
vec2int = Array["2", int]

"""
Methods for manipulating data required for DIALS.
"""


def panel_idx_to_name(idx: int) -> str:
    return str(idx).zfill(2)


def panel_name_to_idx(name: str) -> int:
    return int(name)


def panel_axes_flipped(fast_axis: vec3float, slow_axis: vec3float) -> bool:

    """
    The fast axis is assumed to be in the x direction,
    so if it is more closely aligned to the y direction
    than the slow_axis, the axes are assumed to be flipped.
    """

    y_direction = np.array((0, 1, 0))

    fa_dot = abs(np.dot(fast_axis, y_direction))
    sa_dot = abs(np.dot(slow_axis, y_direction))
    return fa_dot > sa_dot


def panel_axes_direction_flipped(fast_axis: vec3float, slow_axis: vec3float) -> bool:

    """
    The slow axis is assumed to be going in the y direction,
    and if the slow axis is in the negative y direction the panel
    axes directions are assumed to be flipped.
    """

    y_direction = np.array((0, 1, 0))
    return np.dot(slow_axis, y_direction) < 0


def get_panel_axes(
    r_mag: float, gam_in_deg: float, nu_in_deg: float, d_angle: float = 1e-5
) -> Tuple[vec3float, vec3float]:

    """
    Calculates the fast and slow axes by
    perturbing nu and gam by d_angle
    """

    r = coords_converter.spherical_to_vector(
        r_mag=r_mag, nu=radians(nu_in_deg), gam=radians(gam_in_deg)
    )

    ry = coords_converter.spherical_to_vector(
        r_mag=r_mag, nu=radians(nu_in_deg - d_angle), gam=radians(gam_in_deg)
    )
    rx = coords_converter.spherical_to_vector(
        r_mag=r_mag, nu=radians(nu_in_deg), gam=radians(gam_in_deg - d_angle)
    )

    fast_axis = r - rx
    slow_axis = r - ry
    fast_axis /= np.linalg.norm(fast_axis)
    slow_axis /= np.linalg.norm(slow_axis)
    return fast_axis, slow_axis
