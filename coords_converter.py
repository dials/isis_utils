from math import asin, atan2, cos, degrees, sin

import numpy as np

from typing_utils import Array

vec3float = Array["3", float]
vec2float = Array["2", float]
vec3int = Array["3", int]
vec2int = Array["2", int]


"""
Methods to convert between different coordinate systems
"""


def shift_origin_top_left_to_centre(
    top_left_origin_in_m: vec3float,
    fast_axis: vec3float,
    slow_axis: vec3float,
    panel_size_in_m: vec2float,
) -> vec3float:

    # Assumes fast and slow axes are unit vectors

    slow_axis = slow_axis * panel_size_in_m[1] * 0.5
    fast_axis = fast_axis * panel_size_in_m[0] * 0.5
    return top_left_origin_in_m + slow_axis + fast_axis


def shift_origin_centre_to_top_left(
    centre_origin_in_m: vec3float,
    fast_axis: vec3float,
    slow_axis: vec3float,
    panel_size_in_m: vec2float,
) -> vec3float:

    # Assumes fast and slow axes are unit vectors

    slow_axis = slow_axis * panel_size_in_m[1] * 0.5
    fast_axis = fast_axis * panel_size_in_m[0] * 0.5
    return centre_origin_in_m - slow_axis - fast_axis


def vector_to_spherical(r: vec3float) -> vec3float:

    """
    Converts from cartesian coordinates to the vector magnitude
    and the two angles in spherical coordinates,

    Returns (r_mag, theta, phi)
    """

    x, y, z = r
    r_mag = np.linalg.norm(r)
    gam = get_gam_in_deg(x, z)
    nu = get_nu_in_deg(y, r_mag)

    return r_mag, gam, nu


def spherical_to_vector(r_mag: float, gam: float, nu: float) -> vec3float:

    """
    Converts from spherical to a cartesian vector, assuming
    beam direction is along the negative z axis, with gam in
    the zx plane, and nu measured from the plane.
    """

    z = r_mag * cos(gam) * cos(nu)
    x = r_mag * sin(gam) * cos(nu)
    y = r_mag * sin(nu)

    return np.array((x, y, z))


def get_gam_in_deg(x: float, z: float) -> float:
    # Avoid small changes in sign changing around zero
    x = np.round(x, 6)
    z = np.round(z, 6)
    if abs(x) < 1e-3:
        x = 0.0
    if abs(z) < 1e-3:
        z = 0.0
    return degrees(atan2(x, z))


def get_nu_in_deg(y: float, centre_origin_mag_in_mm: float) -> float:
    return degrees(asin(y / centre_origin_mag_in_mm))
