from math import asin, atan2, cos, degrees, radians, sin
from typing import Tuple

import numpy as np

vec3float = Tuple[float, float, float]
vec2float = Tuple[float, float]
vec2int = Tuple[int, int]


"""
Methods to convert between different coordinate systems
"""


def sxd_spherical_to_dials(
    centre_origin_mag_in_mm: Tuple[vec3float, ...],
    gam_in_deg: Tuple[float, ...],
    nu_in_deg: Tuple[float, ...],
    det_orient_x: Tuple[vec2int, ...],
    det_orient_z: Tuple[vec2int, ...],
    panel_size_in_mm: Tuple[vec2float, ...],
) -> Tuple[Tuple[vec3float, ...], Tuple[vec3float, ...], Tuple[vec3float, ...]]:

    """
    Converts spherical coordinates to DIALS cartesian coordinates, with
    the origin in the top left corner of each panel.

    det_orient params set if the panels have flipped or reversed axes.
    e.g. det_orient_x[i] = (1,0), det_orient_z[i] = (0,1) is the default,
    whereas (0,-1), (-1,0) would cause the fast and slow axes to be swapped, with
    the axes directions reversed.

    returns: top_left_origins_in_mm, fast_axes, slow_axes
    """

    assert (
        len(centre_origin_mag_in_mm)
        == len(gam_in_deg)
        == len(nu_in_deg)
        == len(det_orient_x)
        == len(det_orient_z)
        == len(panel_size_in_mm)
    ), "Input must all be the same length."

    top_left_origins_in_mm = []
    fast_axes = []
    slow_axes = []

    for i in range(len(centre_origin_mag_in_mm)):
        top_left_origin_in_mm, fast_axis, slow_axis = spherical_coords_to_dials_panel(
            centre_origin_mag_in_mm=centre_origin_mag_in_mm[i],
            gam_in_deg=gam_in_deg[i],
            nu_in_deg=nu_in_deg[i],
            det_orient_x=det_orient_x[i],
            det_orient_z=det_orient_z[i],
            panel_size_in_mm=panel_size_in_mm[i],
        )

        top_left_origins_in_mm.append(top_left_origin_in_mm)
        fast_axes.append(fast_axis)
        slow_axes.append(slow_axis)

    return top_left_origins_in_mm, fast_axes, slow_axes


def sxd_dials_to_spherical(
    top_left_origin_in_mm: Tuple[vec3float, ...],
    fast_axes: Tuple[vec3float, ...],
    slow_axes: Tuple[vec3float, ...],
    panel_size_in_mm: Tuple[vec2float, ...],
) -> Tuple[Tuple[float, ...], Tuple[float, ...], Tuple[float, ...]]:

    """
    Converts DIALS SXD cartesian coordinates to spherical coordinates.

    Returns r_vals, gam_vals(phi), nu_vals(theta)
    """

    assert len(top_left_origin_in_mm) == len(fast_axes) == len(slow_axes)

    r_vals = []
    gam_vals = []
    nu_vals = []

    for i in range(len(top_left_origin_in_mm)):
        r, nu, gam = dials_panel_to_spherical_coords(
            origin=top_left_origin_in_mm[i],
            fast_axis=fast_axes[i],
            slow_axis=slow_axes[i],
            panel_size_in_mm=panel_size_in_mm[i],
        )

        r_vals.append(r)
        gam_vals.append(gam)
        nu_vals.append(nu)

    return r_vals, gam_vals, nu_vals


def spherical_coords_to_dials_panel(
    centre_origin_mag_in_mm: float,
    gam_in_deg: float,
    nu_in_deg: float,
    panel_size_in_mm: vec2float,
    det_orient_x: vec2int = (1, 0),
    det_orient_z: vec2int = (0, 1),
) -> Tuple[vec3float, vec3float, vec3float]:

    """
    Converts spherical panel coords to Cartesian coords used
    by DIALS
    returns top_left_origin (mm), fast_axis, slow_axis
    """

    def axes_are_flipped():
        abs_x = tuple([abs(i) for i in det_orient_x])
        abs_z = tuple([abs(i) for i in det_orient_z])
        return abs_x == (0, 1) and abs_z == (1, 0)

    def axes_direction_flipped():
        return sum(det_orient_x) == sum(det_orient_z) == -1

    gam = radians(gam_in_deg)
    nu = radians(nu_in_deg)
    origin = spherical_to_cartesian(r_mag=centre_origin_mag_in_mm, theta=nu, phi=gam)
    fast_axis, slow_axis = get_fast_slow_axes(
        r_mag=centre_origin_mag_in_mm, theta=nu, phi=gam
    )

    if axes_are_flipped():
        slow_axis, fast_axis = fast_axis, slow_axis

    if axes_direction_flipped():
        slow_axis = tuple([-i for i in slow_axis])
        fast_axis = tuple([-i for i in fast_axis])

    top_left_origin = shift_origin_centre_to_top_left(
        centre_origin_in_mm=origin,
        fast_axis=fast_axis,
        slow_axis=slow_axis,
        panel_size_in_mm=panel_size_in_mm,
    )

    return top_left_origin, fast_axis, slow_axis


def dials_panel_to_spherical_coords(
    origin: vec3float,
    fast_axis: vec3float,
    slow_axis: vec3float,
    panel_size_in_mm: vec2float,
) -> Tuple[vec3float, float, float, vec2int, vec2int]:

    centre_origin = shift_origin_top_left_to_centre(
        top_left_origin_in_mm=origin,
        fast_axis=fast_axis,
        slow_axis=slow_axis,
        panel_size_in_mm=panel_size_in_mm,
    )

    r_mag, nu, gam = cartesian_to_sperhical(r=centre_origin)

    return centre_origin, nu, gam


def cartesian_to_sperhical(r: vec3float) -> vec3float:

    """
    Converts from cartesian coordinates to the vector magnitude
    and the two angles in spherical coordinates,

    Returns (r_mag, theta, phi)
    """

    x, y, z = r
    r_mag = np.linalg.norm(np.array(r))
    phi = get_gam_in_deg(x, z)
    theta = get_nu_in_deg(y, r_mag)

    return r_mag, theta, phi


def get_gam_in_deg(x: float, z: float) -> float:
    # Avoid small changes in sign changing around zero
    x = np.round(x, 6)
    z = np.round(z, 6)
    if abs(x) < 1e-6:
        x = 0.0
    if abs(z) < 1e-6:
        z = 0.0
    return degrees(atan2(x, z))


def get_nu_in_deg(y: float, centre_origin_mag_in_mm: float) -> float:
    return degrees(asin(y / centre_origin_mag_in_mm))


def spherical_to_cartesian(r_mag: float, theta: float, phi: float) -> vec3float:

    """
    Converts from spherical to cartesian coordinates, assuming
    beam direction is along the negative z axis, with phi in
    the zx plane, and theta measured from the plane.
    """

    z = r_mag * cos(phi) * cos(theta)
    x = r_mag * sin(phi) * cos(theta)
    y = r_mag * sin(theta)

    return (x, y, z)


def get_fast_slow_axes(
    r_mag: float, theta: float, phi: float
) -> Tuple[vec3float, vec3float]:

    """
    Calculates the fast and slow axes by
    perturbing theta and phi by d_angle
    """

    d_angle = 1e-5

    theta = np.array(theta)
    phi = np.array(phi)

    r = spherical_to_cartesian(r_mag=r_mag, theta=theta, phi=phi)

    ry = spherical_to_cartesian(r_mag=r_mag, theta=theta - d_angle, phi=phi)
    rx = spherical_to_cartesian(r_mag=r_mag, theta=theta, phi=phi - d_angle)

    r = np.array(r)
    rx = np.array(rx)
    ry = np.array(ry)

    fast_axis = r - rx
    slow_axis = r - ry
    fast_axis /= np.linalg.norm(fast_axis)
    slow_axis /= np.linalg.norm(slow_axis)
    return tuple(fast_axis), tuple(slow_axis)


def shift_origin_top_left_to_centre(
    top_left_origin_in_mm: vec3float,
    fast_axis: vec3float,
    slow_axis: vec3float,
    panel_size_in_mm: vec2float,
) -> vec3float:

    # Assumes fast and slow axes are unit vectors

    slow_axis = np.array(slow_axis) * panel_size_in_mm[1] * 0.5
    fast_axis = np.array(fast_axis) * panel_size_in_mm[0] * 0.5
    return tuple(np.array(top_left_origin_in_mm) + slow_axis + fast_axis)


def shift_origin_centre_to_top_left(
    centre_origin_in_mm: vec3float,
    fast_axis: vec3float,
    slow_axis: vec3float,
    panel_size_in_mm: vec2float,
) -> vec3float:

    # Assumes fast and slow axes are unit vectors

    slow_axis = np.array(slow_axis) * panel_size_in_mm[1] * 0.5
    fast_axis = np.array(fast_axis) * panel_size_in_mm[0] * 0.5
    return tuple(np.array(centre_origin_in_mm) - slow_axis - fast_axis)
