from math import radians

import numpy as np
import pytest

import coords_converter


def test_get_fast_slow_axes():

    gam_vals = (142.5, 90.0, 37.5, -38.0, -90.0, -142.5, 90.0, 0.0, -90.0, 180.0, 0.0)

    nu_vals = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -44.0, -44.0, -45.0, -44.0, -90.0)

    r_mag_vals = (
        225.0,
        225.0,
        225.0,
        225.0,
        225.0,
        225.0,
        270.0,
        270.0,
        270.0,
        270.0,
        278.0,
    )

    expected_fast_axes = (
        (-0.793, 0.0, -0.609),
        (0.0, 0.0, -1.0),
        (0.793, 0.0, -0.609),
        (0.788, 0.0, 0.616),
        (-0.0, 0.0, 1.0),
        (-0.793, 0.0, 0.609),
        (0.0, 0.0, -1.0),
        (1.0, 0.0, 0.0),
        (-0.0, 0.0, 1.0),
        (-1.0, 0.0, -0.0),
        (1.0, -0.0, -0.0),
    )

    expected_slow_axes = (
        (0.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (-0.0, 1.0, 0.0),
        (-0.0, 1.0, 0.0),
        (-0.0, 1.0, -0.0),
        (0.695, 0.719, 0.0),
        (0.0, 0.719, 0.695),
        (-0.707, 0.707, 0.0),
        (-0.0, 0.719, -0.695),
        (-0.0, 0.0, 1.0),
    )

    for i in range(len(gam_vals)):
        fast_axis, slow_axis = coords_converter.get_fast_slow_axes(
            r_mag=r_mag_vals[i], phi=radians(gam_vals[i]), theta=radians(nu_vals[i])
        )

        assert np.allclose(
            np.array(fast_axis), np.array(expected_fast_axes[i]), atol=1e-3
        )

        assert np.allclose(
            np.array(slow_axis), np.array(expected_slow_axes[i]), atol=1e-3
        )


def test_sxd_spherical_to_dials():

    panel_size_in_mm = tuple([(192, 192) for i in range(11)])

    centre_origin_mag_in_mm = (
        225.0,
        225.0,
        225.0,
        225.0,
        225.0,
        225.0,
        270.0,
        270.0,
        270.0,
        270.0,
        278.0,
    )

    gam_in_deg = (142.5, 90.0, 37.5, -38.0, -90.0, -142.5, 90.0, 0.0, -90.0, 180.0, 0.0)

    nu_in_deg = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -44.0, -44.0, -45.0, -44.0, -90.0)

    det_orient_x = (
        (0, -1),
        (1, 0),
        (1, 0),
        (1, 0),
        (1, 0),
        (1, 0),
        (1, 0),
        (1, 0),
        (1, 0),
        (1, 0),
        (-1, 0),
    )

    det_orient_z = (
        (-1, 0),
        (0, 1),
        (0, 1),
        (0, 1),
        (0, 1),
        (0, 1),
        (0, 1),
        (0, 1),
        (0, 1),
        (0, 1),
        (0, -1),
    )

    expected_top_left_origins_in_mm = (
        (60.81, 96.0, -236.946),
        (224.999, -96.0, 96.0),
        (60.809, -96.0, 236.945),
        (-214.172, -96.0, 118.198),
        (-224.999, -96.0, -96.0),
        (-60.809, -96.0, -236.945),
        (127.534, -256.614, 96.0),
        (-96.0, -256.614, 127.534),
        (-123.036, -258.801, -96.0),
        (96.0, -256.614, -127.534),
        (96.0, -278.0, 96.0),
    )

    expected_fast_axes = (
        (-0.0, -1.0, 0.0),
        (0.0, 0.0, -1.0),
        (0.793, 0.0, -0.609),
        (0.788, 0.0, 0.616),
        (-0.0, 0.0, 1.0),
        (-0.793, 0.0, 0.609),
        (0.0, 0.0, -1.0),
        (1.0, 0.0, 0.0),
        (-0.0, 0.0, 1.0),
        (-1.0, 0.0, -0.0),
        (-1.0, -0.0, -0.0),
    )

    expected_slow_axes = (
        (0.793, -0.0, 0.609),
        (0.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (-0.0, 1.0, 0.0),
        (-0.0, 1.0, 0.0),
        (-0.0, 1.0, -0.0),
        (0.695, 0.719, 0.0),
        (0.0, 0.719, 0.695),
        (-0.707, 0.707, 0.0),
        (-0.0, 0.719, -0.695),
        (-0.0, 0.0, -1.0),
    )

    (
        top_left_origins_in_mm,
        fast_axes,
        slow_axes,
    ) = coords_converter.sxd_spherical_to_dials(
        centre_origin_mag_in_mm=centre_origin_mag_in_mm,
        gam_in_deg=gam_in_deg,
        nu_in_deg=nu_in_deg,
        det_orient_x=det_orient_x,
        det_orient_z=det_orient_z,
        panel_size_in_mm=panel_size_in_mm,
    )

    for i in range(len(expected_top_left_origins_in_mm)):

        assert np.allclose(
            np.array(fast_axes[i]), np.array(expected_fast_axes[i]), atol=1e-3
        )
        assert np.allclose(
            np.array(slow_axes[i]), np.array(expected_slow_axes[i]), atol=1e-3
        )

        assert np.allclose(
            np.array(top_left_origins_in_mm[i]),
            np.array(expected_top_left_origins_in_mm[i]),
            atol=1e-3,
        )


def test_sxd_dials_to_spherical():

    panel_size_in_mm = tuple([(192, 192) for i in range(11)])

    dials_origins = (
        (60.809985, 96.0, -236.94636),
        (224.99904, -96.0, 96.0),
        (60.808816, -96.0, 236.944837),
        (-214.172273, -96.0, 118.198161),
        (-224.99904, -96.0, -96.0),
        (-60.808816, -96.0, -236.944837),
        (127.533717, -256.614047, 96.0),
        (-96.0, -256.614047, 127.533717),
        (-123.035761, -258.800743, -96.0),
        (96.0, -256.614047, -127.533717),
        (96.0, -278.00048, 96.00048),
    )

    dials_fast_axes = (
        (-3e-06, -1.0, 4e-06),
        (5e-06, 0.0, -1.0),
        (0.793356, 0.0, -0.608757),
        (0.788008, 0.0, 0.615665),
        (-5e-06, 0.0, 1.0),
        (-0.793356, 0.0, 0.608757),
        (5e-06, 0.0, -1.0),
        (1.0, 0.0, 5e-06),
        (-5e-06, 0.0, 1.0),
        (-1.0, 0.0, -5e-06),
        (-1.0, -0.0, -5e-06),
    )

    dials_slow_axes = (
        (0.79335, -0.0, 0.608765),
        (5e-06, 1.0, 0.0),
        (3e-06, 1.0, 4e-06),
        (-3e-06, 1.0, 4e-06),
        (-5e-06, 1.0, 0.0),
        (-3e-06, 1.0, -4e-06),
        (0.694662, 0.719336, 0.0),
        (0.0, 0.719336, 0.694662),
        (-0.70711, 0.707103, 0.0),
        (-0.0, 0.719336, -0.694662),
        (-0.0, 5e-06, -1.0),
    )

    expected_r_mag_vals = (
        225.0,
        225.0,
        225.0,
        225.0,
        225.0,
        225.0,
        270.0,
        270.0,
        270.0,
        270.0,
        278.0,
    )

    expected_gam_vals = (
        142.5,
        90.0,
        37.5,
        -38.0,
        -90.0,
        -142.5,
        90.0,
        0.0,
        -90.0,
        180.0,
        0.0,
    )

    expected_nu_vals = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -44.0, -44.0, -45.0, -44.0, -90.0)

    r_mag_vals, gam_vals, nu_vals = coords_converter.sxd_dials_to_spherical(
        top_left_origin_in_mm=dials_origins,
        fast_axes=dials_fast_axes,
        slow_axes=dials_slow_axes,
        panel_size_in_mm=panel_size_in_mm,
    )

    for i in range(len(r_mag_vals)):

        assert r_mag_vals[i] == pytest.approx(expected_r_mag_vals[i], abs=1e-3)
        assert gam_vals[i] == pytest.approx(expected_gam_vals[i], abs=1e-3)
        assert nu_vals[i] == pytest.approx(expected_nu_vals[i], abs=1e-3)
