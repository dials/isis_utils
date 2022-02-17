from math import radians

import numpy as np
import pytest

import coords_converter


def test_shift_origin_top_left_to_centre():

    fast_axis = np.array((0.7931071, -0.0, -0.6090822))
    slow_axis = np.array((0.0, 1.0, 0.0))
    top_left_origin_in_m = np.array((0.060809, -0.096, 0.0236945))

    panel_size_in_m = np.array((0.003, 0.003))

    centre_origin_in_m = coords_converter.shift_origin_top_left_to_centre(
        top_left_origin_in_m=top_left_origin_in_m,
        fast_axis=fast_axis,
        slow_axis=slow_axis,
        panel_size_in_m=panel_size_in_m,
    )
    expected_origin = np.array((0.06199866, -0.0945, 0.02278088))

    assert np.allclose(centre_origin_in_m, expected_origin, atol=1e-3)


def test_shift_origin_centre_to_top_left():

    fast_axis = np.array((0.7931071, -0.0, -0.6090822))
    slow_axis = np.array((0.0, 1.0, 0.0))
    centre_origin_in_m = np.array((0.06199866, -0.0945, 0.02278088))
    panel_size_in_m = np.array((0.003, 0.003))

    top_left_origin_in_m = coords_converter.shift_origin_centre_to_top_left(
        centre_origin_in_m=centre_origin_in_m,
        fast_axis=fast_axis,
        slow_axis=slow_axis,
        panel_size_in_m=panel_size_in_m,
    )
    expected_origin = np.array((0.060809, -0.096, 0.0236945))

    assert np.allclose(top_left_origin_in_m, expected_origin, atol=1e-3)


def test_vector_to_spherical():

    r = np.array((0.1369713215, 0.000, 0.1785045016))

    r_mag, gam, nu = coords_converter.vector_to_spherical(r=r)

    expected_r_mag = 0.225
    expected_gam = 37.4999
    expected_nu = 0

    assert r_mag == pytest.approx(expected_r_mag, abs=1e-4)
    assert gam == pytest.approx(expected_gam, abs=1e-4)
    assert nu == pytest.approx(expected_nu, abs=1e-4)


def test_spherical_to_vector():

    gam = radians(37.5)
    nu = radians(0)
    r_mag = 0.225

    r = coords_converter.spherical_to_vector(r_mag=r_mag, gam=gam, nu=nu)

    expected_r = np.array((0.1369713215, 0.000, 0.1785045016))

    assert np.allclose(r, expected_r)
