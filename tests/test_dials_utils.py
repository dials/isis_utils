import numpy as np

import dials_utils


def test_panel_axes_flipped():

    fast_axis = np.array((0.0, -1.0, 0.0))
    slow_axis = np.array((0.793, 0.0, 0.609))

    assert dials_utils.panel_axes_flipped(fast_axis=fast_axis, slow_axis=slow_axis)

    fast_axis = np.array((-0.0, -0.0, -1.0))
    slow_axis = np.array((0.0, 1.0, 0.0))

    assert not dials_utils.panel_axes_flipped(fast_axis=fast_axis, slow_axis=slow_axis)


def test_panel_axes_direction_flipped():

    slow_axis = np.array((0.0, -1.0, 0.0))
    fast_axis = np.array((0.793, 0.0, 0.609))

    assert dials_utils.panel_axes_direction_flipped(
        fast_axis=fast_axis, slow_axis=slow_axis
    )

    fast_axis = np.array((-0.0, -0.0, -1.0))
    slow_axis = np.array((0.0, 1.0, 0.0))

    assert not dials_utils.panel_axes_direction_flipped(
        fast_axis=fast_axis, slow_axis=slow_axis
    )


def test_get_panel_axes():

    r = np.array((0.1369713215, 0.000, 0.1785045016))

    gam_in_deg = 37.5
    nu_in_deg = 0.0

    fast_axis, slow_axis = dials_utils.get_panel_axes(
        r=r, gam_in_deg=gam_in_deg, nu_in_deg=nu_in_deg
    )

    expected_fast_axis = np.array((0.79335, 0, -0.60876))

    expected_slow_axis = np.array((0, 1, 0))

    assert np.allclose(fast_axis, expected_fast_axis, atol=1e-3)
    assert np.allclose(slow_axis, expected_slow_axis, atol=1e-3)
