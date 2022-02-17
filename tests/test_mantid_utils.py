import pytest

import mantid_utils


def test_rotations_to_spherical_coordinates():

    rotations = ((142.5, (0, 1, 0)), (0, (1, 0, 0)))

    zeroth_pixel_origin = (0.0945, -0.0945)

    gam, nu = mantid_utils.rotations_to_spherical_coordinates(
        zeroth_pixel_origin=zeroth_pixel_origin, rotations=rotations
    )

    expected_gam = 142.5
    expected_nu = 0

    assert gam == pytest.approx(expected_gam)
    assert nu == pytest.approx(expected_nu)

    rotations = ((-90, (0, 1, 0)), (0, (1, 0, 0)))

    zeroth_pixel_origin = (-0.0945, -0.0945)

    gam, nu = mantid_utils.rotations_to_spherical_coordinates(
        zeroth_pixel_origin=zeroth_pixel_origin, rotations=rotations
    )

    expected_gam = 90
    expected_nu = 0

    assert gam == pytest.approx(expected_gam)
    assert nu == pytest.approx(expected_nu)

    rotations = ((90, (0, 1, 0)), (-45, (1, 0, 0)))

    zeroth_pixel_origin = (0.0945, -0.0945)

    gam, nu = mantid_utils.rotations_to_spherical_coordinates(
        zeroth_pixel_origin=zeroth_pixel_origin, rotations=rotations
    )

    expected_gam = -90
    expected_nu = -45

    assert gam == pytest.approx(expected_gam)
    assert nu == pytest.approx(expected_nu)


def test_spherical_coordinates_to_rotations():

    gam = 142.5
    nu = 0

    zeroth_pixel_origin = (0.0945, -0.0945)

    rotations = mantid_utils.spherical_coordinates_to_rotations(
        zeroth_pixel_origin=zeroth_pixel_origin, gam=gam, nu=nu
    )

    expected_rotations = ((142.5, (0, 1, 0)), (0, (1, 0, 0)))

    assert rotations[0][0] == pytest.approx(expected_rotations[0][0])
    assert rotations[0][1] == expected_rotations[0][1]
    assert rotations[1][0] == pytest.approx(expected_rotations[1][0])
    assert rotations[1][1] == expected_rotations[1][1]

    gam = 90
    nu = 0

    zeroth_pixel_origin = (-0.0945, -0.0945)

    rotations = mantid_utils.spherical_coordinates_to_rotations(
        zeroth_pixel_origin=zeroth_pixel_origin, gam=gam, nu=nu
    )

    expected_rotations = ((-90, (0, 1, 0)), (0, (1, 0, 0)))

    assert rotations[0][0] == pytest.approx(expected_rotations[0][0])
    assert rotations[0][1] == expected_rotations[0][1]
    assert rotations[1][0] == pytest.approx(expected_rotations[1][0])
    assert rotations[1][1] == expected_rotations[1][1]

    gam = -90
    nu = -45

    zeroth_pixel_origin = (0.0945, -0.0945)

    rotations = mantid_utils.spherical_coordinates_to_rotations(
        zeroth_pixel_origin=zeroth_pixel_origin, gam=gam, nu=nu
    )

    expected_rotations = ((90, (0, 1, 0)), (-45, (1, 0, 0)))

    assert rotations[0][0] == pytest.approx(expected_rotations[0][0])
    assert rotations[0][1] == expected_rotations[0][1]
    assert rotations[1][0] == pytest.approx(expected_rotations[1][0])
    assert rotations[1][1] == expected_rotations[1][1]
