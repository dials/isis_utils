import numpy as np
import pytest

import coords_converter
import dials_utils
import mantid_utils
from panel import Panel


@pytest.fixture
def dials_sxd_panels():
    return {
        "01": {
            "fast_axis": (0.0, -1.0, 0.0),
            "slow_axis": (0.793, 0.0, 0.609),
            "origin": (60.81, 96.0, -236.946),
        },
        "02": {
            "fast_axis": (-0.0, -0.0, -1.0),
            "slow_axis": (0.0, 1.0, 0.0),
            "origin": (224.999, -96.0, 96.0),
        },
        "03": {
            "fast_axis": (0.793, -0.0, -0.609),
            "slow_axis": (0.0, 1.0, 0.0),
            "origin": (60.809, -96.0, 236.945),
        },
        "04": {
            "fast_axis": (0.788, -0.0, 0.616),
            "slow_axis": (0.0, 1.0, 0.0),
            "origin": (-214.172, -96.0, 118.198),
        },
        "05": {
            "fast_axis": (-0.0, -0.0, 1.0),
            "slow_axis": (0.0, 1.0, 0.0),
            "origin": (-224.999, -96.0, -96.0),
        },
        "06": {
            "fast_axis": (-0.793, -0.0, 0.609),
            "slow_axis": (0.0, 1.0, 0.0),
            "origin": (-60.809, -96.0, -236.945),
        },
        "07": {
            "fast_axis": (0.0, -0.0, -1.0),
            "slow_axis": (0.695, 0.719, -0.0),
            "origin": (127.534, -256.614, 96.0),
        },
        "08": {
            "fast_axis": (1.0, -0.0, -0.0),
            "slow_axis": (0.0, 0.719, 0.695),
            "origin": (-96.0, -256.614, 127.534),
        },
        "09": {
            "fast_axis": (-0.0, -0.0, 1.0),
            "slow_axis": (-0.707, 0.707, -0.0),
            "origin": (-123.036, -258.801, -96.0),
        },
        "10": {
            "fast_axis": (-1.0, -0.0, -0.0),
            "slow_axis": (0.0, 0.719, -0.695),
            "origin": (96.0, -256.614, -127.534),
        },
        "11": {
            "fast_axis": (-1.0, -0.0, -0.0),
            "slow_axis": (-0.0, 0.0, -1.0),
            "origin": (96.0, -278.0, 96.0),
        },
    }


@pytest.fixture
def mantid_sxd_panels():
    return {
        "bank1": {
            "origin": (0.136938, 0.0, -0.17848199999999997),
            "rotations": ((142.5, (0, 1, 0)), (0.0, (1, 0, 0))),
        },
        "bank2": {
            "origin": (0.224999, 0.0, 0.0),
            "rotations": ((-90, (0, 1, 0)), (0.0, (1, 0, 0))),
        },
        "bank3": {
            "origin": (0.136937, 0.0, 0.178481),
            "rotations": ((37.5, (0, 1, 0)), (0.0, (1, 0, 0))),
        },
        "bank4": {
            "origin": (-0.13852399999999998, 0.0, 0.177334),
            "rotations": ((-38.0, (0, 1, 0)), (0.0, (1, 0, 0))),
        },
        "bank5": {
            "origin": (-0.224999, 0.0, 0.0),
            "rotations": ((90, (0, 1, 0)), (0.0, (1, 0, 0))),
        },
        "bank6": {
            "origin": (-0.136937, 0.0, -0.178481),
            "rotations": ((-142.5, (0, 1, 0)), (0.0, (1, 0, 0))),
        },
        "bank7": {
            "origin": (0.194254, -0.18758999999999995, 0.0),
            "rotations": ((90.0, (0, 1, 0)), (44.0, (1, 0, 0))),
        },
        "bank8": {
            "origin": (0.0, -0.18758999999999995, 0.194254),
            "rotations": ((0.0, (0, 1, 0)), (44.0, (1, 0, 0))),
        },
        "bank9": {
            "origin": (-0.19090800000000002, -0.19092900000000002, 0.0),
            "rotations": ((90.0, (0, 1, 0)), (-45.0, (1, 0, 0))),
        },
        "bank10": {
            "origin": (0.0, -0.18758999999999995, -0.194254),
            "rotations": ((0.0, (0, 1, 0)), (-44.0, (1, 0, 0))),
        },
        "bank11": {
            "origin": (0.000, -0.278, 0.000),
            "rotations": ((0.0, (0, 1, 0)), (90.0, (1, 0, 0))),
        },
    }


@pytest.fixture
def mantid_panel_info():
    return {
        "panel": {"xstart": -0.0945, "ystart": -0.0945},
        "panel1346": {"xstart": 0.0975, "ystart": -0.0945},
        "panel11": {"xstart": 0.0975, "ystart": 0.0975},
    }


@pytest.fixture
def mantid_sxd_panel_types():

    return {
        1: "panel1346",
        2: "panel",
        3: "panel1346",
        4: "panel1346",
        5: "panel",
        6: "panel1346",
        7: "panel",
        8: "panel",
        9: "panel1346",
        10: "panel1346",
        11: "panel11",
    }


def test_sxd_panels_dials_to_mantid(
    dials_sxd_panels, mantid_sxd_panels, mantid_panel_info, mantid_sxd_panel_types
):

    dials_panels = dials_sxd_panels
    dials_num_pixels = np.array((64, 64))
    dials_pixel_size_in_m = np.array((0.003, 0.003))

    panels = []
    for i in dials_panels:

        # Extract required properties
        idx = dials_utils.panel_name_to_idx(i)
        top_left_origin_in_m = np.array([i / 1000 for i in dials_panels[i]["origin"]])
        fast_axis = np.array(dials_panels[i]["fast_axis"])
        slow_axis = np.array(dials_panels[i]["slow_axis"])
        num_pixels = dials_num_pixels
        pixel_size_in_m = dials_pixel_size_in_m
        panel_size_in_m = np.multiply(num_pixels, pixel_size_in_m)

        # Get spherical coordinates
        centre_origin_in_m = coords_converter.shift_origin_top_left_to_centre(
            top_left_origin_in_m=top_left_origin_in_m,
            fast_axis=fast_axis,
            slow_axis=slow_axis,
            panel_size_in_m=panel_size_in_m,
        )

        _, gam_in_deg, nu_in_deg = coords_converter.vector_to_spherical(
            r=centre_origin_in_m
        )

        # Get orientation
        x_orientation = (1, 0)
        y_orientation = (0, 1)

        if dials_utils.panel_axes_flipped(fast_axis=fast_axis, slow_axis=slow_axis):
            y_orientation, x_orientation = x_orientation, y_orientation

        if dials_utils.panel_axes_direction_flipped(
            fast_axis=fast_axis, slow_axis=slow_axis
        ):
            y_orientation = tuple([-i for i in y_orientation])
            x_orientation = tuple([-i for i in x_orientation])

        panels.append(
            Panel(
                idx=idx,
                centre_origin_in_m=centre_origin_in_m,
                gam_in_deg=gam_in_deg,
                nu_in_deg=nu_in_deg,
                num_pixels=num_pixels,
                pixel_size_in_m=pixel_size_in_m,
                x_orientation=x_orientation,
                y_orientation=y_orientation,
            )
        )

    mantid_panels = {}
    for i in panels:
        name = mantid_utils.panel_idx_to_name(i.idx)
        panel_type = mantid_sxd_panel_types[i.idx]
        zeroth_pixel_origin = (
            mantid_panel_info[panel_type]["xstart"],
            mantid_panel_info[panel_type]["ystart"],
        )

        if i.idx == 10:
            print("y")

        rotations = mantid_utils.spherical_coordinates_to_rotations(
            gam=i.gam_in_deg, nu=i.nu_in_deg, zeroth_pixel_origin=zeroth_pixel_origin
        )

        mantid_panels[name] = {"origin": i.centre_origin_in_m, "rotations": rotations}

    expected_mantid_panels = mantid_sxd_panels

    for i in expected_mantid_panels:
        assert np.allclose(
            np.array(expected_mantid_panels[i]["origin"]),
            np.array(mantid_panels[i]["origin"]),
            atol=1e-1,
        )

        for j in range(len(expected_mantid_panels[i]["rotations"])):
            assert expected_mantid_panels[i]["rotations"][j][0] == pytest.approx(
                mantid_panels[i]["rotations"][j][0], abs=1e-1
            )
            assert (
                expected_mantid_panels[i]["rotations"][j][1]
                == mantid_panels[i]["rotations"][j][1]
            )


def test_sxd_panels_mantid_to_dials(
    dials_sxd_panels, mantid_sxd_panels, mantid_panel_info, mantid_sxd_panel_types
):

    num_pixels = np.array((64, 64))
    pixel_size_in_m = np.array((0.003, 0.003))

    panels = []
    for i in mantid_sxd_panels:

        idx = mantid_utils.panel_name_to_idx(i)
        panel_info = mantid_panel_info[mantid_sxd_panel_types[idx]]

        zeroth_pixel_origin = (panel_info["xstart"], panel_info["ystart"])
        gam_in_deg, nu_in_deg = mantid_utils.rotations_to_spherical_coordinates(
            zeroth_pixel_origin=zeroth_pixel_origin,
            rotations=mantid_sxd_panels[i]["rotations"],
        )
        if i == "bank1":
            x_or, y_or = np.array((0, -1)), np.array((-1, 0))
        elif i == "bank11":
            x_or, y_or = np.array((-1, 0)), np.array((0, -1))
        else:
            x_or, y_or = np.array((1, 0)), np.array((0, 1))

        panels.append(
            Panel(
                idx=idx,
                centre_origin_in_m=mantid_sxd_panels[i]["origin"],
                gam_in_deg=gam_in_deg,
                nu_in_deg=nu_in_deg,
                num_pixels=num_pixels,
                pixel_size_in_m=pixel_size_in_m,
                x_orientation=x_or,
                y_orientation=y_or,
            )
        )

        dials_panels = {}

    for i in panels:

        fast_axis, slow_axis = dials_utils.get_panel_axes(
            r_mag=np.linalg.norm(i.centre_origin_in_m),
            gam_in_deg=i.gam_in_deg,
            nu_in_deg=i.nu_in_deg,
        )

        # Check orientation
        if i.orientations_flipped():
            slow_axis, fast_axis = fast_axis, slow_axis

        if i.orientation_direction_flipped():
            slow_axis *= -1
            fast_axis *= -1

        top_left_origin_in_m = coords_converter.shift_origin_centre_to_top_left(
            centre_origin_in_m=i.centre_origin_in_m,
            fast_axis=fast_axis,
            slow_axis=slow_axis,
            panel_size_in_m=i.panel_size_in_m(),
        )
        top_left_origin_in_mm = np.array([i * 1000 for i in top_left_origin_in_m])

        name = dials_utils.panel_idx_to_name(i.idx)

        dials_panels[name] = {
            "fast_axis": fast_axis,
            "slow_axis": slow_axis,
            "origin": top_left_origin_in_mm,
        }

    expected_dials_panels = dials_sxd_panels
    for i in expected_dials_panels:
        assert np.allclose(
            np.array(expected_dials_panels[i]["origin"]),
            np.array(dials_panels[i]["origin"]),
            atol=1e-1,
        )
        assert np.allclose(
            np.array(expected_dials_panels[i]["fast_axis"]),
            np.array(dials_panels[i]["fast_axis"]),
            atol=1e-3,
        )
        assert np.allclose(
            np.array(expected_dials_panels[i]["slow_axis"]),
            np.array(dials_panels[i]["slow_axis"]),
            atol=1e-3,
        )
