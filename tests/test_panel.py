import numpy as np

from panel import Panel


def test_orientations_flipped():

    panel = Panel(
        idx=0,
        centre_origin_in_m=(0, 0, 0),
        gam_in_deg=0,
        nu_in_deg=0,
        num_pixels=np.array((0, 0)),
        pixel_size_in_m=np.array((0, 0)),
        x_orientation=np.array((1, 0)),
        y_orientation=np.array((0, 1)),
    )

    assert not panel.orientations_flipped()

    panel.x_orientation = np.array((-1, 0))
    panel.y_orientation = np.array((0, -1))

    assert not panel.orientations_flipped()

    panel.x_orientation = np.array((-1, 0))
    panel.y_orientation = np.array((0, 1))

    assert not panel.orientations_flipped()

    panel.x_orientation = np.array((0, 1))
    panel.y_orientation = np.array((1, 0))

    assert panel.orientations_flipped()

    panel.x_orientation = np.array((0, -1))
    panel.y_orientation = np.array((-1, 0))

    assert panel.orientations_flipped()


def test_orientation_direction_flipped():

    panel = Panel(
        idx=0,
        centre_origin_in_m=(0, 0, 0),
        gam_in_deg=0,
        nu_in_deg=0,
        num_pixels=np.array((0, 0)),
        pixel_size_in_m=np.array((0, 0)),
        x_orientation=np.array((1, 0)),
        y_orientation=np.array((0, 1)),
    )

    assert not panel.orientation_direction_flipped()

    panel.x_orientation = np.array((0, 1))
    panel.y_orientation = np.array((1, 0))

    assert not panel.orientation_direction_flipped()

    panel.x_orientation = np.array((-1, 0))
    panel.y_orientation = np.array((0, -1))

    assert panel.orientation_direction_flipped()
