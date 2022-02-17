from __future__ import annotations

import numpy as np

from typing_utils import Array

vec3float = Array["3", float]
vec2float = Array["2", float]
vec3int = Array["3", int]
vec2int = Array["2", int]


class Panel:

    """
    Class to store panel properties in spherical coordinates
    that different readers can understand.
    """

    def __init__(
        self,
        idx: int,
        centre_origin_in_m: float,
        gam_in_deg: float,
        nu_in_deg: float,
        num_pixels: vec2int,
        pixel_size_in_m: vec2float,
        x_orientation: vec2int = np.array((1, 0)),
        y_orientation: vec2int = np.array((0, 1)),
    ):

        self.idx = idx
        self.centre_origin_in_m = centre_origin_in_m
        self.gam_in_deg = gam_in_deg
        self.nu_in_deg = nu_in_deg
        self.num_pixels = num_pixels
        self.pixel_size_in_m = pixel_size_in_m
        self.x_orientation = x_orientation
        self.y_orientation = y_orientation

    def panel_size_in_m(self) -> vec2float:
        panel_size = []
        for i in range(len(self.num_pixels)):
            panel_size.append(self.num_pixels[i] * self.pixel_size_in_m[i])
        return np.array(panel_size)

    def orientations_flipped(self) -> bool:
        abs_x = tuple([abs(i) for i in self.x_orientation])
        abs_y = tuple([abs(i) for i in self.y_orientation])
        return abs_x == (0, 1) and abs_y == (1, 0)

    def orientation_direction_flipped(self) -> bool:
        return sum(self.x_orientation) == sum(self.y_orientation) == -1

    def __repr__(self) -> None:
        return f"idx: {self.idx} \n \
        centre origin (m): {self.centre_origin_in_m} \n \
        gam (deg): {self.gam_in_deg} \n \
        nu (deg): {self.nu_in_deg} \n \
        num pixels: {self.num_pixels} \n \
        pixel size (m): {self.pixel_size_in_m} \n \
        x_orientation: {self.x_orientation} \n \
        y_orientation: {self.y_orientation}"
