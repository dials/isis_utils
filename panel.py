from __future__ import annotations

from abc import ABCMeta, abstractmethod
from typing import Tuple

from coords_converter import (
    cartesian_to_spherical,
    dials_panel_to_spherical_coords,
    spherical_coords_to_dials_panel,
)

vec3float = Tuple[float, float, float]
vec3int = Tuple[int, int, int]
vec2float = Tuple[float, float]


class Panel(metaclass=ABCMeta):

    """
    Base class to define interface for methods related
    to panel conversion.
    """

    @abstractmethod
    def to_dials(self) -> DIALSPanel:
        """
        Returns DIALSPanel using data in self
        """
        pass

    @abstractmethod
    def to_mantid(self) -> MantidPanel:
        """
        Returns MantidPanel using data in self
        """
        pass

    @abstractmethod
    def show(self) -> None:
        """
        Prints data
        """
        pass


class MantidPanel(Panel):

    """
    Class to store Mantid panel data
    """

    def __init__(
        self,
        name: str,
        origin: vec3float,
        rotations: Tuple[Tuple[float, vec3int], ...],
        panel_size_in_mm: vec2float,
    ) -> None:

        self.name = name
        self.origin = origin
        self.rotations = rotations
        self.panel_size_in_mm = panel_size_in_mm

    def show(self) -> None:
        print(f"name: {self.name}")
        print(f"origin: {self.origin}")
        print("rotations:")
        for i in self.rotations:
            print(i)
        print(f"panel_size_in_mm: {self.panel_size_in_mm}")

    def to_dials(self) -> DIALSPanel:
        # Mantid works in m, DIALS works in mm
        origin_mm = tuple([i * 1000 for i in self.origin])
        r_mag, gam, nu = cartesian_to_spherical(origin_mm)
        origin, fast_axis, slow_axis = spherical_coords_to_dials_panel(
            centre_origin_mag_in_mm=r_mag,
            gam_in_deg=gam,
            nu_in_deg=nu,
            panel_size_in_mm=self.panel_size_in_mm,
        )
        name = MantidPanel.get_dials_panel_name(self.name)
        return DIALSPanel(
            name=name,
            origin=origin,
            fast_axis=fast_axis,
            slow_axis=slow_axis,
            panel_size_in_mm=self.panel_size_in_mm,
        )

    def to_mantid(self) -> MantidPanel:
        return self

    @staticmethod
    def get_dials_panel_name(mantid_panel_name: str) -> str:
        num = mantid_panel_name.split("bank")[1]
        if len(num) == 1:
            return "0" + num
        return num


class DIALSPanel(Panel):
    """
    Class to store DIALS panel data
    """

    def __init__(
        self,
        name: str,
        origin: vec3float,
        fast_axis: vec3float,
        slow_axis: vec3float,
        panel_size_in_mm: vec2float,
    ) -> None:

        self.name = name
        self.origin = origin
        self.fast_axis = fast_axis
        self.slow_axis = slow_axis
        self.panel_size_in_mm = panel_size_in_mm

    def show(self) -> None:
        print(f"name: {self.name}")
        print(f"origin: {self.origin}")
        print(f"fast_axis: {self.fast_axis}")
        print(f"slow_axis: {self.slow_axis}")
        print(f"panel_size_in_mm: {self.panel_size_in_mm}")

    def to_dials(self) -> DIALSPanel:
        return self

    def to_mantid(self) -> MantidPanel:

        # Mantid works in m, DIALS works in mm
        origin_m = tuple([i / 1000 for i in self.origin])
        origin, nu, gam = dials_panel_to_spherical_coords(
            origin=origin_m,
            fast_axis=self.fast_axis,
            slow_axis=self.slow_axis,
            panel_size_in_mm=self.panel_size_in_mm,
        )
        name = DIALSPanel.get_mantid_panel_name(self.name)
        rotations = ((gam, (0, 1, 0)), (nu, (1, 0, 0)))
        return MantidPanel(
            name=name,
            rotations=rotations,
            origin=origin,
            panel_size_in_mm=self.panel_size_in_mm,
        )

    @staticmethod
    def get_mantid_panel_name(dials_panel_name: str) -> str:
        return "bank" + str(int(dials_panel_name))
