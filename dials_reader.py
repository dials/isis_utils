from __future__ import annotations

import json
import logging
from typing import Tuple

import numpy as np

import coords_converter
import dials_utils
from dials_refl_loader import load as load_refl
from experiment_reader import ExperimentReader, Panel, PeakTable

logger = logging.getLogger(__name__)


class DIALSReader(ExperimentReader):

    """
    Class to access and hold data from
    DIALS .expt (experiment) and .refl (reflection table) files.

    See https://dials.github.io/documentation/data_files.html
    """

    def __init__(self, expt_file_path: str, refl_file_path: str = None) -> None:
        self.file_path = expt_file_path
        self.refl_file_path = refl_file_path
        self._expt_file = None
        self._refl_file = None

    def _open(
        self, mode: str = "r", open_refl=False, refl_mode: str = "rb", expt_idx: int = 0
    ) -> None:
        with open(self.file_path, mode) as g:
            self._expt_file = json.load(g)
        if self.refl_file_path is not None and open_refl:
            self._open_refl(mode=refl_mode)

    def _open_refl(self, mode: str = "rb") -> None:
        if self.refl_file_path is not None:
            self._refl_file = load_refl(self.refl_file_path)

    def _close(self, close_refl: bool = True) -> None:
        self._expt_file = None
        if close_refl:
            self._close_refl()

    def _close_refl(self) -> None:
        self._refl_file = None

    def get_panels(self, expt_idx: int = 0) -> Tuple[Panel, ...]:

        """
        Extracts panel information from self._expt_file and returns
        it as a series of Panels.
        """

        self._open(mode="r", expt_idx=expt_idx)
        panels = []
        for panel in self._expt_file["detector"][expt_idx]["panels"]:

            # Extract required properties
            idx = dials_utils.panel_name_to_idx(panel["name"])
            top_left_origin_in_m = np.array([i / 1000 for i in panel["origin"]])
            fast_axis = np.array(panel["fast_axis"])
            slow_axis = np.array(panel["slow_axis"])
            num_pixels = np.array(panel["image_size"])
            pixel_size_in_m = np.array([i / 1000 for i in panel["pixel_size"]])
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
        self._close(close_refl=False)
        logger.debug(f"Extracted {len(panels)} panels.")
        return tuple(panels)

    def replace_panels(
        self,
        new_panels: Tuple[Panel, ...],
        expt_idx: int = 0,
    ) -> None:

        """
        Replaces any panels in self._expt that match those in new_panels with
        data in new_panels.
        """

        self._open(expt_idx=expt_idx)
        new_panel_dict = {dials_utils.panel_idx_to_name(i.idx): i for i in new_panels}
        panel_mod_count = 0
        for panel in self._expt_file["detector"][expt_idx]["panels"]:
            if panel["name"] in new_panel_dict:

                new_panel = new_panel_dict[panel["name"]]

                fast_axis, slow_axis = dials_utils.get_panel_axes(
                    r=new_panel.centre_origin_in_m,
                    gam_in_deg=new_panel.gam_in_deg,
                    nu_in_deg=new_panel.nu_in_deg,
                )

                # Check orientation
                if new_panel.orientations_flipped():
                    slow_axis, fast_axis = fast_axis, slow_axis

                if new_panel.orientation_direction_flipped():
                    slow_axis *= -1
                    fast_axis *= -1

                top_left_origin_in_m = coords_converter.shift_origin_centre_to_top_left(
                    centre_origin_in_m=new_panel.centre_origin_in_m,
                    fast_axis=fast_axis,
                    slow_axis=slow_axis,
                    panel_size_in_m=new_panel.panel_size_in_m(),
                )
                top_left_origin_in_mm = np.array(
                    [i * 1000 for i in top_left_origin_in_m]
                )

                panel["origin"] = tuple(top_left_origin_in_mm)
                panel["fast_axis"] = tuple(fast_axis)
                panel["slow_axis"] = tuple(slow_axis)
                panel["image_size"] = tuple(new_panel.num_pixels)
                panel["pixel_size"] = tuple(new_panel.pixel_size_in_m * 1000)
                panel_mod_count += 1

        logger.debug(f"Replaced {panel_mod_count} panels.")

        with open(self.file_path, "w") as g:
            json.dump(self._expt_file, g, indent=2, separators=(",", ": "))
        self._close(close_refl=False)

    def has_peak_table(self, expt_idx: int = 0) -> bool:
        return self.refl_file_path is not None

    def get_peak_table(self, expt_idx: int = 0) -> PeakTable:

        """
        Reads the .refl file at self.refl_file_path and returns a PeakTable
        """

        if not self.has_peak_table(expt_idx=expt_idx):
            raise ValueError("Tried to get PeakTable but self.refl_file_path is None")

        if self._refl_file is None:
            self._open_refl()

        # DIALS has separate calculated and observed fields
        # Use observed by preference but check for calculated values if observed not available
        # TODO Let user choose which to use
        using_calc_fields = False
        required_fields = [
            "spectra_idx_1D",
            "intensity.sum.value",
            "energy",
            "d_spacing",
        ]

        required_calc_fields = [
            "wavelength_calc",
            "tof_calc",
        ]

        required_obs_fields = ["wavelength", "tof"]

        for i in required_fields:
            if i not in self._refl_file.keys():
                logger.error(f"Cannot get peak table: {i} not found in .refl file")
                raise ValueError
        for i in required_obs_fields:
            if i not in self._refl_file.keys():
                using_calc_fields = True
                logger.warning(
                    f"Cannot find {i} in table, trying calculated values instead.."
                )
                break
        if using_calc_fields:
            for i in required_calc_fields:
                if i not in self._refl_file.keys():
                    logger.error(f"Cannot get peak table: {i} not found in .refl file")
                    raise ValueError

        # Get required values
        idxs = self._refl_file["id"] == expt_idx
        spectra_idx_1D = self._refl_file["spectra_idx_1D"][idxs]
        intensity = self._refl_file["intensity.sum.value"][idxs]
        energy = self._refl_file["energy"][idxs]
        d_spacing = self._refl_file["d_spacing"][idxs]

        # Get calculated or observed values
        if using_calc_fields:
            wavelength = self._refl_file["wavelength_calc"][idxs]
            tof = self._refl_file["tof_calc"][idxs] * 10 ** 6
        else:
            wavelength = self._refl_file["wavelength"][idxs]
            tof = self._refl_file["tof"][idxs] * 10 ** 6

        # Optionally get miller indices
        if "miller_indices" in self._refl_file.keys():
            miller_indices = self._refl_file["miller_indices"][idxs]
        else:
            miller_indices = None

        peak_table = PeakTable(
            spectra_idx_1D=spectra_idx_1D,
            intensity=intensity,
            energy=energy,
            wavelength=wavelength,
            d_spacing=d_spacing,
            tof=tof,
            miller_indices=miller_indices,
        )

        logger.debug(f"Extracted PeakTable of size {len(peak_table)}")

        return peak_table

    def replace_peak_table(self, new_peak_table: PeakTable, expt_idx: int = 0) -> None:
        raise NotImplementedError
