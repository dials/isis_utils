from __future__ import annotations

import json
from typing import Tuple

from dials_refl_loader import load as load_refl
from experiment_reader import ExperimentReader, Panel, PeakTable
from panel import DIALSPanel

vec3float = Tuple[float, float, float]
vec2float = Tuple[float, float]


class DIALSReader(ExperimentReader):

    """
    Class to access and hold data from
    DIALS .expt (experiment) and .refl (reflection table) files

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

    def get_panels(self, expt_idx: int = 0) -> Tuple[DIALSPanel, ...]:

        self._open(mode="r", expt_idx=expt_idx)
        panels = []
        for panel in self._expt_file["detector"][expt_idx]["panels"]:
            panels.append(
                DIALSPanel(
                    name=panel["name"],
                    origin=panel["origin"],
                    fast_axis=panel["fast_axis"],
                    slow_axis=panel["slow_axis"],
                    panel_size_in_mm=(
                        panel["image_size"][0] * panel["pixel_size"][0],
                        panel["image_size"][1] * panel["pixel_size"][1],
                    ),
                )
            )
        self._close(close_refl=False)
        return tuple(panels)

    def replace_panels(
        self,
        new_panels: Tuple[DIALSPanel, ...],
        expt_idx: int = 0,
    ) -> None:

        self._open(expt_idx=expt_idx)
        new_panel_dict = {i.name: i for i in new_panels}
        for panel in self._expt_file["detector"][expt_idx]["panels"]:
            if panel["name"] in new_panel_dict:
                new_panel = new_panel_dict[panel["name"]]
                panel["origin"] = new_panel.origin
                panel["fast_axis"] = new_panel.fast_axis
                panel["slow_axis"] = new_panel.slow_axis

        with open(self.file_path, "w") as g:
            json.dump(self._expt_file, g)
        self._close(close_refl=False)

    def convert_panels(self, panels: Tuple[Panel, ...]) -> Tuple[Panel, ...]:
        return [i.to_dials() for i in panels]

    def get_peak_table(self, expt_idx: int = 0) -> PeakTable:

        """
        Reads the .refl file at self.refl_file_path and returns a PeakTable
        """

        if self.refl_file_path is None:
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
                raise ValueError(f"Cannot get peak table: {i} not found in .refl file")
        for i in required_obs_fields:
            if i not in self._refl_file.keys():
                using_calc_fields = True
                print(f"Cannot find {i} in table, trying calculated values instead..")
                break
        if using_calc_fields:
            for i in required_calc_fields:
                if i not in self._refl_file.keys():
                    raise ValueError(
                        f"Cannot get peak table: {i} not found in .refl file"
                    )

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

        return peak_table

    def replace_peak_table(self, new_peak_table: PeakTable, expt_idx: int = 0) -> None:

        pass
