from __future__ import annotations

import xml.etree.ElementTree as ET
from typing import Tuple

import h5py
import numpy as np

from experiment_reader import ExperimentReader, Panel, PeakTable
from panel import MantidPanel

vec3float = Tuple[float, float, float]
vec2float = Tuple[float, float]


class MantidReader(ExperimentReader):

    """
    Class to access and hold data from Mantid .nxs (NeXus) files

    See https://docs.mantidproject.org/nightly/concepts/NexusFile.html
    """

    _xml_path = "mantid_workspace_1/instrument/instrument_xml/data"
    _peaks_workspace_path = "mantid_workspace_1/peaks_workspace"
    _peak_workspace_columns = {
        "column_1": "detectorID",
        "column_2": "miller_idx_h",
        "column_3": "miller_idx_k",
        "column_4": "miller_idx_l",
        "column_5": "intensity",
        "column_8": "energy",
        "column_9": "energy",
        "column_10": "wavelength",
        "column_12": "d_spacing",
        "column_13": "tof",
    }

    def __init__(self, nxs_file_path: str):
        self.file_path = nxs_file_path
        self._nxs_file = None
        self._xml = None

    def _open(self, mode: str = "r", expt_idx: int = 0, open_xml: bool = True) -> None:
        self._nxs_file = h5py.File(self.file_path, mode)
        if open_xml:
            self._open_xml()

    def _open_xml(self, expt_idx: int = 0) -> None:
        self._xml = ET.fromstring(self._nxs_file[self._xml_path][expt_idx].decode())

    def _close(self) -> None:
        self._close_xml()
        self._nxs_file.close()
        self._nxs_file = None

    def _close_xml(self) -> None:
        self._xml = None

    def get_panels(self, expt_idx: int = 0) -> Tuple[MantidPanel, ...]:
        def get_rotation_vals(rot):
            val = float(rot.attrib["val"])
            x = int(rot.attrib["x-axis"])
            y = int(rot.attrib["y-axis"])
            z = int(rot.attrib["z-axis"])
            return (val, (x, y, z))

        def get_rotations(line, rotations):
            rotations.append(get_rotation_vals(line))
            try:
                return get_rotations(line[0], rotations=rotations)
            except IndexError:
                return rotations

        self._open(mode="r", expt_idx=expt_idx)

        panels = []
        panel_size_in_mm = None
        for child in self._xml:
            if self._is_panel(child):
                panel = child[0]
                name = panel.attrib["name"]
                x = float(panel.attrib["x"])
                y = float(panel.attrib["y"])
                z = float(panel.attrib["z"])
                rotation_start = panel[0]
                rotations = []
                rotations = get_rotations(
                    line=rotation_start,
                    rotations=rotations,
                )
                panels.append(
                    MantidPanel(
                        name=name,
                        origin=(x, y, z),
                        rotations=rotations,
                        panel_size_in_mm=None,
                    )
                )
            elif self._is_panel_settings(child):
                panel_size_in_mm = (
                    float(child.attrib["xpixels"])
                    * float(child.attrib["xstep"])
                    * 1000,
                    float(child.attrib["ypixels"])
                    * float(child.attrib["ystep"])
                    * 1000,
                )

        if panel_size_in_mm is not None:
            for i in panels:
                i.panel_size_in_mm = panel_size_in_mm

        self._close()

        return panels

    def replace_panels(
        self, new_panels: Tuple[MantidPanel, ...], expt_idx: int = 0
    ) -> None:

        self._open(mode="r+", expt_idx=expt_idx)

        def set_rotations(rotation_line, rotations, idx=0):
            if idx < len(rotations):
                try:
                    rotation_line.set("val", str(rotations[idx][0]))
                    rotation_line.set("x-axis", str(rotations[idx][1][0]))
                    rotation_line.set("y-axis", str(rotations[idx][1][1]))
                    rotation_line.set("z-axis", str(rotations[idx][1][2]))
                    idx += 1
                    if idx < len(rotations):
                        set_rotations(rotation_line[0], rotations, idx)
                except IndexError:
                    new_line = ET.SubElement(rotation_line, "rot")
                    new_line.set("val", "")
                    new_line.set("x-axis", "")
                    new_line.set("y-axis", "")
                    new_line.set("z-axis", "")
                    idx += 1
                    set_rotations(new_line, rotations, idx)

        panel_dict = {i.name: i for i in new_panels}
        for child in self._xml:
            if self._is_panel(child):
                panel = child[0]
                name = panel.attrib["name"]
                if name in panel_dict:
                    new_panel = panel_dict[name]
                    x, y, z = new_panel.origin
                    panel.attrib["x"] = str(x)
                    panel.attrib["y"] = str(y)
                    panel.attrib["z"] = str(z)
                    set_rotations(panel[0], rotations=new_panel.rotations)

        self._nxs_file[self._xml_path][...] = ET.tostring(self._xml)

        self._close()

    def _is_panel(self, tree_component):
        if "type" not in tree_component.attrib:
            return False
        return (
            "panel" in tree_component.attrib["type"]
            and "location" in tree_component[0].tag
        )

    def _is_panel_settings(self, tree_element):
        required_fields = ["xpixels", "ypixels", "xstep", "ystep"]
        for i in required_fields:
            if i not in tree_element.attrib:
                return False

        return True

    def convert_panels(self, panels: Tuple[Panel, ...]) -> Tuple[MantidPanel, ...]:
        return [i.to_mantid() for i in panels]

    def replace_peak_table(self, new_peak_table: PeakTable, expt_idx: int = 0) -> None:

        self._open(mode="r+", expt_idx=expt_idx, open_xml=False)

        # peak workspace
        pws = self._nxs_file[self._peaks_workspace_path]

        # map of workspace columns to PeakTable values
        pws_d = self._peak_workspace_columns

        # Mantid stores each miller index in a different column
        miller_idx_h = np.array([i[0] for i in new_peak_table["miller_indices"]])
        miller_idx_k = np.array([i[1] for i in new_peak_table["miller_indices"]])
        miller_idx_l = np.array([i[2] for i in new_peak_table["miller_indices"]])
        miller_idxs = {"h": miller_idx_h, "k": miller_idx_k, "l": miller_idx_l}

        for column in pws.keys():
            # Pad with zeros all unknown columns
            if column not in self._peak_workspace_columns:
                pws[column][:] = np.zeros(new_peak_table.get_size())
            elif "miller_idx" in self._peak_workspace_columns[column]:
                pws[column][:] = miller_idxs[pws_d[column][-1]]
            else:
                pws[column][:] = new_peak_table[pws_d[column]]

        self._close()
