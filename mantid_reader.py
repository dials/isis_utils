from __future__ import annotations

import logging
import xml.etree.ElementTree as ET
from typing import Dict, Tuple

import h5py
import numpy as np

import mantid_utils
from experiment_reader import ExperimentReader, Panel, PeakTable
from typing_utils import Array

vec2int = Array["2", int]

logger = logging.getLogger(__name__)


class MantidReader(ExperimentReader):

    """
    Class to access and hold data from Mantid .nxs (NeXus) files

    See https://docs.mantidproject.org/nightly/concepts/NexusFile.html
    """

    _xml_path = "mantid_workspace_1/instrument/instrument_xml/data"
    _peaks_workspace_path = "mantid_workspace_1/peaks_workspace"
    _process_path = "mantid_workspace_1/process"
    _peak_workspace_columns = {
        "column_1": "spectra_idx_1D",
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
    _peak_workspace_fixed_columns = ("column_11", "column_14", "column_15", "column_16")

    _peak_workspace_types = {
        "column_1": np.dtype("<i4"),
        "column_2": np.dtype("<f8"),
        "column_3": np.dtype("<f8"),
        "column_4": np.dtype("<f8"),
        "column_5": np.dtype("<f8"),
        "column_6": np.dtype("<f8"),
        "column_7": np.dtype("<f8"),
        "column_8": np.dtype("<f8"),
        "column_9": np.dtype("<f8"),
        "column_10": np.dtype("<f8"),
        "column_11": np.dtype("<f8"),
        "column_12": np.dtype("<f8"),
        "column_13": np.dtype("<f8"),
        "column_14": np.dtype("<i4"),
        "column_15": np.dtype("<f8"),
        "column_16": np.dtype("|S16"),
        "column_17": np.dtype("<i4"),
        "column_18": np.dtype("<f8"),
    }

    # The second dimension size for 2D columns
    _peak_workspace_shapes = {"column_15": 9, "column_16": 1}

    def __init__(self, nxs_file_path: str):
        self.file_path = nxs_file_path
        self._nxs_file = None
        self._xml = None

    def _open(self, mode: str = "r", expt_idx: int = 0, open_xml: bool = True) -> None:
        if self._nxs_file is not None:
            return
        self._nxs_file = h5py.File(self.file_path, mode)
        if open_xml:
            self._open_xml(expt_idx=expt_idx)

    def _open_xml(self, expt_idx: int = 0) -> None:
        self._xml = ET.fromstring(self._nxs_file[self._xml_path][expt_idx].decode())

    def _close(self) -> None:
        self._close_xml()
        self._nxs_file.close()
        self._nxs_file = None

    def _close_xml(self) -> None:
        self._xml = None

    def get_panels(self, expt_idx: int = 0) -> Tuple[Panel, ...]:
        def get_rotation_vals(rot):
            val = float(rot.attrib["val"])
            x = int(rot.attrib["axis-x"])
            y = int(rot.attrib["axis-y"])
            z = int(rot.attrib["axis-z"])
            return (val, (x, y, z))

        def get_rotations(line, rotations):
            rotations.append(get_rotation_vals(line))
            try:
                return get_rotations(line[0], rotations=rotations)
            except IndexError:
                return rotations

        def get_panel_orentation(name: str) -> Tuple[vec2int, vec2int]:
            if name == "bank1":
                return np.array((0, -1)), np.array((-1, 0))
            else:
                return np.array((1, 0)), np.array((0, 1))

        self._open(mode="r", expt_idx=expt_idx)

        panels = []
        panel_types = self.get_panel_types(self._xml)

        for child in self._xml:
            if self._is_panel(child):

                panel = child[0]
                name = panel.attrib["name"]
                idx = mantid_utils.panel_name_to_idx(name)
                x = float(panel.attrib["x"])
                y = float(panel.attrib["y"])
                z = float(panel.attrib["z"])
                rotation_start = panel[0]
                rotations = []
                rotations = get_rotations(
                    line=rotation_start,
                    rotations=rotations,
                )

                panel_type = child.attrib["type"]
                panel_info = panel_types[panel_type]

                zeroth_pixel_origin = (panel_info["xstart"], panel_info["ystart"])
                gam_in_deg, nu_in_deg = mantid_utils.rotations_to_spherical_coordinates(
                    zeroth_pixel_origin=zeroth_pixel_origin, rotations=rotations
                )

                num_pixels = (
                    panel_info["xpixels"],
                    panel_info["ypixels"],
                )

                pixel_size_in_m = (panel_info["xpixel_size"], panel_info["ypixel_size"])

                # TODO needs orientation information
                # Mantid does not appear to see SXD panel 1 as upsidedown
                x_or, y_or = get_panel_orentation(name=name)

                panels.append(
                    Panel(
                        idx=idx,
                        centre_origin_in_m=(x, y, z),
                        gam_in_deg=gam_in_deg,
                        nu_in_deg=nu_in_deg,
                        num_pixels=num_pixels,
                        pixel_size_in_m=pixel_size_in_m,
                        x_orientation=x_or,
                        y_orientation=y_or,
                    )
                )

        self._close()
        logger.debug(f"Extracted {len(panels)} panels.")
        return tuple(panels)

    def get_panel_types(self, xml):
        panel_types = {}
        for child in xml:
            if self._is_panel_settings(child):
                key = child.attrib["name"]
                xstart = float(child.attrib["xstart"])
                ystart = float(child.attrib["ystart"])
                xpixels = int(child.attrib["xpixels"])
                ypixels = int(child.attrib["ypixels"])
                xpixel_size = abs(float(child.attrib["xstep"]))
                ypixel_size = abs(float(child.attrib["ystep"]))
                panel_types[key] = {
                    "xstart": xstart,
                    "ystart": ystart,
                    "xpixels": xpixels,
                    "ypixels": ypixels,
                    "xpixel_size": xpixel_size,
                    "ypixel_size": ypixel_size,
                }

        logger.debug(f"Extracted {len(panel_types)} panel types.")
        return panel_types

    def replace_panels(self, new_panels: Tuple[Panel, ...], expt_idx: int = 0) -> None:
        def set_rotations(rotation_line, rotations, idx=0):

            """
            Recursively set rotations in self._xml
            """

            if idx < len(rotations):
                rotation_line.set("val", str(rotations[idx][0]))
                rotation_line.set("axis-x", str(rotations[idx][1][0]))
                rotation_line.set("axis-y", str(rotations[idx][1][1]))
                rotation_line.set("axis-z", str(rotations[idx][1][2]))
                idx += 1
                if idx < len(rotations):
                    # If there is still a rotation to set but it does not exist
                    # create it
                    try:
                        set_rotations(rotation_line[0], rotations, idx)
                    except IndexError:
                        new_line = ET.SubElement(rotation_line, "rot")
                        new_line.set("val", "")
                        new_line.set("axis-x", "")
                        new_line.set("axis-y", "")
                        new_line.set("axis-z", "")
                        set_rotations(new_line, rotations, idx)

        self._open(mode="r+", expt_idx=expt_idx)

        panel_dict = {mantid_utils.panel_idx_to_name(i.idx): i for i in new_panels}
        panel_types = self.get_panel_types(self._xml)

        panel_mod_count = 0
        for child in self._xml:
            if self._is_panel(child):
                panel_type = child.attrib["type"]
                panel = child[0]
                name = panel.attrib["name"]
                if name in panel_dict:
                    new_panel = panel_dict[name]
                    panel_info = panel_types[panel_type]

                    x, y, z = new_panel.centre_origin_in_m
                    panel.attrib["x"] = str(x)
                    panel.attrib["y"] = str(y)
                    panel.attrib["z"] = str(z)

                    zeroth_pixel_origin = (panel_info["xstart"], panel_info["ystart"])
                    rotations = mantid_utils.spherical_coordinates_to_rotations(
                        gam=new_panel.gam_in_deg,
                        nu=new_panel.nu_in_deg,
                        zeroth_pixel_origin=zeroth_pixel_origin,
                    )
                    set_rotations(panel[0], rotations=rotations)
                    panel_mod_count += 1

        ET.register_namespace("xsi", "http://www.w3.org/2001/XMLSchema-instance")
        ET.register_namespace("", "http://www.mantidproject.org/IDF/1.0")
        self._nxs_file[self._xml_path][...] = ET.tostring(
            self._xml, encoding="ASCII", method="xml"
        )
        logger.debug(f"Replaced {panel_mod_count} panels.")

        self._close()

    def _is_panel(self, tree_component):
        if "type" not in tree_component.attrib:
            return False
        return (
            "panel" in tree_component.attrib["type"]
            and "location" in tree_component[0].tag
        )

    def _is_panel_settings(self, tree_element):
        required_fields = ["xstart", "ystart", "xpixels", "ypixels", "xstep", "ystep"]
        for i in required_fields:
            if i not in tree_element.attrib:
                return False

        return True

    def get_peak_table(self, expt_idx: int = 0) -> PeakTable:
        raise NotImplementedError

    def has_peak_table(self, expt_idx: int = 0) -> bool:
        self._open(mode="r+", expt_idx=expt_idx, open_xml=False)
        try:
            _ = self._nxs_file[self._peaks_workspace_path]
            self._close()
            return True
        except KeyError:
            self._close()
            return False

    def replace_peak_table(self, new_peak_table: PeakTable, expt_idx: int = 0) -> None:

        """
        Updates peaks workspace at self._peaks_workspace_path with values in new_peak_table.
        Data that are not in new_peak_table are replaced with zero values,
        unless they are in columns self._peak_workspace_fixed_columns,
        where they are just sliced or padded with the zeroth idx value.
        """

        def get_resized_array(arr: np.array, new_size: int) -> np.array:

            """
            Returns copy of arr resized to new_size using zeroth idx element for padding,
            taking into account the shape of arr.
            """

            if len(arr) > new_size:
                return arr[:new_size]
            if arr.ndim == 1:
                pad_val = arr[0]
                return np.pad(
                    arr, ((0, new_size)), mode="constant", constant_values=pad_val
                )
            elif arr.ndim == 2:
                pad_val = arr[0, 0]
                return np.pad(
                    arr,
                    ((0, new_size), (0, 0)),
                    mode="constant",
                    constant_values=pad_val,
                )
            raise NotImplementedError("Cannot handle columns with dimensions > 2")

        def get_zero_array(column: str, arr_size: int) -> np.array:

            """
            Returns an array of arr_size with all elements as 0,
            taking into account the shape required by column.
            """

            shape_2d = self._peak_workspace_shapes.get(column)
            if shape_2d is not None:
                return np.zeros(arr_size * shape_2d).reshape(arr_size, shape_2d)
            return np.zeros(arr_size)

        def get_column_attributes(column: "h5py.dataset") -> Dict[str, np.bytes_]:

            """
            Gets the attributes for a given peak table
            column and returns them.
            """

            attrib_dict = {}
            for i in column.attrs.items():
                attrib_dict[i[0]] = i[1]
            return attrib_dict

        def set_column_attributes(
            column: "h5py.dataset", attrib_dict: Dict[str, np.bytes_]
        ) -> None:

            """
            Creates attribute columns for column based on values in attrib_dict.
            """

            tid = h5py.h5t.TypeID.copy(h5py.h5t.C_S1)
            tid.set_strpad(h5py.h5t.STR_NULLTERM)
            sid = h5py.h5s.create(0)
            for i in attrib_dict:
                h5py.h5a.create(
                    loc=column.id, name=i.encode("ASCII"), tid=tid, space=sid
                )
                column.attrs.__setitem__(i, np.array(attrib_dict[i], dtype="S"))

        if not self.has_peak_table(expt_idx=expt_idx):
            raise ValueError(
                f"Tried to get PeakTable but not found in {self.file_path}"
            )

        self._open(mode="r+", expt_idx=expt_idx, open_xml=False)

        # peak workspace
        pws = self._nxs_file[self._peaks_workspace_path]
        # map of workspace columns to PeakTable values
        pws_d = self._peak_workspace_columns
        # datatypes for each column
        pws_dt = self._peak_workspace_types

        # Mantid stores each miller index in a different column
        miller_idx_h = np.array([i[0] for i in new_peak_table["miller_indices"]])
        miller_idx_k = np.array([i[1] for i in new_peak_table["miller_indices"]])
        miller_idx_l = np.array([i[2] for i in new_peak_table["miller_indices"]])
        miller_idxs = {"h": miller_idx_h, "k": miller_idx_k, "l": miller_idx_l}
        size = new_peak_table.get_size()

        for column in list(pws.keys()):

            column_attributes = get_column_attributes(column=pws[column])

            # Coordinate system is not modified
            if column == "coordinate_system":
                continue

            # Just resize these columns (e.g. run number)
            if column in self._peak_workspace_fixed_columns:
                data = get_resized_array(arr=pws[column][:], new_size=size)
                del pws[column]
                pws.create_dataset(column, data=np.array(data, dtype=pws_dt[column]))
                continue

            # Pad with zeros all unknown columns
            if column not in self._peak_workspace_columns:
                del pws[column]
                pws.create_dataset(
                    column,
                    data=get_zero_array(column=column, arr_size=size),
                    dtype=pws_dt[column],
                )
            # Miller indices are stored in separate columns in Mantid
            elif "miller_idx" in self._peak_workspace_columns[column]:
                del pws[column]
                pws.create_dataset(
                    column,
                    data=np.array(miller_idxs[pws_d[column][-1]], dtype=pws_dt[column]),
                )
            # Upate with value from PeakTable, ensuring dtype is correct
            else:
                del pws[column]
                pws.create_dataset(
                    column,
                    data=np.array(new_peak_table[pws_d[column]], dtype=pws_dt[column]),
                )

            set_column_attributes(column=pws[column], attrib_dict=column_attributes)

        self._close()
