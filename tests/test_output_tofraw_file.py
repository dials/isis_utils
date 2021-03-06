from os.path import isfile, join

import h5py
import numpy as np
import pytest

"""

The following params are not tested
as an equivalent in .raw was not found:

raw_data_1/instrument/dae/vetos/fermi_chopper3
raw_data_1/instrument/dae/vetos/fifo
raw_data_1/instrument/dae/vetos/msmode
raw_data_1/instrument/dae/vetos/writing_table_file
raw_data_1/instrument/detector_table_file
raw_data_1/instrument/spectra_table_file
raw_data_1/instrument/dae/type
raw_data_1/instrument/dae/vetos/ISIS_50Hz
raw_data_1/instrument/dae/vetos/TS2_pulse
raw_data_1/instrument/dae/vetos/ext2
raw_data_1/instrument/dae/vetos/ext3
raw_data_1/instrument/moderator/distance
raw_data_1/instrument/source/name
raw_data_1/instrument/source/probe
raw_data_1/instrument/source/type
raw_data_1/measurement/first_run
raw_data_1/measurement/id
raw_data_1/measurement/label
raw_data_1/measurement/subid
raw_data_1/measurement/type
raw_data_1/monitor_events_not_saved
raw_data_1/notes
raw_data_1/periods/good_frames_daq
raw_data_1/periods/highest_used
raw_data_1/periods/labels
raw_data_1/periods/output
raw_data_1/periods/sequences
raw_data_1/periods/total_counts
raw_data_1/run_cycle
raw_data_1/run_log/*
raw_data_1/sample/distance
raw_data_1/sample/id
raw_data_1/sample/name
raw_data_1/sample/shape
raw_data_1/script_name
raw_data_1/seci_config
raw_data_1/selog
raw_data_1/total_counts
raw_data_1/total_uncounted_counts


The following params are not tested
as they are hard links to other params:

raw_data_1/instrument/detector_1/counts
raw_data_1/instrument/dae/time_channels_1/time_of_flight
raw_data_1/instrument/detector_1/spectrum_index
raw_data_1/instrument/detector_1/time_of_flight
raw_data_1/instrument/detector_1/time_of_flight_raw
raw_data_1/measurement_first_run
raw_data_1/measurement_id
raw_data_1/measurement_label
raw_data_1/measurement_subid
raw_data_1/measurement_type
raw_data_1/proton_charge
raw_data_1/proton_charge_raw
raw_data_1/raw_frames

The following params are not tested due to
known inconsistencies between file formats:

detector_1/spectrum_index
raw_data_1/framelog
raw_data_1/instrument/detector_1/distance
raw_data_1/instrument/detector_1/polar_angle
raw_data_1/monitor_*/*
raw_data_1/periods/proton_charge
raw_data_1/periods/proton_charge_raw
raw_data_1/periods/type
raw_data_1/run_number
raw_data_1/sample/type
raw_data_1/title

"""


@pytest.fixture(scope="session")
def nacl_gen_nxs(nacl_raw_reader, tmpdir_factory):

    nxs_filename = tmpdir_factory.mktemp("data").join("nacl_gen.nxs")
    nacl_raw_reader.output_tofraw_file(nxs_filename)

    return h5py.File(nxs_filename, "r")


@pytest.fixture(scope="session")
def nacl_nxs(dials_data):

    location = dials_data("isis_sxd_example_data")
    nacl_filename = join(location, "sxd_nacl_run.nxs")
    assert isfile(nacl_filename), f"{nacl_filename} not found"

    return h5py.File(nacl_filename, "r")


@pytest.mark.parametrize(
    "path",
    [
        "raw_data_1/beamline",
        "raw_data_1/collection_time",
        "raw_data_1/definition",
        "raw_data_1/definition_local",
        "raw_data_1/detector_1/period_index",
        "raw_data_1/duration",
        "raw_data_1/end_time",
        "raw_data_1/experiment_identifier",
        "raw_data_1/good_frames",
        "raw_data_1/instrument/dae/vetos/ext0",
        "raw_data_1/instrument/dae/vetos/ext1",
        "raw_data_1/instrument/dae/vetos/fermi_chopper0",
        "raw_data_1/instrument/dae/vetos/fermi_chopper1",
        "raw_data_1/instrument/dae/vetos/fermi_chopper2",
        "raw_data_1/instrument/dae/vetos/smp",
        "raw_data_1/instrument/detector_1/source_detector_distance",
        "raw_data_1/instrument/name",
        "raw_data_1/name",
        "raw_data_1/periods/frames_requested",
        "raw_data_1/periods/good_frames",
        "raw_data_1/periods/number",
        "raw_data_1/periods/raw_frames",
        "raw_data_1/sample/height",
        "raw_data_1/sample/thickness",
        "raw_data_1/sample/width",
        "raw_data_1/start_time",
        "raw_data_1/user_1/affiliation",
        "raw_data_1/user_1/name",
    ],
)
def test_same_value(nacl_nxs, nacl_gen_nxs, path):
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert gen_val == expected_val


@pytest.mark.parametrize(
    "path, atol",
    [
        ("raw_data_1/detector_1/counts", 1e-3),
        ("raw_data_1/detector_1/time_of_flight", 1e-2),
        ("raw_data_1/instrument/dae/time_channels_1/time_of_flight_raw", 100),
        ("raw_data_1/instrument/detector_1/delt", 1e-3),
    ],
)
def test_same_array(nacl_nxs, nacl_gen_nxs, path, atol):
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert gen_val.shape == expected_val.shape
    assert np.allclose(gen_val, expected_val, atol=atol)
