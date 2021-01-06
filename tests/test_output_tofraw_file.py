import pytest
import numpy as np

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
"raw_data_1/instrument/source/name"
"raw_data_1/instrument/source/probe"
"raw_data_1/instrument/source/type"

The following params are not tested
as they are hard links to other params:

raw_data_1/instrument/detector_1/counts
raw_data_1/instrument/dae/time_channels_1/time_of_flight
raw_data_1/instrument/detector_1/spectrum_index
raw_data_1/instrument/detector_1/time_of_flight
raw_data_1/instrument/detector_1/time_of_flight_raw

The following params are not tested due to 
known inconsistencies between file formats:

detector_1/spectrum_index
raw_data_1/framelog
raw_data_1/instrument/detector_1/distance
raw_data_1/instrument/detector_1/polar_angle

"""

@pytest.mark.parametrize("path", 
["raw_data_1/beamline",
"raw_data_1/collection_time",
"raw_data_1/definition",
"raw_data_1/definition_local",
"raw_data_1/detector_1/period_index",
"raw_data_1/duration",
"raw_data_1/end_time",
"raw_data_1/experiment_identifier",
"raw_data_1/good_frames",
"raw_data_1/instrument/good_frames",
"raw_data_1/instrument/dae/vetos/ext0",
"raw_data_1/instrument/dae/vetos/ext1",
"raw_data_1/instrument/dae/vetos/fermi_chopper_0",
"raw_data_1/instrument/dae/vetos/fermi_chopper_1",
"raw_data_1/instrument/dae/vetos/fermi_chopper_2",
"raw_data_1/instrument/dae/vetos/smp",
"raw_data_1/instrument/detector_1/source_detector_distance",
"raw_data_1/instrument/name",

])
def test_same_value(nacl_nxs, nacl_gen_nxs, path):
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val == expected_val)

@pytest.mark.parameterize("path, atol",
[("raw_data_1/detector_1/counts", 1E-3),
("raw_data_1/detector_1/time_channels_1/time_of_flight", 1E-2),
("raw_data_1/instrument/dae/time_channels_1/time_of_flight_raw", 100),
("raw_data_1/instrument/detector_1/delt", 1E-3),

])
def test_same_array(nacl_nxs, nacl_gen_nxs, path, atol):
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val.shape == expected_val.shape)
    assert(np.allclose(gen_val, expected_val, atol=atol))


