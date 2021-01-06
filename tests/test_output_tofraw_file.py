import pytest
import numpy as np

"""

The following params are not tested 
as an equivalent in .raw was not found:

instrument/dae/vetos/fermi_chopper3
instrument/dae/vetos/fifo
instrument/dae/vetos/msmode
instrument/dae/vetos/writing_table_file
instrument/detector_table_file
instrument/spectra_table_file
instrument/dae/type
instrument/dae/vetos/ISIS_50Hz
instrument/dae/vetos/TS2_pulse
instrument/dae/vetos/ext2
instrument/dae/vetos/ext3

The following params are not tested
as they are hard links to other params:

instrument/detector_1/counts
instrument/dae/time_channels_1/time_of_flight

The following params are not tested due to 
known inconsistencies between file formats:

detector_1/spectrum_index
framelog

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
])
def test_same_value(nacl_nxs, nacl_gen_nxs, path):
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val == expected_val)

@pytest.mark.parameterize("path, atol",
[("raw_data_1/detector_1/counts", 1E-3),
("raw_data_1/detector_1/time_of_flight", 1E-2),
("raw_data_1/instrument/dae/time_channels_1/time_of_flight_raw", 100),
])
def test_same_array(nacl_nxs, nacl_gen_nxs, path, atol):
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val.shape == expected_val.shape)
    assert(np.allclose(gen_val, expected_val, atol=atol))


