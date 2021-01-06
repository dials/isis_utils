import pytest
import numpy as np

def test_load_beamline(nacl_nxs, nacl_gen_nxs):
    path = "raw_data_1/beamline"
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val == expected_val)

def test_load_collection_time(nacl_nxs, nacl_gen_nxs):
    path = "raw_data_1/collection_time"
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val == expected_val)

def test_load_definition(nacl_nxs, nacl_gen_nxs):
    path = "raw_data_1/definition"
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val == expected_val)

def test_load_definition_local(nacl_nxs, nacl_gen_nxs):
    path = "raw_data_1/definition_local"
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val == expected_val)

def test_load_detector_1_counts(nacl_nxs, nacl_gen_nxs):
    path = "raw_data_1/detector_1/counts"
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val.shape == expected_val.shape)
    assert(np.allclose(gen_val, expected_val))

def test_load_detector_1_period_index(nacl_nxs, nacl_gen_nxs):
    path = "raw_data_1/detector_1/period_index"
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val == expected_val)

# detector_1 spectrum_index not tested due to known inconsistency
# (see README.md)

def test_load_detector_1_time_of_flight(nacl_nxs, nacl_gen_nxs):
    path = "raw_data_1/detector_1/time_of_flight"
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val.shape == expected_val.shape)
    assert(np.allclose(gen_val, expected_val, atol=1E-2))

def test_load_duration(nacl_nxs, nacl_gen_nxs):
    path = "raw_data_1/duration"
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val == expected_val)

def test_load_end_time(nacl_nxs, nacl_gen_nxs):
    path = "raw_data_1/end_time"
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val == expected_val)

def test_load_end_time(nacl_nxs, nacl_gen_nxs):
    path = "raw_data_1/experiment_identifier"
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path[:]
    assert(gen_val == expected_val)

# framelog is empty in observed TOFRAW files and so is not tested

def test_load_end_time(nacl_nxs, nacl_gen_nxs):
    path = "raw_data_1/good_frames"
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val == expected_val)

# instrument detector_table_file is not tested as equivalent in .raw not found

def test_load_instrument_dae_period_index(nacl_nxs, nacl_gen_nxs):
    path = "raw_data_1/good_frames"
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val == expected_val)

# instrument spectra_table_file is not tested as equivalent in .raw not found

def test_load_instrument_dae_time_of_flight_raw(nacl_nxs, nacl_gen_nxs):
    path = "raw_data_1/instrument/dae/time_of_flight_raw"
    gen_val = nacl_gen_nxs[path][:]
    expected_val = nacl_nxs[path][:]
    assert(gen_val.shape == expected_val.shape)
    assert(np.allclose(gen_val, expected_val, atol=100))
