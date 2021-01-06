import pytest
import numpy as np

def test_load_beamline(nacl_nxs, nacl_gen_nxs):
    gen_val = nacl_gen_nxs['raw_data_1']['beamline'][:]
    expected_val = nacl_nxs['raw_data_1']['beamline'][:]
    assert(gen_val == expected_val)

def test_load_collection_time(nacl_nxs, nacl_gen_nxs):
    gen_val = nacl_gen_nxs['raw_data_1']['collection_time'][:]
    expected_val = nacl_nxs['raw_data_1']['collection_time'][:]
    assert(gen_val == expected_val)

def test_load_definition(nacl_nxs, nacl_gen_nxs):
    gen_val = nacl_gen_nxs['raw_data_1']['definition'][:]
    expected_val = nacl_nxs['raw_data_1']['definition'][:]
    assert(gen_val == expected_val)

def test_load_definition_local(nacl_nxs, nacl_gen_nxs):
    gen_val = nacl_gen_nxs['raw_data_1']['definition_local'][:]
    expected_val = nacl_nxs['raw_data_1']['definition_local'][:]
    assert(gen_val == expected_val)

def test_load_detector_1_counts(nacl_nxs, nacl_gen_nxs):
    gen_val = nacl_gen_nxs['raw_data_1']['detector_1']['counts'][:]
    expected_val = nacl_nxs['raw_data_1']['detector_1']['counts'][:]
    assert(gen_val.shape == expected_val.shape)
    assert(np.allclose(gen_val, expected_val))

def test_load_detector_1_period_index(nacl_nxs, nacl_gen_nxs):
    gen_val = nacl_gen_nxs['raw_data_1']['detector_1']['period_index'][:]
    expected_val = nacl_nxs['raw_data_1']['detector_1']['period_index'][:]
    assert(gen_val == expected_val)

# detector_1 spectrum_index not tested due to known inconsistency
# (see README.md)

def test_load_detector_1_time_of_flight(nacl_nxs, nacl_gen_nxs):
    gen_val = nacl_gen_nxs['raw_data_1']['detector_1']['time_of_flight'][:]
    expected_val = nacl_nxs['raw_data_1']['detector_1']['time_of_flight'][:]
    assert(gen_val.shape == expected_val.shape)
    assert(np.allclose(gen_val, expected_val, atol=1E-2))
