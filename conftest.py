import pytest
import h5py
from isis_raw_reader import IsisRawReader
from os.path import join, isfile

def nacl_raw_reader(dials_data):

    location = dials_data("isis_sxd_example_data")
    nacl_filename = join(location, "sxd_nacl_run.raw")
    assert(isfile(nacl_filename)), f"{nacl_filename} not found"

    raw_reader = IsisRawReader()
    raw_reader.read_raw_file(nacl_filename)

    return raw_reader
