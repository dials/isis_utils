import pytest
import h5py
from isis_raw_reader import IsisRawReader
from os.path import join, isfile

@pytest.fixture(scope="session")
def nacl_gen_nxs(dials_data, tmpdir_factory):

    location = dials_data("isis_sxd_example_data")
    nacl_filename = join(location, "sxd_nacl_run.raw")
    assert(isfile(nacl_filename)), f"{nacl_filename} not found"

    isis_raw_reader = IsisRawReader()
    isis_raw_reader.read_raw_file(nacl_filename)
    nxs_filename = tmpdir_factory.mktemp("data").join("nacl_gen.nxs")
    isis_raw_reader.output_tofraw_file(nxs_filename)

    return h5py.File(nxs_filename, "r")

@pytest.fixture(scope="session")
def nacl_nxs(dials_data):

    location = dials_data("isis_sxd_example_data")
    nacl_filename = join(location, "sxd_nacl_run.nxs")
    assert(isfile(nacl_filename)), f"{nacl_filename} not found"

    return h5py.File(nacl_filename, "r")
