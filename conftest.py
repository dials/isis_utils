from os.path import isfile, join

import pytest

from isis_raw_reader import IsisRawReader


@pytest.fixture(scope="session")
def nacl_raw_reader(dials_data):

    location = dials_data("isis_sxd_example_data")
    nacl_filename = join(location, "sxd_nacl_run.raw")
    assert isfile(nacl_filename), f"{nacl_filename} not found"

    raw_reader = IsisRawReader()
    raw_reader.read_raw_file(nacl_filename)

    return raw_reader
