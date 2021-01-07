import pytest
from isis_raw_reader import IsisRawReader
import numpy as np

def test_read_summary_info(nacl_raw_reader):

    assert(nacl_raw_reader.summary.header.instrument_name == b"SXD")
    assert(nacl_raw_reader.summary.header.run_number == b"33298")
    assert(nacl_raw_reader.summary.header.user_name == b"Team SXD            ")
    assert(nacl_raw_reader.summary.header.title == b"NaCl sphere 6mm diameter")
    assert(nacl_raw_reader.summary.header.start_date == b"17-DEC-2020 ")
    assert(nacl_raw_reader.summary.header.start_time == b"09:37:12")
    assert(nacl_raw_reader.summary.header.run_duration == b"   30.02")
    assert(nacl_raw_reader.summary.header.format_number == 2)

    assert(nacl_raw_reader.info_positions.run == 32)
    assert(nacl_raw_reader.info_positions.instrument == 126)
    assert(nacl_raw_reader.info_positions.sample_environment == 225724)
    assert(nacl_raw_reader.info_positions.dae == 225822)
    assert(nacl_raw_reader.info_positions.time_channels == 451407)
    assert(nacl_raw_reader.info_positions.user == 453517)
    assert(nacl_raw_reader.info_positions.data == 453520)
    assert(nacl_raw_reader.info_positions.log == 21112833)
    assert(nacl_raw_reader.info_positions.file_end == 21112835)

    assert(nacl_raw_reader.summary.format[:] == [0, 1, 33298])

#def test_read_run_info(nacl_raw_reader)
#def test_read_instrument_info(nacl_raw_reader)
#def test_read_sample_info(nacl_raw_reader)
#def test_read_dae_info(nacl_raw_reader)
#def test_read_time_channel_info(nacl_raw_reader)
#def test_read_user_info(nacl_raw_reader)
#def test_read_data_info(nacl_raw_reader)
#def test_read_log_info(nacl_raw_reader)
