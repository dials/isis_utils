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

def test_read_run_info(nacl_raw_reader):

    assert(nacl_raw_reader.run.title[:] == b"NaCl sphere 6mm diameter RT j:14,14")
    
    assert(nacl_raw_reader.run.user.user == b"Team SXD            ")
    assert(nacl_raw_reader.run.user.daytime_phone == b"                    ")
    assert(nacl_raw_reader.run.user.daytime_phone_2 == b"                    ")
    assert(nacl_raw_reader.run.user.nighttime_phone == b"                    ")
    assert(nacl_raw_reader.run.user.institute == b"(Unspecified)       ")

    assert(nacl_raw_reader.run.params.run_duration == 557) 
    assert(nacl_raw_reader.run.params.duration_units == 1) 
    assert(nacl_raw_reader.run.params.duration_freq == 0) 
    assert(nacl_raw_reader.run.params.dump_interval == 0) 
    assert(nacl_raw_reader.run.params.dump_units == 0) 
    assert(nacl_raw_reader.run.params.dump_freq == 0) 
    assert(nacl_raw_reader.run.params.frequency == 1) 
    assert(nacl_raw_reader.run.params.good_frames == 27882) 
    assert(nacl_raw_reader.run.params.raw_frames == 27882) 
    assert(nacl_raw_reader.run.params.requested_duration == 0) 
    assert(nacl_raw_reader.run.params.duration_seconds == 557) 
    assert(nacl_raw_reader.run.params.monitor_sum_1 == 39949211) 
    assert(nacl_raw_reader.run.params.monitor_sum_2 == 0) 
    assert(nacl_raw_reader.run.params.monitor_sum_3 == 0) 
    assert(nacl_raw_reader.run.params.end_date == b"17-DEC-2020 ") 
    assert(nacl_raw_reader.run.params.end_time == b"09:46:29") 
    assert(nacl_raw_reader.run.params.proposal_num == 10000) 
    
def test_read_instrument_info(nacl_raw_reader):

    assert(nacl_raw_reader.instrument.version.value == 2)
    assert(nacl_raw_reader.instrument.name[:] == b"SXD     ")
    
    assert(nacl_raw_reader.instrument.params.chopper_freq_1 == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.chopper_freq_2 == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.chopper_freq_3 == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.chopper_delay_1 == 0)
    assert(nacl_raw_reader.instrument.params.chopper_delay_2 == 0)
    assert(nacl_raw_reader.instrument.params.chopper_delay_3 == 0)
    assert(nacl_raw_reader.instrument.params.chopper_delay_err_1 == 0)
    assert(nacl_raw_reader.instrument.params.chopper_delay_err_2 == 0)
    assert(nacl_raw_reader.instrument.params.chopper_delay_err_3 == 0)
    assert(nacl_raw_reader.instrument.params.i_chopsiz == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.aperture_c2 == 0)
    assert(nacl_raw_reader.instrument.params.aperture_c3 == 0)
    assert(nacl_raw_reader.instrument.params.chopper_status_1 == 0)
    assert(nacl_raw_reader.instrument.params.chopper_status_2 == 0)
    assert(nacl_raw_reader.instrument.params.chopper_status_3 == 0)
    assert(nacl_raw_reader.instrument.params.main_shutter == 0)
    assert(nacl_raw_reader.instrument.params.thermal_shutter == 0)
    assert(nacl_raw_reader.instrument.params.beam_aperture_h == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.beam_aperture_v == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.scattering_pos == 0)
    assert(nacl_raw_reader.instrument.params.moderator_type == 0)
    assert(nacl_raw_reader.instrument.params.detector_tank_vacuum == 0)
    assert(nacl_raw_reader.instrument.params.L1_scattering_length == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.rotor_freq == 0)
    assert(nacl_raw_reader.instrument.params.rotor_energy == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.rotor_phase == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.rotor_slit_package == 0)
    assert(nacl_raw_reader.instrument.params.slow_chop == 0)
    assert(nacl_raw_reader.instrument.params.LOQ_x_centre == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.LOQ_y_centre == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.LOQ_beam_stop == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.LOQ_beam_stop_radius == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.LOQ_source_detector_distance == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.LOQ_foe_angle == pytest.approx(0.0))
    assert(nacl_raw_reader.instrument.params.angle_of_incidence == pytest.approx(0.0))

    assert(nacl_raw_reader.num_detectors == 45104)
    assert(nacl_raw_reader.monitor_detector_number[:] == [45101, 45102, 45103, 45104])
    assert(nacl_raw_reader.monitor_prescale_val[:] == [1, 1, 1, 1])
    assert(len(nacl_raw_reader.instrument.spectrum_number_table[:]) == 45104)
    assert(len(nacl_raw_reader.instrument.hold_off_table[:]) == 45104)
    assert(len(nacl_raw_reader.instrument.L2_table[:]) == 45104)
    assert(len(nacl_raw_reader.instrument.UTn_tables_code[:]) == 45104)
    assert(len(nacl_raw_reader.instrument.two_theta[:]) == 45104)

    

def test_read_sample_info(nacl_raw_reader):



    

def test_read_sample_info(nacl_raw_reader):



#def test_read_dae_info(nacl_raw_reader)
#def test_read_time_channel_info(nacl_raw_reader)
#def test_read_user_info(nacl_raw_reader)
#def test_read_data_info(nacl_raw_reader)
#def test_read_log_info(nacl_raw_reader)
