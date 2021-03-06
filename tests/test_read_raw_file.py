import pytest


def test_read_summary_info(nacl_raw_reader):

    assert nacl_raw_reader.summary.header.instrument_name == b"SXD"
    assert nacl_raw_reader.summary.header.run_number == b"33298"
    assert nacl_raw_reader.summary.header.user_name == b"Team SXD            "
    assert nacl_raw_reader.summary.header.title == b"NaCl sphere 6mm diameter"
    assert nacl_raw_reader.summary.header.start_date == b"17-DEC-2020 "
    assert nacl_raw_reader.summary.header.start_time == b"09:37:12"
    assert nacl_raw_reader.summary.header.run_duration == b"   30.02"
    assert nacl_raw_reader.summary.header.format_number == 2

    assert nacl_raw_reader.info_positions.run == 32
    assert nacl_raw_reader.info_positions.instrument == 126
    assert nacl_raw_reader.info_positions.sample_environment == 225724
    assert nacl_raw_reader.info_positions.dae == 225822
    assert nacl_raw_reader.info_positions.time_channels == 451407
    assert nacl_raw_reader.info_positions.user == 453517
    assert nacl_raw_reader.info_positions.data == 453520
    assert nacl_raw_reader.info_positions.log == 21112833
    assert nacl_raw_reader.info_positions.file_end == 21112835

    assert nacl_raw_reader.summary.format[:] == [0, 1, 33298]


def test_read_run_info(nacl_raw_reader):

    assert nacl_raw_reader.run.title[:] == b"NaCl sphere 6mm diameter RT j:14,14"

    assert nacl_raw_reader.run.user.user == b"Team SXD            "
    assert nacl_raw_reader.run.user.daytime_phone == b"                    "
    assert nacl_raw_reader.run.user.daytime_phone_2 == b"                    "
    assert nacl_raw_reader.run.user.nighttime_phone == b"                    "
    assert nacl_raw_reader.run.user.institute == b"(Unspecified)       "

    assert nacl_raw_reader.run.params.run_duration == 557
    assert nacl_raw_reader.run.params.duration_units == 1
    assert nacl_raw_reader.run.params.duration_freq == 0
    assert nacl_raw_reader.run.params.dump_interval == 0
    assert nacl_raw_reader.run.params.dump_units == 0
    assert nacl_raw_reader.run.params.dump_freq == 0
    assert nacl_raw_reader.run.params.frequency == 1
    assert nacl_raw_reader.run.params.good_frames == 27882
    assert nacl_raw_reader.run.params.raw_frames == 27882
    assert nacl_raw_reader.run.params.requested_duration == 0
    assert nacl_raw_reader.run.params.duration_seconds == 557
    assert nacl_raw_reader.run.params.monitor_sum_1 == 39949211
    assert nacl_raw_reader.run.params.monitor_sum_2 == 0
    assert nacl_raw_reader.run.params.monitor_sum_3 == 0
    assert nacl_raw_reader.run.params.end_date == b"17-DEC-2020 "
    assert nacl_raw_reader.run.params.end_time == b"09:46:29"
    assert nacl_raw_reader.run.params.proposal_num == 10000


def test_read_instrument_info(nacl_raw_reader):

    assert nacl_raw_reader.instrument.version.value == 2
    assert nacl_raw_reader.instrument.name[:] == b"SXD     "

    assert nacl_raw_reader.instrument.params.chopper_freq_1 == pytest.approx(0.0)
    assert nacl_raw_reader.instrument.params.chopper_freq_2 == pytest.approx(0.0)
    assert nacl_raw_reader.instrument.params.chopper_freq_3 == pytest.approx(0.0)
    assert nacl_raw_reader.instrument.params.chopper_delay_1 == 0
    assert nacl_raw_reader.instrument.params.chopper_delay_2 == 0
    assert nacl_raw_reader.instrument.params.chopper_delay_3 == 0
    assert nacl_raw_reader.instrument.params.chopper_delay_err_1 == 0
    assert nacl_raw_reader.instrument.params.chopper_delay_err_2 == 0
    assert nacl_raw_reader.instrument.params.chopper_delay_err_3 == 0
    assert nacl_raw_reader.instrument.params.i_chopsiz == pytest.approx(0.0)
    assert nacl_raw_reader.instrument.params.aperture_c2 == 0
    assert nacl_raw_reader.instrument.params.aperture_c3 == 0
    assert nacl_raw_reader.instrument.params.chopper_status_1 == 0
    assert nacl_raw_reader.instrument.params.chopper_status_2 == 0
    assert nacl_raw_reader.instrument.params.chopper_status_3 == 0
    assert nacl_raw_reader.instrument.params.main_shutter == 0
    assert nacl_raw_reader.instrument.params.thermal_shutter == 0
    assert nacl_raw_reader.instrument.params.beam_aperture_h == pytest.approx(0.0)
    assert nacl_raw_reader.instrument.params.beam_aperture_v == pytest.approx(0.0)
    assert nacl_raw_reader.instrument.params.scattering_pos == 0
    assert nacl_raw_reader.instrument.params.moderator_type == 0
    assert nacl_raw_reader.instrument.params.detector_tank_vacuum == 0
    assert nacl_raw_reader.instrument.params.L1_scattering_length == pytest.approx(0.0)
    assert nacl_raw_reader.instrument.params.rotor_freq == 0
    assert nacl_raw_reader.instrument.params.rotor_energy == pytest.approx(0.0)
    assert nacl_raw_reader.instrument.params.rotor_phase == pytest.approx(0.0)
    assert nacl_raw_reader.instrument.params.rotor_slit_package == 0
    assert nacl_raw_reader.instrument.params.slow_chop == 0
    assert nacl_raw_reader.instrument.params.LOQ_x_centre == pytest.approx(0.0)
    assert nacl_raw_reader.instrument.params.LOQ_y_centre == pytest.approx(0.0)
    assert nacl_raw_reader.instrument.params.LOQ_beam_stop == pytest.approx(0.0)
    assert nacl_raw_reader.instrument.params.LOQ_beam_stop_radius == pytest.approx(0.0)
    assert (
        nacl_raw_reader.instrument.params.LOQ_source_detector_distance
        == pytest.approx(0.0)
    )
    assert nacl_raw_reader.instrument.params.LOQ_foe_angle == pytest.approx(0.0)
    assert nacl_raw_reader.instrument.params.angle_of_incidence == pytest.approx(0.0)

    assert nacl_raw_reader.instrument.num_detectors == 45104
    assert nacl_raw_reader.instrument.monitor_detector_number[:] == [
        45101,
        45102,
        45103,
        45104,
    ]
    assert nacl_raw_reader.instrument.monitor_prescale_val[:] == [1, 1, 1, 1]
    assert len(nacl_raw_reader.instrument.spectrum_number_table[:]) == 45104
    assert len(nacl_raw_reader.instrument.hold_off_table[:]) == 45104
    assert len(nacl_raw_reader.instrument.L2_table[:]) == 45104
    assert len(nacl_raw_reader.instrument.UTn_tables_code[:]) == 45104
    assert len(nacl_raw_reader.instrument.two_theta[:]) == 45104


def test_read_sample_info(nacl_raw_reader):

    assert nacl_raw_reader.sample.env_version.value == 2
    assert nacl_raw_reader.sample.params.sample_changer_pos == 0
    assert nacl_raw_reader.sample.params.sample_type == 1
    assert nacl_raw_reader.sample.params.sample_geom == 1
    assert nacl_raw_reader.sample.params.sample_thickness == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.sample_height == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.sample_width == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.sample_omega_angle == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.sample_chi_angle == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.sample_phi_angle == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.scattering_geom == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.sample_coh_scattering_xsec == pytest.approx(
        0.0
    )
    assert nacl_raw_reader.sample.params.sample_incoh_xsec == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.sample_absorb_xsec == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.sample_number_density == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.can_wall_thickness == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.can_coh_scattering_xsec == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.can_cs_inc == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.can_cs_abs == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.can_number_density == pytest.approx(0.0)
    assert nacl_raw_reader.sample.params.sample_name_chem_form == b""
    assert nacl_raw_reader.sample.params.e_equip == 0
    assert nacl_raw_reader.sample.params.e_eqname == 0

    assert nacl_raw_reader.sample.num_sample_env_params.value == 1
    assert len(nacl_raw_reader.sample.environment[:]) == 1

    assert nacl_raw_reader.sample.environment[0].sample_env_name == b""
    assert nacl_raw_reader.sample.environment[0].sep_value == 0
    assert nacl_raw_reader.sample.environment[0].sep_exponent == 0
    assert nacl_raw_reader.sample.environment[0].sample_env_units == b""
    assert nacl_raw_reader.sample.environment[0].sep_low_trip == 0
    assert nacl_raw_reader.sample.environment[0].sep_high_trip == 0
    assert nacl_raw_reader.sample.environment[0].sep_cur_val == 0
    assert nacl_raw_reader.sample.environment[0].sample_env_status == 0
    assert nacl_raw_reader.sample.environment[0].sample_env_ctrl == 0
    assert nacl_raw_reader.sample.environment[0].sample_env_run_ctrl == 0
    assert nacl_raw_reader.sample.environment[0].sample_env_log == 0
    assert nacl_raw_reader.sample.environment[0].sample_env_stable == pytest.approx(0.0)
    assert nacl_raw_reader.sample.environment[0].monitor_repeat_period == pytest.approx(
        0.0
    )
    assert nacl_raw_reader.sample.environment[0].camac_loc_N == 0
    assert nacl_raw_reader.sample.environment[0].camac_loc_A == 0
    assert nacl_raw_reader.sample.environment[0].camac_offset == 0
    assert nacl_raw_reader.sample.environment[0].camac_register_group == 0
    assert nacl_raw_reader.sample.environment[0].pre_proc_routine_num == 0
    assert nacl_raw_reader.sample.environment[0].camac_values[:] == [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ]


def test_read_dae_info(nacl_raw_reader):

    assert nacl_raw_reader.dae.version.value == 2
    assert nacl_raw_reader.dae.params.word_length == 4
    assert nacl_raw_reader.dae.params.mem_size == 328725240
    assert nacl_raw_reader.dae.params.ppp_min_val == 0
    assert nacl_raw_reader.dae.params.ppp_good_high == 0
    assert nacl_raw_reader.dae.params.ppp_good_low == 17271872
    assert nacl_raw_reader.dae.params.ppp_raw_high == 0
    assert nacl_raw_reader.dae.params.ppp_raw_low == 17271872
    assert nacl_raw_reader.dae.params.neutron_good_high == 0
    assert nacl_raw_reader.dae.params.neutron_good_low == 0
    assert nacl_raw_reader.dae.params.neutron_raw_high == 0
    assert nacl_raw_reader.dae.params.neutron_raw_low == 0
    assert nacl_raw_reader.dae.params.neutron_gate_t1 == 0
    assert nacl_raw_reader.dae.params.neutron_gate_t2 == 0
    assert nacl_raw_reader.dae.params.mon1_detector == 0
    assert nacl_raw_reader.dae.params.mon1_module == 4
    assert nacl_raw_reader.dae.params.mon1_crate == 0
    assert nacl_raw_reader.dae.params.mon1_mask == 16384
    assert nacl_raw_reader.dae.params.mon2_detector == 0
    assert nacl_raw_reader.dae.params.mon2_module == 5
    assert nacl_raw_reader.dae.params.mon2_crate == 0
    assert nacl_raw_reader.dae.params.mon2_mask == 20480
    assert nacl_raw_reader.dae.params.total_good_events_high == 0
    assert nacl_raw_reader.dae.params.total_good_events_low == 0
    assert nacl_raw_reader.dae.params.frame_sync_delay == 124
    assert nacl_raw_reader.dae.params.frame_sync_origin == 0
    assert nacl_raw_reader.dae.params.secondary_master_pulse == 0
    assert nacl_raw_reader.dae.params.external_vetoes[:] == [
        0,
        0,
        0,
    ]
    assert nacl_raw_reader.dae.params.num_shifted_time_regimes == 0
    assert nacl_raw_reader.dae.params.time_regime_shift_value[:] == [
        0,
        0,
        0,
    ]

    assert len(nacl_raw_reader.dae.detector_crate_nums) == 45104
    assert len(nacl_raw_reader.dae.detector_module_nums) == 45104
    assert len(nacl_raw_reader.dae.detector_module_positions) == 45104
    assert len(nacl_raw_reader.dae.detector_time_regimes) == 45104
    assert len(nacl_raw_reader.dae.detector_user_nums) == 45104


def test_read_time_channel_info(nacl_raw_reader):

    assert nacl_raw_reader.time_channel.version == 1
    assert nacl_raw_reader.time_channel.num_time_channels == 1821
    assert len(nacl_raw_reader.time_channel.params[:]) == 4
    assert nacl_raw_reader.time_channel.num_spectra == 45104
    assert len(nacl_raw_reader.time_channel.raw_boundaries[:]) == 1822
    assert len(nacl_raw_reader.time_channel.boundaries[:]) == 1822


def test_read_user_info(nacl_raw_reader):

    assert nacl_raw_reader.user.version.value == 1
    assert nacl_raw_reader.user.length == 1
    assert nacl_raw_reader.user.data[:] == pytest.approx([0.0])


def test_read_data_info(nacl_raw_reader):
    assert nacl_raw_reader.data.version.value == 2
    assert nacl_raw_reader.data.data_header.compression_type == 1
    assert nacl_raw_reader.data.data_header.spectrum_array_offset == 33
    assert nacl_raw_reader.data.data_header.equiv_v1_filesize == 645584
    assert len(nacl_raw_reader.data.ddes[:]) == 45105
    assert len(nacl_raw_reader.data.compressed_data[:]) == 82181310


def test_read_log_info(nacl_raw_reader):
    assert nacl_raw_reader.log.version[:] == [2, 0]
