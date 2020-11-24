from ctypes import Structure, Union, c_char, c_int, c_float, c_long, c_char_p, c_ulong

class IsisRawReader:

    """
    Class to read from ISIS .raw files and covert to other formats (e.g. NeXus).
    (see https://www.isis.stfc.ac.uk/Pages/ISIS-Raw-File-Format.aspx)
    
    .raw files can also be read using 
    - OpenGenie (http://www.opengenie.org/Main_Page)
    - Mantid (https://www.mantidproject.org/)
    
    The purpose of this class is to give a simple standalone Python implementation.
    
    Main source of truth when writing this was the Mantid source code
    (https://github.com/mantidproject/mantid//Framework/DataHandling/src/LoadRaw)

    Raw files broadly contain the sections:

    - Summary
    - Run
    - Instrument
    - Sample
    - DAE
    - Time Channel
    - User
    - Data
    - Log

    These are stored as nested classes of the IsisRawReader, and have specific locations
    within the file, stored in self.info_positions.    

    Usage:
    raw_file_path = "path/to/raw/file.raw"
    isis_raw_reader = IsisRawReader()
    isis_raw_reader.read_raw_file(raw_file_path)

    """

    def __init__(self):
        self.raw_file = None
        self.summary = self.SummaryInfo()
        self.run = self.RunInfo()
        self.instrument = self.InstrumentInfo()
        self.sample = self.SampleInfo()
        self.dae = self.DAEInfo()
        self.time_channel = self.TimeChannelInfo()
        self.user = self.UserInfo()
        self.data = self.DataInfo()
        self.log = self.RawLog()
        self.info_positions = self.RawSectionPositions()

    class SummaryInfo:
        class RawHeader(Structure):
            _fields_ = [ 
                ("instrument_name", c_char * 3),
                ("run_number", c_char * 5),
                ("user_name", c_char * 20),
                ("title", c_char * 24),
                ("start_date", c_char * 12),
                ("start_time", c_char * 8),
                ("run_duration", c_char * 8),
                ("format_number", c_int)
         ]
       
        def __init__(self):
            self.header = self.RawHeader()
            self.format_version = c_int()
            self.format = (c_int*3)()  # 0 = by TC, 1 = by spectrum
            
    class RunInfo:

        class RawUser(Structure):
            _fields_ = [
                ("user", c_char * 20),
                ("daytime_phone", c_char * 20),
                ("daytime_phone_2", c_char * 20),
                ("nighttime_phone", c_char * 20),
                ("institute", c_char * 20),
                ("unused", (c_char*20) * 3)
            ]

        class RawRun(Structure):
            _fields_ = [
                ("run_duration", c_int),
                ("duration_units", c_int),
                ("duration_freq", c_int),
                ("dump_interval", c_int),
                ("dump_units", c_int),
                ("dump_freq", c_int),
                ("frequency", c_int),
                ("good_proton_charge", c_float),
                ("total_proton_charge", c_float),
                ("good_frames", c_int),
                ("raw_frames", c_int),
                ("requested_duration", c_int),
                ("duration_seconds", c_int),
                ("monitor_sum_1", c_int),
                ("monitor_sum_2", c_int),
                ("monitor_sum_3", c_int),
                ("end_date", c_char * 12), 
                ("end_time", c_char * 8),
                ("proposal_num", c_int),
                ("unused", c_int * 10)
            ]
            
        def __init__(self):
            self.title = (c_char * 80)()
            self.user = self.RawUser()
            self.params = self.RawRun()
            self.version = c_int(1)
            self.number = c_int()
            
            
    class InstrumentInfo:

        class RawInstrument(Structure):
            _fields_ = [
                ("chopper_freq_1", c_float),  # (Hz)
                ("chopper_freq_2", c_float),
                ("chopper_freq_3", c_float),
                ("chopper_delay_1", c_int),  # (us)
                ("chopper_delay_2", c_int), 
                ("chopper_delay_3", c_int),
                ("chopper_delay_err_1", c_int),  # (us)
                ("chopper_delay_err_2", c_int),
                ("chopper_delay_err_3", c_int),
                ("i_chopsiz", c_float),  # unknown param
                ("aperture_c2", c_int),  # unknown param
                ("aperture_c3", c_int),  # unknown param
                ("chopper_status_1", c_int),
                ("chopper_status_2", c_int),
                ("chopper_status_3", c_int),
                ("main_shutter", c_int),  # open = 1
                ("thermal_shutter", c_int),  # open = 1
                ("beam_aperture_h", c_float),  # (mm)
                ("beam_aperture_v", c_float),  # (mm)
                ("scattering_pos", c_int),
                ("moderator_type", c_int),
                ("detector_tank_vacuum", c_int),
                ("L1_scattering_length", c_float),
                ("rotor_freq", c_int),
                ("rotor_energy", c_float),
                ("rotor_phase", c_float),
                ("rotor_slit_package", c_int),  # 0="unknown",1="L",2="Med",3="Hi"
                ("slow_chop", c_int),  # on = 1
                ("LOQ_x_centre", c_float),
                ("LOQ_y_centre", c_float),
                ("LOQ_beam_stop", c_float),
                ("LOQ_beam_stop_radius", c_float),
                ("LOQ_souce_detector_distance", c_float),
                ("LOQ_foe_angle", c_float),
                ("angle_of_incidence", c_float),
                ("unused", c_int * 29)
            ]
            
        def __init__(self):
            self.params = self.RawInstrument()
            self.version = c_int()
            self.name = (c_char * 8)()
            self.num_detectors = (c_int * 3)()
            self.num_monitors = c_int()
            self.num_user_tables = c_int()
            self.monitor_detector_number = (c_int * 4)()
            self.monitor_prescale_val = (c_int * 4)()
            self.spectrum_number_table = c_int()
            self.hold_off_table = c_int()
            self.L2_table = c_float()
            self.UTn_tables_code = c_int()
            self.two_theta = c_float()
            self.num_user_tables = c_float()
            
    class SampleInfo:
        
        class RawSample(Structure):
            _fields_ = [
                ("sample_changer_pos", c_int),
                ("sample_type", c_int),  # 1 = sample + can, 2 = empty can
                ("sample_geom", c_int),
                ("sample_thickness", c_float),  # (mm)
                ("sample_height", c_float),  # (mm)
                ("sample_width", c_float),  # (mm)
                ("sample_omega_angle", c_float),  # (degrees)
                ("sample_chi_angle", c_float),  # (degrees)
                ("sample_phi_angle", c_float),  # (degrees)
                ("scattering_geom", c_float),  # 1 = trans, 2 = reflect
                ("sample_coh_scattering_xsec", c_float),  # (barn)
                ("sample_incoh_xsec", c_float), 
                ("sample_absorb_xsec", c_float),
                ("sample_number_density", c_float),  # (atoms A^-3)
                ("can_wall_thickness", c_float),  # (mm)
                ("can_coh_scattering_xsec", c_float),  # (barn)
                ("can_cs_inc", c_float),  # unknown param
                ("can_cs_abs", c_float),  # unknown param
                ("can_number_density", c_float),  # (atoms A^-3)
                ("sample_name_chem_form", c_char * 40),
                ("e_equip", c_int),  # unknown param
                ("e_eqname", c_int),  # unknown param
                ("unused", c_int * 33)
            ]

        class RawSampleEnvironment(Structure):
            _fields_ = [
                ("sample_env_name", c_char * 8),
                ("sep_value", c_int),  # unknown param
                ("sep_exponent", c_int),  # unknown param
                ("sample_env_units", c_char * 8),
                ("sep_low_trip", c_int),  # unknown param
                ("sep_high_trip", c_int),  # unknown param
                ("sep_cur_val", c_int),  # unknown param
                ("sample_env_status", c_int),  # (bool)
                ("sample_env_ctrl", c_int),  # (bool)
                ("sample_env_run_ctrl", c_int),  # (bool)
                ("sample_env_log", c_int),  # (bool)
                ("sample_env_stable", c_float),  # (Hz)
                ("monitor_repeat_period", c_float),
                ("camac_loc_N", c_int),
                ("camac_loc_A", c_int),
                ("camac_offset", c_int),
                ("camac_register_group", c_int),  # 1 or 2
                ("pre_proc_routine_num", c_int),
                ("camac_values", c_int * 12)
            ]
            
        def __init__(self):
            self.params = self.RawSample()
            self.environment = self.RawSampleEnvironment()
            self.env_version = c_int()
            self.num_sample_env_params = c_int()
            
    class DAEInfo:

        class RawDAE(Structure):
            _fields_ = [
                ("word_length", c_int),
                ("mem_size", c_int),  # (bytes) **A
                ("ppp_min_val", c_int),
                ("ppp_good_high", c_int),
                ("ppp_good_low", c_int),
                ("ppp_raw_high", c_int),
                ("ppp_raw_low", c_int),
                ("neutron_good_high", c_int),
                ("neutron_good_low", c_int),
                ("neutron_raw_high", c_int),
                ("neutron_raw_low", c_int),
                ("neutron_gate_t1", c_int),
                ("neutron_gate_t2", c_int),
                ("mon1_detector", c_int),
                ("mon1_module", c_int),
                ("mon1_crate", c_int),
                ("mon1_mask", c_int),
                ("mon2_detector", c_int),
                ("mon2_module", c_int),
                ("mon2_crate", c_int),
                ("mon2_mask", c_int),
                ("total_good_events_high", c_int),
                ("total_good_events_low", c_int),
                ("frame_sync_delay", c_int),
                ("frame_sync_origin", c_int),  # 0:none/1:ext/2:int
                ("secondary_master_pulse", c_int),  # 0:en, 1:dis
                ("external_vetoes", c_int * 3),  # 0,1,2 
                ("num_shifted_time_regimes", c_int),  # for alternative PC DAE
                ("time_regime_shift_value", c_int * 3),  # for alternative PC DAE
                ("unused", c_int * 31)
            ]
            
        def __init__(self):
            self.params = self.RawDAE()
            self.dae_version = c_int()
            self.detector_crate_nums = c_int()
            self.detector_module_nums = c_int()
            self.detector_module_positions = c_int()
            self.detector_time_regimes = c_int()
            self.detector_user_nums = c_int()
            
    class TimeChannelInfo:
        
        def __init__(self):
            self.version = (c_int * 267)()
            self.num_time_regimes = c_int()
            self.num_frames_per_period = c_int()
            self.num_periods = c_int(1)
            self.period_num_per_period = (c_int * 256)()
            self.num_spectra = c_int()
            self.num_time_channels = c_int()
            self.mode = (c_int * 5)()
            self.params = ((c_float * 5) * 4)()
            self.prescale = c_int()
            self.boundaries = c_int()
        
    class UserInfo:
        
        def __init__(self):
            self.version = c_int()
            self.length = c_int()
            self.data = c_float()
            
    class DataInfo:
        
        def __init__(self):
            self.data_version = c_int()
            self.data_header = self.RawDataHeader()
            self.ddes = self.RawDDES()
            self.compressed_data = c_long()

            
        class RawDataHeader(Structure):
            _fields_ = [
                ("compression_type", c_int),  # 0=none, 1=byte relative
                ("unused1", c_int),
                ("spectrum_array_offset", c_int),
                ("compression_ratio", c_float),
                ("file_compression_ratio", c_float),
                ("equiv_v1_filesize", c_int),
                ("unused2", c_int * 26)
            ]

        class RawDDES(Structure):
            _fields_ = [
                ("num_spec_compressed_words", c_int),
                ("compressed_spec_offset", c_int)
            ]

            
            
    class RawSectionPositions(Structure):
        _fields_ = [
            ("run", c_int),
            ("instrument", c_int),
            ("sample_environment", c_int),
            ("dae", c_int),
            ("time_channels", c_int),
            ("user", c_int),
            ("data", c_int),
            ("log", c_int),
            ("file_end", c_int)
        ]

    class RawLog(Structure):
        _fields_ = [
            ("version", c_int * 2),
            ("num_lines", c_int),
        ] 
            

        
    def decompress_data(self, data_in, data_out):
        
        # Byte code used to indicate number cannot be 
        # written +/-127 of previous value
        nan_val = -128  
        
        def number_in_byte(byte_val, nan_val):
            return byte_val != nan_val
        
        def byte_to_int(byte):
            result = 0
            for b in byte:
                result = result * 256 + int(b)
            if result > 127:
                return (256-result)*(-1)
            return result
        
        class BytePacket(Union):
            _fields_ = [
                ("pos", c_int),
                ("val", c_char * 4)
            ]
        
        counter = 0
        byte_packet = BytePacket()
        byte_packet.pos = 0

        for i in range(len(data_out)):
            
            assert(counter < len(data_in)), "channel_counter out of range"
            
            new_byte = byte_to_int(data_in[counter])
            if number_in_byte(new_byte, nan_val):
                byte_packet.pos += new_byte
                counter += 1
            else:
                assert(counter + 4 < len(data_in)),\
                "byte out of range"
                
                byte_packet.val = data_in[counter:counter+4]             
                counter += 5       
            data_out[i] = byte_packet.pos
            
        return data_out
    
    def read_summary_info(self):
        self.read_into_buffer(self.summary.header)
        self.read_into_buffer(self.info_positions)
        self.read_into_buffer(self.summary.format)
        
    def read_run_info(self):
        self.read_into_buffer(self.run.title)
        self.read_into_buffer(self.run.user)
        self.read_into_buffer(self.run.params)
        
    def read_instrument_info(self):
        
        def update_param_sizes(new_size):
            self.instrument.spectrum_number_table = (c_int * new_size)()
            self.instrument.hold_off_table = (c_float * new_size)()
            self.instrument.L2_table = (c_float * new_size)()
            self.instrument.UTn_tables_code = (c_float * new_size)()
            self.instrument.two_theta = (c_float * new_size)()
            
        self.read_into_buffer(self.instrument.version)
        self.read_into_buffer(self.instrument.name)
        self.read_into_buffer(self.instrument.params)
        self.read_into_buffer(self.instrument.num_detectors)
        self.instrument.num_detectors = self.instrument.num_detectors[0]
        self.read_into_buffer(self.instrument.monitor_detector_number)
        self.read_into_buffer(self.instrument.monitor_prescale_val)
        
        update_param_sizes(self.instrument.num_detectors)
        
        self.read_into_buffer(self.instrument.spectrum_number_table)
        self.read_into_buffer(self.instrument.hold_off_table)
        self.read_into_buffer(self.instrument.L2_table)
        self.read_into_buffer(self.instrument.UTn_tables_code)
        self.read_into_buffer(self.instrument.two_theta)
        
    def read_sample_info(self):
        self.read_into_buffer(self.sample.env_version)
        self.read_into_buffer(self.sample.params)
        self.read_into_buffer(self.sample.num_sample_env_params)
        self.sample.environment = \
        (self.sample.RawSampleEnvironment * self.sample.num_sample_env_params.value)()
        self.read_into_buffer(self.sample.environment)
        
    def read_dae_info(self):
        
        def update_param_sizes(new_size):
            self.dae.detector_crate_nums = (c_int * new_size)()
            self.dae.detector_module_nums = (c_int * new_size)()
            self.dae.detector_module_positions = (c_int * new_size)()
            self.dae.detector_time_regimes = (c_int * new_size)()
            self.dae.detector_user_nums = (c_int * new_size)()

        update_param_sizes(self.instrument.num_detectors)
        
        self.read_into_buffer(self.dae.dae_version)
        self.read_into_buffer(self.dae.params)
        self.read_into_buffer(self.dae.detector_crate_nums)
        self.read_into_buffer(self.dae.detector_module_nums)
        self.read_into_buffer(self.dae.detector_module_positions)
        self.read_into_buffer(self.dae.detector_time_regimes)
        self.read_into_buffer(self.dae.detector_user_nums)
        
    def read_time_channel_info(self):
        
        # The version buffer is not currectly defined
        # But in lieu of original documentation of the file structure
        # the hack below is necessary
        self.read_into_buffer(self.time_channel.version)
        self.time_channel.num_time_channels = self.time_channel.version[-6]
        self.time_channel.version = self.time_channel.version[0]
        
        self.read_into_buffer(self.time_channel.params)
        self.time_channel.num_spectra = self.instrument.num_detectors
        self.read_into_buffer(self.time_channel.prescale)
        
        self.time_channel.boundaries = \
        (c_int * (self.time_channel.num_time_channels + 1))()
        self.read_into_buffer(self.time_channel.boundaries)
        
    def read_user_info(self):
        self.read_into_buffer(self.user.version)
        
        # Stray byte here
        dummy_int = c_int()
        self.read_into_buffer(dummy_int)
        
        self.user.length = self.info_positions.data - self.info_positions.user - 2
        self.user.data = (c_float * self.user.length)()
        self.read_into_buffer(self.user.data)
        
    def read_data_info(self):
        
        self.read_into_buffer(self.data.data_version)
        self.read_into_buffer(self.data.data_header) 
        
        compression_type = self.data.data_header.compression_type
        assert(compression_type == 1),\
        "Unable to decompress data (compression_type is {} but must be 1)".format(compression_type)
        
        num_ddes = self.time_channel.num_periods.value * \
                   (self.time_channel.num_spectra + 1) 
        
        self.data.ddes = \
            (self.data.RawDDES * num_ddes)()
        self.read_into_buffer(self.data.ddes)
        
        ddes_size = self.time_channel.num_time_channels + 1
        self.data.compressed_data = \
        (c_ulong * (num_ddes * (ddes_size)))()
        
        for i in range(num_ddes):
            nwords = self.data.ddes[i].num_spec_compressed_words
            read_buffer = (c_char * (4 * nwords))()
            self.read_into_buffer(read_buffer)
            self.data.compressed_data[i * ddes_size : (i+1) * ddes_size] = self.decompress_data(read_buffer, 
                            self.data.compressed_data[i * ddes_size : (i+1) * ddes_size])
            
        
    def read_log_info(self):
        self.read_into_buffer(self.log.version)

    def open_raw_file(self, raw_file_path):
        assert(raw_file_path.split(".")[-1].lower() == "raw"), "raw_file must be a .raw file"
        self.raw_file = open(raw_file_path, "rb")
        
    def read_raw_file(self, raw_file_path):
        
        self.open_raw_file(raw_file_path)
        
        if self.raw_file.tell() != 0:
            print("Reseting raw_file cursor from {} to 0 before reading".format(self.raw_file.tell()))
            self.reset_file_cursor()
            
        self.read_summary_info()
        self.read_run_info()
        self.read_instrument_info()
        self.read_sample_info()
        self.read_dae_info()
        self.read_time_channel_info()
        self.read_user_info()
        self.read_data_info()
        self.read_log_info()
        self.reset_file_cursor()
    
    def read_into_buffer(self, buffer):
        self.raw_file.readinto(buffer)
        
    def move_file_cursor(self, pos):
        self.raw_file.seek(pos)
        
    def reset_file_cursor(self):
        self.move_file_cursor(0)
        
    def print_struct(self, data_struct):
        try:
            for field_name, field_type in data_struct._fields_:
                print(field_name, getattr(data_struct, field_name))
        except AttributeError:
            try:
                for elem in data_struct:
                    for field_name, field_type in elem._fields_:
                        print(field_name, getattr(elem, field_name))
            except AttributeError:
                print("Cannot read data_struct (can only print ctypes Structure instances)")
            