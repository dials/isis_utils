from ctypes import Structure, Union, c_char, c_int, c_float, c_long, c_char_p, c_ulong
import h5py
from os import path
import numpy as np

class IsisRawReader:

    """
    Class to read from ISIS .raw files and covert to other formats (e.g. NeXus).
    (see https://www.isis.stfc.ac.uk/Pages/ISIS-Raw-File-Format.aspx)
    
    .raw files can also be read using 
    - OpenGenie (http://www.opengenie.org/Main_Page)
    - Mantid (https://www.mantidproject.org/)
    
    The purpose of this class is to give a simple standalone Python implementation.
    
    The Main source of truth when writing this was the Mantid source code
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

    # Read .raw file into memory
    raw_file_path = "path/to/raw/file.raw"
    isis_raw_reader = IsisRawReader()
    isis_raw_reader.read_raw_file(raw_file_path)

    # Output .raw file as .nxs file
    output_filename = "SXD.nxs"
    isis_raw_reader.output_tofraw_file(output_filename)

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
            self.format_version = c_int(2)
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
                ("LOQ_source_detector_distance", c_float),
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
            self.monitor_detector_number = (c_int * 4)()
            self.monitor_prescale_val = (c_int * 4)()
            self.spectrum_number_table = c_int()
            self.hold_off_table = c_int()
            self.L2_table = c_float()
            self.UTn_tables_code = c_int()
            self.two_theta = c_float()
            self.num_user_tables = c_int()
            
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
            self.version = c_int()
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
            self.raw_boundaries = c_int()
            self.boundaries = c_float()

        
    class UserInfo:
        
        def __init__(self):
            self.version = c_int()
            self.length = c_int()
            self.data = c_float()
            
    class DataInfo:
        
        def __init__(self):
            self.version = c_int()
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
                assert(counter + 4 < len(data_in)), "byte out of range"
                byte_packet.val = data_in[counter+1:counter+5]             
                counter += 5       
            data_out[i] = byte_packet.pos
            
        return data_out
    
    def read_summary_info(self):
        self.read_into_buffer(self.summary.header)
        self.read_into_buffer(self.info_positions)
        self.read_into_buffer(self.summary.format)
        
    def read_run_info(self):
        self.read_into_buffer(self.run.title)
        self.run.title = self.run.title.rstrip()
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
        num_envs = self.sample.num_sample_env_params.value
        self.sample.environment = (self.sample.RawSampleEnvironment * num_envs)()
        self.read_into_buffer(self.sample.environment)
        
    def read_dae_info(self):
        
        def update_param_sizes(new_size):
            self.dae.detector_crate_nums = (c_int * new_size)()
            self.dae.detector_module_nums = (c_int * new_size)()
            self.dae.detector_module_positions = (c_int * new_size)()
            self.dae.detector_time_regimes = (c_int * new_size)()
            self.dae.detector_user_nums = (c_int * new_size)()

        update_param_sizes(self.instrument.num_detectors)
        
        self.read_into_buffer(self.dae.version)
        self.read_into_buffer(self.dae.params)
        self.read_into_buffer(self.dae.detector_crate_nums)
        self.read_into_buffer(self.dae.detector_module_nums)
        self.read_into_buffer(self.dae.detector_module_positions)
        self.read_into_buffer(self.dae.detector_time_regimes)
        self.read_into_buffer(self.dae.detector_user_nums)
        
    def read_time_channel_info(self):
        
        def get_time_channel_boundaries(raw_boundaries, prescale, s_delay, fmt_ver):

            added_delay = 4. * s_delay
            if fmt_ver.value <= 1:
                added_delay = 0
            
            boundaries = []
            for i in raw_boundaries:
                boundaries.append(i * (prescale.value / 32.) + added_delay)
            return boundaries

        # The version buffer is not currectly defined
        # But in lieu of original documentation of the file structure
        # the hack below is necessary
        self.read_into_buffer(self.time_channel.version)
        self.time_channel.num_time_channels = self.time_channel.version[-6]
        self.time_channel.version = self.time_channel.version[0]
        
        self.read_into_buffer(self.time_channel.params)
        self.time_channel.num_spectra = self.instrument.num_detectors
        self.read_into_buffer(self.time_channel.prescale)
        
        num_boundaries = self.time_channel.num_time_channels + 1
        self.time_channel.raw_boundaries = (c_int * (num_boundaries))()

        self.time_channel.boundaries = (c_float * (num_boundaries))()

        self.read_into_buffer(self.time_channel.raw_boundaries)
        boundaries_lst = get_time_channel_boundaries(self.time_channel.raw_boundaries,
                                            self.time_channel.prescale, 
                                            self.dae.params.frame_sync_delay,
                                            self.summary.format_version)
        self.time_channel.boundaries = (c_float * len(boundaries_lst))(*boundaries_lst)
        
    def get_counts_per_pixel(self, num_detectors, num_channels):

        raw_shape = (1, num_detectors, num_channels)
        raw_data = self.data.compressed_data[:]
        return np.array(raw_data).reshape(raw_shape)

    def read_user_info(self):
        self.read_into_buffer(self.user.version)
        
        # Stray byte here
        dummy_int = c_int()
        self.read_into_buffer(dummy_int)
        
        self.user.length = self.info_positions.data - self.info_positions.user - 2
        self.user.data = (c_float * self.user.length)()
        self.read_into_buffer(self.user.data)
        
    def read_data_info(self):
        
        self.read_into_buffer(self.data.version)
        self.read_into_buffer(self.data.data_header) 
        
        compression_type = self.data.data_header.compression_type
        assert(compression_type == 1), "Unable to decompress data"
        
        num_periods = self.time_channel.num_periods.value
        num_ddes = num_periods * (self.time_channel.num_spectra + 1) 
        
        self.data.ddes = (self.data.RawDDES * num_ddes)()
        self.read_into_buffer(self.data.ddes)
        
        ddes_size = self.time_channel.num_time_channels + 1
        self.data.compressed_data = (c_ulong * (num_ddes * (ddes_size)))()
        
        for i in range(num_ddes):
            nwords = self.data.ddes[i].num_spec_compressed_words
            read_buffer = (c_char * (4 * nwords))()
            self.read_into_buffer(read_buffer)
            dmin = i * ddes_size 
            dmax = (i+1) * ddes_size
            self.data.compressed_data[dmin:dmax] = self.decompress_data(read_buffer, 
                            self.data.compressed_data[dmin:dmax])
            
        
    def read_log_info(self):
        self.read_into_buffer(self.log.version)

    def open_raw_file(self, raw_file_path):
        file_ext = raw_file_path.split(".")[-1].lower()
        assert(file_ext == "raw"), "raw_file must be a .raw file"
        self.raw_file = open(raw_file_path, "rb")
        
    def read_raw_file(self, raw_file_path):
        
        self.open_raw_file(raw_file_path)
        
        if self.raw_file.tell() != 0:
            pos = self.raw_file.tell()
            print(f"Reseting raw_file cursor from {pos} to 0 before reading")
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
    
    def read_into_buffer(self, _buffer):
        self.raw_file.readinto(_buffer)
        
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
                print("Cannot read data_struct (can only print ctypes Structures)")



    def output_tofraw_file(self, output_filename):
        
        """
        Converts info in memory to TOFRaw NeXus format 
        and outputs to output_filename.
        
        The NeXus format used here follows the guide here: 
        https://www.nexusformat.org/TOFRaw.html
        
        """
        
        def load_beamline(nxs_file):
            name = "raw_data_1/beamline"
            beamline = nxs_file.create_dataset(name, (1,), dtype='|S3')
            beamline[0] = self.summary.header.instrument_name
            
        def load_collection_time(nxs_file):
            name = "raw_data_1/collection_time"
            _type = np.dtype("f4")
            collection_time = nxs_file.create_dataset(name, (1,), dtype=_type)
            collection_time[0] = self.run.params.run_duration
            
        def load_definition(nxs_file):
            name = "raw_data_1/definition"
            definition = nxs_file.create_dataset(name, (1,), dtype='|S6')
            definition[0] = b'TOFRAW'
            
        def load_definition_local(nxs_file):
            name = "raw_data_1/definition_local"
            definition_local = nxs_file.create_dataset(name, (1,), dtype='|S10')
            definition_local[0] = b'ISISTOFRAW'
            
        def load_detector_1(nxs_file, num_detectors, num_channels):
            
            grp = nxs_file.create_group("raw_data_1/detector_1")

            name = "counts"
            # 1D raw counts need reshaping per pixel with
            # monitor channels removed to be consistent with .nxs
            raw_data = self.get_counts_per_pixel(num_detectors, num_channels)
            raw_shape = raw_data.shape
            nxs_shape = (raw_shape[0], raw_shape[1] - 5, raw_shape[2] - 1)
            _type = np.dtype("i4")
            counts = grp.create_dataset(name, nxs_shape, dtype=_type)
            counts[0] = raw_data[:, 1:-4, 1:]
            
            name = "period_index"
            _type = np.dtype("i4")
            period_idx = grp.create_dataset(name, (1,), dtype=_type)
            period_idx[0] = self.time_channel.num_periods
            
            """
            NOTE:
            The .raw spectrum_number_table is inconsistent with the .nxs
            detector_1/spectrum_index. 
            The former gives a 1D array from the top right of each detector
            moving down each column (i.e for 64x64 detector [64,63,..]).
            The .nxs index just gives [i for i in range(num pixels)].
            When visualising both give the same laue plots, consistent with
            SXD2001 only when using the indexing e.g. for a 64x64 detector:
            np.arange(min_pixel_idx_range, max_pixel_idx_range).reshape(64,64).T
            """

            name = "spectrum_index"
            shape = (num_detectors - 1,)
            _type = np.dtype("i4")
            spectrum_idx = grp.create_dataset(name, shape, dtype=_type)
            spectrum_idx[:] = self.instrument.spectrum_number_table[:]
            
            name = "time_of_flight"
            shape = (num_channels,)
            _type = np.dtype("f4")
            tof = grp.create_dataset(name, shape, dtype=_type)
            tof[:] = np.array(self.time_channel.boundaries[:])
            
        def load_duration(nxs_file):
            name = "raw_data_1/duration"
            _type = np.dtype("f4")
            duration = nxs_file.create_dataset(name, (1,), dtype=_type)
            duration[0] = self.run.params.run_duration
            
        def convert_to_nxs_date(raw_date):

            date_dict = {"JAN" : "1", "FEB" : "2", "MAR" : "3", "APR": "4",
                        "MAY" : "5", "JUN" : "6", "JUL" : "7", "AUG" : "8",
                        "SEP" : "9", "OCT" : "10", "NOV" : "11", "DEC" : "12"}

            date = raw_date.decode().rstrip()
            day, month, year = date.split("-")
            return "-".join([year, date_dict[month], day])
            
        def load_end_time(nxs_file):
            name = "raw_data_1/end_time"
            end_time = nxs_file.create_dataset(name, (1,), dtype="|S19")
            raw_date = convert_to_nxs_date(self.run.params.end_date)
            raw_time = self.run.params.end_time.decode()
            end_time[0] = (raw_date + "T" + raw_time).encode()
            
        def load_experiment_identifier(nxs_file):
            name = "raw_data_1/experiment_identifier"
            exp_identifier = nxs_file.create_dataset(name, (1,), dtype="|S7")
            exp_identifier[0] = str(self.run.params.proposal_num).encode()
            
        def load_frame_log(nxs_file):
            nxs_file.create_group("raw_data_1/framelog")
            
        def load_good_frames(nxs_file):
            name = "raw_data_1/good_frames"
            _type = np.dtype("i4")
            good_frames = nxs_file.create_dataset(name, (1,), dtype=_type)
            good_frames[0] = self.run.params.good_frames
            
        def load_instrument(nxs_file):
            
            inst_grp = nxs_file.create_group("raw_data_1/instrument")
            
            dae_grp = nxs_file.create_group("raw_data_1/instrument/dae")
            
            # Unable to find equivalent in .raw 
            name = "detector_table_file"
            dae_grp.create_dataset(name, (1,), dtype="|S28")
           
            name =  "period_index"
            _type = np.dtype("i4")
            period_idx = dae_grp.create_dataset(name, (1,), dtype=_type)
            period_idx[0] = self.time_channel.num_periods
            
            # Unable to find equivalent in .raw 
            name = "spectra_table_file"
            _type = np.dtype("S35")
            dae_grp.create_dataset(name, (1,), dtype=_type)
            
            tc1_grp = dae_grp.create_group("time_channels_1")
           
            name = "time_of_flight"
            tc1_grp[name] = nxs_file["raw_data_1/detector_1/time_of_flight"]

            name = "time_of_flight_raw"
            shape = (len(self.time_channel.boundaries[:]),)
            _type =  np.dtype("f4")
            tof_raw = tc1_grp.create_dataset(name, shape, dtype=_type)
            tof_raw[:] = np.array(self.time_channel.raw_boundaries[:]) * 31.25
            
            # Unable to find equivalent in .raw 
            dae_grp.create_dataset("type", (1,), dtype=np.dtype("f4"))
            
            vetos_grp = dae_grp.create_group("vetos")
            
            # Unable to find equivalent in .raw 
            name = "ISIS_50Hz"
            _type = np.dtype("i4")
            vetos_grp.create_dataset(name, (1,), dtype=_type)
            
            # Unable to find equivalent in .raw 
            name = "TS2_pulse"
            _type = np.dtype("i4")
            vetos_grp.create_dataset(name, (1,), dtype=_type)
            
            name = "ext0"
            _type = np.dtype("i4")
            ext0 = vetos_grp.create_dataset(name, (1,), dtype=_type)
            ext0[0] = self.dae.params.external_vetoes[0]
            
            name = "ext1"
            _type = np.dtype("i4")
            ext1 = vetos_grp.create_dataset(name, (1,), dtype=_type)
            ext1[0] = self.dae.params.external_vetoes[1]
            
            name = "ext2"
            _type = np.dtype("i4")
            ext2 = vetos_grp.create_dataset(name, (1,), dtype=_type)
            ext2[0] = self.dae.params.external_vetoes[2]
            
            # Unable to find equivalent in .raw 
            name = "ext3"
            _type = np.dtype("i4")
            vetos_grp.create_dataset(name, (1,), dtype=_type)
            
            name = "fermi_chopper0"
            _type = np.dtype("i4")
            f_chopper_0 = vetos_grp.create_dataset(name, (1,), dtype=_type)
            f_chopper_0[0] = self.instrument.params.chopper_freq_1
            
            name = "fermi_chopper1"
            _type = np.dtype("i4")
            f_chopper_1 = vetos_grp.create_dataset(name, (1,), dtype=_type)
            f_chopper_1[0] = self.instrument.params.chopper_freq_2
            
            name = "fermi_chopper2"
            _type = np.dtype("i4")
            f_chopper_2 = vetos_grp.create_dataset(name, (1,), dtype=_type)
            f_chopper_2[0] = self.instrument.params.chopper_freq_3
            
            # Unable to find equivalent in .raw 
            name = "fermi_chopper3"
            _type = np.dtype("i4")
            vetos_grp.create_dataset(name, (1,), dtype=_type)

            # Unable to find equivalent in .raw 
            name = "fifo"
            _type = np.dtype("i4")
            vetos_grp.create_dataset(name, (1,), dtype=_type)
            
            # Unable to find equivalent in .raw      
            name = "msmode"
            _type = np.dtype("i4")
            vetos_grp.create_dataset(name, (1,), dtype=_type)
            
            name = "smp"
            _type = np.dtype("i4")
            smp = vetos_grp.create_dataset(name, (1,), dtype=_type)
            smp[0] = self.dae.params.secondary_master_pulse
            
            # Unable to find equivalent in .raw 
            name = "writing_table_file"
            _type = np.dtype("S30")
            dae_grp.create_dataset(name, (1,), dtype=_type)
            
            d1_grp = inst_grp.create_group("detector_1")
            d1_grp["counts"] = nxs_file['raw_data_1']['detector_1']['counts']
            
            name = "delt"
            shape = (len(self.instrument.hold_off_table[:-4]),)
            _type = np.dtype("f4")
            d1_delt = d1_grp.create_dataset(name, shape, dtype=_type)
            d1_delt[:] = self.instrument.hold_off_table[:-4]
            
            name = "distance"
            shape = (len(self.instrument.L2_table[:-4]),)
            _type = np.dtype("f4")
            d1_L2 = d1_grp.create_dataset(name, shape,  dtype=_type)
            d1_L2[:] = self.instrument.L2_table[:-4]
            
            # From .nxs examples, looks to be a sum of the longitude angles
            # See https://www.nexusformat.org/TOFRaw.html
            name = "polar_angle"
            shape = (len(self.instrument.two_theta) - 4,)
            _type = np.dtype("f4")
            d1_polar_angle = d1_grp.create_dataset(name, shape, dtype=_type)
            d1_polar_angle[:] = self.instrument.two_theta[:-4]
            
            name = "source_detector_distance"
            _type = np.dtype("f4")
            d1_detector_dist = d1_grp.create_dataset(name, (1,), dtype=_type)
            d1_detector_dist[0] = self.instrument.params.LOQ_souce_detector_distance
            
            name = "spectrum_index"
            d1_grp[name] = nxs_file['raw_data_1']['detector_1']['spectrum_index'][:]
            
            name = "time_of_flight"
            nxs_path = "raw_data_1/detector_1/time_of_flight"
            d1_grp[name] = nxs_file[nxs_path]
            
            name = "time_of_flight_raw"
            nxs_path = "raw_data_1/instrument/dae/time_channels_1/time_of_flight_raw"
            d1_grp[name] = nxs_file[nxs_path]
            
            moderator_grp = inst_grp.create_group("moderator")
            
            # Unable to find equivalent in .raw
            name = "distance"
            _type = np.dtype("f4")
            m_distance = moderator_grp.create_dataset(name, (1,), dtype=_type)
            
            name = "name"
            _type = np.dtype("S3")
            m_name = inst_grp.create_dataset(name, (1,), dtype=_type)
            m_name[0] = self.summary.header.instrument_name
            
            # Unable to find equivalent in .raw
            source_grp = inst_grp.create_group("source")
            source_grp.create_dataset("name", (1,), dtype=np.dtype("S4"))
            source_grp.create_dataset("probe", (1,), dtype=np.dtype("S8"))
            source_grp.create_dataset("type", (1,), dtype=np.dtype("S21"))
            
        def load_isis_vms_compat(nxs_file):
                
            def extract_struct_to_arr(data_struct, data_arr, data_type=None):
                count=0
                for i in data_struct._fields_:
                    try:
                        if data_type is None:
                            data_arr[count] = getattr(data_struct, i[0])
                        else:
                            data_arr[count] = data_type(getattr(data_struct, i[0]))
                        count += 1
                    except TypeError:
                        val = getattr(data_struct, i[0])
                        for j in val: 
                            if data_type is None:
                                data_arr[count] = j
                            else:
                                data_arr[count] = data_type(j)
                            count += 1


            grp = nxs_file.create_group("raw_data_1/isis_vms_compat")
            
            name = "ADD"
            _type = np.dtype("i4")
            add = grp.create_dataset(name, (9,), dtype=_type)
            for c,i in enumerate(self.info_positions._fields_):
                add[c] = getattr(self.info_positions, i[0])
        
            name = "CODE"
            shape = (len(self.instrument.UTn_tables_code),)
            _type = np.dtype("i4")
            code = grp.create_dataset(name, shape, dtype=_type)
            code[:] = self.instrument.UTn_tables_code[:]
            
            name = "CRAT"
            shape = (len(self.dae.detector_crate_nums),) 
            _type = np.dtype("i4")
            crat = grp.create_dataset(name, shape, dtype=_type)
            crat[:] = self.dae.detector_crate_nums[:]
            
            # Unable to find equivalent in .raw
            name = "CRPB"
            shape = (32,1)
            _type = np.dtype("S4")
            grp.create_dataset(name, shape, dtype=_type)
            
            # Unable to find equivalent in .raw
            name = "CSPB"
            shape = (64,1)
            _type = np.dtype("S4")
            grp.create_dataset(name, shape, dtype=_type)
            
            name = "DAEP"
            _type = np.dtype("i4")
            daep = grp.create_dataset(name, (64,), dtype=_type)
            extract_struct_to_arr(self.dae.params, daep)

            name = "DELT"
            shape = (len(self.instrument.hold_off_table),)
            _type = np.dtype("f4")
            delt = grp.create_dataset(name, shape, dtype=_type)
            delt[:] = self.instrument.hold_off_table[:]
            
            name = "FORM"
            _type = np.dtype("i4")
            form = grp.create_dataset(name, (1,), dtype=_type)
            form[0] = self.summary.format[0]
           
            _type =  np.dtype("S80")
            #hdr = grp.create_dataset("HDR", (1,), dtype=_type)
            #extract_struct_to_arr(self.summary.header, hdr)
            
            name = "IRPB"
            _type = np.dtype("i4")
            irpb = grp.create_dataset(name, (32,), dtype=_type)
            extract_struct_to_arr(self.run.params, irpb)
            
            name = "ISPB"
            _type = np.dtype("i4")
            ispb = grp.create_dataset(name, (64,), dtype=_type)
            extract_struct_to_arr(self.sample.params, ispb)
            
            name = "IVPB"
            _type = np.dtype("i4")
            ivpb = grp.create_dataset(name, (64,), dtype=_type)
            extract_struct_to_arr(self.instrument.params, ivpb)
            
            name = "LEN2"
            shape = (len(self.instrument.L2_table),)
            _type = np.dtype("f4")
            len2 = grp.create_dataset(name, shape, dtype=_type)
            len2[:] = self.instrument.L2_table[:]

            name = "MDET"            
            shape = (len(self.instrument.monitor_detector_number),)
            _type = np.dtype("i4")
            mdet = grp.create_dataset(name, shape, dtype=_type)
            mdet[:] = self.instrument.monitor_detector_number[:]
            
            name = "MODN"
            shape = (len(self.dae.detector_module_nums),)
            _type = np.dtype("i4")
            modn = grp.create_dataset(name, shape,  dtype=_type)
            modn[:] = self.dae.detector_module_nums[:]
            
            name = "MONP"
            shape = len(self.instrument.monitor_prescale_val,)
            _type = np.dtype("i4")
            monp = grp.create_dataset(name, shape,  dtype=_type)
            monp[:] = self.instrument.monitor_prescale_val[:]
            
            name = "MPOS"
            shape = (len(self.dae.detector_module_positions),)
            _type = np.dtype("i4")
            mpos = grp.create_dataset(name, shape,  dtype=_type)
            mpos[:] = self.dae.detector_module_positions[:]
            
            name = "NAME"
            _type = np.dtype("S8")
            name = grp.create_dataset(name, (1,), dtype=_type)
            name[:] = self.instrument.name[:]
            
            name = "NDET"
            _type = np.dtype("i4")
            ndet = grp.create_dataset(name, (1,), dtype=_type)
            ndet[0] = self.instrument.num_detectors
            
            name = "NFPP"
            _type = np.dtype("i4")
            nfpp = grp.create_dataset(name, (1,), dtype=_type)
            nfpp[0] = self.time_channel.num_frames_per_period.value
            
            name = "NMON"
            _type = np.dtype("i4")
            nmon = grp.create_dataset(name, (1,), dtype=_type)
            nmon[0] = self.instrument.num_monitors.value

            name = "NOTE"            
            _type = np.dtype("S19")
            note = grp.create_dataset(name, (1,), dtype=_type)
            note[0] = b" No notes were made"
            
            name = "NPER"
            _type = np.dtype("i4")
            nper = grp.create_dataset(name, (1,), dtype=_type)
            nper[0] = self.time_channel.num_periods.value
            
            name = "NSP1"
            _type = np.dtype("i4")
            nsp1 = grp.create_dataset(name, (1,), dtype=_type)
            nsp1[0] = self.time_channel.num_spectra.value
            
            name = "NTC1"
            _type = np.dtype("i4")
            ntc1 = grp.create_dataset(name, (1,), dtype=_type)
            ntc1[0] = self.time_channel.num_time_channels
            
            # Unable to find equivalent in .raw
            name = "NTLL"
            _type = np.dtype("i4")
            grp.create_dataset(name, (1,), dtype=_type)
            
            name = "NTRG"
            _type = np.dtype("i4")
            ntrg = grp.create_dataset(name, (1,), dtype=_type)
            ntrg[0] = self.time_channel.num_time_regimes.value
            
            name = "NUSE"
            _type = np.dtype("i4")
            nuse = grp.create_dataset(name, (1,), dtype=_type)
            nuse[0] = self.instrument.num_user_tables.value
            
            name = "PMAP"
            shape = (1,len(self.time_channel.period_num_per_period))
            _type = dtype=np.dtype("i4")
            pmap = grp.create_dataset(name, shape, _type)
            pmap[0] = self.time_channel.period_num_per_period[:]
            
            name = "PRE1"
            _type = np.dtype("i4")
            pre1 = grp.create_dataset(name, (1,), dtype=_type)
            pre1[0] = self.time_channel.prescale.value
            
            name = "RRPB"
            _type = np.dtype("f4")
            rrpb = grp.create_dataset(name, (32,), dtype=_type)
            extract_struct_to_arr(self.run.params, rrpb)
            
            name = "RSPB"
            _type = np.dtype("f4")
            rspb = grp.create_dataset(name, (64,), dtype=_type)
            extract_struct_to_arr(self.sample.params, rspb)
            
            name = "RUN"
            _type = np.dtype("i4")
            run = grp.create_dataset(name, (1,), dtype=_type)
            run[0] = self.summary.header.run_number
            
            name = "RVPB"
            _type = np.dtype("f4")
            rvpb = grp.create_dataset(name, (64,), dtype=_type)
            extract_struct_to_arr(self.run.params, rvpb)
            
            name = "SPB"
            _type = np.dtype("i4")
            spb = grp.create_dataset(name, (64,), dtype=_type)
            extract_struct_to_arr(self.sample.params, spb)
            
            name = "SPEC"
            shape = (len(self.instrument.spectrum_number_table[:]),)
            _type = np.dtype("i4")
            spec = grp.create_dataset(name, shape, dtype=_type)
            spec[:] = self.instrument.spectrum_number_table[:]
            
            name = "TCM1"
            _type = np.dtype("i4")
            tcm1 = grp.create_dataset(name, (5,), dtype=_type)
            tcm1[:] = self.time_channel.mode[:]
            
            name = "TCP1"
            _type = np.dtype("f4")
            tcp1 = grp.create_dataset(name, (20,), dtype=_type)
            extract_struct_to_arr(self.time_channel.params, tcp1)
            
            name = "TIMR"
            shape = (len(self.dae.detector_time_regimes[:]),)
            _type = np.dtype("f4")
            timr = grp.create_dataset(name, shape, dtype=_type)
            timr[:] = self.dae.detector_time_regimes[:]
            
            name = "TITL"
            _type = np.dtype("S80")
            titl = grp.create_dataset(name, (1,), dtype=_type)
            titl[:] = self.run.title[:]
            
            name = "TTHE"
            shape = (len(self.instrument.two_theta),)
            _type = np.dtype("f4")
            tthe = grp.create_dataset(name, shape, dtype=_type)
            tthe[:] = self.instrument.two_theta[:]
            
            name = "UDET"
            shape = (len(self.dae.detector_user_nums),)
            _type = np.dtype("f4")
            udet = grp.create_dataset(name, shape, dtype=_type)
            udet[:] = self.dae.detector_user_nums[:]
            
            name = "ULEN"
            _type = np.dtype("i4")
            ulen= grp.create_dataset(name, (1,), dtype=_type)
            ulen[0] = self.user.length
            
            name = "USER"
            _type = np.dtype("S160")
            user = grp.create_dataset(name, (1,), dtype=_type)
            user[0] = self.run.user
            
            name = "VER1"
            _type = np.dtype("i4")
            ver1 = grp.create_dataset(name, (1,), dtype=_type)
            ver1[0] = self.summary.format_version.value
            
            name = "VER2"
            _type = np.dtype("i4")
            ver2 = grp.create_dataset(name, (1,), dtype=_type)
            ver2[0] = self.run.version.value
            
            name = "VER3"
            _type = np.dtype("i4")
            ver3 = grp.create_dataset(name, (1,), dtype=_type)
            ver3[0] = self.instrument.version.value

            name = "VER4"
            _type = np.dtype("i4")
            ver4 = grp.create_dataset(name, (1,), dtype=_type)
            ver4[0] = self.sample.env_version.value
            
            name = "VER5"
            _type = np.dtype("i4")
            ver5 = grp.create_dataset(name, (1,), dtype=_type)
            ver5[0] = self.dae.version.value
            
            name = "VER6"
            _type = np.dtype("i4")
            ver6 = grp.create_dataset(name, (1,), dtype=_type)
            ver6[0] = self.summary.format_version.value
            
            name = "VER7"
            _type = np.dtype("i4")
            ver7 = grp.create_dataset(name, (1,), dtype=_type)
            ver7[0] = self.time_channel.version.value
            
            name = "VER8"
            _type = np.dtype("i4")
            ver8 = grp.create_dataset(name, (1,), dtype=_type)
            ver8[0] = self.data.version.value
            
            name = "VER9"
            _type = np.dtype("i4")
            ver9 = grp.create_dataset(name, (1,), dtype=_type)
            ver9[0] = self.log.version[0].value
            
        def load_measurement(nxs_file):
            # Unable to find equivalent in .raw
            grp = nxs_file.create_group("raw_data_1/measurement")
            grp.create_dataset("first_run", (1,), dtype=np.dtype("i4"))
            grp.create_dataset("id", (1,), dtype=np.dtype("S1"))
            grp.create_dataset("label", (1,), dtype=np.dtype("S1"))
            grp.create_dataset("subid", (1,), dtype=np.dtype("S1"))
            grp.create_dataset("type", (1,), dtype=np.dtype("S1"))

        def load_measurement_first_run(nxs_file):
            path = "raw_data_1/measurement_first_run"
            link_path = "raw_data_1/measurement/first_run"
            nxs_file[path] = nxs_file[link_path]
    
        def load_measurement_id(nxs_file):
            path = "raw_data_1/measurement_id"
            link_path = "raw_data_1/measurement/id"
            nxs_file[path] = nxs_file[link_path]
    
        def load_measurement_label(nxs_file):
            path = "raw_data_1/measurement_label"
            link_path = "raw_data_1/measurement/label"
            nxs_file[path] = nxs_file[link_path]
    
        def load_measurement_subid(nxs_file):
            path = "raw_data_1/measurement_subid"
            link_path = "raw_data_1/measurement/subid"
            nxs_file[path] = nxs_file[link_path]

        def load_measurement_type(nxs_file):
            path = "raw_data_1/measurement_type"
            link_path = "raw_data_1/measurement/label"
            nxs_file[path] = nxs_file[link_path]

        def load_monitor(nxs_file, monitor_num, raw_data):
            
            grp = nxs_file.create_group("raw_data_1/monitor_" + str(monitor_num))
            monitor_idx = self.instrument.monitor_detector_number[monitor_num]

            name = "data"            
            shape = (1,1,self.time_channel.num_time_channels+1)
            _type = np.dtype("i4")
            data = grp.create_dataset(name, shape, dtype=_type)
            data[0,0,:] = raw_data[0][monitor_idx][:]
            
            # Unable to find equivalent in .raw
            name = "period_index"
            _type = np.dtype("i4")
            grp.create_dataset(name, (1,), dtype=_type)
            
            name = "spectrum_index"
            _type = np.dtype("i4")
            spectrum_idx = grp.create_dataset(name, (1,), dtype=_type)
            spectrum_idx[0] = monitor_idx 

            name = "time_of_flight"
            shape = (self.time_channel.num_time_channels+1,)
            _type = np.dtype("f4")
            tof = grp.create_dataset(name, shape, dtype=_type)
            tof[:] = self.time_channel.boundaries[:]
            
        def load_unsaved_monitor_events(nxs_file):
            # Unable to find equivalent in .raw
            name = "raw_data_1/monitor_events_not_saved"
            _type = np.dtype("i8")
            nxs_file.create_dataset(name, (1,), dtype=_type)
            
        def load_name(nxs_file):
            _name = "raw_data_1/name"
            _type = np.dtype("|S3")
            name = nxs_file.create_dataset(_name, (1,), dtype=_type)
            name[0] = self.summary.header.instrument_name

        def load_notes(nxs_file):
            name = "raw_data_1/notes"
            _type = np.dtype("S1")
            # Unable to find equivalent in .raw
            nxs_file.create_dataset(name, (1,), dtype=_type)
            
        def load_periods(nxs_file):
            grp = nxs_file.create_group("raw_data_1/periods")
            
            name = "frames_requested"
            _type = np.dtype("i4")
            frames_req = grp.create_dataset(name, (1,), dtype=_type)
            frames_req[0] = self.run.params.requested_duration
            
            name = "good_frames"
            _type = np.dtype("i4")
            good_frames = grp.create_dataset(name, (1,), dtype=_type)
            good_frames[0] = self.run.params.good_frames
            
            # Unable to find equivalent in .raw
            name = "good_frames_daq"
            _type = np.dtype("i4")
            grp.create_dataset(name, (1,), dtype=_type)
        
            # Unable to find equivalent in .raw
            name = "highest_used"
            _type = np.dtype("i4")
            grp.create_dataset(name, (1,), dtype=_type)
            
            # Unable to find equivalent in .raw
            name = "labels"
            _type =  np.dtype("S8")
            grp.create_dataset(name, (1,), dtype=_type)
            
            name = "number"
            _type = np.dtype("i4")
            number = grp.create_dataset(name, (1,), dtype=_type)
            number[0] = self.time_channel.num_periods.value
            
            # Unable to find equivalent in .raw
            name = "output"
            _type = np.dtype("i4")
            grp.create_dataset(name, (1,), dtype=_type)  
            
            name = "proton_charge"
            _type = np.dtype("f4")
            proton_charge = grp.create_dataset(name, (1,), dtype=_type)
            proton_charge[0] = self.run.params.good_proton_charge
            
            name = "proton_charge_raw"
            _type = np.dtype("f4")
            proton_charge_raw = grp.create_dataset(name, (1,), dtype=_type)
            proton_charge_raw[0] = self.run.params.total_proton_charge
            
            name = "raw_frames"
            _type = np.dtype("i4")
            raw_frames = grp.create_dataset(name, (1,), dtype=_type)
            raw_frames[0] = self.run.params.raw_frames
            
            # Unable to find equivalent in .raw
            name = "sequences"
            _type = np.dtype("i4")
            grp.create_dataset(name, (1,), dtype=_type) 
            
            # Unable to find equivalent in .raw
            name = "total_counts"
            _type = np.dtype("f4")
            grp.create_dataset(name, (1,), dtype=_type) 
            
            name = "type"
            _type = np.dtype("i4")
            mod_type = grp.create_dataset(name, (1,), dtype=_type)
            mod_type[0] = self.instrument.params.moderator_type
            
        def load_program_name(nxs_file):
            # Unable to find equivalent in .raw
            name = "raw_data_1/program_name"
            _type = np.dtype("S11")
            nxs_file.create_dataset(name, (1,), dtype=_type)
            
        def load_proton_charge(nxs_file):
            path = "raw_data_1/proton_charge"
            link_path = "raw_data_1/periods/proton_charge"
            nxs_file[path] = nxs_file[link_path]
                                                    
        def load_proton_charge_raw(nxs_file):
            path = "raw_data_1/proton_charge_raw"
            link_path = "raw_data_1/periods/proton_charge_raw"
            nxs_file[path] = nxs_file[link_path]
            
        def load_raw_frames(nxs_file):
            path = "raw_data_1/raw_frames"
            link_path = "raw_data_1/periods/raw_frames"
            nxs_file[path] = nxs_file[link_path]
                                                
        def load_run_cycle(nxs_file):
            # Unable to find equivalent in .raw
            name = "raw_data_1/run_cycle"
            _type = np.dtype("S4")
            nxs_file.create_dataset(name, (1,), dtype=_type)
            
        def load_run_number(nxs_file):
            name = "raw_data_1/run_number"
            _type = np.dtype("i4")
            run_number = nxs_file.create_dataset(name, (1,), dtype=_type)
            run_number[0] = self.run.number.value
            
        def load_runlog(nxs_file):
            # Unable to find equivalent in .raw
            nxs_file.create_group("raw_data_1/runlog")
            
        def load_sample(nxs_file):
            grp = nxs_file.create_group("raw_data_1/sample")
            
            # Unable to find equivalent in .raw
            name = "distance"
            _type = np.dtype("f8")
            grp.create_dataset(name, (1,), dtype=_type) 
            
            name = "height"
            _type = np.dtype("f4")
            height = grp.create_dataset(name, (1,), dtype=_type)
            height[0] = self.sample.params.sample_height
            
            # Unable to find equivalent in .raw
            grp.create_dataset("id", (1,), dtype=np.dtype("S1"))
            
            # Unable to find equivalent in .raw
            grp.create_dataset("name", (1,), dtype=np.dtype("S1"))
            
            # Unable to find equivalent in .raw
            grp.create_dataset("shape", (1,), dtype=np.dtype("S1"))    
            
            name = "thickness"
            _type = np.dtype("f4")
            thickness = grp.create_dataset(name, (1,), dtype=_type)
            thickness[0] = self.sample.params.sample_thickness
            
            name = "type"
            _type = np.dtype("S1")
            sample_type = grp.create_dataset(name, (1,), dtype=_type)
            sample_type[0] = str(self.sample.params.sample_type).encode()
            
            name = "width"
            _type = np.dtype("f4")
            width = grp.create_dataset(name, (1,), dtype=_type)
            width[0] = self.sample.params.sample_width
            
        def load_script_name(nxs_file):
            # Unable to find equivalent in .raw
            name = "raw_data_1/script_name"
            _type = np.dtype("S1")
            nxs_file.create_dataset(name, (1,), dtype=_type)
            
        def load_seci_config(nxs_file):
            # Unable to find equivalent in .raw
            name = "raw_data_1/seci_config"
            _type = np.dtype("S25")
            nxs_file.create_dataset(name, (1,), dtype=_type)
            
        def load_selog(nxs_file):
            # Unable to find equivalent in .raw
            nxs_file.create_group("raw_data_1/selog")
            
        def load_start_time(nxs_file):
            name = "raw_data_1/start_time"
            _type = np.dtype("S19")
            start_time = nxs_file.create_dataset(name, (1,), dtype=_type)
            raw_date = convert_to_nxs_date(self.summary.header.start_date)
            raw_time = self.summary.header.start_time.decode()
            start_time[0] = (raw_date + "T" + raw_time).encode()
            
        def load_title(nxs_file):
            name = "raw_data_1/title"
            _type = np.dtype("S80")
            title = nxs_file.create_dataset(name, (1,), dtype=_type)
            title[:] = self.run.title[:]
            
                
        def load_total_counts(nxs_file):
            # Unable to find equivalent in .raw
            name = "raw_data_1/total_counts"
            _type = np.dtype("f4")
            nxs_file.create_dataset(name, (1,), dtype=_type)
            
        def load_total_uncounted_counts(nxs_file):
            # Unable to find equivalent in .raw
            name = "raw_data_1/total_uncounted_counts"
            _type = np.dtype("i4")
            nxs_file.create_dataset(name, (1,), dtype=_type)

        def load_user_1(nxs_file):
            grp = nxs_file.create_group("raw_data_1/user_1")
            name = "affiliation"
            _type = np.dtype("S13")
            affiliation = grp.create_dataset(name, (1,), dtype=_type)
            affiliation[0] = self.run.user.institute

            name = grp.create_dataset("name", (1,), dtype=np.dtype("S8"))
            name[0] = self.run.user.user
                    
        if path.isfile(output_filename):
            raise IOError("{} already exists.".format(output_filename))
        
        nxs_file = h5py.File(output_filename, "a")
        num_detectors = self.instrument.num_detectors
        num_time_channels = self.time_channel.num_time_channels
        
        load_beamline(nxs_file)
        load_collection_time(nxs_file)
        load_definition(nxs_file)
        load_definition_local(nxs_file)
        load_detector_1(nxs_file, num_detectors + 1, num_time_channels + 1)
        load_duration(nxs_file)
        load_end_time(nxs_file)
        load_experiment_identifier(nxs_file)
        load_frame_log(nxs_file)
        load_good_frames(nxs_file)
        load_instrument(nxs_file)
        #load_isis_vms_compat(nxs_file)
        load_measurement(nxs_file)
        load_measurement_first_run(nxs_file)
        load_measurement_id(nxs_file)
        load_measurement_label(nxs_file)
        load_measurement_subid(nxs_file)
        load_measurement_type(nxs_file)
        
        raw_data = self.get_counts_per_pixel(num_detectors + 1, num_time_channels + 1)
        load_monitor(nxs_file, 1, raw_data)
        load_monitor(nxs_file, 2, raw_data)
        load_monitor(nxs_file, 3, raw_data)
        load_unsaved_monitor_events(nxs_file)
        load_name(nxs_file)
        load_notes(nxs_file)
        load_periods(nxs_file)
        load_program_name(nxs_file)
        load_proton_charge(nxs_file)
        load_proton_charge_raw(nxs_file)
        load_raw_frames(nxs_file)
        load_run_cycle(nxs_file)
        load_run_number(nxs_file)
        load_runlog(nxs_file)
        load_sample(nxs_file)
        load_script_name(nxs_file)
        load_seci_config(nxs_file)
        load_selog(nxs_file)
        load_start_time(nxs_file)
        load_title(nxs_file)
        load_total_counts(nxs_file)
        load_total_uncounted_counts(nxs_file)
        load_user_1(nxs_file)
        
        
    

    
            
