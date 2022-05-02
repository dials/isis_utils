# isis_utils
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![basic_tests](https://github.com/dials/isis_utils/workflows/basic%20tests/badge.svg)
[![Python 3.8](https://img.shields.io/badge/python-3.8-blue.svg)](https://www.python.org/downloads/release/python-380/)

The DIALS project provides an extensible framework to analyse X-ray diffraction data.
Much of this work is agnostic to the method used to obtain diffraction patterns, and can equally be applied to neutron diffraction data.
This repository provides a bridge between data obtained from the ISIS Neutron and Muon source and DIALS, converting raw output to formats that
DIALS can read and, additional post-processing for integration into other packages (e.g. Mantid).

### Installation

```
git clone git@github.com:dials/isis_utils.git
cd isis_utils
pip install -r requirements.txt
```

### Convert ISIS RAW files to TOFRAW

[TOFRAW format](https://www.nexusformat.org/TOFRaw.html) | [ISIS formats](https://www.isis.stfc.ac.uk/Pages/ISIS-Raw-File-Format.aspx) | [Example data](https://doi.org/10.5281/zenodo.4415768)

Raw files can also be read using [OpenGenie](http://www.opengenie.org/Main_Page) and [Mantid](https://www.mantidproject.org/).


**Convert single .raw file**
```
python isis_raw_reader.py my_raw_file.raw --convert_to_tofraw
Reading my_raw_file.raw into memory..
Outputting data to my_raw_file.nxs..
```

**Convert multiple .raw files**
```
python isis_raw_reader.py my_raw_file_1.raw my_raw_file_2.raw --convert_to_tofraw
Reading my_raw_file_1.raw into memory..
Outputting data to my_raw_file_1.nxs..
Reading my_raw_file_2.raw into memory..
Outputting data to my_raw_file_2.nxs..
```
**Convert directory of .raw files**
```
python isis_raw_reader.py my_dir/*.raw --convert_to_tofraw
Reading my_dir/my_raw_file_1.raw into memory..
Outputting data to my_dir/my_raw_file_1.nxs..
etc.
```

**Known Issues**

Although the data in a generated TOFRAW file aims to be as close as possible to the RAW file, there are some known issues with the conversion used here.
This is due to changes/additions/removal of parameters between the two formats, inconsistencies between the two formats, and unknown parameters (the main source of truth for parameter names used when writing this code was [Mantid](https://github.com/mantidproject/mantid/tree/master/Framework/DataHandling/src/LoadRaw)).

*The following TOFRAW params are not found (by the author!) in RAW files*

raw_data_1/instrument/dae/vetos/fermi_chopper3 \
raw_data_1/instrument/dae/vetos/fifo \
raw_data_1/instrument/dae/vetos/msmode \
raw_data_1/instrument/dae/vetos/writing_table_file \
raw_data_1/instrument/detector_table_file \
raw_data_1/instrument/spectra_table_file \
raw_data_1/instrument/dae/type \
raw_data_1/instrument/dae/vetos/ISIS_50Hz \
raw_data_1/instrument/dae/vetos/TS2_pulse \
raw_data_1/instrument/dae/vetos/ext2 \
raw_data_1/instrument/dae/vetos/ext3 \
raw_data_1/instrument/moderator/distance \
raw_data_1/instrument/source/name \
raw_data_1/instrument/source/probe \
raw_data_1/instrument/source/type \
raw_data_1/measurement/first_run \
raw_data_1/measurement/id \
raw_data_1/measurement/label \
raw_data_1/measurement/subid \
raw_data_1/measurement/type \
raw_data_1/monitor_events_not_saved \
raw_data_1/notes \
raw_data_1/periods/good_frames_daq \
raw_data_1/periods/highest_used \
raw_data_1/periods/labels \
raw_data_1/periods/output \
raw_data_1/periods/sequences \
raw_data_1/periods/total_counts \
raw_data_1/run_cycle \
raw_data_1/run_log/* \
raw_data_1/sample/distance \
raw_data_1/sample/id \
raw_data_1/sample/name \
raw_data_1/sample/shape \
raw_data_1/script_name \
raw_data_1/seci_config \
raw_data_1/selog \
raw_data_1/total_counts \
raw_data_1/total_uncounted_counts

*The following TOFRAW params are not consistent (as far as I can tell!) between formats*

detector_1/spectrum_index \
raw_data_1/framelog \
raw_data_1/instrument/detector_1/distance \
raw_data_1/instrument/detector_1/polar_angle \
raw_data_1/monitor_*/* \
raw_data_1/periods/proton_charge \
raw_data_1/periods/proton_charge_raw \
raw_data_1/periods/type \
raw_data_1/run_number \
raw_data_1/sample/type

E.g. the .raw spectrum_number_table is inconsistent with the .nxs detector_1/spectrum_index.
The former gives a 1D array from the top right of each detector moving down each column (i.e for 64x64 detector [64,63,..]).
The .nxs index just gives [i for i in range(num pixels)].
When visualising both give the same laue plots, consistent with SXD2001 only when using the indexing e.g. for a 64x64 detector:
```
np.arange(min_pixel_idx_range, max_pixel_idx_range).reshape(64,64).T
```

### Update SXD Mantid Workspace with DIALS Output

```
python update_experiment.py mantid_workspace.nxs -take_data dials_experiment.expt
INFO: Update Experiment:mantid_workspace.nxs copied to mantid_workspace_updated.nxs
INFO: Update Experiment:Replacing detector panels in mantid_workspace_updated.nxs with those in dials_experiment.expt
```

### Update SXD Mantid Peaks Workspace with DIALS Output
```
python update_experiment.py mantid_peaks_workspace.nxs -take_data dials_experiment.expt dials_reflection_table.refl
INFO: Update Experiment:mantid_peaks_workspace.nxs copied to mantid_peaks_workspace_updated.nxs
INFO: Update Experiment:Replacing detector panels in mantid_peaks_workspace_updated.nxs with those in dials_experiment.expt
INFO: Update Experiment:Replacing peak table in mantid_peaks_workspace_updated.nxs with table from dials_reflection_table.refl
```
