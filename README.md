# isis_utils

The DIALS project provides an extensible framework to analyse X-ray diffraction data.
Much of this work is agnostic to the method used to obtain diffraction patterns, and can equally be applied to neutron diffraction data.
This repository provides a bridge between data obtained from the ISIS Neutron and Muon source and DIALS, converting raw output to formats that
DIALS can read and, additional post-processing for integration into other packages (e.g. Mantid).

### Known Issues

**Spectrum table inconsistencies between TOFRAW and RAW files**

The .raw spectrum_number_table is inconsistent with the .nxs detector_1/spectrum_index.
The former gives a 1D array from the top right of each detector moving down each column (i.e for 64x64 detector [64,63,..]).
The .nxs index just gives [i for i in range(num pixels)].
When visualising both give the same laue plots, consistent with SXD2001 only when using the indexing e.g. for a 64x64 detector:
np.arange(min_pixel_idx_range, max_pixel_idx_range).reshape(64,64).T

