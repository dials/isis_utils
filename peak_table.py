import numpy as np


class PeakTable:
    """
    Class to store peak properties that different readers
    can understand.
    """

    def __init__(
        self,
        spectra_idx_1D: np.array,
        intensity: np.array,
        energy: np.array,  # keV
        wavelength: np.array,  # A
        d_spacing: np.array,
        tof: np.array,  # usec
        miller_indices: np.array = None,
    ) -> None:
        self._size = None
        self._spectra_idx_1D = None
        self._intensity = None
        self._energy = None
        self._wavelength = None
        self._d_spacing = None
        self._tof = None
        self._miller_indices = None
        self.spectra_idx_1D = spectra_idx_1D
        self.intensity = intensity
        self.energy = energy
        self.wavelength = wavelength
        self.d_spacing = d_spacing
        self.tof = tof
        self.miller_indices = miller_indices

        self._dict = {
            "spectra_idx_1D": self.spectra_idx_1D,
            "intensity": self.intensity,
            "energy": self.energy,
            "wavelength": self.wavelength,
            "d_spacing": self.d_spacing,
            "tof": self.tof,
            "miller_indices": self.miller_indices,
        }

    def __len__(self):
        return self._size

    def __getitem__(self, index):
        return self._dict[index]

    def __setitem__(self, index, value):
        self._dict[index] = value

    @property
    def spectra_idx_1D(self):
        return self._spectra_idx_1D

    @spectra_idx_1D.setter
    def spectra_idx_1D(self, value):
        if self._size is None:
            self._size = len(value)
        elif len(value) != self._size:
            raise ValueError("Size is not consistent with current values")
        self._spectra_idx_1D = value

    @property
    def intensity(self):
        return self._intensity

    @intensity.setter
    def intensity(self, value):
        if self._size is None:
            self._size = len(value)
        elif len(value) != self._size:
            raise ValueError("Size is not consistent with current values")
        self._intensity = value

    @property
    def energy(self):
        return self._energy

    @energy.setter
    def energy(self, value):
        if self._size is None:
            self._size = len(value)
        elif len(value) != self._size:
            raise ValueError("Size is not consistent with current values")
        self._energy = value

    @property
    def wavelength(self):
        return self._wavelength

    @wavelength.setter
    def wavelength(self, value):
        if self._size is None:
            self._size = len(value)
        elif len(value) != self._size:
            raise ValueError("Size is not consistent with current values")
        self._wavelength = value

    @property
    def d_spacing(self):
        return self._d_spacing

    @d_spacing.setter
    def d_spacing(self, value):
        if self._size is None:
            self._size = len(value)
        elif len(value) != self._size:
            raise ValueError("Size is not consistent with current values")
        self._d_spacing = value

    @property
    def tof(self):
        return self._tof

    @tof.setter
    def tof(self, value):
        if self._size is None:
            self._size = len(value)
        elif len(value) != self._size:
            raise ValueError("Size is not consistent with current values")
        self._tof = value

    @property
    def miller_indices(self):
        return self._miller_indices

    @miller_indices.setter
    def miller_indices(self, value):
        if self._size is None:
            self._size = len(value)
        elif value is None:
            self._miller_indices = np.zeros(3 * self._size).reshape(self._size, 3)
        elif len(value) != self._size:
            raise ValueError("Size is not consistent with current values")
        else:
            self._miller_indices = value
