from __future__ import annotations

from abc import ABCMeta, abstractmethod
from typing import Tuple

from panel import Panel
from peak_table import PeakTable


class ExperimentReader(metaclass=ABCMeta):

    """
    Base class to define the interface for for reading and writing
    to different input files.
    """

    @abstractmethod
    def __init__(self, file_path: str) -> None:
        self.file_path = file_path

    @abstractmethod
    def _open(self) -> None:
        pass

    @abstractmethod
    def _close(self) -> None:
        pass

    @abstractmethod
    def get_panels(self, expt_idx: int = 0) -> Tuple[Panel, ...]:
        pass

    @abstractmethod
    def replace_panels(
        self,
        new_panels: Tuple[Panel, ...],
        expt_idx: int = 0,
    ) -> None:
        pass

    @abstractmethod
    def get_peak_table(self, expt_idx: int = 0) -> PeakTable:
        pass

    @abstractmethod
    def replace_peak_table(self, new_peak_table: PeakTable, expt_idx: int = 0) -> None:
        pass

    @abstractmethod
    def convert_panels(self, panels: Tuple[Panel, ...]) -> Tuple[Panel, ...]:
        pass

    def replace_from_reader(self, reader: ExperimentReader, expt_idx: int = 0) -> None:

        """
        Updates file at self.file_path with values from reader
        """

        reader_panels = reader.get_panels(expt_idx=expt_idx)
        converted_panels = self.convert_panels(reader_panels)
        self.replace_panels(new_panels=converted_panels, expt_idx=expt_idx)

        reader_peak_table = reader.get_peak_table(expt_idx=expt_idx)
        self.replace_peak_table(new_peak_table=reader_peak_table)
