from dials_reader import DIALSReader
from experiment_reader import ExperimentReader
from mantid_reader import MantidReader


def experiment_reader_factory(
    file_path: str, file_path2: str = None
) -> ExperimentReader:
    def is_dials_file(file_path: str, file_path2: str = None) -> bool:
        if file_path.split(".")[1] == "expt" and file_path2 is None:
            return True
        elif file_path2 is not None:
            if file_path.split(".")[1] == "expt" and file_path2.split(".")[1] == "refl":
                return True
        return False

    def is_mantid_file(file_path: str, file_path2: str = None) -> bool:
        if file_path2 is not None:
            return False
        if file_path.split(".")[1] == "nxs":
            return True
        return False

    if is_dials_file(file_path=file_path, file_path2=file_path2):
        return DIALSReader(expt_file_path=file_path, refl_file_path=file_path2)

    elif is_mantid_file(file_path=file_path, file_path2=file_path2):
        return MantidReader(nxs_file_path=file_path)

    raise Exception("Cannot understand input files.")
