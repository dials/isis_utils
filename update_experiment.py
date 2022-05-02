import logging
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from os.path import splitext
from shutil import copyfile

logging.basicConfig(stream=sys.stdout, level=logging.INFO)

from experiment_reader_factory import experiment_reader_factory

logger = logging.getLogger("Update Experiment")


class ArgParser(ArgumentParser):
    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)


def get_output_path(file_path: str) -> str:
    name, extension = splitext(file_path)
    return name + "_updated" + extension


def main():

    parser = ArgParser(
        description="Updating experiments using experiment files\
        from other packages\n\n"
        "Example usage: \n"
        "python update_experiment.py mantid_example.nxs\
            --take_data_from dials_example.expt\n"
        "output: mantid_example.nxs with fields changed from dials_example.expt",
        formatter_class=RawTextHelpFormatter,
    )

    parser.add_argument(
        "source_files", nargs="+", help="source files to update", type=str
    )

    parser.add_argument(
        "-take_data",
        "--take_data_from",
        nargs="+",
        help="outputs copies of all raw_files as tofraw .nxs files",
        type=str,
    )

    parser.add_argument(
        "-in_place",
        "--edit_in_place",
        action="store_true",
        help="Modifies source_files directly rather than creating a copy",
        default=False,
    )

    args = parser.parse_args()
    if args.take_data_from:
        if len(args.take_data_from) == 2:
            file1, file2 = args.take_data_from
            reader_to_extract_from = experiment_reader_factory(
                file_path=file1, file_path2=file2
            )
        else:
            reader_to_extract_from = experiment_reader_factory(
                file_path=args.take_data_from[0]
            )
        for file_path in args.source_files:
            if args.edit_in_place:
                source_reader = experiment_reader_factory(file_path=file_path)
            else:
                copy_file_path = get_output_path(file_path=file_path)
                copyfile(file_path, copy_file_path)
                logger.info(f"{file_path} copied to {copy_file_path}")
                source_reader = experiment_reader_factory(file_path=copy_file_path)
            source_reader.replace_from_reader(reader=reader_to_extract_from)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
