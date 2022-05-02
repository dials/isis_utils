"""
taken from https://gist.github.com/ndevenish/079eb17365b21d2044ac146f67dbb499

DIALS .refl file loader

This loads msgpack-type DIALS reflection files, without having DIALS or
cctbx in the python environment.

Note: All modern .refl files are at time of writing msgpack-based. Some
much older files might be in pickle format, which this doesn't read.

Usage:

>>> import refl_loader, pathlib
>>> refl_file = Path("/path/to/a/reftable.refl")
>>> print(refl_loader.load(refl_file))

or

>>> print(refl_loader.loads(refl_file.read_bytes()))
"""

import functools
import operator
import os
import struct
import sys

# from collections.abc import Iterable
from io import BytesIO
from pathlib import Path
from typing import IO, Dict, Iterable, List, NamedTuple, Optional, Tuple, Union, cast

import msgpack
import numpy as np


class Shoebox(NamedTuple):
    panel: int
    bbox: Tuple[int]
    data: np.array = None
    mask: np.array = None
    background: np.array = None


def _decode_raw_numpy(dtype, shape: Union[int, Iterable[int]] = 1):
    """
    Decoding a column that maps straight to a numpy array.

    Args:
        dtype: The numpy dtype for the array
        shape:
            The shape of a single item. Either an int, or a collection
            of ints, in C-array order (row major)
    """
    # Convert to a shape tuple
    if isinstance(shape, int):
        shape = (shape,)
    else:
        shape = tuple(shape)

    def _decode_specific(data, copy):
        num_items, raw = data
        array = np.frombuffer(raw, dtype=dtype)

        if shape != (1,):
            item_width = functools.reduce(operator.mul, shape)
            assert len(raw) % item_width == 0
            assert num_items * item_width == len(array)
            array = array.reshape(num_items, *shape)
        if copy:
            return np.copy(array)
        return array

    return _decode_specific


def _decode_shoeboxes(data: List, copy) -> List[Optional[Shoebox]]:
    # Shoebox is float
    num_items, raw = data
    shoeboxes: List[Optional[Shoebox]] = []
    pos = 0
    while pos < len(raw):
        sbox_header_fmt = "<IiiiiiiB"
        sb_info = struct.unpack_from(sbox_header_fmt, raw, pos)
        pos += struct.calcsize(sbox_header_fmt)
        panel = sb_info[0]
        bbox = sb_info[1:7]
        data_present = sb_info[7]
        shoebox = {"panel": panel, "bbox": bbox}
        if data_present:
            bbox_size = (bbox[5] - bbox[4], bbox[3] - bbox[2], bbox[1] - bbox[0])
            data_size = (bbox_size[0] * bbox_size[1] * bbox_size[2]) * 4
            # Read three sets of data: data, mask and background
            shoebox["data"] = np.frombuffer(
                raw[pos : pos + data_size], dtype=np.float32
            ).reshape(bbox_size)
            pos += data_size
            shoebox["mask"] = np.frombuffer(
                raw[pos : pos + data_size], dtype=np.int32
            ).reshape(bbox_size)
            pos += data_size
            shoebox["background"] = np.frombuffer(
                raw[pos : pos + data_size], dtype=np.float32
            ).reshape(bbox_size)
            pos += data_size
            if copy:
                shoebox["data"] = np.copy(shoebox["data"])
                shoebox["mask"] = np.copy(shoebox["mask"])
                shoebox["background"] = np.copy(shoebox["background"])

        # Although this is technically a divergence, return None instead of an empty shoebox
        if not data_present and all(x == 0 for x in bbox) and panel == 0:
            shoeboxes.append(None)
        else:
            shoeboxes.append(Shoebox(**shoebox))
    assert num_items == len(shoeboxes)
    return np.array(shoeboxes, dtype=np.object_)


# Mapping from type name to decoder function
_reftable_decoders = {
    "bool": _decode_raw_numpy(bool),
    "int": _decode_raw_numpy(np.int32),
    "double": _decode_raw_numpy(np.double),
    "int6": _decode_raw_numpy(np.int32, shape=6),
    "std::size_t": _decode_raw_numpy(np.uint64),
    "vec3<double>": _decode_raw_numpy(np.double, shape=3),
    "cctbx::miller::index<>": _decode_raw_numpy(np.int32, shape=3),
    "Shoebox<>": _decode_shoeboxes,
    "vec2<double>": _decode_raw_numpy(np.double, shape=2),
    "mat3<double>": _decode_raw_numpy(np.double, shape=(3, 3)),
    # "std::string": _decode_wip, # - string writing broken; dials/dials#1858
}


def decode_column(column_entry, copy):
    """Decode a single column value"""
    datatype, data = column_entry

    converter = _reftable_decoders.get(datatype)
    if not converter:
        print(f"Warning: Data type '{datatype}' does not have a converter; cannot read")
        return None
    return converter(data, copy=copy)


def _get_unpacked(stream_or_path: Union[str, IO, bytes, os.PathLike]):
    """Works out the logic to pass a stream/pathlike to msgpack"""
    try:
        path = os.fspath(cast(str, stream_or_path))
        is_fspathlike = True
    except (TypeError, ValueError):
        is_fspathlike = isinstance(stream_or_path, str)

    if is_fspathlike:
        with open(path, "rb") as f:
            un = msgpack.Unpacker(f, strict_map_key=False)
            return un.unpack()
    else:
        un = msgpack.Unpacker(stream_or_path, strict_map_key=False)
        return un.unpack()


def loads(data: bytes, copy=False):
    """
    Load a DIALS msgpack-encoded .refl file.

    Args:
        data: bytes data, already read from the file.
        copy: Should the data be copied into writable numpy arrays.

    Returns: See .load(stream_or_path)
    """
    return load(BytesIO(data), copy)


def load(stream_or_path: Union[IO, os.PathLike], copy=False) -> Dict:
    """
    Load a DIALS msgpack-encoded .refl file

    Args:
        stream_or_path: The filename or data to load
        copy:
            Should the data be copied. This will cause more memory usage
            whilst loading the raw data.

    Returns:

        A dictionary with each column in the reflection table. If there
        is an identifier mapping as part of the reflection table, then
        this is returned as an extra 'experiment_identifier' column.
        All columns except Shoeboxes are returned as numpy arrays,
        except Shoebox columns, which are returned as NamedTuple objects
        which contains the portions of data from the file.

        With copy=False, all numpy arrays are pointing against the raw
        memory returned by msgpack, which means they are read-only.
        With copy=True, an immediate copy is done. This causes memory
        usage to double while loading, but the created numpy arrays own
        their own memory.
    """
    root_data = _get_unpacked(stream_or_path)

    if not root_data[0] == "dials::af::reflection_table":
        raise ValueError("Does not appear to be a dials reflection table file")
    if not root_data[1] == 1:
        raise ValueError(
            f"reflection_table data is version {root_data[1]}. Only Version 1 is understood"
        )
    refdata = root_data[2]

    rows = refdata["nrows"]
    data = refdata["data"]

    decoded_data = {
        name: decode_column(value, copy=copy) for name, value in data.items()
    }

    # Filter out empty (unknown) columns
    decoded_data = {k: v for k, v in decoded_data.items() if v is not None}

    # Cross-check the columns are the expected lengths
    for name, column in decoded_data.items():
        if len(column) != rows:
            print(
                f"Warning: Mismatch of column lengths: {name} is {len(column)} instead of expected {rows}"
            )

    return decoded_data


# Everything under here is optional stuff for demoing capabilities or
# generating and running test data
if __name__ == "__main__":
    import argparse
    import pprint

    def _write_test_file():
        import dials.array_family.flex as flex

        ref = flex.reflection_table()

        ref["bool"] = flex.bool([True, False] * 5)
        ref["int"] = flex.int(range(10))
        ref["std::size_t"] = flex.size_t(range(10))
        ref["double"] = flex.double(range(10))
        ref["vec2<double>"] = flex.vec2_double([(x + 1, x + 2) for x in range(10)])
        ref["vec3<double>"] = flex.vec3_double(
            [(x + 1, x + 2, x + 3) for x in range(10)]
        )
        ref["int6"] = flex.int6(
            [(x + 1, x + 2, x + 3, x + 4, x + 5, x + 6) for x in range(10)]
        )
        ref["cctbx::miller::index<>"] = flex.miller_index(
            [(x + 1, x + 2, x + 3) for x in range(10)]
        )
        ref["Shoebox<>"] = flex.shoebox(10)
        ref["mat3<double>"] = flex.mat3_double(
            [[x + y for y in range(9)] for x in range(10)]
        )

        ref.as_msgpack_file("test.refl")
        print(f"Written test reflection file {Path.cwd()/'test.refl'}")

    parser = argparse.ArgumentParser(
        description="Read a DIALS reflection table with only numpy"
    )
    parse_group = parser.add_mutually_exclusive_group(required=True)
    parse_group.add_argument(
        "--write-test",
        help="Write a test .refl file. Must be run inside cctbx environment.",
        action="store_true",
    )
    parse_group.add_argument(
        "FILE", help="Reflection filename to read", type=Path, nargs="?"
    )
    args = parser.parse_args()

    if args.write_test:
        try:
            _write_test_file()
        except ModuleNotFoundError:
            sys.exit(
                "Error: Could not import flex. Please run --write-test inside a cctbx environment"
            )
    else:
        pprint.pprint(load(args.FILE))


def test_reading():
    test_file = Path("test.refl")
    assert (
        test_file.is_file()
    ), "Please generate test file inside cctbx environment with 'libtbx.python refl_loader.py --write-test'"

    ref = load(test_file)
    expected_data = {
        "bool": np.array([True, False] * 5, dtype=bool),
        "int": np.array(range(10), dtype=np.int32),
        "std::size_t": np.array(range(10), dtype=np.uint64),
        "double": np.array(range(10), dtype=np.double),
        "vec2<double>": np.array([(x + 1, x + 2) for x in range(10)], dtype=np.double),
        "vec3<double>": np.array(
            [(x + 1, x + 2, x + 3) for x in range(10)], dtype=np.double
        ),
        "int6": np.array(
            [(x + 1, x + 2, x + 3, x + 4, x + 5, x + 6) for x in range(10)],
            dtype=np.int32,
        ),
        "cctbx::miller::index<>": np.array(
            [(x + 1, x + 2, x + 3) for x in range(10)], dtype=np.double
        ),
        "Shoebox<>": np.array([None] * 10, dtype=np.object_),
        "mat3<double>": np.array(
            [np.array([x + y for y in range(9)]).reshape((3, 3)) for x in range(10)],
            dtype=np.double,
        ),
    }
    unexpected_columns = set(ref.keys()) - set(expected_data.keys())
    assert not unexpected_columns

    # Go through each column and compare the value we got with expected
    for column, value in ref.items():
        assert (value == expected_data[column]).all()
