# -*- coding: utf-8 -*-
"""
Created on Mon May 11 17:38:04 2020

@author: Connor
"""
import warnings
import os.path
from Series import DataSet


def load_file(path):
    """Return loader function depending on extension"""
    # Split parameters
    if "?" in path:
        filename, param_string = path.split("?", 1)

        # Use "&" as separator of parameters
        if "?" in param_string:
            warnings.warn("The use of '?' as separator of dataset parameters is deprecated.", DeprecationWarning)
            parts = param_string.split("?")
        else:
            parts = param_string.split("&")
    else:
        filename = path
        parts = []
    # Use extension as type
    type = os.path.splitext(filename)[1].lower()
    if type:
        type = type[1:]
    # Parse parameters
    param = {}
    for part in parts:
        if not part:
            continue
        subparts = part.split('=', 1)
        key = subparts[0]
        if len(subparts) != 2:
            value = ''
        else:
            value = subparts[1]
        param[key] = value

    # Override file type
    type = param.pop("type", type)
    if type == "dm3":
        return DataSet.load_dm3(filename)
    elif type == "hdf5" or type == "h5":
        if "dataset" in param:
            dataset = param["dataset"]
        elif (len(parts) == 1) and not ("=" in parts[0]):
            dataset = parts[0]
            warnings.warn("Passing the dataset name after question mark for HDF5 files is deprecated.", DeprecationWarning)
        else:
            raise ValueError("Parameter 'dataset' missing for HDF5 file: %s" % path)
        return DataSet.load_hdf5(filename, dataset)
    elif type == "raw":
        shape = int(param["ysize"]), int(param["xsize"])
        dtype = param["dtype"]
        offset = int(param.get("offset", 0))
        swap_bytes = int(param.get("swap_bytes", 0))
        return DataSet.load_raw(filename, shape, dtype, offset=offset, swap_bytes=swap_bytes)
    else:
        raise ValueError("Unrecognized image extension.")