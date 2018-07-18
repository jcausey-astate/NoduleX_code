"""
Creates an HDF5 for CNN training/testing.
Outputs the rectangular solid around each point of interest, 
along with metadata:

images            - voxel intensity tensor representing the image volume
classes           - the nodule class
ids               - unique id for this image
patients          - patient id for the patient associated with this image
rescale_slope     - rescale slope value for this nodule (to scale into HU)
rescale_intercept - rescale intercept for this nodule (to scale into HU)

Other metadata may be added by including a header line in the `candidate_file`.
The column labels there will be used to add additional metadata if the --use-header
option is specified.

The candidate file header line must begin with a '#' character.

Required values in the candidate file:

    patientID   Z(sliceNo)   X(col)  Y(row)   ([class] |  class  [optional metadata if specified in a header line])

    Patient IDs must correspond to directory names below `dicom_dir`.


"""

from __future__ import print_function
import sys, os, re, h5py
import numpy as np
import matplotlib.pyplot as plt
import crop_3d as c3d
import volumetric_ct_image as ct3d
from get_tumor_centroids import *

__h5str = h5py.special_dtype(vlen=bytes)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def usage():
    print(__doc__.format(sys.argv[0]))


def get_args():
    import argparse

    # construct the argument parser and parse the arguments
    ap = argparse.ArgumentParser(
        prog="{0}".format(os.path.basename(sys.argv[0])),
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    ap.add_argument(
        "-s",
        "--size",
        dest="size",
        required=False,
        help="Size of each slice WxH.",
        default=47,
    )
    ap.add_argument(
        "-d",
        "--depth",
        dest="depth",
        required=False,
        help="Number of slices per cube.",
        default=5,
    )
    ap.add_argument(
        "-H",
        "--use-header",
        dest="use_header",
        action="store_true",
        help="Use header row to specify additional metadata columns.",
    )
    ap.add_argument(
        "-m",
        "--mm-per-voxel",
        dest="mm_per_voxel",
        required=False,
        default=1,
        help='The size of each output voxel in mm.  Use "-1" to keep\n'
        "original voxel resolution. Use (x,y,z) format to specify\n"
        "non-isometric voxels.",
    )
    ap.add_argument(
        "dicom_dir",
        metavar="dicom-dir",
        help="Directory containing DICOM source images",
    )
    ap.add_argument(
        "candidate_file",
        metavar="candidate-file",
        help="Candidate nodule coordinate file list.",
    )
    ap.add_argument(
        "-o",
        "--output",
        dest="output_file",
        metavar="output-file",
        required=False,
        help="Name of output file.",
    )

    args = vars(ap.parse_args())
    try:
        args["mm_per_voxel"] = float(args["mm_per_voxel"])
    except ValueError:
        try:
            # In case of a tuple, parse each part and convert to float:
            args["mm_per_voxel"] = [
                float(d) for d in args["mm_per_voxel"].strip()[1:-1].split(",")
            ]
        except:
            ap.error(
                "argument -m/--mm-per-voxel: invalid value: '{}'".format(
                    args["mm_per_voxel"]
                )
            )
    return args


def get_dicomdir_files(dicom_dir):
    files = []
    dicom_dir = os.path.normpath(dicom_dir)
    if os.path.isdir(dicom_dir):
        for rt, dr, fl in os.walk(dicom_dir):
            files.extend(fl)
            break
    else:
        raise RuntimeError("dicom-dir must be a directory")

    dicomfiles = [f for f in files if f[-4:] == ".dcm"]
    return dicomfiles


def read_candidates(dicom_dir, candidate_file, use_header=False):
    """
    Reads candidate locations and metadata from `candidate_file`.  The format must match 
    the following:

    patientID   Z(sliceNo)   X(col)  Y(row)   ([class] |  class  [optional metadata if specified in a header line])

    Patient IDs must correspond to directory names below `dicom_dir`.
    
    If a header line is specified, it must contain one column heading (not containing spaces) for each
    column of data (including the required ones AND the class -- the first five will be ignored).  
    If you need to specify the datatype, you can use a colon followed by one of {str, float, int}; example: 
        probability:float   stringClass:str   anIntValue:int
    If you do not specify the type, `float` will be used.
    """
    candidates = []
    file_list = []
    header_line = None
    data_started = False
    counts_by_patient = {}
    header_schema = None
    col_names = []

    with open(candidate_file, "rU") as fin:
        for idx, line in enumerate(fin):
            if line.strip() == "":
                continue
            line = line.strip()
            line_parts = line.split()
            patient_id = line_parts[0]

            if patient_id[0] == "#":
                if not data_started and header_line is None:
                    header_schema = {}
                    line = line[1:].strip()
                    line_parts = line.split()
                    col_names, header_schema = sanitize_headers(line_parts[5:])
                continue

            data_started = True
            try:
                counts_by_patient[patient_id] += 1
            except:
                counts_by_patient[patient_id] = 1

            cz, cx, cy = line_parts[1:4]
            class_rating = -1
            if len(line_parts) > 4:
                class_rating = line_parts[4]
            center = {}
            center["sliceNo"] = int(round(float(cz)))
            center["x"] = float(cx)
            center["y"] = float(cy)
            center["noduleID"] = "{}-{}".format(
                patient_id, counts_by_patient[patient_id]
            )
            center["class"] = int(class_rating) if int(class_rating) >= 0 else 0
            center["patient"] = patient_id
            center["dicom_directory"] = os.path.join(dicom_dir, center["patient"])

            # If there are additional metadata items, add them to the centroid dict:
            if use_header and header_schema is not None:
                line_parts = line_parts[5:]
                for idx, title in enumerate(col_names):
                    center[title] = line_parts[idx]

            candidates.append(center)
    if not use_header:
        header_schema = None
    return candidates, header_schema


def write_3D_solids(
    centroids, output_file, s=35, depth=None, additional_metadata=None, mm_per_voxel=1
):
    global __h5str
    s = int(s)
    depth = int(depth) if depth is not None else s
    labels_written = {}
    if os.path.splitext(output_file)[1] not in [".hdf5", ".hd5"]:
        output_file = os.path.splitext(output_file)[0] + ".hd5"

    # Sort centroids by by Z offset
    centroids.sort(key=lambda x: int(x["sliceNo"]))
    # Then by dicom_dir
    centroids.sort(key=lambda x: x["dicom_directory"])

    with h5py.File(output_file, "w") as fout:
        shape1d = (len(centroids), 1)
        output_images = fout.create_dataset(
            "image",
            (len(centroids), depth, s, s),
            maxshape=(len(centroids), depth, s, s),
            dtype="float32",
        )
        output_classes = fout.create_dataset(
            "class", shape1d, maxshape=shape1d, dtype="int"
        )
        output_patient = fout.create_dataset(
            "patient", shape1d, maxshape=shape1d, dtype=__h5str
        )
        output_nodule_id = fout.create_dataset(
            "id", shape1d, maxshape=shape1d, dtype=__h5str
        )
        output_pixel_min = fout.create_dataset(
            "pixel_min", shape1d, maxshape=shape1d, dtype="float32"
        )
        output_pixel_max = fout.create_dataset(
            "pixel_max", shape1d, maxshape=shape1d, dtype="float32"
        )
        output_rescale_slope = fout.create_dataset(
            "rescale_slope", shape1d, maxshape=shape1d, dtype="float32"
        )
        output_rescale_intercept = fout.create_dataset(
            "rescale_intercept", shape1d, maxshape=shape1d, dtype="float32"
        )
        if additional_metadata is not None:
            for key in additional_metadata:
                additional_metadata[key] = {
                    "dataset": fout.create_dataset(
                        key, shape1d, maxshape=shape1d, dtype=additional_metadata[key]
                    ),
                    "dtype": additional_metadata[key],
                }

        cached_scan = {
            "dicom_directory": None,
            "min": 0,
            "max": 0,
            "scan": None,
            "fail_flag": False,
        }
        n_idx = (
            0
        )  # Some samples may fail; don't use enumerate, so we can skip bad nodules
        for nodule in centroids:  # without leaving "blank" records in the hd5 file.
            cx = int(round(nodule["x"]))
            cy = int(round(nodule["y"]))
            cslice_no = int(nodule["sliceNo"])
            nodule_class = int(nodule["class"])
            nodule_patient = nodule["patient"]
            nodule_dir = nodule["dicom_directory"]
            if nodule_dir == "":
                print(
                    "Warning: missing nodule directory for nodule {}".format(nodule),
                    file=sys.stderr,
                )
                continue

            nodule_id = "{}~{}-{}-{}~{}".format(
                nodule_patient, cx, cy, cslice_no, nodule_class
            )
            label = nodule_class

            # Get the standardized image cube in slice-major order [slice,r,c]
            full_scan = cached_scan["scan"]
            if nodule_dir != cached_scan["dicom_directory"]:
                try:
                    full_scan = ct3d.StandardizedScan(
                        nodule_dir, mm_per_voxel=mm_per_voxel
                    )
                    cached_scan["fail_flag"] = False
                except Exception as e:
                    eprint(
                        "FAILED: {} failed with error: {}".format(
                            nodule_patient, str(e)
                        )
                    )
                    cached_scan["dicom_directory"] = nodule_dir
                    cached_scan["fail_flag"] = True
                    cached_scan["scan"] = None
                    continue

            if cached_scan["fail_flag"]:
                continue

            mapped_center = full_scan.map_original_coordinates((cslice_no, cy, cx))
            result_cube = c3d.crop_centered_at_point(
                full_scan.img, mapped_center, (depth, s, s)
            )

            pix_min, pix_max = None, None
            if nodule_dir == cached_scan["dicom_directory"]:
                pix_min, pix_max = cached_scan["pix_min"], cached_scan["pix_max"]
            else:
                pix_min, pix_max = full_scan.img.min(), full_scan.img.max()
                cached_scan["pix_min"], cached_scan["pix_max"] = pix_min, pix_max
                cached_scan["dicom_directory"] = nodule_dir
                cached_scan["scan"] = full_scan

            del full_scan

            output_images[n_idx] = result_cube
            output_patient[n_idx] = nodule_patient
            output_nodule_id[n_idx] = nodule_id
            output_classes[n_idx] = nodule_class
            output_pixel_min[n_idx] = pix_min
            output_pixel_max[n_idx] = pix_max
            output_rescale_slope[n_idx] = 1.0
            output_rescale_intercept[n_idx] = 0
            if additional_metadata is not None:
                for key in additional_metadata:
                    additional_metadata[key]["dataset"][n_idx] = astype(
                        nodule[key], additional_metadata[key]["dtype"]
                    )

            if label in labels_written:
                labels_written[label] += 1
            else:
                labels_written[label] = 1
            n_idx += (
                1
            )  # Now we are certain that the nodule was written OK; increment the index.
        # Now, resize the datasets to "fit" the number of values actually written if necessary:
        if n_idx != len(centroids):
            for key in fout.keys():
                shape = list(fout[key].shape)
                shape[0] = n_idx
                shape = tuple(shape)
                fout[key].resize(shape)

    print(
        "Wrote {0} cubes; labels: {1}".format(
            sum([labels_written[l] for l in labels_written]),
            ["(" + str(l) + ":" + str(labels_written[l]) + ")" for l in labels_written],
        )
    )


def astype(value, to_type):
    """
    Cast the `value` to the selected `to_type`, which must be one of 
    {'str', 'int', 'float'} (technically, any value other than 'str' or 'int'
    will be cast to float32 by default).
    """
    global __h5str

    if to_type in ["str", __h5str]:
        return str(value)
    elif to_type == "int":
        return int(value)
    else:
        return np.float32(value)


def sanitize_headers(header_line):
    """
    Remove any non-alphanumeric characters from elements of the list 
    `header_line`, and make sure each element begins with a non-digit.
    Also, look for type specifications of the form ':type' as a prefix
    to header titles.  Returns a list of header titles and a dict of
    titles => types.
    """
    global __h5str

    header_schema = {}
    titles = []
    for idx, text in enumerate(header_line):
        header_parts = text.split(":")
        text = header_parts[0]
        field_type = "float32"

        if len(header_parts) > 1:
            field_type = header_parts[1].lower()
            if field_type not in ["int", "float", "str"]:
                field_type = "float32"
            elif (
                field_type == "float"
            ):  # 'int' is OK already; float=>float32 and 'str'=>__h5str
                field_type = "float32"
            elif field_type == "str":
                field_type = __h5str

        text = re.sub("[^0-9a-zA-Z]+", "_", text.strip())
        if re.match("^[0-9]", text) is not None:
            text = "_{}".format(text)
        titles.append(text)
        header_schema[text] = field_type
    return titles, header_schema


def run_standalone():
    args = get_args()
    dicom_dir = args["dicom_dir"]
    centroids, file_list = [], []
    centroids, header_schema = read_candidates(
        dicom_dir, args["candidate_file"], use_header=args["use_header"]
    )
    write_3D_solids(
        centroids,
        output_file=args["output_file"]
        if args["output_file"] is not None
        else args["candidate_file"],
        s=args["size"],
        depth=args["depth"],
        additional_metadata=header_schema,
        mm_per_voxel=args["mm_per_voxel"],
    )


if __name__ == "__main__":
    run_standalone()
