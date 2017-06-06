"""
Creates an HDF5 for CNN training/testing.
Outputs the rectangular solid around each point of interest, 
along with metadata:

nodule_images            - voxel intensity tensor representing the image volume
nodule_classes           - the nodule class
nodule_malignancy        - malignancy rating
nodule_inclusion         - 1 if this is an "inclusion' 0 if it is non-nodule, nodule < 3mm, or is an "exclusion"        
nodule_is_nodule         - 1 if this is a nodule >= 3mm       
nodule_patients          - patient id for this nodule
nodule_ids               - unique id for this nodule
nodule_pixel_min         - minimum intensity for the scan this nodule was selected from
nodule_pixel_max         - maximum intensity for the scan this nodule was selected from
nodule_rescale_slope     - rescale slope value for this nodule (to scale into HU)
nodule_rescale_intercept - rescale intercept for this nodule (to scale into HU)

Run with --help to see command options.
"""

from __future__ import print_function
import sys, dicom, re, os, pylab, math, json, h5py
import pickle
import numpy as np, scipy.signal as sp
from scipy.signal import argrelextrema
from scipy.spatial import distance
from PIL import Image
from get_tumor_centroids import *
import matplotlib.pyplot as plt

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def usage():
    print(__doc__.format(sys.argv[0]))

def get_args():
    import argparse
    # construct the argument parser and parse the arguments
    ap = argparse.ArgumentParser(prog='{0}'.format(os.path.basename(sys.argv[0])), description=__doc__)
    ap.add_argument("-s", dest = 'size', required = False, help = "Size of each slice WxH.", default=47)
    ap.add_argument("-d", dest = 'depth', required = False, help = "Number of slices per cube.", default=5)
    ap.add_argument("-c", dest = 'n_classes', required = False, help = "Number of classes.", default=2)
    ap.add_argument("-r", "--raw-intensity", dest = 'raw_intensity', action='store_true', 
        help = "Do not re-scale intensity according to DICOM header before extracting cube.")
    ap.add_argument("dicom_dir", metavar='dicom-dir' , help = "Directory containing DICOM source images")
    ap.add_argument("candidate_file", metavar='candidate-file', help = "Candidate nodule coordinate file list.")
    ap.add_argument("-o", dest = 'output_file', metavar='output-file', required = False, help = "Name of output file.")
    
    args = vars(ap.parse_args())
    return args

def get_dicomdir_files(dicom_dir):
    files       = []
    dicom_dir = os.path.normpath(dicom_dir)
    if os.path.isdir(dicom_dir):
        for rt, dr, fl in os.walk(dicom_dir):
            files.extend(fl)
            break
    else:
        raise RuntimeError("dicom-dir must be a directory")

    dicomfiles = [f for f in files if f[-4:] == '.dcm']
    return dicomfiles

def read_candidates(dicom_dir, candidate_file, return_file_list=True):
    candidates  = []
    file_list   = []
    multiple_patients = False

    with open(candidate_file , 'rU') as fin:
        for idx, line in enumerate(fin):
            line_parts   = line.split()
            patient_name = None
            if len(line_parts) > 6:
                multiple_patients = True
                patient_name      = line_parts[0]
                line_parts        = line_parts[1:]
            cz, cx, cy, is_nodule, is_included, malignancy = line_parts
            center = {}
            center['sliceNo']    = int(cz)
            center['x']          = float(cx)
            center['y']          = float(cy)
            center['noduleID']   = 'candidate-{0}'.format(idx)
            center['isNodule']   = int(is_nodule)    if int(is_nodule)   >= 0 else 0
            center['inclusion']  = int(is_included)  if int(is_included) >= 0 else 0
            center['malignancy'] = int(malignancy)   if int(malignancy)  >= 0 else 0
            center['is_auto']    = 0 if int(is_included) >= 0 else 1
            if patient_name is not None:
                center['patient']    = patient_name
                center['dicom_directory'] = os.path.join(dicom_dir, center['patient'])
            candidates.append(center)
    patients    = {}
    dicom_files = []
    if multiple_patients:
        for candidate in candidates:
            new_dicom_files = get_dicomdir_files(os.path.join(dicom_dir, candidate['patient']))
            for idx, fname in enumerate(new_dicom_files):
                new_dicom_files[idx] = os.path.join(candidate['patient'], fname)
            dicom_files.extend(new_dicom_files)
    else:
        dicom_files = get_dicomdir_files(dicom_dir)
    files_by_slice = {}
    for df in dicom_files:
        slice_key = ""
        if multiple_patients:
            parts = os.path.split(df)
            slice_key = "{}-".format(parts[0])
        ddat      = dicom.read_file(os.path.join(dicom_dir, df))
        try:
            zCoord    = ddat.SliceLocation
        except:
            zCoord    = ddat.InstanceNumber
        try:
            sopUID    = ddat.SOPInstanceUID
        except:
            import uuid
            sopUID    = os.path.splitext(os.path.basename(df))[0] + str(uuid.uuid4().hex)
        sliceNo   = ddat.InstanceNumber
        slope     = 1.0
        intercept = 0.0
        try:
            slope     = np.float32(ddat.RescaleSlope)
            intercept = np.float32(ddat.RescaleIntercept)
        except AttributeError:  # In case the fields aren't filled in the dicom header
            pass
        fileName  = df
        slice_key = "{}{}".format(slice_key, sliceNo)
        files_by_slice[slice_key] = {'z': zCoord, 'fileName': df, 'sopUID': sopUID, 'filePath': os.path.join(dicom_dir, df), 'slope': slope, 'intercept': intercept}
        file_list.append({'sliceNo':sliceNo, 'fileName':os.path.join(dicom_dir,df)})
    for idx, center in enumerate(candidates):
        slice_key = "{}-{}".format(center['patient'], center['sliceNo']) if center['patient'] is not None else center['sliceNo']
        file_info = files_by_slice[slice_key] if slice_key in files_by_slice else None
        center['z']                = float(file_info['z']) if file_info is not None else float(cz)
        center['sopUID']           = file_info['sopUID'] if file_info is not None else 'Unknown'
        center['fileName']         = file_info['fileName'] if file_info is not None else 'Unknown'
        center['filePath']         = file_info['filePath'] if file_info is not None else ''
        center['rescaleSlope']     = file_info['slope'] if file_info is not None else 1.0
        center['rescaleIntercept'] = file_info['intercept'] if file_info is not None else 0.0
        center['dicom_directory']  = os.path.split(center['filePath'])[0]
        if center['dicom_directory'] == '':
            center['dicom_directory'] = os.path.join(dicom_dir, center['patient'])
        candidates[idx]            = center
    return (candidates, file_list) if return_file_list else candidates

def write_3D_solids(centroids, file_list, s=47, output_file=None, classes=2, depth=5, rescale_intensity=False):
    s       = int(s)
    classes = int(classes)
    depth   = int(depth)
    labels_written = {}
    if output_file.split('.')[-1] not in ['hdf5', 'hd5']:
        output_file += '.hd5'

    # Sort centroids by by Z offset
    centroids.sort(key=lambda x: int(x['z']))
    # Then by dicom_dir
    centroids.sort(key=lambda x: x['dicom_directory'])

    with h5py.File(output_file, "w") as fout:
        h5str    = h5py.special_dtype(vlen=bytes)
        output_X = fout.create_dataset("nodule_images",  (len(centroids), depth, s, s), dtype='float32')
        output_Y = fout.create_dataset("nodule_classes",  (len(centroids), 1,), dtype='int')
        output_malignancy = fout.create_dataset("nodule_malignancy",  (len(centroids), 1,), dtype='int')
        output_inclusion  = fout.create_dataset("nodule_inclusion",  (len(centroids), 1,), dtype='int')
        output_is_nodule  = fout.create_dataset("nodule_is_nodule",  (len(centroids), 1,), dtype='int')
        output_patient    = fout.create_dataset("nodule_patients",  (len(centroids), 1,), dtype=h5str)
        output_nodule_id  = fout.create_dataset("nodule_ids",  (len(centroids), 1,), dtype=h5str)
        output_pixel_min  = fout.create_dataset("nodule_pixel_min",  (len(centroids), 1,), dtype='float32')
        output_pixel_max  = fout.create_dataset("nodule_pixel_max",  (len(centroids), 1,), dtype='float32')
        output_rescale_slope      = fout.create_dataset("nodule_rescale_slope",      (len(centroids), 1,), dtype='float32')
        output_rescale_intercept  = fout.create_dataset("nodule_rescale_intercept",  (len(centroids), 1,), dtype='float32')
        cached_min_max            = {'dicom_directory':None, 'min':0, 'max':0}
        for n_idx, nodule in enumerate(centroids):
            cx = int(round(nodule['x']))
            cy = int(round(nodule['y']))
            cslice_no   = int(nodule['sliceNo'])
            nodule_is_nodule  = int(nodule['isNodule'])
            nodule_inclusion  = int(nodule['inclusion'])
            nodule_malignancy = int(nodule['malignancy'])
            nodule_patient    = nodule['patient']
            nodule_dir        = nodule['dicom_directory']
            if nodule_dir == '':
                print("Warning: empty nodule directory for nodule {}".format(nodule), file=sys.stderr)
                continue
            nodule_rescale_slope     = nodule['rescaleSlope']
            nodule_rescale_intercept = nodule['rescaleIntercept']
            nodule_id                = "{}~{}-{}-{}~{}-{}-{}".format(nodule_patient, cx, cy, cslice_no, nodule_is_nodule, nodule_inclusion, nodule_malignancy)
            label                    = nodule_malignancy
            if nodule_is_nodule and nodule_inclusion == 0:
                continue      # inclusion False (for a nodule) indicates a subtraction from another nodule - not useful here
            if classes == 2:
                label = 0 if nodule_is_nodule == 0 else ( 1 if nodule_malignancy >= 1 else -1)  # 1 for nodules >= 3mm, 0 for non-nodules, -1 for nodules < 3mm
            elif classes == 3:
                label = 0 if label in [0,1] else 1 if label in [2,3] else 2
            nodule['label']   = label
            dicom_files       = get_dicomdir_files(nodule_dir)
            # print("nodule dir {}: {}".format(n_idx + 1, nodule_dir))
            files_by_slice = {}
            for df_name in dicom_files:
                df_name  = os.path.join(nodule_dir,df_name)
                df = dicom.read_file(df_name)
                slice_no = df.InstanceNumber
                files_by_slice[slice_no] = df_name
            hd           = int(depth / 2)
            hs           = int(s/2)
            result_cube  = np.zeros((depth,s,s))
            mid_slice    = cslice_no
            first_slice  = mid_slice - hd
            last_slice   = first_slice + depth - 1
            # Get all slices between z_min, z_max
            included_slices = range(first_slice, last_slice + 1)
            # print("{0}: {1}".format(cslice_no, included_slices))
            for idx, slice_no in enumerate(included_slices):
                if slice_no in files_by_slice:
                    pixels              = dicom.read_file(files_by_slice[slice_no]).pixel_array.astype('float32')
                    ysize,xsize         = pixels.shape
                    result_cube[idx]    = np.zeros((s,s)).astype('float32')
                    left_x              = cx - hs
                    top_y               = cy - hs
                    left_x_offset       = abs(left_x) if left_x < 0 else 0
                    top_y_offset        = abs(top_y)  if top_y  < 0 else 0
                    right_x             = left_x + s
                    bottom_y            = top_y  + s
                    left_x              = max(0, left_x)
                    top_y               = max(0, top_y)
                    right_x_offset      = xsize - right_x  if right_x  > xsize else s
                    bottom_y_offset     = ysize - bottom_y if bottom_y > ysize else s
                    right_x             = min(right_x, xsize)
                    bottom_y            = min(bottom_y, ysize)
                    # print('Extracting cube from {} : cx{} cy{} cz{} left_x{} right_x{} top_y{} bot_y{} o_left_x{}, o_right_x{}, o_top_y{} o_bot_y{}'.format(
                    #         nodule['patient'], cx, cy, cslice_no, left_x, right_x, top_y, bottom_y, 
                    #         left_x_offset, right_x_offset, top_y_offset, bottom_y_offset
                    #     ))
                    result_cube[idx][top_y_offset:bottom_y_offset, left_x_offset:right_x_offset] \
                        = pixels[top_y:bottom_y, left_x:right_x].copy()
                    r_size = result_cube[idx].shape
                else:
                    pixels = np.zeros((s,s)).astype('float32')
                    result_cube[idx] = pixels
            # Write image cube in slice-major order [slice,r,c]
            pix_min, pix_max = None, None
            if nodule_dir == cached_min_max['dicom_directory']:
                pix_min, pix_max = cached_min_max['pix_min'], cached_min_max['pix_max']
            else:
                pix_min, pix_max = get_dicom_min_max(nodule_dir)
                cached_min_max['pix_min'], cached_min_max['pix_max'] = pix_min, pix_max
                cached_min_max['dicom_directory'] = nodule_dir
            if rescale_intensity:  # If user requests intensity rescaling here, do so and adjust the parameters to reflect the change before writing to the file
                result_cube = (result_cube * nodule_rescale_slope) + nodule_rescale_intercept
                pix_min     = (pix_min * nodule_rescale_slope)     + nodule_rescale_intercept
                pix_max     = (pix_max * nodule_rescale_slope)     + nodule_rescale_intercept
                nodule_rescale_slope     = 1.0
                nodule_rescale_intercept = 0.0

            output_Y[n_idx]          = [int(nodule['label'])]
            output_X[n_idx]          = result_cube
            output_patient[n_idx]    = nodule_patient
            output_nodule_id[n_idx]  = nodule_id
            output_is_nodule[n_idx]  = nodule_is_nodule
            output_inclusion[n_idx]  = nodule_inclusion
            output_malignancy[n_idx] = nodule_malignancy
            output_pixel_min[n_idx]  = pix_min
            output_pixel_max[n_idx]  = pix_max
            output_rescale_slope[n_idx]     = nodule_rescale_slope
            output_rescale_intercept[n_idx] = nodule_rescale_intercept
            if nodule['label'] in labels_written:
                labels_written[nodule['label']] += 1
            else:
                labels_written[nodule['label']] = 1
            #plt.imshow(sbs, cmap='bone')
            #plt.show()     
    if rescale_intensity:
        print("Intensity rescaled for all output cubes.")       
    print("Wrote {0} cubes; labels: {1}".format(sum([labels_written[l] for l in labels_written]), ["("+str(l)+":"+str(labels_written[l])+")" for l in labels_written]))

def get_dicom_min_max(dicom_dir):
    files = get_dicomdir_files(dicom_dir)
    if len(files) == 0:
        raise RuntimeError("dicom_dir '{}' seems to contain no DICOM files".format(dicom_dir))
    min_intensity = []
    max_intensity = []
    for fname in files:
        fname   = os.path.join(dicom_dir, fname)
        fp      = dicom.read_file(fname)
        pixels  = fp.pixel_array
        pix_min = pixels.min()
        pix_max = pixels.max()
        min_intensity.append(pix_min)
        max_intensity.append(pix_max)
    return min(min_intensity), max(max_intensity)

def run_standalone():
    args = get_args()
    dicom_dir = args['dicom_dir']
    centroids, file_list = [],[]
    centroids, file_list = read_candidates(dicom_dir, args['candidate_file'], return_file_list=True)
    rescale_intensity    = not args['raw_intensity']
    write_3D_solids(centroids, file_list, args['size'], args['output_file'], args['n_classes'], depth=args['depth'], rescale_intensity=rescale_intensity)

if __name__ == "__main__":
    run_standalone()
