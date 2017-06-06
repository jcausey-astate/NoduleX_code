"""
Creates a list of unique tumor centroids given one or more LIDC-IDRI 
DICOM directories.  

Run with --help to see usage info.

Centroid file is in the form:,

patient Z X Y is_nodule is_included malignancy

unless the --short option is specified; in that case the format excluded the
'patient' field, and output files will be saved in the individual dicom directories.

The --combine option overrides --short, and requires the long form, producing a single
file with all centroids listed.

Uniqueness of nodules is determined by plotting all nodules with segmentations, then clustering
such that connected regions become part of the same resulting nodule.
Non-nodules (with no segmentation) are plotted as "spheres" of diameter DIAMETER-mm, but will never
add to the size of a segmentation (they only "silently" group if they intersect, removing them from
the output set.
If --nodules-only is specified, non-nodules will not be listed at all.
"""

from __future__ import print_function
import sys, dicom, re, os
import numpy as np
from scipy.ndimage import measurements, morphology
from PIL import Image, ImageDraw
from get_tumor_centroids import *
from candidate_dict_utilities import *

import matplotlib.pyplot as plt

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_args():
    import argparse
    # construct the argument parser and parse the arguments
    ap = argparse.ArgumentParser(prog='{0}'.format(os.path.basename(sys.argv[0])))
    ap.add_argument("dicom_dir", metavar='dicom-dir', nargs="+", help = "Directory (or) containing DICOM source images")
    ap.add_argument("--short", dest = 'short', action='store_true', help = "Use short output form (no patient name)")
    ap.add_argument("--combine", dest = 'combine', action='store_true', help = "Combine output into a single file")
    ap.add_argument("--nodules-only", dest = 'nodules_only', action='store_true', help = "Only output centroids for nodules with segmentation")
    ap.add_argument("--diameter", dest = 'diameter_mm', metavar="diameter-mm", required = False, help = "Diameter (mm) for spheres representing non-nodules; default=5", default=5)
    ap.add_argument("--output", dest = 'output_file', metavar="output-filename", required=False, help = "Filename for output file")
    
    args = vars(ap.parse_args())
    return args

def coords_to_polygon(coords):
    # it = iter(coords)
    # return izip(it,it)
    return zip(coords[::2], coords[1::2])

def draw_polygon_from_coords(pixels, coords, outline=1):
    if outline == True:
        outline  = 1
    pixels       = np.array(pixels)
    boundary     = coords_to_polygon(coords)
    (rows, cols) = pixels.shape
    img          = Image.new('L', (cols, rows), 0)
    ImageDraw.Draw(img).polygon(boundary, outline=outline, fill=1)
    pixels[np.array(img) == 1] = 1
    return pixels

def draw_circle(pixels, center, radius=3, outline=1):
    if outline == True:
        outline = 1
    pixels      = np.array(pixels)
    x,y,r       = center[0], center[1], radius - 0.5
    if r > 0:
        (rows, cols) = pixels.shape
        img          = Image.new('L', (cols, rows), 0)
        ImageDraw.Draw(img).ellipse((x-r, y-r, x+r, y+r), outline=outline, fill=1)
        pixels[np.array(img) == 1] = 1
    return pixels

def collect_centroids_from_xml(dicom_dirs):
    if type(dicom_dirs) != type([]):
        dicom_dirs = [dicom_dirs]
    centroids_by_scan = {}
    for dicom_dir in dicom_dirs:
        centroids, file_list = get_tumor_annotations(dicom_dir, return_file_list=True)
        centroids_by_scan[dicom_dir] = {'centroids':centroids, 'file_list':file_list}
    return centroids_by_scan

def remove_non_nodules(centroids_by_scan):
    nodules_by_scan = {}
    for scan in centroids_by_scan:
        for candidate in centroids_by_scan[scan]['centroids']:
            if 'coords' in candidate and len(candidate['coords']) > 2:
                if scan not in nodules_by_scan:
                    nodules_by_scan[scan] = {'centroids': [], 'file_list': centroids_by_scan[scan]['file_list']}
                nodules_by_scan[scan]['centroids'].append(candidate)
    return nodules_by_scan

def get_slice_range(dicom_files):
    if type(dicom_files) == type(""):
        dicom_dir   = dicom_files
        dicom_files = get_dicomdir_files(dicom_dir)
        dicom_files = [os.path.join(dicom_dir, x) for x in dicom_files]
    min_idx = 9e9
    max_idx = 0
    for fname in dicom_files:
        sliceNo = 0
        if type(fname) == type({}):
            fname = fname['fileName']
        df = dicom.read_file(fname)
        sliceNo = int(df.InstanceNumber)
        df = None
        if sliceNo < min_idx:
            min_idx = sliceNo
        if sliceNo > max_idx:
            max_idx = sliceNo
    return min_idx, max_idx

def merge_centroids(centroids_by_scan, non_nodule_diameter=5):
    merged_by_scan = {}
    for scan in centroids_by_scan:
        # scan is a dicom dir; load the file list:
        dicom_dir        = scan
        candidates       = centroids_by_scan[dicom_dir]['centroids']
        file_list        = centroids_by_scan[dicom_dir]['file_list']
        candidates.sort(key=lambda c: c['sliceNo'])
        pixels           = None
        cpixels          = None
        sphere_radii     = []
        sphere_radius_mm = non_nodule_diameter / 2.0
        n_slices         = None
        min_slice        = candidates[0]['sliceNo']
        segmentations    = []
        centroids        = []
        min_slice_idx, max_slice_idx = get_slice_range(file_list)
        # If the scan doesn't start at slice 1, pad the beginning - inefficient, 
        # but this is rare and it is less error prone than calculating offset indexes everywhere.
        padding_at_start = min_slice_idx - 1
        candidates.sort(key=lambda x: x['sliceNo'])
        candidates.sort(key=lambda x: 'coords' in x and len(x['coords']) > 2, reverse=True)
        candidates = fill_file_paths(candidates, dicom_dir)
        for candidate in candidates:
            # print("candidate slice {} x{} y{}".format(candidate['sliceNo'], candidate['x'], candidate['y']))
            filePath   = candidate['filePath']
            df         = dicom.read_file(filePath)
            sliceIndex = int(candidate['sliceNo']) - 1
            if pixels is None:  # If we are initializing this scan...
                n_slices   = len(file_list)
                img_width  = int(df.Columns)
                img_height = int(df.Rows)
                pixels     = np.zeros((n_slices + padding_at_start, img_height, img_width), dtype=int) # holds real segmentations
                cpixels    = np.zeros((n_slices + padding_at_start, img_height, img_width), dtype=int) # holds pseuo-segmentations from non-nodules
                # Calculate radii-per-slice (in pixels) for the "dummy" sphere used with non-nodules:
                y_scale, x_scale   = [float(v) for v in df.PixelSpacing]
                z_scale            = float(df.SliceThickness)
                xy_scale           = (x_scale + y_scale) / 2.0
                if xy_scale == 0:  # Information is sometimes missing; assume unit scale if this happens
                    xy_scale = 1
                if z_scale  == 0:
                    z_scale  = 1
                sphere_radius_px   = int(round(sphere_radius_mm / xy_scale))
                sphere_z_radius_px = int(round(sphere_radius_mm / z_scale))
                for i in range(-sphere_z_radius_px, sphere_z_radius_px+1, 1):
                    offset_mm   = np.abs(i * z_scale)
                    if offset_mm < sphere_radius_mm:
                        r_offset_mm = np.sqrt(sphere_radius_mm**2 - offset_mm**2)
                        r_offset_px = int(round(r_offset_mm / xy_scale))
                        sphere_radii.append(r_offset_px)
            # draw the polygon if we have one:
            if 'coords' in candidate and len(candidate['coords']) > 2:
                # print("plot polygon {} on slice {} ".format(nodule_unique_id(candidate), sliceIndex))
                pix_slice  = pixels[sliceIndex]
                pix_slice  = pix_slice + draw_polygon_from_coords(pix_slice, candidate['coords'])
                pix_slice[pix_slice > 0] = 1
                pixels[sliceIndex] = pix_slice
                candidate['area']  = pix_slice.sum()
                segmentations.append(candidate)
            else:
                # Sphere radius in other slices =  r_offset = sqrt(R^2 - d^2) 
                # where d is the offset distance from
                # the sphere's central slice and R is the sphere's radius.
                candidate['area']  = 0
                centroids.append(candidate)
                first_slice = int(sliceIndex - (len(sphere_radii) / 2))
                first_i     = 0
                last_i      = len(sphere_radii)
                last_slice  = first_slice + len(sphere_radii)
                if first_slice < 0:
                    first_i += np.abs(first_slice)
                    first_slice = 0
                if last_slice >= max_slice_idx:
                    last_i     -= (last_slice - max_slice_idx)
                    last_slice  = max_slice_idx
                for i in range(last_i - first_i):
                    slice_i  = first_slice + i
                    r_i      = first_i + i
                    # print("plot candidate {} on slice {} from center slice {}".format(nodule_unique_id(candidate), slice_i, sliceIndex))
                    s_radius = sphere_radii[r_i]
                    center   = (candidate['x'], candidate['y'])
                    cpix_slice = cpixels[slice_i]
                    cpix_slice = cpix_slice + draw_circle(cpix_slice, center, s_radius)
                    cpix_slice[cpix_slice > 0] = 1
                    cpixels[slice_i] = cpix_slice
        # Now we have plotted all polygons and spheres; create groups.
        structure             = morphology.generate_binary_structure(3,2)
        n_centroid_groups     = 0
        n_segmentation_groups = 0
        if cpixels is not None:
            n_centroid_groups     = measurements.label(cpixels, output=cpixels, structure=structure)
        if pixels is not None:
            n_segmentation_groups = measurements.label(pixels, output=pixels, structure=structure)
                
        # slices = [x for x in range(len(pixels)) if pixels[x].sum() > 0]
        # for s in slices:
        #     print("groups: {} - in slice {} : {} : {}".format(n_segmentation_groups, s, len(np.unique(pixels[s])) - 1, np.unique(pixels[s])))
        #     plt.imshow(pixels[s])
        #     plt.show()

        # cslices = [x for x in range(len(cpixels)) if cpixels[x].sum() > 0]
        # for s in cslices:
        #     print("groups: {} - in slice {} : {} : {}".format(n_centroid_groups, s, len(np.unique(cpixels[s])) - 1, np.unique(cpixels[s])))
        #     plt.imshow(cpixels[s])
        #     plt.show()
        
        # print("n_centroid_groups: {}, n_segmentation_groups: {}".format(n_centroid_groups, n_segmentation_groups))
        centroid_groups       = [[] for i in range(n_centroid_groups)]
        centroid_labels       = range(1,n_centroid_groups+1)
        ignore_cgroups        = []
        segmentation_groups   = [[] for i in range(n_segmentation_groups)]
        segmentation_labels   = range(1,n_segmentation_groups+1)
        # Make note of any centroid groups that have overlap with a segmentation group:
        for i in range(n_centroid_groups):
            group = i + 1
            if pixels[cpixels == group].sum() > 0:
                # print("Adding ignore cgroup {}".format(group))
                ignore_cgroups.append(group)
        # Now group candidate listings:      
        for centroid in centroids:
            r, c   = int(round(centroid['y'])), int(round(centroid['x']))
            s      = int(centroid['sliceNo']) - 1
            try:
                cgroup = int(cpixels[s,r,c])
            except Exception as e:
                if r >= 512 or c >= 512 or r < 0 or c < 0:
                    print("ILLEGAL ANNOTATION on {} r:{} c:{} s:{} \n{}".format(centroid['patient'], r,c,s, centroid))
                    continue  # A coordinate out-of-bounds is an error in the XML; ignore the annotation and move on.
                else: # Something unexpected went wrong.
                    print("FAIL on {} r:{} c:{} s:{} \n{}\n{}".format(centroid['patient'], r,c,s, e, centroid))
                    sys.exit(1)
            if cgroup <= 0:
                # This can happen for "concave" polygon shapes, where center falls into a "cutout" in the middle.
                objects = measurements.find_objects(cpixels)
                was_placed = False
                for g_idx, box in enumerate(objects):
                    s_start, s_end = box[0].start, box[0].stop # objects are given by slices, which have start, stop attributes
                    ulr, lrr = box[1].start, box[1].stop       # one slice per dimension, so we have three: one for slices,
                    ulc, lrc = box[2].start, box[2].stop       # one for rows and one for cols
                    if s >= s_start and s < s_end and r >= ulr and r < lrr and c >= ulc and c < lrc:
                        was_placed = True
                        cgroup = g_idx + 1
                        break
                # plt.imshow(cpixels[s])
                # plt.show()
                if not was_placed:
                    print("WARN: centroid did not map to a group! p:{} x:{} y:{} s:{}".format(scan,c,r,s))
                    print(centroid)
                    continue
            if not cgroup in ignore_cgroups: # Don't add to groups that overlap segmentations
                centroid_groups[cgroup-1].append(centroid)
            else:
                # print("Not adding centroid group {} due to overlap with a segmentation.".format(cgroup+1))
                pass
        for segmentation in segmentations:
            r, c  = int(round(segmentation['y'])), int(round(segmentation['x']))
            s     = int(segmentation['sliceNo']) - 1
            group = int(pixels[s,r,c])
            if group <= 0:
                # This can happen for "concave" polygon shapes, where center falls into a "cutout" in the middle.
                objects = measurements.find_objects(pixels)
                was_placed = False
                for g_idx, box in enumerate(objects):
                    s_start, s_end = box[0].start, box[0].stop # objects are given by slices, which have start, stop attributes
                    ulr, lrr = box[1].start, box[1].stop       # one slice per dimension, so we have three: one for slices,
                    ulc, lrc = box[2].start, box[2].stop       # one for rows and one for cols
                    if s >= s_start and s < s_end and r >= ulr and r < lrr and c >= ulc and c < lrc:
                        was_placed = True
                        group = g_idx + 1
                        break
                # plt.imshow(pixels[s])
                # plt.show()
                if not was_placed:
                    print("WARN: polygon did not map to a group! p:{} x:{} y:{} s:{}".format(scan,c,r,s))
                    continue
            segmentation_groups[group-1].append(segmentation)
        # Now create candidate groups for all:
        candidate_groups = []
        for idx, segmentation in enumerate(segmentation_groups):
            if len(segmentation) > 0:
                new_group = {'patient':scan, 'members':segmentation}
                candidate_groups.append(set_group_stats(new_group))
            else:
                print("WARN: no segmentations in group {} for {}".format(idx+1, scan))
        for idx, centroid in enumerate(centroid_groups):
            if len(centroid) > 0:
                new_group = {'patient':scan, 'members':centroid}
                candidate_groups.append(set_group_stats(new_group))
        # Now add to the full "merged by scan" group, keyed by scan (dicomdir):
        merged_by_scan[scan] = candidate_groups
    return merged_by_scan

def write_candidates(centroid_groups_by_scan, output_name=None, combined=False, short=False):
    output_text = {}
    scans = centroid_groups_by_scan.keys()
    scans.sort()
    for scan in scans:
        dicom_dir = scan
        output_file_name = os.path.join(dicom_dir, output_name) if not combined else output_name
        print("output_file_name: {}".format(output_file_name))
        centroid_groups  = centroid_groups_by_scan[scan]
        centroid_groups.sort(key=lambda x: x['sliceNo'])
        print("centroids: {}".format(centroid_groups))
        for group in centroid_groups:
            # patient Z X Y is_nodule is_included malignancy
            output = [group['patient'], int(round(group['sliceNo'])), round(group['x'],1), round(group['y'],1), group['isNodule'], group['inclusion'], group['malignancy']]    
            if short:
                output = output[1:]
            output = '\t'.join([str(x) for x in output])
            if output_file_name in output_text:
                output_text[output_file_name].append(output)
            else:
                output_text[output_file_name] = [output]
    for output_file in output_text:
        print("writing output_file_name: {}".format(output_file))
        with open(output_file, 'wb') as fout:
            fout.write('\n'.join(output_text[output_file]) + '\n')

def run_standalone():
    args = get_args()
    if args['short'] and args['combine']:
        print("ERROR: Cannot use short output format with combined output.")
        sys.exit(1)

    centroids_by_scan       = collect_centroids_from_xml(args['dicom_dir'])
    
    if args['nodules_only']:
        centroids_by_scan   = remove_non_nodules(centroids_by_scan)

    centroid_groups_by_scan = merge_centroids(centroids_by_scan)
    write_candidates(centroid_groups_by_scan, output_name=args['output_file'], combined=args['combine'], short=args['short'])

if __name__ == "__main__":
    run_standalone()
