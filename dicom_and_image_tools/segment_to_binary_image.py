"""
Segment a selected nodule or nodules from annotation file(s) to 
produce binary segmentation image(s).

Run with --help to see usage info.
"""

from __future__ import print_function
import sys, pydicom as dicom, os, errno
import numpy as np
from PIL import Image, ImageDraw
import matplotlib.path as mplPath
from get_tumor_centroids import *
from get_tumor_polygons import *
from segment_tumors import *
from candidate_dict_utilities import *
from scipy.ndimage import measurements
from itertools import izip
import SimpleITK as sitk
import create_centroid_lists_from_LIDC_xml as xmlReader
from dicom_dir_utilities import *

import matplotlib.pyplot as plt


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def usage():
    eprint(__doc__.format(sys.argv[0]))


def construct_path(path):
    """
    construct a directory by creating necessary directories
    http://stackoverflow.com/a/600612
    """
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def is_polygon(candidate):
    has_coords = (
        "has_segmentation" in candidate and candidate["has_segmentation"] == True
    )
    if (not has_coords) and ("members" in candidate):
        for member in candidate["members"]:
            if "coords" in member and len(member["coords"]) > 0:
                has_coords = True
                break
    return has_coords or "coords" in candidate and len(candidate["coords"]) > 0


def contour_to_coords_list(contour):
    coord_list = []
    for pair in contour:
        coord_list.extend([pair[0], pair[1]])
    return coord_list


def contour_to_boundary(contour):
    return [{"x": pair[0], "y": pair[1]} for pair in contour]


def boundary_to_contour(boundary):
    return [(v["x"], v["y"]) for v in boundary]


def coords_to_contour(coords):
    return coords_to_polygon(coords)


def coords_to_polygon(coords):
    # it = iter(coords)
    # return izip(it,it)
    return zip(coords[::2], coords[1::2])


def point_is_inside_coords(point, coords):
    rx, ry = point
    return mplPath.Path(coords_to_polygon(coords)).contains_point((rx, ry))


def centroid_for_coords(coords):
    return centroid_for_polygon(boundary_to_polygon(coords_to_polygon(coords)))


def centroid_to_point(centroid):
    return (centroid["x"], centroid["y"])


def centroid_is_inside_polygon(centroid, polygon):
    return point_is_inside_coords(centroid_to_point(centroid), polygon["coords"])


def get_nifti_from_dicom_dir(dicom_dir, tmp_dir=None):
    from tempfile import mkstemp
    import dicom2nifti

    tmp_fd, nii_file = mkstemp(prefix="from_dicom_", suffix=".nii", dir=tmp_dir)
    os.close(tmp_fd)
    nifti = dicom2nifti.dicom_series_to_nifti(dicom_dir, nii_file, reorient_nifti=True)
    return (nii_file, sitk.ReadImage(nifti["NII_FILE"]))


def read_dicom_pixels(filename):
    df = dicom.read_file(filename)
    return get_dicom_pixels(df)


def read_dicom_blank(filename):
    df = dicom.read_file(filename)
    df.pixel_array[:] = np.zeros(df.pixel_array.shape, df.pixel_array.dtype)
    df.PixelData = df.pixel_array.tostring()
    return df


def save_dicom_image(filename, dicom_obj, reset_intercept=False):
    intercept = dicom_obj.RescaleIntercept
    if reset_intercept is not False:
        intercept = reset_intercept if type(reset_intercept) == type(0) else 0
    dicom_obj.RescaleIntercept = intercept
    dicom_obj.save_as(filename)


def get_dicom_pixels(dicom_obj, pixel_type="float"):
    pixels = dicom_obj.pixel_array
    if pixel_type != None:
        pixels = pixels.astype(pixel_type)
    return pixels


def set_dicom_pixels(dicom_obj, new_pixels):
    new_pixels = np.array(new_pixels, dtype=dicom_obj.pixel_array.dtype)
    dicom_obj.pixel_array[:] = new_pixels[:]
    dicom_obj.Rows, dicom_obj.Cols = dicom_obj.pixel_array.shape
    dicom_obj.PixelData = dicom_obj.pixel_array.tostring()
    return dicom_obj


def build_choices_list(polygons, centroids=[]):
    choices = []
    for poly in polygons:
        choices.append({"type": "polygon", "value": poly})
    for center in centroids:
        is_inside = False
        for pidx, poly in enumerate(polygons):
            rx, ry, rz = center["x"], center["y"], center["sliceNo"]
            if poly["sliceNo"] == rz and point_is_inside_coords(
                (rx, ry), poly["coords"]
            ):
                is_inside = True
                eprint(
                    "{0} is inside poly {1}".format(
                        [
                            round(p)
                            for p in (center["x"], center["y"], center["sliceNo"])
                        ],
                        pidx + 1,
                    )
                )
                break
        if not is_inside:
            center = centroid_to_polygon(center)
            choices.append({"type": "centroid", "value": center})
    choices.sort(key=lambda x: x["value"]["sliceNo"])
    return choices


def choose_candidate(polygons, centroids=[], auto_choose=None):
    if auto_choose is None:
        return manually_choose_candidate(polygons, centroids)
    strategy_parts = auto_choose.split(":")
    strategy_name = strategy_parts[0]
    strategy_args = strategy_parts[1:] if len(strategy_parts) > 1 else None
    print("Filter Strategy: {}".format(strategy_name))
    if strategy_name == "most-malignant":
        print(" - filtering by Most malignant")
        return auto_choose_candidate_by_most_malignant(polygons, centroids)
    elif strategy_name == "least-malignant":
        print(" - filtering by Least malignant")
        return auto_choose_candidate_by_most_malignant(polygons, centroids, -1)
    elif strategy_name == "malignancy-filter":
        threshold, comparison = strategy_args
        threshold = int(threshold)
        print(
            " - filtering by malignancy filter (threshold: {}, comparision: {})".format(
                threshold, comparison
            )
        )
        return auto_choose_candidate_filter_by_malignancy(
            polygons, centroids, threshold, comparison
        )
    elif strategy_name == "malignancy-range":
        low = None
        high = None
        sort_order = None
        seed = None
        if len(strategy_args) == 3:
            low, high, sort_order = strategy_args
        else:
            low, high, sort_order, seed = strategy_args
        low = int(low)
        high = int(high)
        if seed is not None:
            seed = int(seed)
        direction = -1 if sort_order == "ascending" else 1
        randomize = sort_order == "random"
        print(
            " - filtering by Malignancy range {} to {}, randomize is {} seed is {}.".format(
                low, high, randomize, seed
            )
        )
        return auto_choose_candidate_filter_by_malignancy_range(
            polygons, centroids, low, high, direction, randomize, seed
        )
    elif strategy_name == "largest":
        print(" - filtering by Largest")
        return auto_choose_candidate_by_area(polygons, centroids)
    elif strategy_name == "is-nodule":
        is_nodule_choice = 1
        if strategy_args is not None:
            if strategy_args[0].lower() == "false" or int(strategy_args[0]) == 0:
                is_nodule_choice = 0
        return auto_choose_candidates_by_is_nodule(
            polygons, centroids, is_nodule_choice
        )
    elif strategy_name == "random":
        seed = None
        if strategy_args is not None:
            seed = int(strategy_args)
        print(" - filtering by RANDOM (seed: {})".format(seed))
        return auto_choose_candidate_by_random(polygons, centroids, seed)
    elif strategy_name == "all":  # all candidates selected
        print(" - filtering by ALL")
        return auto_choose_candidate_all(polygons, centroids)
    else:
        raise RuntimeError("Illegal filter strategy: {}".format(strategy_name))


def auto_choose_candidate_all(polygons, centroids):
    # prefer polygons; combine poly and centroids so that any
    # centroids that are included inside a polygon are ignored
    chosen_list = copy.deepcopy(polygons)
    for centroid in centroids:
        # see if this centroid is inside any polygon:
        is_contained = False
        for polygon in polygons:
            if centroid_is_inside_polygon(centroid, polygon):
                is_contained = True
                break
        if not is_contained:
            chosen_list.append(centroid)
    chosen_list = sort_candidates_by_is_nodule(chosen_list)
    chosen_list = sort_candidates_by_inclusion(chosen_list)
    return chosen_list


def auto_choose_candidate_by_random(polygons, centroids, seed=None):
    # Same as "all", but shuffle the result:
    import random

    chosen_list = auto_choose_candidate_all(polygons, centroids)
    if seed is not None:
        random.seed(int(seed))
    random.shuffle(chosen_list)
    return chosen_list


def auto_choose_candidate_by_most_malignant(polygons, centroids, direction=1):
    combined_list = auto_choose_candidate_all(polygons, centroids)
    direction = True if direction == 1 else False
    # Sort list by malignancy:
    combined_list = sort_candidates_by_area(combined_list)
    combined_list = sort_candidates_by_malignancy(combined_list, direction)
    return [combined_list[0]]


def auto_choose_candidate_filter_by_malignancy(
    polygons, centroids, threshold=0, comparison=">"
):
    combined_list = auto_choose_candidate_all(polygons, centroids)
    combined_list = sort_candidates_by_area(combined_list)
    combined_list = sort_candidates_by_malignancy(combined_list, direction)
    return filter_candidates_by_malignancy(combined_list, threshold, comparison)


def auto_choose_candidate_filter_by_malignancy_range(
    polygons, centroids, low=1, high=5, direction=1, randomize=False, seed=None
):
    print(
        "Choosing by malignancy range {}-{} (random: {}, seed {})".format(
            low, high, randomize, seed
        )
    )
    if low > high:
        low, high = (high, low)
    if low > 0:
        centroids = [
            c for c in centroids if "malignancy" in c and c["malignancy"] > 0
        ]  # Only use centroids if they are marked for malignancy (original ones aren't)
    combined_list = auto_choose_candidate_all(polygons, centroids)
    combined_list = sort_candidates_by_area(combined_list)
    combined_list = sort_candidates_by_malignancy(combined_list, direction)
    combined_list = filter_candidates_by_malignancy_range(combined_list, low, high)
    if randomize:
        import random

        if seed is not None:
            random.seed(seed)
        random.shuffle(combined_list)
    return combined_list


def auto_choose_candidates_by_is_nodule(polygons, centroids, is_nodule_choice):
    combined_list = auto_choose_candidate_all(polygons, centroids)
    return filter_candidates_by_is_nodule(combined_list, is_nodule_choice)


def auto_choose_candidate_by_area(polygons, centroids):
    # only polygons matter here:
    sorted_list = sort_candidates_by_malignancy(polygons)
    sorted_list = sort_candidates_by_area(sorted_list)
    return [sorted_list[0]]


def sort_candidates_by_malignancy(candidates, descending=True):
    return sorted(candidates, key=lambda x: x["malignancy"], reverse=descending)


def sort_candidates_by_inclusion(candidates, descending=True):
    return sorted(candidates, key=lambda x: x["inclusion"], reverse=descending)


def sort_candidates_by_is_nodule(candidates, descending=True):
    return sorted(candidates, key=lambda x: x["isNodule"], reverse=descending)


def sort_candidates_by_area(candidates, descending=True):
    candidates = set_area_for_all_candidates(candidates)
    return sorted(candidates, key=lambda x: x["area"], reverse=descending)


def filter_candidates_by_malignancy(candidates, threshold, comparison=">"):
    """ 
        only return candidates with malignancy VS threshold 
        where VS is given by `comparison` in {>, < , >=, <=, ==} 
    """
    filter_func = (
        lambda x: x["malignancy"] > threshold
        if comparison == ">"
        else lambda x: x["malignancy"] < threshold
        if comparison == "<"
        else lambda x: x["malignancy"] >= threshold
        if comparison == ">="
        else lambda x: x["malignancy"] <= threshold
        if comparison == "<="
        else lambda x: x["malignancy"] == threshold
    )
    return filter(filter_func, candidates)


def filter_candidates_by_malignancy_range(candidates, low=1, high=5):
    """ 
        only return candidates with malignancy in the range [low, high]
        (inclusive on both ends)
    """
    print("Filtering {} \n{} {}".format(candidates, low, high))
    filter_func = lambda x: x["malignancy"] >= low and x["malignancy"] <= high
    return filter(filter_func, candidates)


def filter_candidates_by_area(candidates, threshold, comparison=1):
    """ 
        only return candidates with area VS threshold 
        where VS is given by `comparison` in {>, < , >=, <=, ==} 
    """
    candidates = set_area_for_all_candidates(candidates)
    filter_func = (
        lambda x: x["area"] > threshold
        if comparison == ">"
        else lambda x: x["area"] < threshold
        if comparison == "<"
        else lambda x: x["area"] >= threshold
        if comparison == ">="
        else lambda x: x["area"] <= threshold
        if comparison == "<="
        else lambda x: x["area"] == threshold
    )
    return filter(filter_func, candidates)


def filter_candidates_by_volume(candidates, threshold, comparison=1):
    """ 
        only return candidates with volume VS threshold 
        where VS is given by `comparison` in {>, < , >=, <=, ==} 
    """
    candidates = set_area_for_all_candidates(candidates)
    filter_func = (
        lambda x: x["volume"] > threshold
        if comparison == ">"
        else lambda x: x["volume"] < threshold
        if comparison == "<"
        else lambda x: x["volume"] >= threshold
        if comparison == ">="
        else lambda x: x["volume"] <= threshold
        if comparison == "<="
        else lambda x: x["volume"] == threshold
    )
    return filter(filter_func, candidates)


def filter_candidates_by_is_nodule(candidates, is_nodule):
    return filter(lambda x: int(x["isNodule"]) == is_nodule, candidates)


def filter_groups(candidate_groups, auto_choose=None):
    if auto_choose is None or auto_choose == False:
        auto_choose = "all"
    strategy_parts = auto_choose.split(":")
    strategy_name = strategy_parts[0]
    strategy_args = strategy_parts[1:] if len(strategy_parts) > 1 else None
    eprint("Group Filter Strategy: {}".format(strategy_name))
    candidate_groups = sort_candidates_by_is_nodule(candidate_groups)
    candidate_groups = sort_candidates_by_inclusion(candidate_groups)

    if strategy_name == "most-malignant":
        print(" - filtering groups by Most malignant")
        return filter_groups_by_most_malignant(candidate_groups)
    elif strategy_name == "least-malignant":
        print(" - filtering groups by Least malignant")
        return filter_groups_by_most_malignant(candidate_groups, -1)
    elif strategy_name == "malignancy-filter":
        threshold, comparison = strategy_args
        threshold = int(threshold)
        print(
            " - filtering groups by malignancy filter (threshold: {}, comparision: {})".format(
                threshold, comparison
            )
        )
        return filter_groups_by_malignancy(candidate_groups, threshold, comparison)
    elif strategy_name == "malignancy-range":
        low = None
        high = None
        sort_order = None
        seed = None
        if len(strategy_args) == 3:
            low, high, sort_order = strategy_args
        else:
            low, high, sort_order, seed = strategy_args
        low = int(low)
        high = int(high)
        if seed is not None:
            seed = int(seed)
        direction = -1 if sort_order == "ascending" else 1
        randomize = sort_order == "random"
        print(
            " - filtering groups by Malignancy range {} to {}, randomize is {} seed is {}.".format(
                low, high, randomize, seed
            )
        )
        return filter_groups_filter_by_malignancy_range(
            candidate_groups, low, high, direction, randomize, seed
        )
    elif strategy_name == "largest":
        print(" - filtering groups by Largest")
        return filter_groups_by_size(candidate_groups)
    elif strategy_name == "is-nodule":
        is_nodule_choice = 1
        if strategy_args is not None:
            if strategy_args[0].lower() == "false" or int(strategy_args[0]) == 0:
                is_nodule_choice = 0
        print(" - filtering groups by is-nodule")
        return filter_groups_by_is_nodule(candidate_groups, is_nodule_choice)
    elif strategy_name == "random":
        seed = None
        if strategy_args is not None:
            seed = int(strategy_args)
        print(" - filtering groups by RANDOM (seed: {})".format(seed))
        return filter_groups_by_random(candidate_groups, seed)
    elif strategy_name == "all":  # all candidates selected
        print(" - filtering groups by ALL")
        return candidate_groups
    else:
        raise RuntimeError("Illegal group filter strategy: {}".format(strategy_name))


def filter_groups_by_most_malignant(candidate_groups, direction=1):
    direction = True if direction == 1 else False
    # Sort list by malignancy:
    combined_list = sort_groups_by_size(candidate_groups)
    combined_list = sort_groups_by_malignancy(combined_list, direction)
    return [combined_list[0]]


def filter_groups_filter_by_malignancy(candidate_groups, threshold=0, comparison=">"):
    return filter_candidates_by_malignancy(candidate_groups, threshold, comparison)


def filter_groups_filter_by_malignancy_range(
    candidate_groups, low=1, high=5, direction=1, randomize=False, seed=None
):
    print(
        "Choosing by malignancy range {}-{} (random: {}, seed {})".format(
            low, high, randomize, seed
        )
    )
    if low > high:
        low, high = (high, low)
    if low > 0:
        candidate_groups = [
            c for c in candidate_groups if "malignancy" in c and c["malignancy"] > 0
        ]  # Only use centroids if they are marked for malignancy (original ones aren't)
    candidate_groups = sort_candidates_by_area(candidate_groups)
    candidate_groups = sort_candidates_by_malignancy(candidate_groups, direction)
    print("About to filter!!!!")
    filtered_list = filter_candidates_by_malignancy_range(candidate_groups, low, high)
    if randomize:
        import random

        if seed is not None:
            random.seed(seed)
        random.shuffle(filtered_list)
    return filtered_list


def filter_groups_by_area(candidate_groups):
    sorted_list = sort_groups_by_malignancy(candidate_groups)
    sorted_list = sort_groups_by_size(sorted_list)
    return [sorted_list[0]]


def filter_groups_by_random(candidate_groups, seed=None):
    # Same as "all", but shuffle the result:
    import random

    if seed is not None:
        random.seed(int(seed))
    random.shuffle(candidate_groups)
    return candidate_groups


def filter_groups_by_is_nodule(candidate_groups, is_nodule):
    return filter_candidates_by_is_nodule(candidate_groups, is_nodule)


def sort_groups_by_malignancy(candidate_groups, descending=True):
    return sorted(candidate_groups, key=lambda x: x["malignancy"], reverse=descending)


def sort_groups_by_size(candidate_groups, descending=True):
    return sorted(candidate_groups, key=lambda x: x["volume"], reverse=descending)


def manually_choose_candidate(polygons, centroids=[]):
    choices = build_choices_list(polygons, centroids)
    done = len(choices) == 0
    refresh = True
    chosen = None
    while not done:
        if refresh:
            print("File: {0}".format(choices[0]["value"]["fileName"]))
            for idx, choice in enumerate(choices):
                value = choice["value"]
                description = "mal: {0}, x: {1}, y: {2}, slice: {3}, area: {5}, id: {4}".format(
                    value["malignancy"],
                    round(value["x"]),
                    round(value["y"]),
                    value["sliceNo"],
                    value["noduleID"],
                    value["area"],
                )
                print("[{0}]: {1} {2}".format(idx + 1, choice["type"], description))
            refresh = False
        v = raw_input(
            "Type # to choose, -# to remove, R to reload list, or X to exit: "
        )
        if v.lower() == "r":
            refresh = True
        elif v.lower() == "x":
            sys.exit(0)
        else:
            v = int(v)
            if v > 0:
                v = v - 1
                chosen = v
                done = True
            else:
                v = (-v) - 1
                del choices[v]
    return choices[chosen]


def set_area_for_all_candidates(candidates):
    for idx, polygon in enumerate(candidates):
        if not "area" in polygon:
            polygon["area"] = (
                area_for_polygon(polygon["coords"]) if "coords" in polygon else 0
            )
        candidates[idx] = polygon
    return candidates


def centroid_to_polygon(centroid, max_area=None):
    polygon = centroid
    if "coords" in polygon and len(polygon["coords"]) > 0:
        return polygon
    (x, y) = (polygon["x"], polygon["y"])
    image = read_dicom_pixels(centroid["fileName"])
    shape = image.shape
    if max_area == None:
        max_area = (shape[0] / 4) * (shape[1] / 4)
    max_diameter = int(round(np.sqrt(max_area)))
    candidates = find_candidate_nodules(
        image, seed_point=(x, y), seed_area=True, max_diameter=max_diameter
    )
    # boundary = trace_contour(region_grow(threshold(image), (x,y), mask=True, max_area=max_area), (x,y))
    print("boundary len: {0}".format(len(candidates)))
    if len(candidates) == 0:
        candidates.append(centroid)
        candidates[0]["boundary"] = [{x: x, y: y}]
        candidates[0]["area"] = 0
    boundary = candidates[0]["boundary"]
    coords = []
    for (c, r) in boundary:
        coords.append({"x": float(c), "y": float(r)})
    polygon["area"] = area_for_polygon(coords)
    polygon["coords"] = coords_to_list(coords)
    return polygon


# prefer polygons; combine poly and centroids so that any
# centroids that are included inside a polygon are ignored
# and overlapping polygons are handled by `strategy`
# (vote, seconded, union, intersect)


def get_3d_boundaries(chosen, polygons, centroids=[], area_tolerance=1):
    boundaries = []
    choices = build_choices_list(polygons, centroids)

    # If chosen is a centroid, we need to grow polygons for all centroids in the choices set first.
    for entry in choices:
        if entry["type"] == "centroid":
            entry["value"] = centroid_to_polygon(entry["value"])
            # entry['type']  = 'polygon'

    chosen_boundary_coords = chosen["value"]["coords"]

    # Any candidates whose centroid is inside the chosen boundary and whose
    # slice numbers are adjacent will be included.
    chosen_slice = chosen["value"]["sliceNo"]
    eprint(chosen_slice)
    by_slice = [(c["value"]["sliceNo"], idx) for idx, c in enumerate(choices)]
    by_slice.sort()
    adjacent_slices = [chosen]
    mid_slice = 0
    for idx, s in enumerate(by_slice):
        if chosen_slice == s[0]:
            mid_slice = idx
            break
    eprint(mid_slice)
    j = mid_slice - 1
    last_slice = chosen_slice
    last_boundary = chosen_boundary_coords
    area_epsilon = 1e-20
    last_area = chosen["value"]["area"] + area_epsilon
    while (
        j >= 0 and (last_slice - by_slice[j][0]) <= 1
    ):  # moving to smaller-numbered slices
        eprint("dn", j, by_slice[j][0], "last", last_slice)
        candidate = choices[by_slice[j][1]]
        center = (candidate["value"]["x"], candidate["value"]["y"])
        area = candidate["value"]["area"] + area_epsilon
        eprint(
            chosen["value"]["y"],
            chosen["value"]["x"],
            center,
            (abs(area - last_area) / float(last_area)) <= area_tolerance,
            (abs(area - last_area) / float(last_area)),
        )
        if (
            point_is_inside_coords(center, last_boundary)
            and (abs(area - last_area) / float(last_area)) <= area_tolerance
        ):
            eprint("appending")
            adjacent_slices.append(candidate)
            last_slice = by_slice[j][0]
            last_boundary = candidate["value"]["coords"]
            last_area = candidate["value"]["area"] + area_epsilon
        j -= 1
    j = mid_slice + 1
    last_slice = chosen_slice
    last_boundary = chosen_boundary_coords
    last_area = chosen["value"]["area"] + area_epsilon
    while (
        j < len(by_slice) and (by_slice[j][0] - last_slice) <= 1
    ):  # moving to larger-numbered slices
        eprint("up", j, by_slice[j][0], "last", last_slice)
        candidate = choices[by_slice[j][1]]
        center = (candidate["value"]["x"], candidate["value"]["y"])
        area = candidate["value"]["area"] + area_epsilon
        eprint(
            chosen["value"]["y"],
            chosen["value"]["x"],
            center,
            (abs(area - last_area) / float(last_area)) <= area_tolerance,
            (abs(area - last_area) / float(last_area)),
        )
        if (
            point_is_inside_coords(center, last_boundary)
            and (abs(area - last_area) / float(last_area)) <= area_tolerance
        ):
            eprint("appending")
            adjacent_slices.append(candidate)
            last_slice = by_slice[j][0]
            last_boundary = candidate["value"]["coords"]
            last_area = candidate["value"]["area"] + area_epsilon
        j += 1
    adjacent_slices.sort(
        key=lambda x: x["value"]["sliceNo"] + 0.000001 * x["value"]["malignancy"]
    )
    for slice in adjacent_slices:
        eprint(
            "Selected slice {0} x: {1} y: {2}".format(
                slice["value"]["sliceNo"], slice["value"]["x"], slice["value"]["y"]
            )
        )
    return adjacent_slices


def produce_3d_binary_image(dicomdir, outputdir, slice_boundaries, shape=(512, 512)):
    files = []
    if os.path.abspath(outputdir) == os.path.abspath(dicomdir):
        raise IOError(
            "Cannot output segmented images to same directory as original images."
        )
    if os.path.isdir(dicomdir):
        for rt, dr, fl in os.walk(dicomdir):
            files.extend(fl)
            break
    else:
        eprint("dicomdir = {0} is not a directory.".format(dicomdir))
        raise IOError("dicomdir = {0} is not a directory.".format(dicomdir))

    dicomfiles = [f for f in files if f[-4:] == ".dcm"]
    construct_path(outputdir)

    included_slices_by_name = {}

    center_x = center_y = center_z = center_slice = None
    for idx, slice in enumerate(slice_boundaries):
        fname = slice["value"]["fileName"]
        center_x = (
            slice["value"]["x"] if center_x == None else center_x + slice["value"]["x"]
        )
        center_y = (
            slice["value"]["y"] if center_y == None else center_y + slice["value"]["y"]
        )
        center_z = (
            slice["value"]["z"] if center_z == None else center_z + slice["value"]["z"]
        )
        center_slice = (
            slice["value"]["sliceNo"]
            if center_slice == None
            else center_slice + slice["value"]["sliceNo"]
        )
        if fname in included_slices_by_name:
            included_slices_by_name[fname].append(slice)
        else:
            included_slices_by_name[fname] = [slice]
    center_x, center_y, center_z, center_slice = [
        int(round(float(v) / len(slice_boundaries)))
        for v in [center_x, center_y, center_z, center_slice]
    ]
    slice_boundaries = None

    for file in dicomfiles:
        filepath = os.path.join(dicomdir, file)
        out_filepath = os.path.join(outputdir, os.path.basename(file))
        df = read_dicom_blank(filepath)
        if filepath in included_slices_by_name:
            pixels = get_dicom_pixels(df)
            # Union all segmentations
            for slice in included_slices_by_name[filepath]:
                pixels = draw_polygon_from_coords(pixels, slice["value"]["coords"])
            df = set_dicom_pixels(df, pixels)
        save_dicom_image(out_filepath, df)
    info_filename = os.path.join(outputdir, "segment_info.txt")
    with open(info_filename, "w") as info_fp:
        info_fp.write("Segmentation Centroid x, y, z, sliceNo, totalSlices\n")
        info_fp.write(
            "{0}, {1}, {2}, {3}, {4}\n".format(
                center_x, center_y, center_z, center_slice, len(dicomfiles)
            )
        )


def get_vote_threshold(consensus="vote"):
    """
    Valid consensus strategies are: 
    {'vote' (default), 'vote:PCT', 'union', 'intersect'}
        'vote'      : any pixel within 50% or more of boundaries is included
        'vote:PCT'  : any pixel within PCT (integer in range [0,100]) of boundaries is included
        'union'     : same as vote:0 - any pixel within at least one boundary is included
        'intersect' : same as 'vote:100' - only pixels within all boundaries are included
    """
    vote_threshold = 0.50
    if consensus[0:4] == "vote" and len(consensus) > 4:
        vote_threshold = int(consensus[5:]) / 100.0
        if vote_threshold == 0:
            vote_threshold += 1e-10
            # epsilon prevents perfect zero to require at least one vote
    elif consensus == "union":
        vote_threshold = 1e-10  # epsilon above zero to require at least one vote
    elif consensus == "intersect":
        vote_threshold = 1.0
    return vote_threshold


def write_3d_binary_image(
    dicomdir, outputdir, segmentation_group, shape=(512, 512), consensus="vote"
):
    """
    Valid consensus strategies are: 
    {'vote' (default), 'vote:PCT', 'union', 'intersect'}
        'vote'      : any pixel within 50% or more of boundaries is included
        'vote:PCT'  : any pixel within PCT (integer in range [0,100]) of boundaries is included
        'union'     : same as vote:0 - any pixel within at least one boundary is included
        'intersect' : same as 'vote:100' - only pixels within all boundaries are included
    """
    vote_threshold = get_vote_threshold(consensus)

    files = []
    if os.path.abspath(outputdir) == os.path.abspath(dicomdir):
        raise IOError(
            "Cannot output segmented images to same directory as original images."
        )
    if os.path.isdir(dicomdir):
        for rt, dr, fl in os.walk(dicomdir):
            files.extend(fl)
            break
    else:
        eprint("dicomdir = {0} is not a directory.".format(dicomdir))
        raise IOError("dicomdir = {0} is not a directory.".format(dicomdir))

    dicomfiles = [f for f in files if f[-4:] == ".dcm"]
    construct_path(outputdir)

    included_slices_by_name = {}

    center_x = segmentation_group["x"]
    center_y = segmentation_group["y"]
    center_z = segmentation_group["z"]
    center_slice = segmentation_group["sliceNo"]
    center_malignancy = segmentation_group["malignancy"]

    dicomfiles_by_slice = {}
    for filename in dicomfiles:
        df = dicom.read_file(os.path.join(dicomdir, filename))
        slice_no = df.InstanceNumber
        dicomfiles_by_slice[slice_no] = filename

    for idx, slice in enumerate(segmentation_group["members"]):
        if "coords" not in slice:
            print("Missing coords for slice; dropping: {}".format(slice))
            continue
        if "fileName" not in slice:
            slice["fileName"] = os.path.join(
                dicomdir, dicomfiles_by_slice[slice["sliceNo"]]
            )
            print(
                "Filling in file name {} for slice {}.".format(
                    slice["fileName"], slice["sliceNo"]
                )
            )
        fname = slice["fileName"]

        if fname in included_slices_by_name:
            included_slices_by_name[fname].append(slice)
        else:
            included_slices_by_name[fname] = [slice]
    slice_boundaries = None

    for file in dicomfiles:
        filepath = os.path.join(dicomdir, file)
        out_filepath = os.path.join(outputdir, os.path.basename(file))
        df = read_dicom_blank(filepath)
        pixels = get_dicom_pixels(df)
        if filepath in included_slices_by_name:
            # Union all segmentations that overlap within the same group
            # (this is rare, but could happen if you are slicing through lobes of
            #  a larger 3D object)
            n_objects_on_slice = len(included_slices_by_name[filepath])
            for s_idx, slice in enumerate(included_slices_by_name[filepath]):
                print(
                    "Drawing pixels for object {} on slice {}".format(
                        s_idx + 1, slice["sliceNo"]
                    )
                )
                pixels = np.array(pixels) + draw_polygon_from_coords(
                    pixels, slice["coords"]
                )
            # Determine absolute (#-of-votes) threshold from vote_threshold (percent):
            # NOTE: If the number of nodules merged is greater than the maximum of the
            #       segmentation matrix, it means that at least one of the merges is
            #       disjoint from one or more of the others (this can happen when muliple lobes of
            #       a large 3D nodule are sliced so that the portions in the current slice are
            #       disjoint).  In these cases, the voting should be based on the maximum
            #       agreement seen.
            if n_objects_on_slice > pixels.max():
                eprint(
                    "WARNING: Multiple segmentations on same slice detected; slice {}".format(
                        included_slices_by_name[filepath][0]["sliceNo"]
                    )
                )
            threshold = vote_threshold * min(n_objects_on_slice, pixels.max())
            pixels[pixels < threshold] = 0
            pixels[pixels != 0] = 1

            df = set_dicom_pixels(df, pixels)
        save_dicom_image(out_filepath, df, 0)
    info_filename = os.path.join(outputdir, "segment_info.txt")
    with open(info_filename, "w") as info_fp:
        info_fp.write(
            "Segmentation Centroid: malignancy, x, y, z, sliceNo, totalSlices\n"
        )
        info_fp.write(
            "{0}, {1}, {2}, {3}, {4}, {5}\n".format(
                center_malignancy,
                center_x,
                center_y,
                center_z,
                center_slice,
                len(dicomfiles),
            )
        )


def draw_polygon_from_coords(pixels, coords, outline=1):
    if outline == True:
        outline = 1
    pixels = np.array(pixels)
    print("Drawing polygon {0}".format(coords))
    boundary = coords_to_polygon(coords)
    (rows, cols) = pixels.shape
    img = Image.new("L", (cols, rows), 0)
    ImageDraw.Draw(img).polygon(boundary, outline=outline, fill=1)
    pixels[np.array(img) == 1] = 1
    return pixels


def draw_circle(pixels, center, radius=3, outline=1):
    if outline == True:
        outline = 1
    pixels = np.array(pixels)
    x, y, r = center[0], center[1], radius - 0.5
    print("Drawing circle {0} {1} {2}".format(x, y, r))
    (rows, cols) = pixels.shape
    img = Image.new("L", (cols, rows), 0)
    ImageDraw.Draw(img).ellipse((x - r, y - r, x + r, y + r), outline=outline, fill=1)
    pixels[np.array(img) == 1] = 1
    print("Circle area: {0}".format(pixels.sum()))
    return pixels


def read_coords_file(coords_file):
    centroids = []
    with open(coords_file, "rU") as fin:
        for line in fin:
            if line.strip() != "":
                parts = line.split()
                if len(parts) == 7:
                    patient, slice_no, x, y, is_nodule, inclusion, malignancy = parts
                    centroids.append(
                        {
                            "patient": patient,
                            "sliceNo": int(round(float(slice_no))),
                            "x": float(x),
                            "y": float(y),
                            "z": float(slice_no),
                            "isNodule": int(is_nodule),
                            "inclusion": int(inclusion),
                            "malignancy": int(malignancy),
                            "noduleID": make_unique_id_string(
                                int(round(float(x))),
                                int(round(float(y))),
                                int(slice_no),
                                is_nodule,
                                inclusion,
                                malignancy,
                                patient,
                            ),
                        }
                    )
                else:
                    slice_no, x, y, is_nodule, inclusion, malignancy = parts
                    centroids.append(
                        {
                            "sliceNo": int(round(float(slice_no))),
                            "x": float(x),
                            "y": float(y),
                            "z": float(slice_no),
                            "isNodule": int(is_nodule),
                            "inclusion": int(inclusion),
                            "malignancy": int(malignancy),
                            "noduleID": make_unique_id_string(
                                int(round(float(x))),
                                int(round(float(y))),
                                int(slice_no),
                                is_nodule,
                                inclusion,
                                malignancy,
                            ),
                        }
                    )
    return centroids


def get_polygon_from_binary_img(img):
    nz = np.nonzero(img)
    seed_y = nz[0][0]
    seed_x = nz[1][0]
    boundary_contour = trace_contour(img, (seed_x, seed_y))
    boundary = contour_to_boundary(boundary_contour)
    coords = contour_to_coords_list(boundary_contour)
    return boundary, coords


def get_binary_segmentation_2d(candidate, img_shape=(512, 512), min_radius=3):
    c = candidate
    slice_bin = np.zeros(img_shape)
    slice_bin = (
        draw_polygon_from_coords(slice_bin, c["coords"])
        if "coords" in c
        else draw_circle(
            slice_bin, (int(round(c["x"])), int(round(c["y"]))), min_radius
        )
    )
    return slice_bin


def make_polygon_dict_given_boundary(boundary=None, coords=None):
    if boundary is None and coords is None:
        raise RuntimeError("either boundary or coords must be specified")
    if coords is None:
        coords = contour_to_coords_list(boundary_to_contour(boundary))
    if boundary is None:
        boundary = contour_to_boundary(coords_to_contour(coords))
    centroid = centroid_for_coords(coords)
    candidate = {
        "x": centroid["x"],
        "y": centroid["y"],
        "coords": coords,
        "boundary": boundary,
        "type": "polygon",
        "has_segmentation": True,
        "area": area_for_polygon(boundary),
        # we can't set the following meaningfully, but create the keys:
        "z": 0,
        "sliceNo": 0,
        "malignancy": -1,
        "inclusion": -1,
        "isNodule": -1,
        "filepath": None,
    }
    candidate["noduleID"] = nodule_unique_id(candidate)
    return candidate


def consolidate_candidates(
    candidates,
    consensus="vote",
    min_radius=3,
    img_shape=(512, 512),
    keep_bin_image=False,
):
    """
    Consolidate any candidates (from 2d slices) that are part of the same 
    "nodule-like-object" into a single list representing a 3-d candidate.
    
    In 2D:
    For candidates that have boundaries defined, the consolidation means that all other
    candidates whose centroid falls within the boundary are part of the same "object" and
    should be combined.  When multiple candidates exist that all have boundaries and which
    overlap, the new candidate's boundary will be determined by a consensus strategy according
    to the `consensus` argument.  Valid consensus strategies are: 
    {'vote' (default), 'vote:PCT', 'union', 'intersect'}
        'vote'      : any pixel within 50% or more of boundaries is included
        'vote:PCT'  : any pixel within PCT (integer in range [0,100]) of boundaries is included
        'union'     : same as vote:0 - any pixel within at least one boundary is included
        'intersect' : same as 'vote:100' - only pixels within all boundaries are included
    For candidates without boundaries defined 2D consolidation occurs by first plotting a circular
    area with radius=`radius`, then applying the consolidation according to the boundaries. 
    Note: If a nodule without a boundary is consolidated into a nodule with a boundary, the 
          boundary is considered the 'truth' and the sans-boundary nodule will not contribute to
          votes or the resulting boundary shape.

    In 3D:
    After 2D consolidation, candidates with overlapping boundaries in consecutive slices are 
    merged to create 3D candidate objects.
    """
    candidates = copy.deepcopy(candidates)
    candidates = consolidate_candidates_2d(
        candidates, consensus, min_radius, img_shape, keep_bin_image=True
    )
    candidates = consolidate_candidates_3d(
        candidates, consensus, min_radius, img_shape, keep_bin_image=keep_bin_image
    )
    return candidates


def consolidate_candidates_2d(
    candidates,
    consensus="vote",
    min_radius=10,
    img_shape=(512, 512),
    keep_bin_image=False,
):
    """
    Consolidate any candidates (from 2d slices) that are part of the same 
    "nodule-like-object" into a single list representing a 3-d candidate.
    
    In 2D:
    For candidates that have boundaries defined, the consolidation means that all other
    candidates whose centroid falls within the boundary are part of the same "object" and
    should be combined.  When multiple candidates exist that all have boundaries and which
    overlap, the new candidate's boundary will be determined by a consensus strategy according
    to the `consensus` argument.  Valid consensus strategies are: 
    {'vote' (default), 'vote:PCT', 'union', 'intersect'}
        'vote'      : any pixel within 50% or more of boundaries is included
        'vote:PCT'  : any pixel within PCT (integer in range [0,100]) of boundaries is included
        'union'     : same as vote:0 - any pixel within at least one boundary is included
        'intersect' : same as 'vote:100' - only pixels within all boundaries are included
    For candidates without boundaries defined 2D consolidation occurs by first plotting a circular
    area with radius=`radius`, then applying the consolidation according to the boundaries. 
    Note: If a nodule without a boundary is consolidated into a nodule with a boundary, the 
          boundary is considered the 'truth' and the sans-boundary nodule will not contribute to
          votes or the resulting boundary shape.
    """
    vote_threshold = 0.50
    if consensus[0:4] == "vote" and len(consensus) > 4:
        vote_threshold = int(consensus[5:]) / 100.0
        if vote_threshold == 0:
            vote_threshold += 1e-10
            # epsilon prevents perfect zero to require at least one vote
    elif consensus == "union":
        vote_threshold = 1e-10  # epsilon above zero to require at least one vote
    elif consensus == "intersect":
        vote_threshold = 1.0
    min_window = max(min_radius ** 2, 64)
    print("Consensus method: {} ; vote threshold {}".format(consensus, vote_threshold))

    # First order all the candidates by slice:
    candidates = sorted(candidates, key=lambda x: x["sliceNo"])
    candidates_by_slice = {}
    merged_by_slice = {}
    for candidate in candidates:
        if candidate["sliceNo"] in candidates_by_slice:
            candidates_by_slice[candidate["sliceNo"]].append(candidate)
        else:
            candidates_by_slice[candidate["sliceNo"]] = [candidate]
    for s in candidates_by_slice:
        # print(s) # TODO: remove
        slice_candidates = candidates_by_slice[s]
        # print(slice_candidates) # TODO: remove
        slice_bin = np.zeros(img_shape)
        # For efficiency, sort by area (ascending, we will pop from end):
        slice_candidates = sorted(
            slice_candidates, key=lambda c: c["area"] if "area" in c else 0
        )
        unmatched = copy.deepcopy(slice_candidates)
        merged = []
        while len(unmatched) > 0:
            c = unmatched.pop()
            found_match = False
            slice_bin = get_binary_segmentation_2d(c, img_shape, min_radius)
            print("slice bin: \n{0}".format(slice_bin))
            print(slice_bin[slice_bin != 0])
            for idx, m in enumerate(merged):
                merged_bin = np.array(m["seg_map"])
                if (
                    (merged_bin != 0) & (slice_bin != 0)
                ).any():  # If any pixel overlaps...
                    # To be the "same nodule", the centroid for the new one must lie inside the
                    # boundary of the merged nodule.
                    center_c, center_r = centroid_to_point(c)
                    center_r = int(round(center_r))
                    center_c = int(round(center_c))
                    print(
                        "Check {0} {1} inside? {2}".format(
                            center_c, center_r, merged_bin[center_r, center_c]
                        )
                    )
                    if merged_bin[center_r, center_c] != 0:
                        print("merging...")
                        print(merged_bin[merged_bin != 0])
                        # Ok, the new one is part of the existing merged one... It will contribute to the
                        # vote if it is a bounded polygon or if no merged point is a bounded polygon.
                        # Otherwise if it is a point and the merged polygon has one or more segmentations, we
                        # just quietly "merge" it (really, just ignored from this point on).
                        found_match = True
                        if is_polygon(c) and (not merged[idx]["has_segmentation"]):
                            print("New point is a polygon; we didn't have one before.")
                            merged[idx]["seg_map"][
                                :
                            ] = 0  # This will be the first "real" segmentation info.
                            merged[idx]["merged_count"] = 0
                        if not merged[idx]["has_segmentation"]:
                            merged[idx]["has_segmentation"] = is_polygon(c)
                        # merge the binary images:
                        if is_polygon(c) or (not merged[idx]["has_segmentation"]):
                            merged[idx]["seg_map"] = merged[idx]["seg_map"] + slice_bin
                            merged[idx]["merged_count"] += 1
                            print(
                                "Adding new segmentation to existing... seg idx {} count now {}".format(
                                    idx, merged[idx]["merged_count"]
                                )
                            )

            if not found_match:
                # If we couldn't merge the new candidate, it becomes a new 'merged' nodule on its own:
                c[
                    "seg_map"
                ] = (
                    slice_bin
                )  # Note: new candidates created from points only will contain the default radius.
                c["merged_count"] = 1
                c["has_segmentation"] = is_polygon(c)
                merged.append(c)

        # At this point we have merged everything; now determine the final segmentation shapes by consensus
        merged_result = []
        for idx, m in enumerate(merged):
            slice_bin = m["seg_map"]
            # Determine absolute (#-of-votes) threshold from vote_threshold (percent):
            # NOTE: If the number of nodules merged is greater than the maximum of the
            #       segmentation matrix, it means that at least one of the merges is
            #       disjoint from one or more of the others (this can happen when muliple lobes of
            #       a large 3D nodule are sliced so that the portions in the current slice are
            #       disjoint).  In these cases, the voting should be based on the maximum
            #       agreement seen.
            threshold = vote_threshold * min(m["merged_count"], slice_bin.max())
            print(
                "Need {} votes to keep voxel. (Merged {} candidates.)".format(
                    threshold, m["merged_count"]
                )
            )
            new_slice_bin = np.array(slice_bin)
            new_slice_bin[new_slice_bin < threshold] = 0
            new_slice_bin[new_slice_bin != 0] = 1
            c_of_m = measurements.center_of_mass(new_slice_bin)
            if np.isnan(c_of_m).any():  # TODO: REMOVE
                print(
                    "NaN detected. Slice {}\n{}".format(
                        m["sliceNo"], slice_bin[slice_bin != 0]
                    )
                )
                import scipy.misc

                scipy.misc.imsave("troublesome_segmentation.png", slice_bin)
                sys.exit(0)
            print("CENTER OF MASS: {}".format(c_of_m))
            merged[idx]["x"] = round(c_of_m[0], 1)
            merged[idx]["y"] = round(c_of_m[1], 1)
            merged[idx]["area"] = new_slice_bin.sum() if m["has_segmentation"] else 0
            merged[idx]["noduleID"] = nodule_unique_id(merged[idx])
            if m["has_segmentation"]:
                cx, cy = int(round(c_of_m[1])), int(round(c_of_m[0]))
                # NOTE: In cases where there were disjoint objects, we have to produce separate
                #       segmentation polygons for each one.
                img_labels, num_objects = measurements.label(new_slice_bin)
                if num_objects > 1:
                    merged[idx]["members"] = []  # This will actually become a group
                    merged[idx]["boundary"] = []
                    merged[idx]["coords"] = []
                for o_i in range(num_objects):
                    print(
                        "Creating polygon for object {} of {}...".format(
                            o_i + 1, num_objects
                        )
                    )
                    isolated_object = np.array(new_slice_bin)
                    isolated_object[img_labels != (o_i + 1)] = 0
                    isolated_object[isolated_object != 0] = 1
                    boundary, coords = get_polygon_from_binary_img(isolated_object)
                    if num_objects == 1:
                        merged[idx]["boundary"] = boundary
                        merged[idx]["coords"] = coords
                        merged_result.append(merged[idx])
                    else:
                        split_segment = copy.deepcopy(merged[idx])
                        split_segment["boundary"] = boundary
                        split_segment["coords"] = coords
                        split_segment["area"] = isolated_object.sum()
                        merged_result.append(split_segment)
                    print("NEW contour: {}".format(boundary))
                    print("NEW coords:  {}".format(coords))
            else:
                for key in ["boundary", "coords"]:
                    try:
                        del merged[idx][key]
                    except:
                        pass
                merged_result.append(merged[idx])
        merged = merged_result
        # One more pass for cleanup if requested:
        if not keep_bin_image:
            for idx, m in enumerate(merged):
                for key in [
                    "seg_map",
                    "merged_count",
                    "has_segmentation",
                    "seg_map_ul_x",
                    "seg_map_ul_y",
                ]:
                    try:
                        del merged[idx][key]
                    except:
                        pass
        merged_by_slice[s] = merged
    merged = []
    for s in merged_by_slice:
        merged.extend(merged_by_slice[s])
    return merged


def consolidate_candidates_3d(
    candidates,
    consensus="vote",
    min_radius=3,
    img_shape=(512, 512),
    keep_bin_image=False,
):
    """
    Consolidate any candidates (from 2d slices) that are part of the same 
    "nodule-like-object" into a single list representing a 3-d candidate.

    In 3D:
    After 2D consolidation, candidates in contiguous slices whose segmentation area 
    overlaps are merged to create 3D candidate objects.
    For candidates with no segmentation (seed points), a circular projection with 
    radius=`radius` is used as the segmentation.
    """
    merged = []
    merged_by_slice = {}

    print(candidates)
    for candidate in candidates:
        if not "seg_map" in candidate:
            candidate["seg_map"] = get_binary_segmentation_2d(
                candidate, img_shape, min_radius
            )
        candidate["seg_hash"] = img_hash(candidate["seg_map"])

    candidates = sorted(candidates, key=lambda x: x["sliceNo"])
    candidates_by_slice = {}
    for candidate in candidates:
        if candidate["sliceNo"] in candidates_by_slice:
            candidates_by_slice[candidate["sliceNo"]].append(candidate)
        else:
            candidates_by_slice[candidate["sliceNo"]] = [candidate]
    slices_with_candidates = sorted([s for s in candidates_by_slice])
    first = int(round(slices_with_candidates[0]))
    last = int(round(slices_with_candidates[-1]))
    print(first, last)

    merged = []  # a 3d merged group consists of a segmentation image and list of (2d) nodules included
    # in the form {'seg_map': <most-recent-binary-image>, 'members': <list-of-candidates>}
    for candidate in candidates_by_slice[first]:
        merged.append(
            {
                "seg_map": np.array(candidate["seg_map"]),
                "last_slice": first,
                "members": [candidate],
                "member_hashes": [candidate["seg_hash"]],
                "has_segmentation": is_polygon(candidate),
            }
        )
    # For each remaining slice add/merge candidates:
    for s in range(first + 1, last + 1):
        if s in candidates_by_slice:
            candidates = candidates_by_slice[s]
            was_merged = [False for candidate in candidates]
            for group in merged:
                if group["last_slice"] == (s - 1):
                    g_map = group["seg_map"]
                    new_map = np.zeros(img_shape)
                    for idx, candidate in enumerate(candidates):
                        c_map = candidate["seg_map"]
                        if ((c_map != 0) & (g_map != 0)).any():
                            # There is overlap.  Merge; true merge if group is poly and so is candidate or if group is all points
                            # and so is candidate; weak merge (drop candidate) if group is poly and candidate is not; restart
                            # group if candidate is poly and group is all points.
                            was_merged[idx] = True
                            if (group["has_segmentation"] == is_polygon(candidate)) or (
                                not group["has_segmentation"]
                            ):
                                new_map = new_map.astype(bool) | c_map.astype(bool)
                                group["members"].append(candidate)
                                group["member_hashes"].append(candidate["seg_hash"])
                                group["last_slice"] = s
                            if (not group["has_segmentation"]) and is_polygon(
                                candidate
                            ):  # group was all points, but a segmentation has been found; prefer that.
                                group["members"] = [candidate]
                                group["member_hashes"] = [candidate["seg_hash"]]
                                new_map = np.array(c_map)
                            group["has_segmentation"] = group[
                                "has_segmentation"
                            ] or is_polygon(candidate)
                    group["seg_map"] = new_map
            for idx, flag in enumerate(was_merged):
                if not flag:
                    # Candidate didn't merge: generate a new group:
                    merged.append(
                        {
                            "seg_map": np.array(candidates[idx]["seg_map"]),
                            "last_slice": s,
                            "members": [candidates[idx]],
                            "member_hashes": [candidates[idx]["seg_hash"]],
                            "has_segmentation": is_polygon(candidates[idx]),
                        }
                    )
    # Now we've run one direction, but in cases of branching nodules, it might be necessary to run the other direction as well:
    # When running backward, every candidate already has a group membership: we don't generate any new groups, we just (maybe) add members.
    # But, along the way we will "revive" existing groups when we get to their last active slice.
    # for s in range(last-1, first-1, -1):
    #     if s in candidates_by_slice:
    #         candidates = candidates_by_slice[s]
    #         was_merged = [False for candidate in candidates]
    #         for group in merged:
    #             if group['last_slice'] == (s+1):
    #                 g_map   = group['seg_map']
    #                 new_map = np.zeros(img_shape)     # Generating the "new map" is trickier here; we need
    #                 for member in group['members']:   # the union of all of the current
    #                     if member['sliceNo'] == s:    # members from this slice as an initial "new_map"
    #                         new_map = new_map.astype(bool) | np.array(member['seg_map'], dtype=bool)
    #                 for idx, candidate in enumerate(candidates):
    #                     if not candidate['seg_hash'] in group['member_hashes']:
    #                         c_map = candidate['seg_map']
    #                         if ((c_map != 0) & (g_map != 0)).any():
    #                             # There is overlap.  Merge; true merge if group is poly and so is candidate or if group is all points
    #                             # and so is candidate; weak merge (drop candidate) if group is poly and candidate is not; restart
    #                             # group if candidate is poly and group is all points.
    #                             was_merged[idx] = True
    #                             if (group['has_segmentation'] == is_polygon(candidate)) or (not group['has_segmentation']):
    #                                 new_map = (new_map.astype(bool) | c_map.astype(bool))
    #                                 group['members'].append(candidate)
    #                                 group['last_slice'] = s
    #                             if (not group['has_segmentation']) and is_polygon(candidate):  # group was all points, but a segmentation has been found; prefer that.
    #                                 group['members'] = [candidate]
    #                                 new_map = np.array(c_map)
    #                             group['has_segmentation'] = (group['has_segmentation'] or is_polygon(candidate))
    #                 group['seg_map'] = new_map
    #             else:
    #                 group['seg_map'] = None
    # Some groups may have gained shared members now (which means they can be further merged).  Do so:
    final_merged = {}
    final_hash_map = {}
    fm_key = 0
    for group in merged:
        common_members_with = []
        for member in group["members"]:
            m_hash = member["seg_hash"]
            if m_hash in final_hash_map:
                group_idx = final_hash_map[m_hash]
                common_members_with.append(group_idx)
        supergroup = {"members": [m for m in group["members"]]}
        remove_keys = []
        for group_idx in common_members_with:
            supergroup["members"].extend(final_merged[group_idx]["members"])
            remove_keys.append(group_idx)
        for key in list(set(remove_keys)):
            del final_merged[key]
        remove_keys = None
        for member in supergroup["members"]:
            final_hash_map[member["seg_hash"]] = fm_key
        final_merged[fm_key] = supergroup
        fm_key += 1
    # Rebuild as a simple list
    merged = []
    for fm_key in final_merged:
        if len(final_merged[fm_key]["members"]) > 0:
            merged.append(final_merged[fm_key]["members"])
    final_merged = final_hash_map = None

    # Now post-process to set group details
    for group_idx, group in enumerate(merged):
        group_malignancy = (
            group_inclusion
        ) = group_is_nodule = group_has_segmentation = 0
        group_x = group_y = group_z = group_slice = group_volume = 0
        for member_idx, member in enumerate(group):
            if not "has_segmentation" in member:
                member["has_segmentation"] = is_polygon(member)
            group_has_segmentation = max(
                group_has_segmentation, int(member["has_segmentation"])
            )
            group_malignancy = max(group_malignancy, member["malignancy"])
            group_inclusion = max(group_inclusion, member["inclusion"])
            group_is_nodule = max(group_is_nodule, member["isNodule"])
            group_x += member["x"]
            group_y += member["y"]
            group_z += member["z"]
            group_slice += member["sliceNo"]
            group_volume += member["area"] if "area" in member else 0
            if not keep_bin_image:
                for key in ["seg_map", "merged_count", "has_segmentation"]:
                    try:
                        del merged[group_idx][member_idx][key]
                    except:
                        pass
        group_size = len(group) if len(group) > 0 else 1
        group_volume = float(group_volume) if group_volume > 0 else 1
        group_x = float(group_x) / group_size
        group_y = float(group_y) / group_size
        group_z = float(group_z) / group_size
        group_slice = int(round(float(group_slice) / group_size))
        merged[group_idx] = {
            "x": group_x,
            "y": group_y,
            "z": group_z,
            "sliceNo": group_slice,
            "isNodule": group_is_nodule,
            "inclusion": group_inclusion,
            "malignancy": group_malignancy,
            "volume": group_volume,
            "members": group,
            "has_segmentation": group_has_segmentation,
            "noduleID": nodule_unique_id(group),
        }
    return merged


def verify_candidate_segmentations(
    candidates,
    dicom_dir,
    tmp_dir=None,
    candidate_limit=None,
    volume_max_threshold=100000,
    volume_min_threshold=4,
):
    import auto_segment_around_seed
    from tempfile import mkdtemp
    from shutil import rmtree

    tmp_dir = mkdtemp(
        dir=tmp_dir if (tmp_dir is not None) and (os.path.isdir(tmp_dir)) else None
    )
    nii_file = None
    nii_img = None
    n_slices = 0
    height = 0
    orig_img = None
    candidates_out = []
    candidate_limit = int(candidate_limit) if candidate_limit is not None else None
    print("candidate limit {}".format(candidate_limit))
    for candidate in candidates:
        if is_polygon(candidate):
            print("candidate is polygon")
            candidates_out.append(candidate)
        else:
            # print("\n\nCandidate Seed Point:")
            # print(candidate)
            print("candidate is seed point; create segmentation")
            # For any seed-point-only candidate, we need a segmentation.
            if nii_file is None:
                nii_file, nii_img = get_nifti_from_dicom_dir(dicom_dir, tmp_dir)
                nii_img_size = nii_img.GetSize()
                n_slices = nii_img_size[2]
                height = nii_img_size[1]
                print("nii image shape: {}".format(nii_img.GetSize()))

            # Seed is given as (x, h-y, n_slices - sliceNo)
            # seed     = (candidate['y']+1, height - candidate['x'] - 1, n_slices - candidate['sliceNo'])
            seed = (
                candidate["x"] + 1,
                height - candidate["y"] - 1,
                n_slices - candidate["sliceNo"],
            )
            seed = tuple([int(round(v)) for v in list(seed)])
            print(
                "candidate at {} {} {} : seed: {}".format(
                    candidate["x"], candidate["y"], candidate["sliceNo"], seed
                )
            )
            run_info = None
            consensus_img = None
            try:
                consensus_img, run_info = auto_segment_around_seed.get_consensus_for_seed(
                    nii_img,
                    candidate["noduleID"],
                    tmp_dir,
                    seed,
                    min_size=volume_min_threshold,
                    max_size=volume_max_threshold,
                )
            except Exception as exc:  # pylint: disable=W0703
                eprint("Encountered critical exception:\n{}".format(exc))
                pass  # Just drop this candidate from consideration

            if (
                run_info is not None
                and consensus_img is not None
                and type(consensus_img) != type(np.array([]))
            ):
                # Now we have the binary segmentation in a temporary nifti file;
                # extract it back into the candidate as a polygon per slice:
                # This seems silly since we have the binary image now, but
                # it is memory conservative in case we run a strategy that
                # involves a large number of candidates in memory at once.
                # NOTE: This can be troublesome if there are disjoint areas in the same
                #       slice (lobes that are part of a larger 3D nodule, for example).
                #       For that reason, we need to label all filled areas and create a
                #       separate segmentation polygon for each of them.
                try:
                    consensus_img = sitk.GetArrayFromImage(consensus_img)
                except Exception as exc:
                    print(
                        "Exception when trying to get array from segmentation image: {}".format(
                            exc
                        )
                    )
                    consensus_img = None
                    candidate = None
                    continue
                pix_volume = consensus_img[consensus_img > 0].sum()
                if pix_volume > 0 and pix_volume <= volume_max_threshold:
                    seg_group = copy.deepcopy(candidate)
                    seg_group["seriesUID"] = get_candidate_seriesUID(candidate)
                    is_nodule = seg_group["isNodule"]
                    inclusion = seg_group["inclusion"]
                    malignancy = seg_group["malignancy"]
                    seg_group["members"] = []
                    consensus_img = consensus_img[
                        ::-1, ::-1
                    ]  # put the image back in our usual orientation
                    img_shape = consensus_img.shape
                    for slice_idx in range(img_shape[0]):
                        if consensus_img[slice_idx].sum() > 0:
                            # This slice is included.  Get the polygon(s).
                            img_labels, num_objects = measurements.label(
                                consensus_img[slice_idx]
                            )
                            for o_i in range(num_objects):
                                print(
                                    "Creating polygon for object {} of {}...".format(
                                        o_i + 1, num_objects
                                    )
                                )
                                isolated_object = np.array(consensus_img[slice_idx])
                                isolated_object[img_labels != (o_i + 1)] = 0
                                isolated_object[isolated_object != 0] = 1
                                boundary, coords = get_polygon_from_binary_img(
                                    isolated_object
                                )
                                slice_polygon_dict = make_polygon_dict_given_boundary(
                                    boundary, coords
                                )
                                slice_polygon_dict["z"] = slice_idx + 1
                                slice_polygon_dict["sliceNo"] = slice_idx + 1
                                # TODO[KLUDGE]: The fix for x,y transposition works, but should probably be fixed
                                #              in a deeper way.  The segmentation pixels end up in the right places in the
                                #              DICOM output, and match the fixed x,y location produced by the fix below
                                #              but this seems like a "band aid".
                                # Due to x,y vs row,col handling here, the axes need to be flipped in the ID:
                                slice_polygon_dict["x"], slice_polygon_dict["y"] = (
                                    slice_polygon_dict["y"],
                                    slice_polygon_dict["x"],
                                )
                                # [END]
                                slice_polygon_dict["noduleID"] = nodule_unique_id(
                                    slice_polygon_dict
                                )
                                slice_polygon_dict["isNodule"] = is_nodule
                                slice_polygon_dict["inclusion"] = inclusion
                                slice_polygon_dict["malignancy"] = malignancy
                                seg_group["members"].append(slice_polygon_dict)
                    seg_group = set_group_stats(seg_group)

                    print("Generated candidate group: {}".format(seg_group))
                    candidates_out.append(seg_group)
                else:
                    print(
                        "Generated segmentation outside allowed volume range: {} vs {}".format(
                            pix_volume, volume_max_threshold
                        )
                    )
                consensus_img = None
            else:
                candidate = None
                consensus_img = None
        if candidate_limit is not None and len(candidates_out) >= candidate_limit:
            break
    # clean the temporary files:
    rmtree(tmp_dir)
    return candidates_out


def filter_by_seed_points(
    candidates, seed_points, consensus="vote", image_shape=(512, 512)
):
    """
    Filter candidates and select only the ones that contain at least one of `seed_points`, which should be a dict
    containing x, y, z, sliceNo
    """
    vote_threshold = get_vote_threshold(consensus)
    selected = []
    for candidate_group in candidates:
        group_by_slice = {}
        if not "members" in candidate_group:
            candidate_group = {"members": candidate_group}
        for candidate in candidate_group["members"]:
            slc = candidate["sliceNo"]
            if not slc in group_by_slice:
                group_by_slice[slc] = []
            group_by_slice[slc].append(candidate)
        for slc in group_by_slice:
            pixels = np.zeros(image_shape)
            # Union all segmentations that overlap within the same group
            # (this is rare, but could happen if you are slicing through lobes of
            #  a larger 3D object)
            n_drawn_polygons = 0
            for s_idx, slice in enumerate(group_by_slice[slc]):
                pixels = np.array(pixels) + draw_polygon_from_coords(
                    pixels, slice["coords"]
                )
                n_drawn_polygons += 1
            # Figure out how many discrete objects appear here:
            object_labels, n_objects_on_slice = measurements.label(
                np.array(pixels, dtype=bool)
            )
            # For each distinct polygon, threshold it.
            for i in range(n_objects_on_slice):
                obj_mask = np.where(object_labels == (i + 1))
                obj_max = pixels[obj_mask].max()
                threshold = vote_threshold * obj_max
                # Apply the threshold to the masked region only
                pixels[obj_mask] = [0 if v < threshold else 1 for v in pixels[obj_mask]]

            # Now if there is any seed point inside the segmented region, we keep it.
            keep_candidate = False
            for seed in seed_points:
                x, y, sliceNo = (
                    int(round(seed["x"])),
                    int(round(seed["y"])),
                    int(round(seed["sliceNo"])),
                )
                if sliceNo == slc:
                    if pixels[y][x] > 0:
                        keep_candidate = True
                        break  # no need to look for other seeds; one that hits the segmentation was found
            if keep_candidate:
                selected.append(candidate_group)
                break  # No need to examine more slices of this group; a seed was found

    return selected


def get_args():
    import argparse

    # construct the argument parser and parse the arguments
    ap = argparse.ArgumentParser(prog="{0}".format(os.path.basename(sys.argv[0])))
    ap.add_argument(
        "dicom_dir",
        metavar="dicom-dir",
        help="Directory containing DICOM source images",
    )
    ap.add_argument(
        "output_dir", metavar="output-dir", help="Name of destination directory."
    )
    ap.add_argument(
        "--tmp",
        dest="tmp_dir",
        required=False,
        default="/tmp",
        help="Directory for storing temporary files. (default=/tmp)",
    )
    ap.add_argument(
        "--candidates",
        dest="candidate_file",
        required=False,
        help="File containing seed coordinates; these will override any XML candidates present in the DICOM directory.",
    )
    ap.add_argument(
        "--seed",
        dest="seed_coords",
        metavar=("x", "y", "z"),
        nargs=3,
        required=False,
        help="Coordinates ((x y z) where z is slice #) of a custom seed point; will override any XML candidates present in the DICOM directory.",
    )
    ap.add_argument(
        "--consensus",
        dest="consensus_strategy",
        required=False,
        default="vote",
        help="Strategy for consensus among multiple segmentations for the same candidate, choose one of {union, intersect, vote, vote:PCT} where PCT is in the range [0,100] and represents the amount of consensus required in the vote.",
    )
    ap.add_argument(
        "--limit",
        dest="candidate_limit",
        metavar="N",
        required=False,
        help="If more than one candidate passes filters, limit to first N candidates. (default=all)",
    )
    ap.add_argument(
        "--segmented-only",
        dest="segmented_only",
        action="store_true",
        help="Do not consider any candidate that doesn't have a pre-defined segmentation boundary."
        " When used in conjunction with --candidates or --seed, it will act to select which segmentations to produce.",
    )
    ap.add_argument(
        "--output-prefix",
        dest="nodule_id_prefix",
        metavar="PREFIX",
        required=False,
        help="Prepend PREFIX to all output directory names.",
    )
    ap.add_argument(
        "--auto",
        dest="auto_choose",
        metavar="STRATEGY",
        required=False,
        help="Auto select candidate according to STRATEGY in {all, malignancy-filter:LEVEL:DIRECTION, malignancy-range:LOW:HIGH:SORT, most-malignant, least-malignant, largest, random[:SEED]}"
        + " where SEED is a seed for the random selection so the choice can be repeated if desired (omit to get a different sequence each time), "
        " valid values for DIRECTION are {>, <, >=, <=, =}; LOW is lowest malignancy to include, HIGH is highest malignancy to include,"
        " SORT must be one of {ascending, descending, random[:SEED]} where SEED behaves the same as previously described.",
    )

    args = vars(ap.parse_args())
    return args


def run_standalone():
    args = get_args()
    dicom_dir = args["dicom_dir"]
    output_dir = args["output_dir"]
    tmp_dir = args["tmp_dir"]
    candidate_limit = (
        int(args["candidate_limit"]) if args["candidate_limit"] is not None else None
    )
    candidate_file = args["candidate_file"]
    seed_coords = args["seed_coords"]
    auto_choose = False if args["auto_choose"] is None else args["auto_choose"]
    nodule_id_prefix = args["nodule_id_prefix"]
    consensus_strategy = args["consensus_strategy"]
    segmented_only = args["segmented_only"]
    selected_segmented_only = (
        candidate_file is not None or seed_coords is not None
    ) and segmented_only == True
    selected_seed_coords = None

    polygons = []
    centroids = []
    if candidate_file is None and seed_coords is None:
        polygons = get_tumor_polygons(dicom_dir)
        centroids = get_tumor_centroids(dicom_dir) if (segmented_only == False) else []
    if candidate_file is not None:
        centroids = read_coords_file(candidate_file)
    if selected_segmented_only:
        polygons = get_tumor_polygons(dicom_dir)

    if candidate_file is not None:
        included_candidates = []
        for centroid in centroids:
            if "patient" in centroid:
                if centroid["patient"] == os.path.basename(dicom_dir):
                    included_candidates.append(centroid)
            else:
                included_candidates.append(centroid)
        centroids = included_candidates

    if seed_coords is not None:
        if candidate_file is None:
            centroids = []
        x, y, slice_no = [int(round(float(x))) for x in seed_coords]
        centroids.append(
            {
                "sliceNo": int(round(float(slice_no))),
                "x": float(x),
                "y": float(y),
                "z": float(slice_no),
                "isNodule": -1,
                "inclusion": -1,
                "malignancy": -1,
                "noduleID": make_unique_id_string(
                    x, y, slice_no, -1, -1, -1, nodule_id_prefix
                ),
            }
        )

    if selected_segmented_only:
        selected_seed_coords = [
            {"x": c["x"], "y": c["y"], "z": c["z"], "sliceNo": c["sliceNo"]}
            for c in centroids
        ]
        centroids = []

    centroids_by_scan = None
    chosen_candidates = None

    if len(polygons) == 0 and len(centroids) == 0:
        centroids_by_scan = xmlReader.collect_centroids_from_xml(dicom_dir)
        centroids_by_scan = xmlReader.merge_centroids(centroids_by_scan)
        chosen_candidates = centroids_by_scan[dicom_dir]
    else:
        file_list = get_dicomdir_files(dicom_dir, full_path=True)
        centroids_by_scan = {
            dicom_dir: {"centroids": centroids, "file_list": file_list}
        }
        centroids_by_scan[dicom_dir]["centroids"].extend(polygons)
        centroids_by_scan = xmlReader.merge_centroids(centroids_by_scan)
        chosen_candidates = centroids_by_scan[dicom_dir]

    filtered_candidates = None
    verified_candidates = None
    retry = True
    while retry and len(chosen_candidates) > 0:
        retry = False
        filtered_candidates = filter_groups(chosen_candidates, auto_choose)
        if len(filtered_candidates) == 0:
            verified_candidates = []
            print("NO CANDIDATES passed the filter for {}".format(dicom_dir))
            break
        # print("\nfiltered: {}".format(filtered_candidates))
        verified_candidates = verify_candidate_segmentations(
            filtered_candidates,
            dicom_dir=dicom_dir,
            tmp_dir=tmp_dir,
            candidate_limit=candidate_limit,
        )
        # print("\nVERIFIED: {}".format(verified_candidates))
        retry = len(verified_candidates) == 0
        if retry:
            filt_ids = [c["noduleID"] for c in filtered_candidates]
            chosen_minus_filtered = []
            for candidate in chosen_candidates:
                if not (candidate["noduleID"] in filt_ids):
                    chosen_minus_filtered.append(candidate)
            chosen_candidates = chosen_minus_filtered

    chosen_candidates = verified_candidates

    if selected_segmented_only:
        chosen_candidates = filter_by_seed_points(
            chosen_candidates, selected_seed_coords, consensus=consensus_strategy
        )

    # print("\n\nFinal grouped candidates: \n")
    # for member in chosen_candidates:
    #     print("\nGroup: {}\n".format(member))

    for chosen in chosen_candidates:
        output_dir_suffix = nodule_unique_id(chosen, nodule_id_prefix)
        output_path = (
            os.path.join(output_dir, output_dir_suffix)
            if len(chosen_candidates) > 0
            else output_dir
        )
        write_3d_binary_image(
            dicom_dir, output_path, chosen, consensus=consensus_strategy
        )

    if len(chosen_candidates) == 0:
        sys.exit(1)


if __name__ == "__main__":
    run_standalone()
