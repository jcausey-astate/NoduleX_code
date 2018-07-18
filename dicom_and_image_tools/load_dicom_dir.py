from __future__ import print_function
import sys, os
import pydicom as dicom, numpy as np


def set_dicom_pixels(dicom_obj, new_pixels):
    new_pixels = np.array(new_pixels, dtype=dicom_obj.pixel_array.dtype)
    dicom_obj.Rows, dicom_obj.Cols = dicom_obj.pixel_array.shape
    dicom_obj.pixel_array[:] = new_pixels[:]
    dicom_obj.PixelData = new_pixels.tostring()
    return dicom_obj


def dicom_file_list(indir):
    """Returns a list of the name of all dicom files in a dicom series given
    the directory in which they are stored.  The list is not in any particuar
    order."""
    files = []
    if os.path.isdir(indir):
        for rt, dr, fl in os.walk(indir):
            files.extend(fl)
            break
    else:
        raise IOError("{0} was not a directory.".format(indir))

    dicomfiles = [os.path.join(indir, f) for f in files if f[-4:] == ".dcm"]
    return dicomfiles


def dicom_files(indir):
    """Returns a list of the names of all files in a dicom series given the
    directory in which they're stored. The list is sorted based upon the
    InstanceNumber DICOM parameter"""
    import dicom

    dicomfiles = dicom_file_list(indir)

    slices = [(f, dicom.read_file(os.path.join(indir, f))) for f in dicomfiles]

    slices.sort(key=lambda x: x[1].InstanceNumber)

    return [os.path.join(indir, x[0]) for x in slices]


def dicom_dir_hash(dicomdir):
    """Compute an the sha1 hash of the contents of a dicom stack. Returns the
    hexadecimal representation as a string."""
    import hashlib

    sha = hashlib.sha1()

    for dcm in dicom_files(dicom):
        with open(dcm, "rb") as f:
            sha.update(f.read())

    return sha.hexdigest()


def dicom_img_hash(dicom_object):
    """Compute the sha1 hash of the pixel values for a dicom object."""
    import hashlib

    sha = hashlib.sha1()
    sha.update(dicom_object.pixel_array)
    return sha.hexdigest()


def load_dicom_dir_back_compat(
    dicomdir, all_headers=False, include_file_list=False, rescale=False
):
    """Load the images from dicomdir into a numpy array"""

    slices = filter_duplicates_from_dicomdir(dicomdir)

    slices.sort(key=lambda x: x[1].InstanceNumber)
    slice_shape = slices[0][1].pixel_array.shape
    img = np.zeros((len(slices), slice_shape[0], slice_shape[1]))

    for idx, slc in enumerate(slices):
        img[idx] = np.array(slc[1].pixel_array)

    hdr = slices[0][1]

    if rescale:
        img = rescale_image(img, hdr)
        try:
            hdr.RescaleIntercept = 0.0
        except:
            pass
        try:
            hdr.RescaleSlope = 1.0
        except:
            pass

    if all_headers or include_file_list:
        hdr = [s[1] for s in slices] if include_file_list == False else slices
        if rescale:
            for h in hdr:
                h.RescaleIntercept = 0.0
                h.RescaleSlope = 1.0

    return img, hdr


def load_dicom_dir(dicomdir, all_headers=False, include_file_list=False, rescale=None):
    """
    Load the images from dicomdir into a numpy array, returned in "view-matrix" dimension ordering
    (Z, Y, X), where Z==0 is the superior ("head") end, increasing Z is toward the inferior ("feet").
    Y increases toward the posterior, X increases toward the patient's left.
    """
    if rescale is not None and rescale == False:
        print(
            "Warning: load_dicom_dir() with `rescale=False` is deprecated; please update dependant code.",
            file=sys.stderr,
        )
        return load_dicom_dir_back_compat(
            dicomdir, all_headers, include_file_list, rescale=False
        )
    if rescale is None:
        print(
            "Notice: load_dicom_dir() now rescales to HU by default; please update dependent code if necessary.",
            file=sys.stderr,
        )

    img, hdr = load_dicom_dir_dicom_orientation(
        dicomdir, all_headers=all_headers, include_file_list=include_file_list
    )
    affine = (
        hdr.Affine
        if not (all_headers or include_file_list)
        else hdr[0].Affine
        if not include_file_list
        else hdr[0][1].Affine
    )

    img = np.transpose(img, (2, 1, 0))  # Move from (X,Y,Z) to (Z,Y,X)
    img = img[::-1]  # Arrange slices so "head" end is at index 0.
    affine[0, 0], affine[2, 2],  # Re-arrange the affine to match the new layout of img
    affine[0, 3], affine[2, 3] = affine[2, 2], affine[0, 0]
    affine[2, 3], affine[0, 3]
    # The re-arrangement above only works if all non-diagonal components are zeros, so check for that.
    af_non_diagonal = affine[0:3, 0:3].copy()
    for i in range(3):
        af_non_diagonal[i, i] = 0
    if af_non_diagonal.sum() != 0:
        raise RuntimeError(
            "DICOM orientation non-standard, leading to affine transform with non-diagonal components."
        )
    del af_non_diagonal

    if all_headers or include_file_list:
        hdr = [s[1] for s in slices] if include_file_list == False else slices
        for h in hdr:
            if include_file_list:
                h = h[1]
            h.Affine = affine
    else:
        hdr.Affine = affine

    return img, hdr


def load_dicom_dir_dicom_orientation(
    dicomdir, all_headers=False, include_file_list=False
):
    """
    Load the images from dicomdir into a numpy array, returned in "DICOM" standard dimension ordering
    (X, Y, Z), where Z==0 is the inferior ("feet") end, increasing Z is toward the superior ("head").
    Y increases toward the posterior, X increases toward the patient's left.
    """
    import logging

    logging.basicConfig()
    import dicom_numpy as dnp  # https://github.com/innolitics/dicom-numpy

    slices = filter_duplicates_from_dicomdir(dicomdir)

    img, affine = dnp.combine_slices([s[1] for s in slices])

    hdr = slices[0][1]
    hdr.Affine = affine
    hdr.RescaleIntercept = 0
    hdr.RescaleSlope = 1.0

    if all_headers or include_file_list:
        hdr = [s[1] for s in slices] if include_file_list == False else slices
        for h in hdr:
            if include_file_list:
                h = h[1]
            h.RescaleIntercept = 0
            h.RescaleSlope = 1.0
            h.Affine = affine

    return img, hdr


def save_dicom_dir(dicomdir, img, header, prefix=None):
    """Save DICOM image series as a dicom directory."""
    header_list = False
    if type(header) == type([]):
        header_list = True
    if prefix is None:
        prefix = os.path.splitext(os.path.basename(dicomdir))[0]
    if not os.path.isdir(dicomdir):
        os.makedirs(dicomdir)

    for idx, slc in enumerate(img):
        fname = os.path.join(dicomdir, "{}_{}.dcm".format(prefix, idx + 1))
        dcm = header if not header_list else header[idx]
        dcm.InstanceNumber = str(idx + 1)
        dcm = set_dicom_pixels(dcm, slc)
        dcm.save_as(fname)


def get_dicom_shape_from_dicomdir(dicomdir):
    """Just gets the (numpy array) shape of the dicom 3D image"""
    img, _ = load_dicom_dir(dicomdir)
    return img.shape


def rescale_image(img, hdr):
    """Apply rescale formula from DICOM header, if that information is available."""
    if type(hdr) == type([]):
        hdr = hdr[0]
    img = np.array(img)
    img_type = img.dtype
    rescale_slope = 1.0
    rescale_intercept = 0

    try:
        rescale_intercept = hdr.RescaleIntercept
    except:
        print(
            "Warning: no rescale intercept available in image header.", file=sys.stderr
        )

    try:
        rescale_slope = hdr.RescaleSlope
    except:
        print("Warning: no rescale slope available in image header.", file=sys.stderr)

    img = float(rescale_slope) * img.astype(np.float64)
    img = img + np.int16(rescale_intercept)
    img = img.astype(img_type)
    return img


def filter_duplicates(dicom_series, return_indices=False):
    """Remove duplicate slices from a series (list) of dicom objects.
    Use return_indices = True to return the list of indices to keep instead
    of returning the modified dicom series."""
    locations_seen = set([])
    images_seen = set([])
    result_series = []
    for idx, dcm in enumerate(dicom_series):
        img_hash = dicom_img_hash(dcm)
        location = ",".join([str(v) for v in dcm.ImagePositionPatient])
        if location in locations_seen:
            if img_hash not in images_seen:
                raise RuntimeError(
                    "DICOM series error: Duplicate image location, but different pixel data."
                )
        else:
            result_series.append(dcm if not return_indices else idx)
            locations_seen.add(location)
            images_seen.add(img_hash)
    return result_series


def filter_duplicates_from_paths(dicom_files, return_indices=False):
    """Remove duplicate slices from a list of DICOM file paths."""
    files = [dicom.read_file(f) for f in dicom_files]
    mask = filter_duplicates(files, return_indices=True)
    return [dicom_files[i] for i in mask] if not return_indices else mask


def filter_duplicates_from_dicomdir(dicomdir, return_paths_only=False):
    """Returns a list of dicom object or just file paths where no duplicate 
    slices are included, given a DICOM directory."""
    files = filter_duplicates_from_paths(dicom_file_list(dicomdir))
    result = files
    if not return_paths_only:
        result = [(f, dicom.read_file(f)) for f in files]
    return result


if __name__ == "__main__":
    print("This file is not meant to run stand-alone", file=sys.stderr)
