'''
Library for loading a 3D image from a direcotry containing multiple 2D DICOM images.
'''

from __future__ import print_function
import sys, os
import dicom, numpy as np

def set_dicom_pixels(dicom_obj, new_pixels):
    '''Set the pixels in a dicom object, given a Numpy array `new_pixels`'''
    new_pixels = np.array(new_pixels, dtype=dicom_obj.pixel_array.dtype)
    dicom_obj.Rows, dicom_obj.Cols = dicom_obj.pixel_array.shape
    dicom_obj.pixel_array[:]       = new_pixels[:]
    dicom_obj.PixelData            = new_pixels.tostring()
    return dicom_obj

def dicom_file_list(indir):
    '''Returns a list of the name of all dicom files in a dicom series given
    the directory in which they are stored.  The list is not in any particuar
    order.'''
    files       = []
    if os.path.isdir(indir):
        for rt, dr, fl in os.walk(indir):
            files.extend(fl)
            break
    else:
        raise IOError('{0} was not a directory.'.format(indir))
    
    dicomfiles = [f for f in files if f[-4:] == '.dcm']
    return dicomfiles    

def dicom_files(indir):
    '''Returns a list of the names of all files in a dicom series given the
    directory in which they're stored. The list is sorted based upon the
    InstanceNumber DICOM parameter'''
    import dicom

    dicomfiles = dicom_file_list(indir)

    slices = [(f, dicom.read_file(os.path.join(indir, f)))
              for f in dicomfiles]

    slices.sort(key=lambda x: x[1].InstanceNumber)

    return [os.path.join(indir, x[0]) for x in slices]

def dicom_dir_hash(dicomdir):
    '''Compute an the sha1 hash of the contents of a dicom stack. Returns the
    hexadecimal representation as a string.'''
    import hashlib

    sha = hashlib.sha1()

    for dcm in dicom_files(dicom):
        with open(dcm, 'rb') as f:
            sha.update(f.read())

    return sha.hexdigest()

def load_dicom_dir(dicomdir, all_headers=False, include_file_list=False, rescale=False):
    '''Load the images from dicomdir into a numpy array'''

    dicomfiles = dicom_file_list(dicomdir)
    slices     = [(f, dicom.read_file(os.path.join(dicomdir, f)))
                    for f in dicomfiles]
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
            hdr.RescaleSlope     = 1.0
        except:
            pass

    if all_headers or include_file_list:
        hdr = [s[1] for s in slices] if include_file_list == False else slices
        if rescale:
            for h in hdr:
                h.RescaleIntercept = 0.0
                h.RescaleSlope     = 1.0
    
    return img, hdr

def save_dicom_dir(dicomdir, img, header, prefix=None):
    '''Save DICOM image series as a dicom directory.'''
    header_list = False
    if type(header) == type([]):
        header_list = True
    if prefix is None:
        prefix = os.path.splitext(os.path.basename(dicomdir))[0]
    if not os.path.isdir(dicomdir):
        os.makedirs(dicomdir)

    for idx, slc in enumerate(img):
        fname = os.path.join(dicomdir, '{}_{}.dcm'.format(prefix,idx+1))
        dcm   = header if not header_list else header[idx]
        dcm.InstanceNumber  = str(idx + 1)
        dcm   = set_dicom_pixels(dcm, slc)
        dcm.save_as(fname)

def get_dicom_shape_from_dicomdir(dicomdir):
    '''Just gets the (numpy array) shape of the dicom 3D image'''
    img, _ = load_dicom_dir(dicomdir)
    return img.shape

def rescale_image(img, hdr):
    '''Apply rescale formula from DICOM header, if that information is available.'''
    if type(hdr) == type([]):
        hdr = hdr[0]
    img = np.array(img)
    img_type = img.dtype
    rescale_slope     = 1.0
    rescale_intercept = 0
    
    try:
        rescale_intercept = hdr.RescaleIntercept
    except:
        print("Warning: no rescale intercept available in image header.", file=sys.stderr)
    
    try:            
        rescale_slope = hdr.RescaleSlope
    except:
        print("Warning: no rescale slope available in image header.", file=sys.stderr)

    img = float(rescale_slope) * img.astype(np.float64)
    img = img + np.int16(rescale_intercept)
    img = img.astype(img_type)
    return img


if __name__ == "__main__":
    print("This file is not meant to run stand-alone", file=sys.stderr)
