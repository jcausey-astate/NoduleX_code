'''Produce a segmentation of the lungs, and produce a set of seeds for that
segmentation.'''

from __future__ import print_function
import sys
import argparse
import SimpleITK as sitk  # pylint: disable=F0401
import os
import numpy as np

def showmid(img):
    import matplotlib.pyplot as plt
    array = sitk.GetArrayFromImage(img)
    mid_i = len(array) // 2
    plt.imshow(array[mid_i], cmap='gray')
    plt.show()

def otsu(img):
    '''Use an 'otsu' thresholding to segment out high- and low-attenuation
    regions. In every chest CT case I've seen, it produces an "air" region and
    a "soft tissue + bone" region, but a constrast CT might produce different
    results.'''

    array = sitk.GetArrayFromImage(img)
    minval = np.min(array)
    frac_minval = np.count_nonzero(array == minval) / float(array.size)

    filt = sitk.OtsuThresholdImageFilter()
    if frac_minval > .1:
        mask = np.logical_not(array == minval)
        mask = mask.astype('uint8')

        filt.SetMaskValue(1)
        filt.SetMaskOutput(False)

        mask = sitk.GetImageFromArray(mask)
        mask.CopyInformation(img)

        return filt.Execute(img, mask)
    else:
        return filt.Execute(img)

def find_hist_peaks(img, drop_first_max=False):
    '''Finds the first two interior peaks in the histogram; usually a good binary
    threshold lies halfway between these two peaks.'''
    import scipy.signal as sig
    array        = sitk.GetArrayFromImage(img)
    
    hist         = np.histogram(array, 16)
    hist_counts  = np.array(hist[0])
    hist_centers = np.array(hist[1])
    peaks        = sig.argrelextrema(hist_counts, np.greater)[0]
    
    if hist_counts[peaks[0]] < (0.01 * (hist_counts[peaks].mean())):  # Don't allow a "noise" peak at the beginning.
        peaks = peaks[1:]
    
    if len(peaks) < 2: # Only one peak is found if the other peak is really at the edge.
        peaks = np.array([0, peaks[0]])
    elif hist_counts[0] > hist_counts[peaks[0]]:
        if not drop_first_max:
            peaks = np.array([0, peaks[0]])
    elif drop_first_max:
        peaks = np.array([peaks[1], peaks[2]])
    
    i_L    = hist_centers[peaks[0]]
    i_FM   = hist_centers[peaks[1]]
    return i_L, i_FM

def hist_threshold(img, drop_first_max=False):
    '''Use histogram-based approach to threshold lung image; this is a fallback
    for the rare cases when otsu thresholding fails'''
    i_L, i_FM = find_hist_peaks(img, drop_first_max=drop_first_max)
    thresh   = (float(i_FM + i_L) / 2.0)
    
    array    = sitk.GetArrayFromImage(img)
    shape    = array.shape
    mask     = np.zeros(array.shape, dtype=np.int64)
    mask[array >= thresh] = 1
    array    = mask
    
    array    = sitk.GetImageFromArray(array)
    array.CopyInformation(img)
    return array

def dialate(img, probe_size):
    '''Once lungs are segmented out specifically, there's a tendency to get
    little islands in the lung fields. This dialates the selection to remove
    islands and to generate a smoother segmentation.'''
    filt = sitk.BinaryDilateImageFilter()
    filt.SetKernelType(filt.Ball)
    filt.SetKernelRadius(probe_size)

    return filt.Execute(img, 0, 1, False)


def find_components(img):
    '''Produce a separate label for each region in the binary image. Takes the
    type at the edge of the image (specifically at 0,0,0) to be 1 and zero for
    all other labels.'''

    array = sitk.GetArrayFromImage(img)

    bg_fixed = array == array[0, 0, 0]
    bg_fixed = bg_fixed.astype(array.dtype)

    new_img = sitk.GetImageFromArray(bg_fixed)
    new_img.CopyInformation(img)

    filt = sitk.ConnectedComponentImageFilter()
    return filt.Execute(new_img)


def dump(img, name):
    '''Dump several slices for the given numpy array.'''
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas

    figure = Figure()
    canvas = FigCanvas(figure)

    nplot = 9
    for i in range(1, nplot+1):
        ax = figure.add_subplot(3, 3, i)  # pylint: disable=C0103
        ax.imshow(sitk.GetArrayFromImage(img)[i*img.GetDepth()/(nplot+1)])

    canvas.print_figure(name)


def isolate_lung_field(img):
    '''Isolate the lung field only by taking the largest object or two objects
    that is/are not the chest wall (identified as 0 due to Otsu filtering) or 
    outside air (identified by appearing at the border).'''

    array = sitk.GetArrayFromImage(img)

    counts = np.bincount(np.ravel(array))
    outside    = array[0, 0, 0]
    chest_wall = 0

    counts[outside]    = 0
    counts[chest_wall] = 0
    ordered_indices    = np.argsort(counts)[::-1]
    top_two            = ordered_indices[0:2]
    top_two_log_counts = [int(round(lc)) if not np.isinf(lc) else -1 for lc in np.log10(counts[top_two])]
    # If lungs are disconnected, there will be two objects with similar counts
    keep_objects = top_two  
    # But if the lungs are connected there will be just one
    if top_two_log_counts[0] != top_two_log_counts[1]: 
        keep_objects = [top_two[0]]

    lung_only = np.array(array == keep_objects[0], dtype=array.dtype)
    if len(keep_objects) > 1:
        lung_only[array == keep_objects[1]] = 1
    lung_only = sitk.GetImageFromArray(lung_only)
    lung_only.CopyInformation(img)

    return lung_only

# Original version:
# def isolate_lung_field(img):
#     '''Isolate the lung field only by taking the largest object that is not
#     the chest wall (identified as 0 due to Otsu filtering) or outside air
#     (identified by appearing at the border).'''

#     array = sitk.GetArrayFromImage(img)

#     counts = np.bincount(np.ravel(array))
#     print("bin counts: {}".format(counts))

#     outside = array[0, 0, 0]
#     chest_wall = 0

#     themax = (0, 0)
#     for (obj_index, count) in enumerate(counts):
#         if obj_index in [outside, chest_wall]:
#             continue
#         elif count > themax[1]:
#             themax = (obj_index, count)

#     lung_only = np.array(array == themax[0], dtype=array.dtype)
#     lung_only = sitk.GetImageFromArray(lung_only)
#     lung_only.CopyInformation(img)

#     return lung_only



def isolate_not_biggest(img):
    '''Takes an sitk image with labels for many regions and produces a binary
    mask with zero for the largest region (by number of voxels) and one
    everywhere else.'''
    array = sitk.GetArrayFromImage(img)

    counts = np.bincount(np.ravel(array))

    big = np.argmax(counts)

    not_big = np.array(array != big, dtype=array.dtype)
    not_big = sitk.GetImageFromArray(not_big)
    not_big.CopyInformation(img)

    return not_big

def checkdist(seeds):
    '''UNDER CONSTRUCTION'''

    dists = {}

    for (i, seed) in enumerate(seeds):
        for oseed in seeds[i+1:]:
            dist = sum([(seed[k] - oseed[k])**2
                        for k in range(len(seed))])**0.5

            dists[(seed, oseed)] = dist

    raise NotImplementedError("Checkdist is under construction.")

def ensure_img_border(img):
    '''Make sure that at least one voxel in the X-Y plane is the same intensity as the 
    (0,0,0) voxel so that background detection will not fail if the body touches the
    left and right sides of the scan.'''
    array        = sitk.GetArrayFromImage(img)
    shape        = array.shape
    bg_color     = array[0,0,0]
    array[:,:,0] = bg_color
    array[:,0,:] = bg_color
    array[:,:,shape[2]-1] = bg_color
    array[:,shape[2]-1,:] = bg_color
    array = sitk.GetImageFromArray(array)
    array.CopyInformation(img)
    return array

def segment_lung(img, options):
    '''Segment lung.'''
    orig   = img
    img    = ensure_img_border(img)
    img    = otsu(img)
    img    = find_components(img)
    counts = np.bincount(np.ravel(sitk.GetArrayFromImage(img)))
    
    if len(counts) < 4:  # otsu (rarely) fails to threshold well; try histogram thresholding if it does
        img = orig
        img = ensure_img_border(img)
        img = hist_threshold(img)
        img = find_components(img)
    counts = np.bincount(np.ravel(sitk.GetArrayFromImage(img)))
    
    if len(counts) < 4: # histogram thresholding can also fail for some images if there is a 
        img = orig      # "black ring" minimum very far below the image's true background.  Try to compensate...
        img = ensure_img_border(img)
        img = hist_threshold(img, True)
        img = find_components(img)
    counts = np.bincount(np.ravel(sitk.GetArrayFromImage(img)))
    
    if len(counts) < 4: # At this point, bail out.
        raise RuntimeError("Unable to find a threshold that allows segmenting the lungs.")
    
    img = isolate_lung_field(img)
    img = dialate(img, options['probe_size'])
    img = find_components(img)
    img = isolate_not_biggest(img)
    
    img = sitk.BinaryErode(img, int(round(options['probe_size'] * (5.0/7.0))),  # matches sitkstrats (7-2) but also scales
                           sitk.BinaryErodeImageFilter.Ball)
    return img

def lungseg(img, options):
    '''Segment lung (deprecated: use `segment_lung` instead)'''
    return segment_lung(img, options)