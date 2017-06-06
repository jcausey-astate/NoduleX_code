'''A collection of segmentation strategies from sitk instead of itk.'''
import SimpleITK as sitk  # pylint: disable=F0401
import numpy as np
import datetime
import os

from functools import wraps
import logging

import lungseg


def write(img, fname, compression=True):
    '''
    A wrapper for sitk WriteImage that defaults to using compression and
    creates directories where required.
    '''

    # ImageFileWriter fails if the directory doesn't exist. Create if req'd
    try:
        os.makedirs(os.path.dirname(fname))
    except OSError:
        pass

    sitk.WriteImage(img, fname, compression)


def read(fname):
    img = sitk.ReadImage(fname)
    return img


def hash_img(img, provenance=""):
    '''
    Calculate the hash of an image and the options that would be used to
    process it using sha512. This is used most frequently to cache images for
    later reuse.
    '''
    import hashlib

    sha = hashlib.sha512()
    sha.update(sitk.GetArrayFromImage(img))
    sha.update(provenance)

    return sha.hexdigest()


def log_size(func):
    '''A decorator that calculates the size of a segmentation.'''
    @wraps(func)
    def exec_func(img, opts=None):
        '''Execute func from outer context and compute the size of the image
        func produces.'''
        if opts is None:
            opts = {}

        (img_out, opts_out) = func(img, opts)

        opts['size'] = np.count_nonzero(  # pylint: disable=E1101
            sitk.GetArrayFromImage(img_out))

        return (img_out, opts_out)

    return exec_func


def options_log(func):
    '''A decorator that will modify the incoming options object to also include
    information about runtime and algorithm choice.'''
    @wraps(func)
    def exec_func(*args, **kwargs):
        '''The inner function for options_log'''
        start = datetime.datetime.now()

        (img, out_opts) = func(*args, **kwargs)

        out_opts['algorithm'] = func.__name__
        out_opts['time'] = datetime.datetime.now() - start

        return (img, out_opts)

    return exec_func


def cached(relevant_opts, max_cache_size=1):
    '''A decorator that uses options and input image to cache an image for
    possible later reuse.'''

    def cached_decorator(func):
        '''
        Produce a decorator configured with the options defined in 'cached'
        above.
        '''
        from collections import OrderedDict

        # build a cache once for each method?
        cache = OrderedDict()

        @wraps(func)
        def exec_func(img_in, opts):  # pylint: disable=C0111

            # SHA should be built based only on the options defined in
            # relevant opts. This set is not passed to the internal function,
            # however, since the output is decorated by other functions to
            # report on outcomes and whatnot
            limited_opts = dict([(i, opts[i]) for i in opts
                                 if i in relevant_opts])
            limited_opts['func'] = func.__name__
            sha = hash_img(img_in,
                           provenance=str(limited_opts))

            if sha in cache:
                logging.info("Loading '" + sha + "' from " + func.__name__ +
                             " cache (size=" + str(len(cache)) + ")")

                img = cache[sha]

                # refresh this image in the order, so we re-add it to the cache
                del cache[sha]
            else:
                logging.info("Cache miss for '" + sha + "' from "
                             + func.__name__ + "(size="+str(len(cache)) + ")")
                (img, opts) = func(img_in, opts)
                cache[sha] = img

            # if the cache size exceeds the max, ditch the oldest entry
            # this shouldn't ever actually get executed more than once, but
            # you never know--recursion, for example, could cause many pushes
            # without appropriate pops.
            while len(cache) > max_cache_size:
                # the last item in the list (FIFO) popped only when last=True
                cache.popitem(last=True)

            return (img, opts)
        return exec_func
    return cached_decorator

def com_calc(img, max_size, min_size, lung_img):
    '''
    Calculate the center of mass of each of the labeled regions in img,
    excluding regions that are outside the lung given in lung_img, or the
    range size [min_size, max_size], which are reported treated as fractions of
    the input lung.
    '''
    from scipy.ndimage.measurements import center_of_mass as com
    # pylint: disable=E1101

    arr = sitk.GetArrayFromImage(img)
    lung_arr = sitk.GetArrayFromImage(lung_img)

    # Take elements from arr only when lung_arr is not zero, i.e. take only
    # regions in the lung.
    arr = np.where(lung_arr != 0, arr,
                   np.zeros(arr.shape,
                            dtype=arr.dtype))

    counts = np.bincount(np.ravel(arr))

    # volume per voxel is encoded in img spacing, with units mm^3
    vox_vol = reduce(lambda x, y: x * y, img.GetSpacing())

    # the size of the lung is the size of a voxel times the number of voxels
    lung_size = np.count_nonzero(lung_arr)*vox_vol

    # print sorted(counts)[-10:], sorted([c for c in counts if c != 0])[0:10]
    # print np.count_nonzero(lung_arr)*vox_vol

    logging.debug("Availiable labels and sizes: %s",
                  [p for p in enumerate(counts) if p[1] > 0])

    # We gate the deterministic seeds for their regions being of a reasonable
    # size.
    labels = [label for (label, n_vox) in enumerate(counts)
              if min_size*lung_size < n_vox*vox_vol < max_size*lung_size]

    # we don't want to calculate the COM of the background, no matter what
    # our min/max sizes are
    try:
        labels.remove(0)
    except ValueError:
        #Zero wasn't in the list
        pass

    # compute the center of mass for each of the elements up to 100 elements
    com_list = com(np.where(arr != 0, np.ones(arr.shape), np.zeros(arr.shape)),
                   labels=arr, index=labels)

    logging.debug("Label-COM correspondence: %s", dict(zip(labels, com_list)))

    # these are array-indexed and we take our seeds to be image-indexed
    # plus, they're floats and need to be cast back to integers
    seeds = [[int(k) for k in reversed(s)] for s in com_list
             if lung_arr[s] == 1]

    info = {'nseeds': len(seeds),
            'max_size': max_size,
            'min_size': min_size,
            'seeds': [s for s in seeds]}  # deep (enough) copy

    logging.info("%s in-lung of %s seeds from %s watershed labels",
                 len(seeds), len(labels), len(counts))

    return (seeds, info)

@options_log
def crop_to_segmentation(img, seg_img, padding_ratio=None, padding_px=None):
    '''
    Given an image and a segmentation of that image, crop the image to include
    only the portions of the image present in that segmentation. Only one of
    padding_px and padding_ratio can be specified.
    '''
    from bounding import bounding_cube

    # it's not allowed to specify both padding_ratio and padding_px
    assert not ((padding_px is not None) and (padding_ratio is not None))

    assert img.GetSize() == seg_img.GetSize()
    assert img.GetSpacing() == seg_img.GetSpacing()

    # calculate the bounding cube of the lung segmentation in numpy coordinates
    lims = bounding_cube(sitk.GetArrayFromImage(seg_img))

    # calculate the amount of each image to be removed (in itk indexing)
    lower_remove = [l[0] for l in reversed(lims)]
    upper_remove = [img.GetSize()[i] - l[1] for (i, l)
                    in enumerate(reversed(lims))]

    if padding_px is None:
        padding_ratio = 0.0
        padding = [int(padding_ratio*img.GetSize()[i])
                   for i in range(len(img.GetSize()))]
    elif padding_ratio is None:
        assert padding_px is not None
        padding = [padding_px]*len(img.GetSize())

    padded_removes = [[r[i] - padding[i] for i in range(len(padding))]
                      for r in [lower_remove, upper_remove]]

    logging.info("Cropped img to %s", lims)

    if np.any(np.array(padded_removes) < 0):
        raise ValueError("Padded removes " + str(padded_removes) +
                         " are invalid.")

    return (sitk.Crop(img, *padded_removes),
            {"origin": (lims[0][0], lims[1][0], lims[2][0]),
             "padding": padding})


def distribute_seeds(img, n_pts=100):
    '''Randomly distribute n seeds amongst all points where img != 0'''
    import random

    array = sitk.GetArrayFromImage(img)

    seeds = list()
    while len(seeds) < n_pts:
        (z, y, x) = [random.randrange(0, i) for i in array.shape]

        if array[z, y, x] != 0 and (z, y, x) not in seeds:
            seeds.append((x, y, z))

    return seeds


@cached(relevant_opts=["anisodiff", "gauss"])
def aniso_gauss(img_in, options):
    '''CurvatureAnisotropicDiffusion + GradientMagnitudeRecursiveGaussian is a
    a common featurization strategy. Compute these for consumption by other
    sitkstrat functions.'''

    img = sitk.Cast(img_in, sitk.sitkFloat32)

    img = sitk.CurvatureAnisotropicDiffusion(
        img,
        timeStep=options['anisodiff']['timestep'],
        conductanceParameter=options['anisodiff']['conductance'],
        # options['anisodiff'].setdefault('scaling_interval', 1),
        numberOfIterations=options['anisodiff']['iterations'])

    img = sitk.GradientMagnitudeRecursiveGaussian(
        img,
        options['gauss']['sigma'])

    return (img, options)


@log_size
@options_log
def segmentation_union(imgs, options):
    '''Compute a consensus segmentation amongst a small set of segmentations'''
    # pylint: disable=E1101
    # Sadly, images and arrays have a different coordinate system (z, y, x) vs
    # (x, y, z) so it's safest just to convert here. Don't worry, it's fast.
    incl_count = np.zeros(sitk.GetArrayFromImage(imgs[0]).shape)
    n_img = 0

    # This would be memory-expensive for large numbers of images, but it's only
    # a couple so ith's hopefully ok.
    for img in [sitk.GetArrayFromImage(i) for i in imgs]:
        assert img.shape == incl_count.shape

        img_size = np.count_nonzero(img)
        if img_size < options['max_size'] and \
           img_size > options['min_size']:
            n_img += 1
            incl_count += (img != 0)

    # store the number of images that passed QC
    options['n_imgs'] = n_img

    if n_img == 0:
        raise RuntimeWarning("No images satisifed image size thresholds" +
                             str((options['min_size'], options['max_size'])))

    # sitk is pretty particular about the datatypes that come in to
    # GetImageFromArray, and the default output from the following (bool?)
    # isn't acceptable apparently
    consensus = np.array(incl_count >= options['threshold'] * n_img,
                         dtype='uint8')

    consensus_size = np.count_nonzero(consensus)
    if consensus_size < options['min_size']:
        raise RuntimeWarning("Consensus image failed size threshold.  " +
                             "Image too small at " + str(consensus_size))

    consensus = sitk.GetImageFromArray(consensus)
    consensus.CopyInformation(imgs[0])
    consensus = sitk.Cast(consensus, sitk.sitkUInt8)

    return (consensus, options)


@log_size
@options_log
def segment_lung(img, options):
    '''Produce a lung segmentation from an input image.'''
    orig = img
    img  = lungseg.otsu(img)
    img  = lungseg.ensure_img_border(img)
    img  = lungseg.find_components(img)
    img  = lungseg.isolate_lung_field(img)
    img  = lungseg.dialate(img, options['probe_size'])
    img  = lungseg.find_components(img)
    img  = lungseg.isolate_not_biggest(img)
    img  = sitk.BinaryErode(img, options['probe_size']-2,
                           sitk.BinaryErodeImageFilter.Ball)
    
    return (img, options)


@options_log
def curvature_flow(img_in, options):

    img = sitk.Cast(img_in, sitk.sitkFloat32)
    img = sitk.CurvatureFlow(
        img,
        options['curvature_flow']['timestep'],
        options['curvature_flow']['iterations'])

    return (img, options)


@log_size
@options_log
def confidence_connected(img_in, options):
    img = sitk.ConfidenceConnected(
        img_in,
        [options['seed']],
        options['conf_connect']['iterations'],
        options['conf_connect']['multiplier'],
        options['conf_connect']['neighborhood'])

    img = sitk.BinaryDilate(
        img,
        options['dialate']['radius'],
        sitk.BinaryDilateImageFilter.Ball)

    return (img, options)


@options_log
def aniso_gauss_watershed(img_in, options_in):
    '''Compute CurvatureAnisotropicDiffusion +
    GradientMagnitudeRecursiveGaussian + Sigmoid featurization of the image.'''

    (img, options) = aniso_gauss(img_in, options_in)

    img = sitk.MorphologicalWatershed(
        img,
        level=options['watershed']['level'],
        markWatershedLine=True,
        fullyConnected=False)

    return (img, options)


@log_size
@options_log
def isolate_watershed(img_in, options):
    '''Isolate a particular one of the watershed segmentations.'''
    seed = options['seed']

    arr = sitk.GetArrayFromImage(img_in)

    label = arr[seed[2], seed[1], seed[0]]
    options['label'] = int(label)

    lab_arr = np.array(arr == label, dtype='float32')  # pylint: disable=E1101

    out_img = sitk.GetImageFromArray(lab_arr)
    out_img.CopyInformation(img_in)

    return (out_img, options)


@options_log
def aniso_gauss_sigmo(img_in, opts_in):
    '''Compute CurvatureAnisotropicDiffusion +
    GradientMagnitudeRecursiveGaussian + Sigmoid featurization of the image.'''

    (img, options) = aniso_gauss(img_in, opts_in)

    img = sitk.Sigmoid(
        img,
        options['sigmoid']['alpha'],
        options['sigmoid']['beta'])

    return (img, options)


@log_size
@options_log
def fastmarch_seeded_geocontour(img_in, options):
    '''Segment img_in using a GeodesicActiveContourLevelSetImageFilter with an
    inital level set built using FastMarchingImageFilter at options['seed']'''

    # The speed of wave propagation should be one everywhere, so we produce
    # an appropriately sized np array of all ones and convert it into an img
    ones_img = sitk.GetImageFromArray(np.ones(  # pylint: disable=E1101
        sitk.GetArrayFromImage(img_in).shape))
    ones_img.CopyInformation(img_in)

    fastmarch = sitk.FastMarchingImageFilter()

    # to save time, we limit the distances we calculate to a quarter of the
    # image size away (i.e. a region no more than half the image in diameter).
    fastmarch.SetStoppingValue(max(img_in.GetSize())*0.25)
    seeds = sitk.VectorUIntList()
    seeds.append(options['seed'])
    fastmarch.SetTrialPoints(seeds)

    seed_img = fastmarch.Execute(ones_img)

    # FastMarch won't output the right PixelType, so we have to cast.
    seed_img = sitk.Cast(seed_img, img_in.GetPixelID())

    # Generally speaking, you're supposed to subtract an amount from the
    # input level set, so that growing algorithm doesn't need to go as far
    img_shifted = sitk.GetImageFromArray(
        sitk.GetArrayFromImage(seed_img) - options['seed_shift'])
    img_shifted.CopyInformation(seed_img)
    seed_img = img_shifted

    geodesic = sitk.GeodesicActiveContourLevelSetImageFilter()
    geodesic.SetPropagationScaling(
        options['geodesic']['propagation_scaling'])
    geodesic.SetNumberOfIterations(
        options['geodesic']['iterations'])
    geodesic.SetCurvatureScaling(
        options['geodesic']['curvature_scaling'])
    geodesic.SetMaximumRMSError(
        options['geodesic']['max_rms_change'])

    out = geodesic.Execute(seed_img, img_in)

    options['geodesic']['elapsed_iterations'] = geodesic.GetElapsedIterations()
    options['geodesic']['rms_change'] = geodesic.GetRMSChange()

    out = sitk.BinaryThreshold(out, insideValue=0, outsideValue=1)

    return (out, options)
