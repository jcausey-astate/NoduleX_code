"""
Adapted from Justin Porter's (https://github.com/justinrporter/nodule-seg)
'masterseg.py' to focus specifically on segmenting the lung
volume from the rest of the CT image.
Takes and returns images as numpy 3D arrays.
"""
from segment import lungseg
import SimpleITK as sitk, numpy as np

def binary_erode(img, probe_size):
    '''Erodes a binary image.'''
    filt = sitk.BinaryErodeImageFilter()
    filt.SetKernelType(filt.Ball)
    filt.SetKernelRadius(probe_size)

    return filt.Execute(img, 0, 1, False)

def segment_lung_volume(img, probe_size = 7, erode=0):  # pylint: disable=C0111
    '''Segment the image to isolate the lung volume.  Returns binary segmentation
    image where 0 = not lung and 1 = lung volume.
    ''' 
    lung_img = lungseg.segment_lung(sitk.GetImageFromArray(img), {'probe_size': probe_size})
    if erode is not None and erode > 0:
        lung_img = binary_erode(lung_img, erode)
    lung_img = sitk.GetArrayFromImage(lung_img)
    return lung_img

def apply_segmentation_to_image(img, seg_image, reverse_mask=False, bg_value=0):
    '''Apply segmentation mask to image, returning the portion not covered by the
    segmentation (or the opposite, if `reverse_mask` is True)'''
    img  = np.array(img)
    mask = 0 if not reverse_mask else 1
    img[seg_image == mask] = bg_value
    return img;

def isolate_lungs(img, probe_size = 7, remove_lungs=False, erode=0, bg_value=0):
    '''Subtract all non-lung areas from the image (unless `remove lungs` is 
    True, then subtract the lungs from the image.'''
    seg_mask = segment_lung_volume(img, probe_size=probe_size, erode=erode)
    return apply_segmentation_to_image(img, seg_mask, reverse_mask=remove_lungs, bg_value=bg_value)

def find_lung_bounds(img, probe_size = 7, ignore_trachea=True):
    '''Returns the bounding cube as three tuples representing the min and max
    bounds along each axis in `img`.'''
    import segment.bounding as bnd
    from   scipy.ndimage.measurements import label
    
    if not is_binary_image(img):
        img = segment_lung_volume(img, probe_size)
    
    bounds = bnd.bounding_cube(img)
    
    if ignore_trachea:
        # The trachea causes trouble at the "top" - ignore slices that only contain
        # one object at the top.  Note: This will remove some upper lung slices in some cases.
        first  = None
        slices = bounds[0][1] - bounds[0][0]
        for i in range(bounds[0][0], bounds[0][1]):
            labeled, count = label(img[i])
            if first is None and count > 1:  # If we have at least 2 objects, one must be a lung.
                first = i
                break;
        
        bounds = [list(ax) for ax in bounds]
        bounds[0][0] = first if first is not None else bounds[0][0]
        bounds = tuple(tuple(t) for t in bounds)

    return bounds

def is_binary_image(img):
    '''Determine if an image (numpy array form) is binary (only 0 or 1) or not.'''
    return np.array_equal(img, img.astype(bool))

if __name__ == '__main__':
    print("This script is not designed to run in standalone mode.")
    exit(1)