'''
Given an HDF5 file with nodule cubes, output the corresponding images
in a single image mosaic in e.g. PNG format.

Run with --help for usage info.
'''

from __future__ import print_function
import sys, os, h5py, numpy as np
from PIL import Image
from medpy.io import save

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_args(do_usage=False):
    import argparse
    # construct the argument parser and parse the arguments
    ap = argparse.ArgumentParser(prog='{0}'.format(os.path.basename(sys.argv[0])), description=__doc__)
    ap.add_argument('-n', '--normalize', dest='normalize', action='store_true', help='Normalize cubes before export.')
    ap.add_argument('-W', '--window', dest='window_normalize', action='store_true', help='HU Window normalize cubes before export.')
    ap.add_argument('-w', '--width', dest='width', default=20, type=int, help='Number of images in each output image row; sets the width of the mosaic. Set to 0 for best square fit.')
    ap.add_argument('-l', '--limit', dest='limit', default=400, type=int, help='Maximum number of cubes to include.  If fewer than available number, results will be randomized. Set to 0 for no limit.')
    ap.add_argument('data_file', metavar='data-filename' , help = 'name of hd5 data file')
    ap.add_argument('output_file', metavar='output-filename' , help = 'name of output file')
    if(do_usage):
        ap.print_help()
        return None
    args = vars(ap.parse_args())
    return args

def run_standalone():
    args = get_args()
    hd5file = args['data_file']
    outfile = args['output_file']

    if args['normalize']:
        print("Normalizing output images.")
    if args['window_normalize']:
        print("Window Normalizing output images.")
    if outfile[-4:] != '.png':
        outfile = os.path.splitext(outfile)[0] + '.png'
    img_out = None
    with h5py.File(hd5file, 'r') as fp:
        Xmin  = fp['nodule_pixel_min'].value
        Xmax  = fp['nodule_pixel_max'].value
        available_count = len(Xmin)
        chosen_nodules  = range(available_count)
        if args['limit'] > 0 and available_count > args['limit']:
            chosen_nodules = np.random.choice(chosen_nodules, args['limit'], replace=False)
        cube_shape = fp['nodule_images'][0].shape
        n_images   = len(chosen_nodules)
        if args['width'] == 0:
            args['width'] = int(np.ceil(np.sqrt(n_images)))
        n_rows     = int(np.ceil(float(n_images) / args['width']))
        cube_w     = cube_shape[1]
        cube_h     = cube_shape[2]
        cube_mid_z = int(cube_shape[0] / 2)
        img_out  = np.zeros((cube_h * n_rows, cube_w * args['width']), dtype='float32')
        for count_idx, idx in enumerate(chosen_nodules):
            image     = fp['nodule_images'][idx]
            image     = np.array(image, dtype='float32')
            image     = image[cube_mid_z]
            if args['normalize'] or args['window_normalize']:
                nXmin = Xmin[idx] if not args['window_normalize'] else -1000.0
                nXmax = Xmax[idx] if not args['window_normalize'] else  4096.0
                image = (image - nXmin) / (nXmax - nXmin)
                image[image < 0] = 0
                image[image > 1] = 1
                image = image * 255  # Scale to 0-255 range
            row = count_idx / args['width']
            col = count_idx % args['width']
            row_px = row * cube_h
            col_px = col * cube_w
            img_out[row_px : row_px + cube_h, col_px : col_px + cube_w] = image
    img_out = img_out - img_out.min()
    img_out = img_out / img_out.max()
    img_out = img_out * 255
    img_out = Image.fromarray(img_out)
    img_out = img_out.convert('RGB')
    img_out.save(outfile, optimize=False, compress_level=0)

if __name__ == '__main__':
    run_standalone()
