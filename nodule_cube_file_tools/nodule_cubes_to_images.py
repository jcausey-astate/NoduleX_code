'''
Given an HDF5 file with nodule cubes, output the corresponding images
in Nifti or PNG format.

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
    ap.add_argument('--output-dir', dest = 'output_dir', metavar='DIR', required = False, help = 'Image output directory')
    ap.add_argument('--format', dest = 'output_type', metavar='FMT', required=False, default='nifti',
        help='Image format for output files, options are [nifti, dicom, analyze, png]; default is nifti.')
    ap.add_argument('--single', dest = 'single_slice', action='store_true', help="Output central slice only in PNG format.")
    ap.add_argument('data_file', metavar='data-filename' , help = 'name of hd5 data file')
    if(do_usage):
        ap.print_help()
        return None
    args = vars(ap.parse_args())
    return args

def run_standalone():
    args = get_args()
    hd5file = args['data_file']
    output_dir = args['output_dir']
    if args['normalize']:
        print("Normalizing output images.")
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), 'cube_images')
    out_ext = 'nii'
    if args['output_type'].lower()[0] == 'd':
        out_ext = 'dcm'
    elif args['output_type'].lower()[0] == 'a':
        out_ext = 'img'
    elif args['output_type'].lower()[0] == 'p':
        out_ext = 'png'

    if args['single_slice']:
        if out_ext != 'png':
            print("Output type changed to PNG for single slice output.")
        out_ext = 'png'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    elif not os.path.isdir(output_dir):
        print("Failed: '{}'' exists and is not a directory.".format(output_dir), file=sys.stderr)
        sys.exit(1)

    with h5py.File(hd5file, 'r') as fp:
        Xmin  = fp['nodule_pixel_min'].value
        Xmax  = fp['nodule_pixel_max'].value
        for idx, cube in enumerate(fp['nodule_images']):
            image     = np.array(cube, dtype='float32')
            image     = np.moveaxis(image, 0, -1)
            if args['normalize']:        
                # print("min: {} max: {}".format(Xmin[idx], Xmax[idx]))
                image = (image - Xmin[idx]) / (Xmax[idx] - Xmin[idx])
                # Scale back into a common range 0 - 4095
                image = image * 4095
                # print("  - min {}  - max {}".format(image.min(), image.max()))

            nodule_id = str(fp['nodule_ids'][idx][0])
            nodule_id = '_'.join(nodule_id.split())  # remove any whitespace from noduleID
            out_file  = '{0}-{1}.{2}'.format(os.path.splitext(hd5file)[0], nodule_id, out_ext)
            if args['output_type'] != 'png':
                # Make sure the type of the image file is compatible with the matlab analyze75read function
                # if not convert to int16
                if image.dtype != np.int16 and image.dtype != np.int32:
                    image = image.astype('int16')
                save(image, os.path.join(output_dir, out_file), force=True)
            else:
                image = np.moveaxis(image, -1, 0)
                if args['single_slice']:
                    image = np.array([image[image.shape[0] // 2]])
                for idx, img_out in enumerate(image):
                    img_out = img_out - img_out.min()
                    img_out = img_out / img_out.max()
                    img_out = img_out * 255
                    img_out = Image.fromarray(img_out)
                    img_out = img_out.convert('RGB')
                    parts   = os.path.splitext(out_file)
                    slice_name = "{}{}".format(parts[0], parts[1])
                    if not args['single_slice']:
                        slice_name = "{}_{}{}".format(parts[0], (idx+1), parts[1])
                    img_out.save(os.path.join(output_dir, slice_name), optimize=False, compress_level=0)

if __name__ == '__main__':
    run_standalone()
