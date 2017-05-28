'''
List (text) info for nodules in hd5 data file.

Run with --help for usage info.
'''

from __future__ import print_function
import sys, os, h5py

allowed_attributes = [   
                        'classes', 'malignancy', 
                        'inclusion', 'is_nodule', 'patients', 
                        'ids', 'pixel_min', 'pixel_max',
                        'rescale_intercept', 'rescale_slope', 'shape'
                     ]

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_args(do_usage=False):
    import argparse
    # construct the argument parser and parse the arguments
    ap = argparse.ArgumentParser(prog='{0}'.format(os.path.basename(sys.argv[0])), description=__doc__)
    ap.add_argument('-a', dest = 'attributes', metavar='ATTRIBUTES', nargs='*', default=['ids', 'classes', 'malignancy'], 
        help="List of attributes to include. (default: ids, classes, malignancy)".format(str(allowed_attributes)))
    ap.add_argument('data_file', metavar='data-filename' , help = 'name of hd5 input data file')
    if(do_usage):
        ap.print_help()
        return None
    args = vars(ap.parse_args())
    return args

def attr_key(attribute):
    return 'nodule_{}'.format(attribute)

def print_nodule(nodule, attributes=['ids', 'classes', 'malignancy']):
    nodule[attr_key('shape')] = [nodule[attr_key('images')].shape]
    available  = list(set(nodule.keys()).intersection(set([attr_key(k) for k in allowed_attributes])))
    fmt_string = '\t'.join(['{}' for a in attributes if attr_key(a) in available])
    print(fmt_string.format(*[nodule[attr_key(k)][0] for k in attributes if attr_key(k) in available]))

def run_standalone():
    args  = get_args()
    hfile = h5py.File(args['data_file'], 'r')
    for idx in range(len(hfile[attr_key('ids')])):
        nodule = {k: hfile[k][idx] for k in hfile.keys()}
        print_nodule(nodule, args['attributes'])

if __name__ == '__main__':
    run_standalone()
