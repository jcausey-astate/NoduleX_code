'''
Sort nodules from hd5 file by an attribute (stable sort).

Run with --help for usage info.
'''

from __future__ import print_function
import sys, os, h5py

attributes        = [   
                        'images', 'classes', 'malignancy', 
                        'inclusion', 'is_nodule', 'patients', 
                        'ids', 'pixel_min', 'pixel_max'
                    ]

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_args(do_usage=False):
    import argparse
    # construct the argument parser and parse the arguments
    ap = argparse.ArgumentParser(prog='{0}'.format(os.path.basename(sys.argv[0])), description=__doc__)
    ap.add_argument('-a', dest = 'attribute', metavar='ATTRIBUTE', default='ids', 
        help="Sort on attribute from {}. (default: ids)".format(str(attributes)))
    ap.add_argument('-r', '--reverse', dest = 'reverse', action='store_true', help = 'reverse the sort ordering')
    ap.add_argument('-v', '--verbose', dest = 'verbose', action='store_true', help = 'show nodules while sorting')
    ap.add_argument('data_file', metavar='data-filename' , help = 'name of hd5 input data file')
    ap.add_argument('output_file', metavar='output-filename', nargs='?', help = 'name of hd5 output data file (if omitted, data will be changed in-place)')
    if(do_usage):
        ap.print_help()
        return None
    args = vars(ap.parse_args())
    return args

def attr_key(attribute):
    return 'nodule_{}'.format(attribute)

def print_nodule(nodule, attribute='ids'):
    print('{}\t{}\t{}\t{}'.format(
        nodule[attr_key('ids')], 
        nodule[attr_key('classes')], 
        nodule[attr_key('malignancy')], 
        nodule[attr_key(attribute)] if attribute not in ['ids', 'classes', 'malignancy'] else ''
    ))

def sort_nodules(hd5file, outfile, attribute, reverse=False, verbose=False):
    hfile        = h5py.File(hd5file, 'r')
    types        = {}
    shapes       = {}
    in_place     = False
    if outfile is None:
        outfile = hd5file + '.tmp'
        in_place = True
    
    output_count = len(hfile[attr_key(attribute)])

    for key in hfile.keys():
        types[key]  = type(hfile[key][0])
        try:
            types[key] = hfile[key][0].dtype
        except:
            pass
        shapes[key] = (1,)
        try:
            shapes[key] = hfile[key].shape
        except:
            pass

    values_to_sort = []
    for idx, value in enumerate(hfile[attr_key(attribute)]):
        values_to_sort.append((value, idx))
    sorted_indices = [x[1] for x in sorted(values_to_sort, reverse=reverse, key=lambda x: x[0])]
    values_to_sort = None

    fout = h5py.File(outfile, 'w')
    for k in hfile.keys():
        shape = tuple(list((output_count,)) + list(shapes[k][1:]))
        fout.create_dataset(k, shape, dtype=types[k])
    if verbose:
        print('ID                               \tclass\tmalignancy\t{}'.format(attribute if attribute not in ['ids', 'classes', 'malignancy'] else ''))
    added = 0
    for count_idx in range(len(hfile[attr_key(attribute)])):
        idx = sorted_indices[count_idx]
        nodule = {k: hfile[k][idx] for k in hfile.keys()}
        if verbose:
            print_nodule(nodule, attribute)
        for k in nodule.keys():
            fout[k][added] = nodule[k]
        added += 1
    fout.close()
    hfile.close()
    if verbose:
        print('{} nodules'.format(added))
    if in_place:
        try: # Rename the temp file over the original if we are sorting "in-place".
            os.rename(outfile, hd5file)
        except WindowsError:
            os.remove(hd5file)
            os.rename(outfile, hd5file)

def run_standalone():
    args       = get_args()
    sort_nodules(args['data_file'], args['output_file'], args['attribute'], reverse=args['reverse'], verbose=args['verbose'])
    
if __name__ == '__main__':
    run_standalone()
