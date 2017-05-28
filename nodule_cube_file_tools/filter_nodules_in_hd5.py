'''
Extract a subset of nodules from a hd5 container, by attribute and 
conditional expression.

Run with --help for usage info.
'''

from __future__ import print_function
import sys, os, h5py

attributes        = [   
                        'images', 'classes', 'malignancy', 
                        'inclusion', 'is_nodule', 'patients', 
                        'ids', 'pixel_min', 'pixel_max'
                    ]

allowed_relations = ['>', '<', '>=', '<=', '==', '!=', 'in', 'not-in', 'in-file', 'not-in-file']

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_args(do_usage=False):
    import argparse
    # construct the argument parser and parse the arguments
    ap = argparse.ArgumentParser(prog='{0}'.format(os.path.basename(sys.argv[0])), description=__doc__)
    ap.add_argument('-a', dest = 'attribute', metavar='ATTRIBUTE', default='classes', 
        help="Select attribute from {}.".format(str(attributes)))
    ap.add_argument('-c', dest = 'condition', metavar='CONDITION', required = False, 
        help = 'Extract nodules whose attribute satisfies this condition.  CONDITION must be of the form "OP VALUE" where'\
               'OP is one of {} and VALUE is a valid value for that attribute'.format(
                    str(allowed_relations)))
    ap.add_argument('--dry-run', dest = 'dry_run', action='store_true', help = 'show extracted nodule set, but do not perform extraction')
    ap.add_argument('-v', '--verbose', dest = 'verbose', action='store_true', help = 'show nodules while extracting')
    ap.add_argument('data_file', metavar='data-filename' , help = 'name of hd5 input data file')
    ap.add_argument('output_file', metavar='output-filename', nargs='?', help = 'name of hd5 output data file (optional only in --dry-run mode)')
    if(do_usage):
        ap.print_help()
        return None
    args = vars(ap.parse_args())
    return args

def attr_key(attribute):
    return 'nodule_{}'.format(attribute)

def check_condition(stored_value, condition, cond_value):
    if condition is not None and not (condition in allowed_relations):
        raise(RuntimeError('Condition {} not allowed; valid conditions are: {}'.format(condition, allowed_relations)))    
    condition_str = 'True'
    if condition is not None:
        condition  = ' '.join(condition.split('-'))
        cond_value = eval(cond_value,{},{})  # Eval with no global or local context!
        if type(cond_value) == type([]):
            cond_value = [type(stored_value)(x) for x in cond_value]
        else:
            cond_value = type(stored_value)(cond_value)
        if type(stored_value) == type('string'):
            stored_value = "'{}'".format(stored_value)
        condition_str = '{} {} {}'.format(stored_value, condition, cond_value)
    return eval(condition_str,{},{}) # Eval with no global or local context!

def print_nodule(nodule, attribute='ids'):
    print('{}\t{}\t{}\t{}'.format(
        nodule[attr_key('ids')], 
        nodule[attr_key('classes')], 
        nodule[attr_key('malignancy')], 
        nodule[attr_key(attribute)] if attribute not in ['ids', 'classes', 'malignancy'] else ''
    ))

def print_nodule_header(attribute, condition=None, cond_value=None):
    print('ID                               \tclass\tmalignancy\t{}'.format(attribute if attribute not in ['ids', 'classes', 'malignancy'] else ''))

def print_nodules(hd5file, attribute, condition=None, cond_value=None):
    hfile        = h5py.File(hd5file, 'r')
    types        = {}
    shapes       = {}
    in_place     = False
    
    print_nodule_header(attribute, condition, cond_value)
    added = 0
    for idx in range(len(hfile[attr_key(attribute)])):
        stored_value  = hfile[attr_key(attribute)][idx][0]
        nodule        = {k: hfile[k][idx] for k in hfile.keys()}
        if check_condition(stored_value, condition, cond_value):
            print_nodule(nodule, attribute)
            added += 1
    hfile.close()
    print('{} nodules'.format(added))

def extract_nodules(hd5file, outfile, attribute, condition=None, cond_value=None, verbose=False):
    hfile        = h5py.File(hd5file, 'r')
    types        = {}
    shapes       = {}
    in_place     = False
    if outfile is None:
        outfile = hd5file + '.tmp'
        in_place = True
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

    fout = h5py.File(outfile, 'w')
    for k in hfile.keys():
        max_shape   = tuple(list((None,)) + list(shapes[k][1:]))
        basic_shape = tuple(list((0,))    + list(shapes[k][1:]))
        fout.create_dataset(k, basic_shape, maxshape=max_shape, dtype=types[k])
    if verbose:
        print_nodule_header(attribute, condition, cond_value)
    added = 0
    for idx in range(len(hfile[attr_key(attribute)])):
        stored_value = hfile[attr_key(attribute)][idx][0]
        nodule       = {k: hfile[k][idx] for k in hfile.keys()}
        if check_condition(stored_value, condition, cond_value):
            if verbose:
                print_nodule(nodule, attribute)
            for k in nodule.keys():
                fout[k].resize(len(fout[k])+1, axis=0)
                fout[k][added] = nodule[k]
            added += 1
    fout.close()
    hfile.close()
    if verbose:
        print('{} nodules'.format(added))
    if in_place:
        try: # Rename the temp file over the original if we are working "in-place".
            os.rename(outfile, hd5file)
        except WindowsError:
            os.remove(hd5file)
            os.rename(outfile, hd5file)

def run_standalone():
    args       = get_args()
    hd5file    = args['data_file']
    condition  = args['condition'].split() if args['condition'] is not None else None
    cond_value = ' '.join(condition[1:])   if condition is not None else None
    condition  = condition[0].strip() if condition is not None else None

    if condition in ['in-file', 'not-in-file']:   # The "in-file" options must supply a file containing the value list
        with open(cond_value, 'rU') as list_file: # format must be whitespace-separated values
            list_values = list_file.read().split()
            try:                                  # Handle both numeric and string list values (float, int, str)
                list_values     = [float(x) for x in list_values]
                list_values_int = [int(x) for x in list_values]
                if all([abs(a-b) < 1e-10 for a,b in zip(list_values, list_values_int)]):
                    list_values = list_values_int
            except Exception as e:
                list_values = ['{}'.format(x) for x in list_values]
            cond_value = str(list_values)
        condition = condition[:-5]

    if args['output_file'] is None and not args['dry_run']:
        print("You must provide an output file name unless you use the --dry-run flag.", file=sys.stderr);
        get_args(True)
        sys.exit(1)
    # try:
    if args['dry_run']:
        print_nodules(hd5file, args['attribute'], condition, cond_value)
    else:
        extract_nodules(hd5file, args['output_file'], args['attribute'], condition, cond_value, verbose=args['verbose'])
    # except Exception as e:
    #     print("FAILED: \n{}".format(str(e)), file=sys.stderr)
    #     sys.exit(1)

if __name__ == '__main__':
    run_standalone()
