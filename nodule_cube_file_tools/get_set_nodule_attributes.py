'''
Gets or sets an attribute for nodules in a hd5 data archive.

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
    ap.add_argument('-s', dest = 'set_value', metavar='SET-TO-VALUE', required = False, help = 'Value to set all matching attributes to.')
    ap.add_argument('-c', dest = 'condition', metavar='CONDITION', required = False, 
        help = 'Change value only if attribute satisfies this condition.  CONDITION must be of the form "OP VALUE" where'\
               'OP is one of {} and VALUE is a valid value for that attribute'.format(
                    str(allowed_relations)))
    ap.add_argument('data_file', metavar='data-filename' , help = 'name of hd5 data file')
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

def change_attribute(hd5file, attribute, to_value, condition=None, cond_value=None):
    hfile      = h5py.File(hd5file)
    for idx in range(len(hfile[attr_key(attribute)])):
        stored_value  = hfile[attr_key(attribute)][idx][0]
        if check_condition(stored_value, condition, cond_value):
            hfile[attr_key(attribute)][idx] = [type(hfile[attr_key(attribute)][idx][0])(to_value)]    
    hfile.close()

def print_attribute(hd5file, attribute, condition=None, cond_value=None):
    hfile = h5py.File(hd5file, 'r')
    for idx in range(len(hfile[attr_key(attribute)])):
        stored_value  = hfile[attr_key(attribute)][idx][0]
        if check_condition(stored_value, condition, cond_value):
            print('{}'.format(hfile[attr_key(attribute)][idx][0]))
    hfile.close()

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

    try:
        if args['set_value'] is not None:
            change_attribute(hd5file, args['attribute'], args['set_value'], condition, cond_value)
        else:
            print_attribute(hd5file, args['attribute'], condition, cond_value)
    except Exception as e:
        print("FAILED: \n{}".format(str(e)), file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    run_standalone()
