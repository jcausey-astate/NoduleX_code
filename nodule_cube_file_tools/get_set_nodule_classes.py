'''
Gets or sets the classes for nodules in a hd5 data archive.
Can set all classes to the same value, or change one class value
to another, or map malignancy to class value.

Run with --help for usage info.
'''

from __future__ import print_function
import sys, os, h5py

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_args(do_usage=False):
    import argparse
    # construct the argument parser and parse the arguments
    ap = argparse.ArgumentParser(prog='{0}'.format(os.path.basename(sys.argv[0])), description=__doc__)
    ap.add_argument('-s', dest = 'set_class', metavar='CLASS', required = False, help = 'Value to set all classes to CLASS')
    ap.add_argument('-m', dest = 'malignancy_map', required = False, 
        nargs=6, metavar=('M0', 'M1', 'M2', 'M3', 'M4', 'M5'),
        help = 'Map malignancy score to class - provide class for all malignancy scores including 0 (not rated).')
    ap.add_argument('-c', dest = 'change_class', required = False, 
        nargs=2, metavar=('FROM', 'TO'),
        help = 'Change all classes that are currently equal to FROM to the class TO.')
    ap.add_argument('data_file', metavar='data-filename' , help = 'name of hd5 data file')
    if(do_usage):
        ap.print_help()
        return None
    args = vars(ap.parse_args())
    return args

def map_malignancies(hd5file, malignancy_map):
    malignancy_map = [[int(v)] for v in malignancy_map]
    hfile = h5py.File(hd5file)
    for idx in range(len(hfile['nodule_classes'])):
        hfile['nodule_classes'][idx] = malignancy_map[hfile['nodule_malignancy'][idx][0]]
    hfile.close()

def change_class(hd5file, from_class, to_class):
    from_class = int(from_class)
    to_class   = int(to_class)
    hfile      = h5py.File(hd5file)
    for idx in range(len(hfile['nodule_classes'])):
        if hfile['nodule_classes'][idx] == [from_class]:
            hfile['nodule_classes'][idx] = [to_class]    
    hfile.close()

def set_class_to(hd5file, to_class):
    to_class = int(to_class)
    hfile    = h5py.File(hd5file)
    for idx in range(len(hfile['nodule_classes'])):
        hfile['nodule_classes'][idx] = [to_class]    
    hfile.close()    

def print_class(hd5file):
    hfile = h5py.File(hd5file, 'r')
    for c in hfile['nodule_classes']:
        print('{}'.format(c[0]))
    hfile.close()        

def run_standalone():
    args = get_args()
    hd5file = args['data_file']
    if args['set_class'] is not None:
        set_class_to(hd5file, args['set_class'])
    elif args['malignancy_map'] is not None:
        map_malignancies(hd5file, args['malignancy_map'])
    elif args['change_class'] is not None:
        change_class(hd5file, from_class=args['change_class'][0], to_class=args['change_class'][1])
    else:
        print_class(hd5file)

if __name__ == '__main__':
    run_standalone()
