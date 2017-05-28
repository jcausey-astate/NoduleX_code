"""
Usage:  {0} output_file input_file1 input_file2 [... input_fileN]

    Concatenates images from all input pickle or hdf5 files into a single output pickle or hdf5 file.
        Note: If the output file is not hdf5, all input files must be pickle-files.  
              Any combination of pickle and hdf5 input files may be used when the output file is hdf5.
"""

import pickle, numpy as np
import sys

def usage(exit_code=None):
    print(__doc__.format(sys.argv[0].split('/')[-1]))
    if exit_code is not None:
        sys.exit(exit_code)


def file_ext(filename):
    return filename.split('.')[-1] if filename is not None else None

def is_hdf5(filename):
    hdf_ext_list  = ['hd5', 'hdf5', 'hdf', 'h5']
    return file_ext(filename) in hdf_ext_list

def run_standalone():
    if len(sys.argv) < 2:
        print(__doc__.format(sys.argv[0].split('/')[-1]))
        usage(1)

    out_file_name = sys.argv[1]
    
    if out_file_name in ['-h', '--help']:
        usage(0)

    if not is_hdf5(out_file_name):
        print("Preparing pickle output.")
        output_list   = [[],[]]

        for in_file_name in sys.argv[2:]:
            with open(in_file_name, 'rb') as fin:
                X,y = pickle.load(fin)
                output_list[0].extend(X)
                output_list[1].extend(y)

        with open(out_file_name, 'wb') as fout:
            pickle.dump(output_list, fout, pickle.HIGHEST_PROTOCOL)
    else:
        import h5py
        with h5py.File(out_file_name, 'w') as fout:
            print("Preparing HDF5 output.")
            images_initialized = False
            for in_file_name in sys.argv[2:]:
                hd5_datasets = {}
                if not is_hdf5(in_file_name):
                    with open(in_file_name, 'rb') as fin:
                        X,y = pickle.load(fin)
                        hd5_datasets['nodule_images']  = X
                        hd5_datasets['nodule_classes'] = y
                else:
                    with h5py.File(in_file_name, 'r') as fin:
                        file_keys = fin.keys()
                        for key in file_keys:
                            hd5_datasets[key]  = fin[key].value
                hd5_datasets['nodule_images']  = np.array(hd5_datasets['nodule_images'], 'float32')
                hd5_datasets['nodule_classes'] = np.array(hd5_datasets['nodule_classes'], 'int')
                length = len(hd5_datasets['nodule_classes'])
                if not images_initialized:
                    h5str    = h5py.special_dtype(vlen=bytes)
                    shape    = list(hd5_datasets['nodule_images'].shape)    
                    shape[0] = None # Let the array grow in the first dimension
                    fout.create_dataset('nodule_images', dtype='float32', maxshape=tuple(shape), data=hd5_datasets['nodule_images'])
                    fout.create_dataset('nodule_classes', dtype='int', maxshape=(None,1), data=hd5_datasets['nodule_classes'])
                    fout.create_dataset('nodule_malignancy',  maxshape=(None,1), dtype='int', shape=(0,1))
                    fout.create_dataset('nodule_inclusion',  maxshape=(None,1), dtype='int', shape=(0,1))
                    fout.create_dataset('nodule_patients',  maxshape=(None,1), dtype=h5str, shape=(0,1))
                    fout.create_dataset('nodule_ids',  maxshape=(None,1), dtype=h5str, shape=(0,1))
                    fout.create_dataset("nodule_pixel_min",  maxshape=(None, 1,), dtype='float32', shape=(0,1))
                    fout.create_dataset("nodule_pixel_max",  maxshape=(None, 1,), dtype='float32', shape=(0,1))
                    for key in ['nodule_malignancy', 'nodule_inclusion', 'nodule_patients', 'nodule_ids', 'nodule_pixel_min', 'nodule_pixel_max']:
                        try:
                            if key in hd5_datasets:
                                fout[key].resize(length, 0)
                                fout[key][:] = hd5_datasets[key]
                            else:
                                fout[key].resize(length, 0)
                                fout[key][:] = np.zeros(length, dtype=int) if fout[key].dtype == 'int' \
                                                                                else np.zeros(length, dtype='float32') if fout[key].dtype == 'float32' \
                                                                                    else ['Unknown' for i in range(length)]
                        except Exception as e:
                            print("Failed on key {} -- {}".format(key, e))
                            sys.exit(0)
                    images_initialized = True
                else: # Extend the dataset by resizing and then adding the new data at the end
                    current_length = len(fout['nodule_images'])
                    fout['nodule_images'].resize(current_length  + length, 0)
                    fout['nodule_classes'].resize(current_length + length, 0)
                    fout['nodule_images'][current_length:, ...]  = hd5_datasets['nodule_images']
                    fout['nodule_classes'][current_length:, ...] = hd5_datasets['nodule_classes']                                        
                    for key in ['nodule_malignancy', 'nodule_inclusion', 'nodule_patients', 'nodule_ids', 'nodule_pixel_min', 'nodule_pixel_max']:
                        if key in hd5_datasets:
                            fout[key].resize(current_length + length, 0)
                            fout[key][current_length:, ...] = hd5_datasets[key]
                        else:
                            fout[key].resize(current_length + length)
                            fout[key][:] = np.zeros(current_length + length, dtype=int) if fout[key].dtype == 'int' \
                                                                                else np.zeros(current_length + length, dtype='float32') if fout[key].dtype == 'float32' \
                                                                                    else ['Unknown' for i in range(current_length + length)]

if __name__ == '__main__':
    run_standalone()