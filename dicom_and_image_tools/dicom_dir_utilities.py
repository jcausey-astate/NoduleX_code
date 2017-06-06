"""
Some common functions on DICOM directories.
"""
import sys, os, numpy as np

def get_dicomdir_files(dicom_dir, include_xml = False, xmldir=None, full_path=False):
    files       = []
    dicomfiles  = []
    xmlfiles    = []
    dicom_dir = os.path.normpath(dicom_dir)
    if os.path.isdir(dicom_dir):
        for rt, dr, fl in os.walk(dicom_dir):
            files.extend(fl)
            break
    else:
        files     = [dicom_dir]
        dicom_dir = os.path.abspath(os.path.join(dicom_dir, os.pardir))
        if xmldir == None:
            for rt, dr, fl in os.walk(dicom_dir):
                for fname in fl:
                    if fname[-4:] == '.xml':
                        if len(files) == 2:
                            eprint("Ambiguous xml files found.")
                        files.append(os.path.join(dicom_dir, fname))
                break

    dicomfiles = [f for f in files if f[-4:] == '.dcm']
    xmlfiles   = [f for f in files if f[-4:] == '.xml']
    if full_path:
        dicomfiles = [os.path.join(dicom_dir, f) for f in dicomfiles]
        xmlfiles   = [os.path.join(dicom_dir, f) for f in xmlfiles]
    return (dicomfiles, xmlfiles) if include_xml else dicomfiles


if __name__ == '__main__':
    print("This file is designed to be included as a library only.")
    sys.exit(1)