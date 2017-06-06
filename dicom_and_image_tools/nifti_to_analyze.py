"""
Converts a file in Nifti format (single .nii file) to the 
Analyze 7.5 (.img / .hdr) format.
"""
import sys, os
_usage = """
Usage:
    {0} nifti-filename [analyze-filename]

        nifti-filename      name of Nifti file to convert
        analyze-filename    name of the analyze file (without extension)
                            will produce a .img/.hdr pair with this name
                            prefix  (default - same name as nift file)
""".format(os.path.basename(sys.argv[0]))
from medpy.io import load, save
import numpy as np

def usage():
    print(_usage)

def convert_to_analyze(nifti_file, analyze_file=None):
    if analyze_file is None:
        analyze_file = nifti_file
    analyze_file = '{0}.img'.format(os.path.splitext(analyze_file)[0])
    image, header = load(nifti_file)
    # Make sure the type of the image file is compatible with the matlab analyze75read function
    # if not convert to int16
    if image.dtype != np.int16 and image.dtype != np.int32:
        image = image.astype('int16')
    save(image, analyze_file, header, force=True)

def main():
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)
    nifti_file    = sys.argv[1]
    analyze_file  = sys.argv[2] if len(sys.argv) > 2 else None
    convert_to_analyze(nifti_file, analyze_file)

if __name__ == '__main__':
    main()

