"""
Converts a set of images in DICOM format (directory of .dcm files) to the 
Analyze 7.5 (.img / .hdr) format.
"""
import sys, os
_usage = """
Usage:
    {0} dicom-directory [analyze-filename]

        dicom-directory     name of directory containing DICOM files
        analyze-filename    name of the analyze file (without extension)
                            will produce a .img/.hdr pair with this name
                            prefix  (default - same prefix as DICOM direcotry)
""".format(os.path.basename(sys.argv[0]))
from nifti_to_analyze import convert_to_analyze
from dicom_to_nifti import dicom_dir_to_nifti_file
from tempfile import mkdtemp
from shutil import rmtree

def usage():
    print(_usage)

def convert_dicom_to_analyze(dicom_dir, analyze_file=None):
    if analyze_file is None:
        analyze_file = dicom_dir
        if analyze_file[-1] == os.path.sep:
            analyze_file = analyze_file[:-1]
    if analyze_file[-1] == os.path.sep:
        analyze_file += dicom_dir if dicom_dir[-1] != os.path.sep else dicom_dir[:-1]
    analyze_file = '{0}.img'.format(analyze_file)
    tmp_dir      = mkdtemp()
    nifti_file,_ = dicom_dir_to_nifti_file(dicom_dir, tmp_dir)
    # Create directory for analyze file if necessary
    try:
        os.makedirs(os.path.dirname(analyze_file))
    except OSError:
        # no need to do anything if the directory already exists.
        pass
    convert_to_analyze(nifti_file, analyze_file)
    rmtree(tmp_dir)

def main():
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)
    if sys.argv[1] in ['-h', '--help']:
        usage()
        sys.exit(0)
    dicom_dir    = sys.argv[1]
    analyze_file = sys.argv[2] if len(sys.argv) > 2 else None
    convert_dicom_to_analyze(dicom_dir, analyze_file)

if __name__ == '__main__':
    main()
