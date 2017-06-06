#!/usr/local/bin/python
# ---------------------------------------------------------------------------------------
# For every file in the directory, rename it so that the name contains the
# Instance Number prior to the extension.  This puts files in the order of 
# the scan position.
#    example.dcm  becomes  example_001.dcm   and similar
# Optionally, you can also rename by Patient ID:
#    example.dcm  becomes  Example-Patient-ID_001.dcm   and similar
#    
# Required library:  
#   pydicom (https://pypi.python.org/pypi/pydicom)
#       If you use Pip, you can install with `pip install pydicom`
# ---------------------------------------------------------------------------------------
# License: MIT (https://opensource.org/licenses/MIT)
# 
# The MIT License (MIT)
# Copyright (c) 2015 Jason L Causey, Arkansas State University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ---------------------------------------------------------------------------------------


from __future__ import print_function

try:
    import dicom
except ImportError, e:
    print("This script requires the pydicom library: https://pypi.python.org/pypi/pydicom")
    exit(1)
import argparse, os, sys

VERBOSE = False

def vprint(msg):
    if VERBOSE:
        print(msg)

parser = argparse.ArgumentParser(description='Add scan ordering numbers to DICOM image filenames.')
parser.add_argument('dir', type=str,
                   help='the directory containing one or more DICOM images to process')
parser.add_argument('--use-patient-id', dest='use_patient_id', action='store_true',
                   help='use Patient ID for the main file name instead of current name')
parser.add_argument('--verbose', action='store_true',
                   help='give more verbose output')
parser.set_defaults(verbose=False, use_patient_id=False)
args = parser.parse_args()

if not os.path.isdir(args.dir):
    print("You must supply a valid directory name.")
    parser.print_help()
    exit(1)
if args.dir[-1] != os.sep:
    args.dir = args.dir + os.sep;

VERBOSE = args.verbose

# for every DICOM (".dcm" or ".dicom") file in 'DIR', rename it by adding 
# the Instance Number before the extension
for filename in os.listdir(args.dir):
    file_base, file_extension = os.path.splitext(filename)
    if not file_extension.lower() in ['.dcm', '.dicom']:
        continue
    try:
        dcm = dicom.read_file(os.path.join(args.dir, filename))
        new_base =  str(dcm.PatientID) if args.use_patient_id else file_base
        newname  = new_base + "_" + str(dcm.InstanceNumber) + file_extension
        vprint("Renaming " + filename + " to " + newname)
        os.rename(os.path.join(args.dir, filename), os.path.join(args.dir, newname))
    except:
        print("Failed renaming " + filename + " to " + newname, file=sys.stderr)
        vprint("Error: " + str(sys.exc_info()[0]))
