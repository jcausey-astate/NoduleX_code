"""
Just reads the file produced by segment_to_binary_image and produces an output containing
the coordinates for the seed point in the orientation that can be used as input to 
"masterseg.py".

Usage: binary_seg_info_to_seed.py segment-info-file [output-file]
"""

from __future__ import print_function
import sys, os

def usage():
    print("Usage\n\t{0} segment-info-file [output-file]".format(os.path.basename(sys.argv[0])))

def main():
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)

    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    input_file  = sys.argv[1]

    with open(input_file, 'rU') as fin:
        fin.readline() # throw away the first line
        info = fin.readline()
        coords = info.split(',')
        coords = [c.strip() for c in coords]
        output_string = "{0} {1} {2}\n".format(coords[0], coords[1], int(coords[4]) - int(coords[3]))
        if not output_file is None:
            with open(output_file, 'w') as fout:
                fout.write(output_string)
        else:
            print(output_string)

if __name__ == '__main__':
    main()