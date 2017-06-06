'''Run a python-based segmentation using a geodesic active contour level set
algorithm.'''

import sys
import argparse
import segstrats
import os.path
import json


def process_command_line(argv):
    '''Parse the command line and do a first-pass on processing them into a
    format appropriate for the rest of the script.'''

    parser = argparse.ArgumentParser(formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)

    parser.add_argument("images", nargs="+",
                        help="The image that should be segmented.")
    parser.add_argument('--seeds', type=str,
                        help="A list of files initial points in JSON format.")
    parser.add_argument('--intermediate_images', action="store_true",
                        default=False, help="Produce pipeline intermediate " +
                        "images (i.e. after each filter stage.")

    segstrats.register_options(segstrats.aniso_gauss_watershed, parser)

    args = parser.parse_args(argv[1:])

    with open(args.seeds) as f:
        args.seeds = json.loads(f.read())

    return args


def input2output(fname, label, path=None):
    '''Build output filename from input filename.'''

    ext = fname[fname.rfind('.'):]
    new_fname = fname.rstrip(ext) + "-" + label + ext

    if path:
        new_fname = os.path.join(path, os.path.basename(new_fname))

    return new_fname


class DateTimeEncoder(json.JSONEncoder):  # pylint: disable=C0111

    def default(self, obj):  # pylint: disable=E0202
        import datetime

        if isinstance(obj, datetime.datetime):
            return str(obj)
        elif isinstance(obj, datetime.timedelta):
            return str(obj)

        return json.JSONEncoder.default(self, obj)


def main(argv=None):
    '''Run the driver script for this module. This code only runs if we're
    being run as a script. Otherwise, it's silent and just exposes methods.'''
    args = process_command_line(argv)

    def get_seed(fname):
        basefname = os.path.basename(fname)
        return args.seeds[basefname]

    seg_opts = {"watershed":  {"level": args.watershed_level,
                               "threshold": args.watershed_threshold},
                "gauss": {"sigma": args.sigma}
                }

    result = segstrats.batch_segment(segstrats.aniso_gauss_watershed,
                                     "waterseg", args.path, args.images,
                                     input2output, seg_opts, get_seed)

    with open("waterseg-opts.json", 'w') as f:
        json_out = json.dumps(result, sort_keys=True,
                              indent=4, separators=(',', ': '),
                              cls=DateTimeEncoder)
        f.write(json_out)

    return 1

if __name__ == "__main__":
    sys.exit(main(sys.argv))
