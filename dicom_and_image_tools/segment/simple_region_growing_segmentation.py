'''Run a python-based version of the Slicer "Simple Region Growing Segmentation
strategy, which uses ConfidenceConnectedImageFilter followed by
CurvatureFlowImageFilter. I found the defaults in both ITK and Slicer to be
silly, so they're set more sensibly here. More complete documentation is
availiable through ITK.'''

import sys
import argparse
import os.path

from segstrats import flow_confidence


def process_command_line(argv):
    '''Parse the command line and do a first-pass on processing them into a
    format appropriate for the rest of the script.'''

    parser = argparse.ArgumentParser(formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)

    parser.add_argument("images", nargs="+",
                        help="The image that should be segmented.")
    parser.add_argument("--connect_iterations", type=int, default=2,
                        help="The number of iterations to run" +
                             " ConfidenceConnectedImageFilter")
    parser.add_argument('--connect_stddevs', type=float, default=2.0,
                        help="The number of voxel property standard devs " +
                        "to consider as connected.")
    parser.add_argument('--connect_neighborhood', type=int, default=1,
                        help='the number of local pixels around the seed to ' +
                        'use as the start of the calculation.')
    parser.add_argument("--smooth_iterations", default=25, type=long,
                        help="The number of iterations to run" +
                             " CurvatureFlowImageFilter")
    parser.add_argument('--smooth_timestep', default=0.01, type=float,
                        help="The step size used by CurvatureFlowImageFilter")
    parser.add_argument('--seeds', type=str,
                        help="A list of files initial points in JSON format.")
    parser.add_argument('--sigma', default=1.0, type=float,
                        help="The stddev in units of image spacing for the " +
                             "GradientMagnitudeRecursiveGaussianImageFilter.")
    parser.add_argument('-p', '--path', default=None,
                        help="The segmented file to output")
    parser.add_argument('--label', default="confseg",
                        help="The label to add to each file.")
    parser.add_argument('--intermediate_images', action="store_true",
                        default=False, help="Produce pipeline intermediate " +
                        "images (i.e. after each filter stage.")
    parser.add_argument('--no_overwrite', default=False, action="store_true",
                        help="Do not overwrite existing files.")

    args = parser.parse_args(argv[1:])

    configs = argparse.Namespace()

    import json
    with open(args.seeds) as f:
        seeds = json.loads(f.read())
    configs.seeds = seeds

    configs.images = args.images
    configs.label = args.label
    configs.path = args.path
    configs.no_overwrite = args.no_overwrite

    configs.intermediate_images = args.intermediate_images

    configs.connect = {}
    configs.connect['neighborhood'] = args.connect_neighborhood
    configs.connect['stddevs'] = args.connect_stddevs
    configs.connect['iterations'] = args.connect_iterations

    configs.smooth = {}
    configs.smooth['iterations'] = args.smooth_iterations
    configs.smooth['timestep'] = args.smooth_timestep

    configs.gauss = {'sigma': args.sigma}

    return configs


def input2output(fname, label, path=None):
    '''Build output filename from input filename.'''

    ext = fname[fname.rfind('.'):]
    new_fname = fname.rstrip(ext) + "-" + label + ext

    if path:
        new_fname = os.path.join(path, os.path.basename(new_fname))

    return new_fname


def main(argv=None):
    '''Run the driver script for this module. This code only runs if we're
    being run as a script. Otherwise, it's silent and just exposes methods.'''
    import datetime

    allstart = datetime.datetime.now()
    configs = process_command_line(argv)

    print "Segmenting", len(configs.images), "images"
    times = []
    skipped = []

    for fname in configs.images:
        basefname = os.path.basename(fname)
        sys.stdout.write("Segmenting " + basefname + "... ")
        sys.stdout.flush()
        start = datetime.datetime.now()

        try:
            configs.connect['seeds'] = configs.seeds[basefname]
        except KeyError as exc:
            skipped.append(fname)
            print "skipped (" + str(exc) + ")"
            continue

        try:
            outname = input2output(fname, configs.label, configs.path)

            if os.path.isfile(outname) and configs.no_overwrite:
                print "outfile exists, skipping"
                continue

            flow_confidence(
                fname, outname,
                smooth=configs.smooth,
                gauss=configs.gauss,
                connect=configs.connect,
                intermediate_images=configs.intermediate_images
                )
        except Exception:  # pylint: disable=W0703
            if len(configs.images) > 1:
                skipped.append(fname)
                continue
            else:
                raise

        times.append(datetime.datetime.now() - start)
        print "took", times[-1]

    if len(times) > 0:
        print "min/avg/max", min(times), \
              sum(times, datetime.timedelta())/len(times), max(times)
        print "total time:", datetime.datetime.now() - allstart
    print "skipped", len(skipped), "files:"
    print "\n".join(skipped)

    return 1

if __name__ == "__main__":
    sys.exit(main(sys.argv))
