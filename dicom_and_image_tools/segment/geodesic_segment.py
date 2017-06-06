'''Run a python-based segmentation using a geodesic active contour level set
algorithm.'''

import sys
import argparse
from segstrats import aniso_gauss_sigmo_geocontour
import os.path


def process_command_line(argv):
    '''Parse the command line and do a first-pass on processing them into a
    format appropriate for the rest of the script.'''

    parser = argparse.ArgumentParser(formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)

    parser.add_argument("images", nargs="+",
                        help="The image that should be segmented.")
    parser.add_argument('--seeds', type=str,
                        help="A list of files initial points in JSON format.")
    parser.add_argument('-p', '--path', default=None,
                        help="The segmented file to output")
    parser.add_argument('--sigma', default=1.0, type=float,
                        help="The stddev in units of image spacing for the " +
                             "GradientMagnitudeRecursiveGaussianImageFilter.")
    parser.add_argument('--alpha', default=-15, type=float,
                        help="Alpha ('A') parameter in sigmoid filter.  " +
                        "Transition width.")
    parser.add_argument('--beta', default=150, type=float,
                        help="Beta ('B') parameter in sigmoid filter. Obeys" +
                        " the expression exp((-x+B)/A). Zero adjustment.")
    parser.add_argument('--propagation_scaling', default=7.0, type=float,
                        help="The weight on propagation force in level set " +
                        "segmentation.")
    parser.add_argument('--curvature_scaling', default=1.0, type=float,
                        help="The weight on propagation force in level set " +
                        "segmentation.")
    parser.add_argument('--geodesic_iterations', default=300, type=int,
                        help="The number of iterations by the " +
                        "GeodesicActiveContourLevelSetImageFilter")
    parser.add_argument('--seed_distance', default=1, type=int,
                        help="The expected distance from the seed to the" +
                        "first level set.")
    parser.add_argument('--intermediate_images', action="store_true",
                        default=False, help="Produce pipeline intermediate " +
                        "images (i.e. after each filter stage.")

    args = parser.parse_args(argv[1:])

    import json
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


def main(argv=None):
    '''Run the driver script for this module. This code only runs if we're
    being run as a script. Otherwise, it's silent and just exposes methods.'''
    import datetime

    allstart = datetime.datetime.now()
    args = process_command_line(argv)

    print "Segmenting", len(args.images), "images"
    times = []
    skipped = []

    with open("geoseg-opts.json", 'w') as f:
        import json
        arg_dict = vars(args)
        arg_dict['begin_time'] = str(datetime.datetime.now())
        json_out = json.dumps(arg_dict, sort_keys=True,
                              indent=4, separators=(',', ': '))
        f.write(json_out)

    for fname in args.images:
        basefname = os.path.basename(fname)
        sys.stdout.write("Segmenting " + basefname + "... ")
        sys.stdout.flush()
        start = datetime.datetime.now()

        try:
            seed = args.seeds[basefname]
        except KeyError as exc:
            skipped.append(fname)
            print "skipped (" + str(exc) + ")"
            continue

        try:
            outname = input2output(fname, "geoseg", args.path)

            res = aniso_gauss_sigmo_geocontour(
                fname, outname,
                gauss={'sigma': args.sigma},
                sigmo={'alpha': args.alpha, 'beta': args.beta},
                seed=[seed],
                seed_distance=args.seed_distance,
                geodesic={"iterations": args.geodesic_iterations,
                          "propagation_scaling": args.propagation_scaling,
                          "curvature_scaling": args.curvature_scaling},
                intermediate_images=args.intermediate_images)

            sys.stdout.write(str(res))

        except Exception as e:  # pylint: disable=W0703,C0103
            if len(args.images) > 1:
                skipped.append(fname)
                sys.stdout.write("skipped (threw " + str(e) + " )")
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
