import sys
import argparse
import os
import datetime
import json
import logging

import SimpleITK as sitk
import numpy as np

import sitkstrats
import bounding

# a flag to run the script in debug mode. ONLY SET in process_command_line.
global DEBUG  # pylint: disable=W0604
DEBUG = False

def process_command_line(argv):
    '''Parse the command line and do a first-pass on processing them into a
    format appropriate for the rest of the script.'''

    parser = argparse.ArgumentParser(formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "image",
        help="The image to process on.")
    parser.add_argument(
        "--nseeds", type=int, default=10,
        help="The number of randomly placed seeds to produce.")
    parser.add_argument(
        '--media_root', default="media_root/",
        help="The directory to store temporary and intermediate media output")
    parser.add_argument(
        '--profile', default=False, action='store_true',
        help="Run cProfile on script execution.")
    parser.add_argument(
        '--log', default="logs/",
        help="The directory to place logs in.")
    parser.add_argument(
        '--seed', default=None, nargs=3, type=int, metavar=('X', 'Y', 'Z'),
        help="Add an additional, manually determined seed to the " +
        "calculation. Seed should be image-indexed (x, y, z not z, y, x).")
    parser.add_argument(
        '--debug', default=False, action="store_true",
        help="Set the script to run in debug mode, where it produces FAR " +
        "fewer intermediate files.")

    args = parser.parse_args(argv[1:])
    args.media_root = os.path.abspath(args.media_root)
    args.image = os.path.abspath(args.image)
    args.log = os.path.abspath(args.log)

    global DEBUG #pylint: disable=W0603
    DEBUG = args.debug

    return args


def set_label(fname, label, labsep='-'):
    '''Set the label (a string addition of labsep + label) for this filename.
    '''
    ext = fname[fname.rfind('.'):]
    fname = fname[:fname.rfind('.')]+labsep+label+ext

    return fname


def opthash(options):
    '''Produce a short hash of the input options.'''
    import hashlib

    sha = hashlib.sha1()
    sha.update(str(options))

    return sha.hexdigest()[0:8]


def debug_log(func, arg, *args, **kwargs):
    '''
    Wrapper for mediadir_log that writes to disk only if DEBUG is set to true.
    '''
    if DEBUG:
        return mediadir_log(func, arg, *args, **kwargs)
    else:
        return func(*arg)


def mediadir_log(func, (in_img, in_opts), mediadir, sha, subdir=None):
    '''
    Invoke some image processing step in the pipeline and write the resulting
    file to the directory appropriate to the algorithm/step
    that generated it using its sha and the function. Also decorates the info
    object with information about the file's location.
    '''
    optha = opthash(in_opts)
    label = func.__name__

    (img, opts) = func(in_img, in_opts)

    if subdir is None:
        subdir = label

    out_fname = os.path.join(mediadir, subdir, sha+"-"+optha+'.nii')

    sitkstrats.write(img, out_fname)

    opts['file'] = out_fname

    return (img, opts)


def configure_strats():
    '''
    Construct a dictionary that represents the configuration of all
    segmentation strategies to be used in the script using command line
    arguments. CURRENTLY ACCEPTS NO INPUT.
    '''

    strats = {
        'confidence_connected': {
            'seed-independent': {
                'strategy': sitkstrats.curvature_flow,
                'opts': {'curvature_flow': {'timestep': 0.01,
                                            'iterations': 25}}
                },
            'seed-dependent': {
                'strategy': sitkstrats.confidence_connected,
                'opts': {'conf_connect': {'iterations': 2,
                                          'multiplier': 1.5,
                                          'neighborhood': 1},
                         'dialate': {'radius': 1}}
                },
        },
        'geodesic': {
            'seed-independent': {
                'strategy': sitkstrats.aniso_gauss_sigmo,
                'opts': {"anisodiff": {'timestep': 0.01,
                                       'conductance': 9.0,
                                       'iterations': 50},
                         "gauss": {'sigma': 1.5},
                         "sigmoid": {'alpha': -20,
                                     'beta': 50}},
            },
            'seed-dependent': {
                'strategy': sitkstrats.fastmarch_seeded_geocontour,
                'opts': {"geodesic": {"propagation_scaling": 2.0,
                                      "iterations": 300,
                                      "curvature_scaling": 1.0,
                                      "max_rms_change": 1e-7},
                         "seed_shift": 3}
            }
        },
        'watershed': {
            'seed-independent': {
                'strategy': sitkstrats.aniso_gauss_watershed,
                'opts': {"anisodiff": {'timestep': 0.01,
                                       'conductance': 9.0,
                                       'iterations': 50},
                         "gauss": {'sigma': 1.5},
                         "watershed": {"level": 20}}
            },
            'seed-dependent': {
                'strategy': sitkstrats.isolate_watershed,
                'opts': {}
            }
        }
    }

    return strats


def seeddep(imgs, seeds, root_dir, sha, segstrats, lung_size, img_in):

    # pick an image, basically at random, from imgs to initialize an array
    # that tracks which areas of the image have already been segmented out
    segmented = np.zeros(sitk.GetArrayFromImage(  # pylint: disable=E1101
        imgs.values()[0]).shape)

    out_info = {}

    for seed in seeds:
        try:
            if segmented[seed[2], seed[1], seed[0]] >= 2:
                logging.info(
                    "Tried to segment %s but it was already segmented", seed)
                continue
        except IndexError as err:
            sys.stderr.write("Tried to access " + str(seed) + " as " +
                             str(list(reversed(seed))) + " in img of size " +
                             str(segmented.shape)+"\n")
            logging.error(" ".join([str(seed), str(segmented.shape)]))
            logging.error(str(err))
            raise

        # We want to hold onto images and info dicts for each segmentation,
        # and we want to automagically store the info we put in seed_info into
        # out_info for returning later => use setdefault
        out_imgs = {}
        seed_info = out_info.setdefault("-".join([str(k) for k in seed]), {})

        # for each strategy we want to segment with, get its name and the
        # function that executes it.
        for (sname, strat) in [(strnam, segstrats[strnam]['seed-dependent'])
                               for strnam in segstrats]:

            img_in = imgs[sname]

            opts = dict(strat['opts'])
            opts['seed'] = seed

            (tmp_img, tmp_info) = debug_log(strat['strategy'],
                                            (img_in, opts),
                                            root_dir,
                                            sha)

            out_imgs[sname] = tmp_img
            seed_info[sname] = tmp_info

            logging.info("Segmented %s with %s", seed, sname)

        # we need the names of the input files so that our options hash is
        # dependent on the input images.
        seed_indep_hashes = [sitkstrats.hash_img(i)[0:8]
                             for i in out_imgs.values()]

        try:
            # First, we compute the segmentation union
            (consensus, consensus_info) = mediadir_log(
                sitkstrats.segmentation_union,
                (out_imgs.values(),
                 {'threshold': 2.0/3.0,
                  'max_size': lung_size * 0.5,
                  'min_size': lung_size * 1e-5,
                  'indep_img_hashes': seed_indep_hashes}),
                root_dir,
                sha)

            # Then we crop down both the initial image ("img_in") AND the
            # segmentation based upon the size of the segmentation.
            (crop_seg, crop_seg_info) = sitkstrats.crop_to_segmentation(
                img=consensus, seg_img=consensus, padding_px=5)
            (crop_img, crop_img_info) = sitkstrats.crop_to_segmentation(
                img=img_in, seg_img=consensus, padding_px=5)

            # Because we used the image "consensus" for both croppings,
            # the images should be cropped the same way
            assert crop_seg.GetSize() == crop_img.GetSize()
            assert crop_seg_info['origin'] == crop_img_info['origin']
            assert crop_seg_info['padding'] == crop_img_info['padding']

            logging.info("Cropped %s to %s.", seed, crop_img.GetSize())

            consensus_info.update(crop_seg_info)

            mediadir_log(
                lambda x, y: (x, y),
                (crop_seg, crop_seg_info),
                root_dir,
                sha,
                subdir="consensus-label")
            mediadir_log(
                lambda x, y: (x, y),
                (crop_img, crop_img_info),
                root_dir,
                sha,
                subdir="consensus-grey")

        except RuntimeWarning as war:
            logging.info("Failed %s during consensus: %s", seed, war)
            seed_info['consensus'] = "failure"
            continue
        except ValueError as err:
            # this ocurrs when the segmentation runs to the edge of the image.
            logging.info("Failed %s during cropping: %s", seed, err)
            seed_info['consensus'] = "failure"
            continue

        logging.info("Finished segmenting %s", seed)

        segmented += sitk.GetArrayFromImage(consensus)

        seed_info['consensus'] = consensus_info

    return out_info


def run_img(img, sha, nseeds, root_dir, addl_seed):  # pylint: disable=C0111
    '''Run the entire protocol on a particular image starting with sha hash'''
    img_info = {}

    lung_img, lung_info = debug_log(sitkstrats.segment_lung,
                                    (img, {'probe_size': 7}),
                                    root_dir, sha)
    img_info['lungseg'] = lung_info

    # (img, tmp_info) = debug_log(sitkstrats.crop_to_segmentation,
    #                             (img_in, lung_img),
    #                             root_dir, sha, subdir="crop")
    # lung_img = sitkstrats.crop_to_segmentation(lung_img, lung_img)[0]
    # img_info['crop'] = tmp_info

    segstrats = configure_strats()
    seed_indep_imgs = {}
    seed_indep_info = {}
    seeds           = []

    for (sname, strat) in [(strnam, segstrats[strnam]['seed-independent'])
                           for strnam in segstrats]:

        try:
            optha = opthash(strat['opts'])
            fname = os.path.join(root_dir, strat['strategy'].__name__,
                                 sha + "-" + optha + ".nii")
            tmp_img = sitkstrats.read(fname)
            tmp_info = strat['opts']
            tmp_info['file'] = os.path.join(fname)
            logging.info(
                "Loaded seed-independent image for '%s', '%s' from file",
                sname, fname)
        except RuntimeError:
            logging.debug(
                "Building seed-independent image '%s', '%s'.", sname, fname)
            (tmp_img, tmp_info) = debug_log(strat['strategy'],
                                            (img, strat['opts']),
                                            root_dir,
                                            sha)
            logging.info(
                "Built seed-independent image '%s', '%s' in %s",
                sname, fname, tmp_info['time'])

        seed_indep_imgs[sname] = tmp_img
        seed_indep_info[sname] = tmp_info

    # compute seeds, first by taking the centers of mass of a bunch of the
    # watershed segemented regions, then by adding a bunch of random ones that
    # are inside the lung field.
    if nseeds > 1 or addl_seed is None:
        (seeds, tmp_info) = sitkstrats.com_calc(img=seed_indep_imgs['watershed'],
                                                max_size=0.05, min_size=1e-5,
                                                lung_img=lung_img)
        img_info['deterministic-seeds'] = tmp_info
        seeds.extend(sitkstrats.distribute_seeds(lung_img, nseeds-len(seeds)))
    else:
        print("Provided seed only.")
        img_info['deterministic-seeds'] = None

    if addl_seed is not None:
        # additional seed is given in terms of the input image, and must
        # be corrected for cropping.
        if 'crop' in img_info:
            origin = img_info['crop']['origin']
            assert len(addl_seed) == len(origin)
            addl_seed = [addl_seed[i] - origin[i]
                         for i in range(len(addl_seed))]
        seeds.insert(0, addl_seed)

    if len(seeds) > nseeds:
        logging.warning("The number of seeds generated in the deterministic " +
                        "phase (%s) is greater than the allowed number of " +
                        "seeds (%s). The list of seeds is being truncated.",
                        len(seeds), nseeds)

    # with many deterministic seeds, this list can be longer than nseeds.
    seeds = seeds[0:nseeds]

    seg_info = seeddep(seed_indep_imgs, seeds,
                       root_dir, sha, segstrats, img_info['lungseg']['size'],
                       img)

    img_info['noduleseg'] = {}
    for seed in seg_info:
        for segstrat in seg_info[seed]:
            combined_info = {'seed-dependent': seg_info[seed][segstrat]}

            if segstrat in seed_indep_info:
                combined_info['seed-independent'] = seed_indep_info[segstrat]

            img_info['noduleseg'].setdefault(
                seed, {})[segstrat] = combined_info

    return img_info


class DateTimeEncoder(json.JSONEncoder):  # pylint: disable=C0111

    def default(self, obj):  # pylint: disable=E0202
        if isinstance(obj, datetime.datetime):
            return str(obj)
        elif isinstance(obj, datetime.timedelta):
            return str(obj)

        return json.JSONEncoder.default(self, obj)


def write_info(info, filename="masterseg-run.json"):
    with open(filename, 'w') as f:
        try:
            json_out = json.dumps(info, sort_keys=True,
                                  indent=2, separators=(',', ': '),
                                  cls=DateTimeEncoder)
        except TypeError as err:
            logging.error("Error encountered serializing for JSON, dumping " +
                          "dict here:\n"+str(info))
            raise err

        f.write(json_out)


def log_name_gen(sha, log_dir):
    '''
    Build a the absolute path name of the log file for an image, given its sha
    and a log file directory.
    '''
    logfilename = "-".join([sha, str(datetime.datetime.now())])+".log"
    logfilename = logfilename.replace(" ", "-")
    logfilename = os.path.join(log_dir, logfilename)

    return logfilename


def main(argv=None):
    '''Run the driver script for this module. This code only runs if we're
    being run as a script. Otherwise, it's silent and just exposes methods.'''
    args = process_command_line(argv)

    basename = os.path.basename(args.image)
    sha = basename[:basename.rfind('.')]

    # this only gets run the first time around, so if >1 image is
    # specified, they will ALL be in the file named by the sha of the first
    # one. But grepping around isn't too hard.
    logging.basicConfig(filename=log_name_gen(sha, args.log),
                        level=logging.DEBUG,
                        format='%(asctime)s %(message)s')

    logging.info("Beginning image %s", args.image)

    try:
        run_info = run_img(sitkstrats.read(args.image), sha,
                           args.nseeds, args.media_root, args.seed)
    except Exception as exc:  # pylint: disable=W0703
        logging.critical("Encountered critical exception:\n%s", exc)
        raise

    write_info(run_info, filename=os.path.join(args.log, sha+"-seg.json"))

    return 0


if __name__ == "__main__":
    if "--profile" in sys.argv:
        import cProfile
        sys.exit(cProfile.runctx("main(sys.argv)", globals(),
                                 {"sys.argv": sys.argv}))
    else:
        sys.exit(main(sys.argv))
