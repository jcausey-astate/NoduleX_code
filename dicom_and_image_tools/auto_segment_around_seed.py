"""
Adapted from Justin Porter's (https://github.com/justinrporter/nodule-seg)
'masterseg.py' to focus specifically on segmentations that involve
growing a region from a seed point inside the lung volume.
"""
import os, SimpleITK as sitk, numpy as np
from segment import masterseg, sitkstrats, segstrats

def get_consensus_for_seed(img, nodule_id, output_dir, seed, proximity=10, min_size=None, max_size=None):  # pylint: disable=C0111
    '''Run the entire protocol on a particular image''' 
    img_info = {}

    lung_img, lung_info = sitkstrats.segment_lung(img, {'probe_size': 7})
    print("type of lung img: {}".format(type(lung_img)))
    img_info['lungseg'] = lung_info
    print("running consensus segmentation for seed {}".format(seed))
    lung_img = sitk.Image(lung_img)
        
    img_arr = sitk.GetArrayFromImage(lung_img)
    shape   = img_arr.shape
    if img_arr[seed[2], seed[1], seed[0]] != 1:
        zprox = int(round(proximity / 2.0))
        proximity = np.array(img_arr[
            max(seed[2]-proximity,0):min(seed[2]+proximity+1, shape[0]), 
            max(seed[1]-proximity,0):min(seed[1]+proximity+1, shape[1]),
            max(seed[0]-zprox,0):min(seed[1]+zprox+1, shape[2]) 
        ])
        if proximity.sum() == 0:
            print("The seed is not inside the lung portion of the image.")
            return None, None

    print("Inside lung OK")
    # (img, tmp_info) = sitkstrats.crop_to_segmentation((img_in, lung_img),
    #                             output_dir, nodule_id, subdir="crop")
    # lung_img = sitkstrats.crop_to_segmentation(lung_img, lung_img)[0]
    # img_info['crop'] = tmp_info

    segstrats = configure_seed_driven_strats()
    seed_indep_imgs = {}
    seed_indep_info = {}
    seeds           = []
    print("Running strategies...")
    for (sname, strat) in [(strnam, segstrats[strnam]['seed-independent'])
                           for strnam in segstrats]:
        print("Running strategy {}".format(sname))

        try:
            optha = masterseg.opthash(strat['opts'])
            fname = os.path.join(output_dir, strat['strategy'].__name__,
                                 nodule_id + "-" + optha + ".nii")
            tmp_img  = sitkstrats.read(fname)
            tmp_info = strat['opts']
            tmp_info['file'] = os.path.join(fname)
            # logging.info(
                # "Loaded seed-independent image for '%s', '%s' from file",
                # sname, fname)
        except RuntimeError:
            # logging.debug(
            #     "Building seed-independent image '%s', '%s'.", sname, fname)
            (tmp_img, tmp_info) = strat['strategy'](img, strat['opts'])
            # logging.info(
            #     "Built seed-independent image '%s', '%s' in %s",
            #     sname, fname, tmp_info['time'])

        seed_indep_imgs[sname] = tmp_img
        seed_indep_info[sname] = tmp_info
    
    print("Provided seed only.")
    img_info['deterministic-seeds'] = None
    seeds = []
    if seed is not None:
        # additional seed is given in terms of the input image, and must
        # be corrected for cropping.
        if 'crop' in img_info:
            origin = img_info['crop']['origin']
            assert len(seed) == len(origin)
            seed = [seed[i] - origin[i]
                         for i in range(len(seed))]
        seeds = [seed]
    
    consensus_img, seg_info = get_seed_dep_union(seed_indep_imgs, seeds,
                       output_dir, nodule_id, segstrats, img_info['lungseg']['size'],
                       img, min_size, max_size)

    img_info['noduleseg'] = {}
    for seed in seg_info:
        for segstrat in seg_info[seed]:
            combined_info = {'seed-dependent': seg_info[seed][segstrat]}

            if segstrat in seed_indep_info:
                combined_info['seed-independent'] = seed_indep_info[segstrat]

            img_info['noduleseg'].setdefault(
                seed, {})[segstrat] = combined_info

    return consensus_img, img_info

def get_seed_dep_union(imgs, seeds, root_dir, nodule_id, segstrats, lung_size, img_in, min_size=None, max_size=None):
    # pick an image, basically at random, from imgs to initialize an array
    # that tracks which areas of the image have already been segmented out
    segmented = np.zeros(sitk.GetArrayFromImage(  # pylint: disable=E1101
        imgs.values()[0]).shape)

    out_info = {}

    for seed in seeds:
        try:
            if segmented[seed[2], seed[1], seed[0]] >= 2:
                # logging.info(
                #     "Tried to segment %s but it was already segmented", seed)
                print("Tried to segment %s but it was already segmented", seed)
                continue
        except IndexError as err:
            sys.stderr.write("Tried to access " + str(seed) + " as " +
                             str(list(reversed(seed))) + " in img of size " +
                             str(segmented.shape)+"\n")
            # logging.error(" ".join([str(seed), str(segmented.shape)]))
            # logging.error(str(err))
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

            (tmp_img, tmp_info) = strat['strategy'](img_in, opts)

            out_imgs[sname]  = tmp_img
            seed_info[sname] = tmp_info

            # logging.info("Segmented %s with %s", seed, sname)

        # we need the names of the input files so that our options hash is
        # dependent on the input images.
        seed_indep_hashes = [sitkstrats.hash_img(i)[0:8]
                             for i in out_imgs.values()]

        consensus = np.zeros(sitk.GetArrayFromImage(  # pylint: disable=E1101
            imgs.values()[0]).shape)
        try:
            # First, we compute the segmentation union
            (consensus, consensus_info) = sitkstrats.segmentation_union(out_imgs.values(),
                 {'threshold': 2.0/3.0,
                  'max_size': lung_size * 0.33 if max_size is None else max_size,
                  'min_size': lung_size * 5e-6 if min_size is None else min_size,
                  'indep_img_hashes': seed_indep_hashes})

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

            # logging.info("Cropped %s to %s.", seed, crop_img.GetSize())

            consensus_info.update(crop_seg_info)

        except RuntimeWarning as war:
            # logging.info("Failed %s during consensus: %s", seed, war)
            print("Failed {} during consensus: {}".format(seed, war))
            seed_info['consensus'] = "failure"
            continue
        except ValueError as err:
            # this ocurrs when the segmentation runs to the edge of the image.
            # logging.info("Failed %s during cropping: %s", seed, err)
            print("Failed {} during cropping: {}".format(seed, err))
            seed_info['consensus'] = "failure"
            continue

        # logging.info("Finished segmenting %s", seed)
        print("Finished segmenting {}".format(seed))

        segmented += sitk.GetArrayFromImage(consensus)
        segmented[segmented != 0] = 1
        consensus  = sitk.GetImageFromArray(segmented)
        if 'consensus' in seed_info and seed_info['consensus'] == 'failure':
            consensus_info = None
            consensus      = None
        out_info['consensus'] = consensus_info

    return consensus, out_info


def configure_seed_driven_strats_aggressive():
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

def configure_seed_driven_strats_conservative():
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
                                            'iterations': 20}}
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
                                       'iterations': 20},
                         "gauss": {'sigma': 1.5},
                         "sigmoid": {'alpha': -20,
                                     'beta': 50}},
            },
            'seed-dependent': {
                'strategy': sitkstrats.fastmarch_seeded_geocontour,
                'opts': {"geodesic": {"propagation_scaling": 2.0,
                                      "iterations": 50,
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
                                       'iterations': 30},
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

def configure_seed_driven_strats():
    return configure_seed_driven_strats_conservative()
    return strats

def main():
    print("This script is not designed to run in standalone mode.")
    exit(1)

if __name__ == '__main__':
    main()