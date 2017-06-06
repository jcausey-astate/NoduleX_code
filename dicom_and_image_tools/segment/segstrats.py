'''A collection of strategies for segmenting an image using python and itk.'''


def single_segment(fname, outname, seg_alg, seg_opts):
    '''Run a single segmentation run of an algorithm and process its output
    into a JSON-able dict.'''
    import datetime

    start = datetime.datetime.now()
    stats = {"start": start}

    res = seg_alg(fname, outname, **seg_opts)  # pylint: disable=W0142

    for key in res:
        stats[key] = res[key]

    stats['time'] = datetime.datetime.now() - start

    return stats


def batch_segment(seg_alg, seg_label, outpath,
                  files, input2output, seg_opts, get_seed):
    import datetime
    import os.path
    import sys

    allstart = datetime.datetime.now()

    stats = {'run': seg_opts}
    stats['run']['begin_time'] = datetime.datetime.now()
    stats['run']['seg_alg'] = seg_alg.__name__
    stats['run']['path'] = outpath

    for fname in files:
        basefname = os.path.basename(fname)
        sys.stdout.write("Segment "+basefname+"... ")
        sys.stdout.flush()

        try:
            seg_opts['seed'] = get_seed(fname)
        except KeyError as exc:
            stats[basefname] = 'skipped'
            print "skipped (" + str(exc) + ")"
            continue

        try:
            outname = input2output(fname, seg_label, outpath)

            res = single_segment(fname, outname, seg_alg, seg_opts)

        except Exception as e:  # pylint: disable=W0703,C0103
            if len(files) > 1:
                sys.stdout.write("skipped (threw " + str(e) + " )")
                stats[basefname] = 'skipped'
                continue
            else:
                raise

        stats[basefname] = res
        stats[basefname]['fullpath'] = fname
        stats[basefname]['outpath'] = outname

        print " took", stats[basefname]['time']

    times = [stats[f]['time'] for f in stats if 'time' in stats[f]]
    skipped = [f for f in stats if stats[f] == 'skipped']

    stats['run']['total_time'] = datetime.datetime.now() - allstart

    return stats


def register_options(segfunc, parser):
    '''Register appropriate options with an options parser.'''

    parser.add_argument('-p', '--path', default=None,
                        help="The segmented file to output")

    if segfunc is aniso_gauss_watershed:
        parser.add_argument(
            '--sigma', default=1.0, type=float,
            help="The stddev in units of image spacing for  the " +
                 "GradientMagnitudeRecursiveGaussian ImageFilter.")
        parser.add_argument(
            '--watershed_level', default=0.01, type=float,
            help="The weight on propagation force in level set segmentation.")
        parser.add_argument(
            '--watershed_threshold', default=.1, type=float,
            help="The number of iterations by the " +
                 "GeodesicActiveContourLevelSetImageFilter")


def aniso_gauss_watershed(in_image, out_image, **kwargs):
    '''Implements a basic watershed-based strategy for image segmentation'''

    import itk_attach

    gauss = kwargs['gauss']
    watershed = kwargs['watershed']

    pipe = itk_attach.FileReader(in_image)
    pipe = itk_attach.AnisoDiffStage(pipe, iterations=25)
    pipe = itk_attach.GradMagStage(pipe)

    pipe = itk_attach.WatershedStage(pipe,
                                     watershed['level'],
                                     watershed['threshold'])

    pipe = itk_attach.ConverterStage(pipe, "UC")

    pipe = itk_attach.FileWriter(pipe, out_image)

    pipe.execute()

    # A hacky solution that writes the file out using ITK and reads it back as
    # a numpy array to choose a segmentation based on a seed.
    import medpy.io
    import numpy as np

    (img, hdr) = medpy.io.load(out_image)

    seed = kwargs['seed']
    chosen_seg = img[seed[0], seed[1], seed[2]]

    img = np.array(img == chosen_seg, dtype='uint8')

    medpy.io.save(img, out_image, hdr)

    return {}


def aniso_gauss_sigmo_geocontour(in_image, out_image, **kwargs):
    '''Implements a basic strategy that relies upon a gradient magnitude
    geodesic level set strategy described in the ITK docs.'''

    import itk_attach

    gauss = kwargs['gauss']
    sigmo = kwargs['sigmo']
    geodesic = kwargs['geodesic']
    binary = kwargs.get('binary', {'threshold': (0.1, 1.5)})

    pipe = itk_attach.FileReader(in_image)
    aniso = itk_attach.AnisoDiffStage(pipe)
    gauss = itk_attach.GradMagRecGaussStage(aniso, gauss['sigma'])
    feature = itk_attach.SigmoidStage(gauss, sigmo['alpha'], sigmo['beta'])

    fastmarch = itk_attach.FastMarchingStage(
        pipe,
        imageless=True,
        seeds=kwargs['seed'],
        seed_value=kwargs['seed_distance'])

    geo = itk_attach.GeoContourLSetStage(
        fastmarch,
        feature,
        propagation_scaling=geodesic['propagation_scaling'],
        curvature_scaling=geodesic['curvature_scaling'],
        iterations=geodesic['iterations'])

    if kwargs.get('intermediate_images', False):
        itk_attach.FileWriter(aniso, 'out-aniso.nii').execute()
        itk_attach.FileWriter(gauss, 'out-gauss.nii').execute()
        itk_attach.FileWriter(feature, 'out-sigmo.nii').execute()
        itk_attach.FileWriter(fastmarch, 'out-march.nii').execute()

    pipe = itk_attach.BinaryThreshStage(geo, binary['threshold'])

    pipe = itk_attach.FileWriter(pipe, out_image)

    # run the pipeline
    pipe.execute()

    return {'geodesic_iterations': geo.instance.GetElapsedIterations()}


def aniso_gauss_confidence(in_image, out_image, **kwargs):
    '''Perform an aniso + gauss + confidence connected segmentation strategy.
    '''

    import itk_attach

    smooth_param = kwargs['smooth']
    gauss_param = kwargs['gauss']
    connect_param = kwargs['connect']
    intermediate_images = kwargs.get('intermediate_images', False)

    pipe = itk_attach.FileReader(in_image)

    pipe = itk_attach.AnisoDiffStage(pipe,
                                     smooth_param['timestep'],
                                     smooth_param['iterations'])
    if intermediate_images:
        itk_attach.FileWriter(pipe, "aniso.nii").execute()

    pipe = itk_attach.GradMagRecGaussStage(pipe, gauss_param['sigma'])
    if intermediate_images:
        itk_attach.FileWriter(pipe, "gauss.nii").execute()

    pipe = itk_attach.ConfidenceConnectStage(pipe,
                                             connect_param['seeds'],
                                             connect_param['iterations'],
                                             connect_param['stddevs'],
                                             connect_param['neighborhood'])

    pipe = itk_attach.FileWriter(pipe, out_image)

    pipe.execute()


def flow_confidence(in_image, out_image, **kwargs):
    '''Perform a curvatureflow + confidence connected segmentation strategy.'''

    import itk_attach

    smooth = kwargs['smooth']
    connect = kwargs['connect']
    intermed_img = kwargs.get('intermediate_images', False)

    pipe = itk_attach.FileReader(in_image)

    pipe = itk_attach.CurvatureFlowStage(pipe, smooth['timestep'],
                                         smooth['iterations'])

    if intermed_img:
        itk_attach.FileWriter(pipe, 'curvature.nii').execute()

    pipe = itk_attach.ConfidenceConnectStage(pipe, connect['seeds'],
                                             connect['iterations'],
                                             connect['stddevs'],
                                             connect['neighborhood'])

    if intermed_img:
        itk_attach.FileWriter(pipe, "confidence.nii").execute()

    binary = kwargs.get('binary', None)
    if binary:
        pipe = itk_attach.VotingIterativeBinaryFillholeStage(
            pipe,
            binary['threshold'],
            binary['iterations'])

        if intermed_img:
            itk_attach.FileWriter(pipe, "binvote.nii").execute()

    pipe = itk_attach.BinaryFillholeStage(pipe)

    pipe = itk_attach.FileWriter(pipe, out_image)

    pipe.execute()
