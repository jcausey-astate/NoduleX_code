'''A library of itk-attach functions, to be used to build out itk pipelines.'''


def IMG_UC(dim=3):  # pylint: disable=invalid-name
    '''dynamically load unsigned character 3d image type to preven super long
    loads on import.'''
    from itk import UC, Image  # pylint: disable=no-name-in-module
    return Image[UC, dim]


def IMG_F(dim=3):  # pylint: disable=invalid-name
    '''dynamically load 3d float image type to preven super long loads when
    on import.'''
    from itk import F, Image  # pylint: disable=no-name-in-module
    return Image[F, dim]


def extract_image_type(img_type):
    '''Awful hack to get around the fact that there's apparently no way to ask
    an itk Image its PixelType'''
    import itk

    type_abbrev = img_type.__name__[img_type.__name__.rfind("Image")+5:-1]
    dim = int(img_type.__name__[-1])

    return (getattr(itk, type_abbrev), dim)


class PipeStage(object):
    '''A stub itk pipeline stage, to be inherited from by other classes.'''

    def __init__(self, template, previous_stage, params=None):
        self.prev = previous_stage
        self.template = template
        self.instance = self._instantiate(template)

        if params is not None:
            try:
                for param in params:
                    set_method = getattr(self.instance, param)
                    set_method(params[param])
            except TypeError:
                print "Failed to set the parameter", param, "on", \
                      type(self.instance)
                raise

    def in_type(self):
        '''Get the itk type that is input for this pipe stage. Default
        behavior is to draw automatically from the output of the previous
        stage.'''
        return self.prev.out_type()

    def out_type(self):
        '''Get the itk type that is output for this pipe stage. By default the
        behavior is to simply output with the same type as provided for input.
        '''
        return self.in_type()

    def _instantiate(self, template):
        '''Instantiate an instance of the wrapped class template for use.
        Useful to override in cases where there is unusual templating.'''
        instance = template[self.in_type(), self.out_type()].New()

        return instance

    def _finished(self, instance):
        '''A hook for asking questions about the instance after instance.
        Update() has been run. Useful for debugging or printing informative
        output.'''
        pass

    def _bind_input(self):
        '''Bind the input of the previous pipeline stage to an instance of the
        class template. Will result in a call to Update() for upstream stages.
        '''
        self.instance.SetInput(self.prev.execute())

    def execute(self):
        '''Execute this and all previous stages recursively to build output
        from this pipeline stage. Returns the result of a GetOutput call to
        the wrapped itk object.'''

        # this can't happen in the constructor since it requires a call to
        # execute()
        self._bind_input()

        self.instance.Update()

        return self.instance.GetOutput()


class StatsStage(PipeStage):
    '''An itk pipestage to get statistics on an imput image.'''

    def __init__(self, previous_stage):
        # pylint: disable=no-name-in-module,no-member
        from itk import StatisticsImageFilter as stats

        super(StatsStage, self).__init__(stats, previous_stage)

    def _instantiate(self, template):
        return template[self.in_type()].New()

    def max(self):
        '''Get the maximum value of the input pipestage. Induces a call to
        Update().'''
        self.execute()
        return self.instance.GetMaximum()

    def min(self):
        '''Get the minimum value of the input pipestage. Induces a call to
        Update().'''
        self.execute()
        return self.instance.GetMinimum()

    def stddev(self):
        '''Get the standard deviation of the input pipestage image. Induces a
        call to Update().'''
        self.execute()
        return self.instance.GetSigma()

    def sum(self):
        '''Get the sum of the input pipestage's image. Results in a call to
        Update().'''
        self.execute()
        return self.instance.GetSum()

    def mean(self):
        '''Get the average of the input pipestage's image. Results in a call
        to Update().'''
        self.execute()
        return self.instance.GetMean()


class BinaryStage(PipeStage):
    '''A generic superclass to manage one-templated binary image filters.'''

    def _instantiate(self, template):
        return template[self.in_type()].New()

    def _bind_input(self):
        # at the last possible second, determine the appropriate forground
        # value for the biary image by looking for the maximum value in the
        # input.
        stats = StatsStage(self.prev)
        self.instance.SetForegroundValue(stats.max())

        super(BinaryStage, self)._bind_input()


class VotingIterativeBinaryFillholeStage(BinaryStage):
    '''An itk PipeStage that implements a
    VotingBinaryIterativeHoleFillingImageFilter'''

    def __init__(self, previous_stage, **kwargs):
        # pylint: disable=no-name-in-module,no-member
        from itk import VotingBinaryIterativeHoleFillingImageFilter as fillhole

        super(VotingIterativeBinaryFillholeStage, self).__init__(
            fillhole,
            previous_stage)

        self.instance.SetMaximumNumberOfIterations(
            kwargs.get('iterations', 10))
        self.instance.SetMajorityThreshold(kwargs.get('threshold', 3))


class BinaryFillholeStage(BinaryStage):
    '''An itk PipeStage that implements BinaryFillholeImageFilter.'''

    def __init__(self, previous_stage):
        # pylint: disable=no-name-in-module,no-member
        from itk import BinaryFillholeImageFilter as fillhole

        super(BinaryFillholeStage, self).__init__(fillhole, previous_stage)

        self.instance.SetForegroundValue(1)


class CurvatureFlowStage(PipeStage):
    '''An itk PipeStage that implements CurvatureFlowImageFilter.'''

    def __init__(self, previous_stage, timestep, iterations):
        # pylint: disable=no-name-in-module,no-member
        from itk import CurvatureFlowImageFilter

        params = {"SetNumberOfIterations": iterations,
                  "SetTimeStep": timestep}

        template = CurvatureFlowImageFilter
        super(CurvatureFlowStage, self).__init__(template,
                                                 previous_stage,
                                                 params)


class ConfidenceConnectStage(PipeStage):
    '''An itk PipeStage that implements ConfidenceConnectedImageFilter.
    Default values for parameters drawn from ITKExamples
    SegmentWithGeodesicActiveContourLevelSet.'''

    # pylint: disable=too-many-arguments
    def __init__(self, previous_stage, seed, iterations=5, stddevs=3.0,
                 neighborhood=1):
        # pylint: disable=no-name-in-module
        from itk import ConfidenceConnectedImageFilter

        params = {"AddSeed": seed,
                  "SetMultiplier": stddevs,
                  "SetNumberOfIterations": iterations,
                  "SetInitialNeighborhoodRadius": neighborhood}

        template = ConfidenceConnectedImageFilter
        super(ConfidenceConnectStage, self).__init__(template,
                                                     previous_stage,
                                                     params)

    def out_type(self):
        '''ConfidenceConnectedImageFilter is only able to output as unsigned
        characters.'''
        availiable_out_types = [pair[1] for pair in self.template
                                if pair[0] == self.in_type()]

        if len(availiable_out_types) == 0:
            s = "".join(["ConfidenceConnectPipeStage could not find an",
                         "acceptable output type based upon output type",
                         str(self.in_type()), ". Options were:",
                         str(self.template.GetTypes())])
            raise TypeError(s)

        # temporary hack
        return availiable_out_types[-1]


class FileReader(object):
    '''A PipeStage that can initiate a pipeline using an itk ImageFileReader.
    '''

    def __init__(self, fname, img_type=IMG_F()):
        self.fname = fname
        self.img_type = img_type

    def out_type(self):
        '''Get type of image read by the wrapped ImageFileReader. This is
        determined based upon user choice at construction.'''
        return self.img_type

    def execute(self):
        '''Execute this pipeline stage--that is, read the image from file and
        build the data into the appropriate itk Image object.'''
        from itk import ImageFileReader  # pylint: disable=no-name-in-module

        reader = ImageFileReader[self.out_type()].New()

        reader.SetFileName(self.fname)

        reader.Update()

        return reader.GetOutput()


class FileWriter(object):
    '''A PipeStage that can close a pipeline by writing to file with an itk
    ImageFileWriter.'''

    def __init__(self, previous_stage, fname):
        self.fname = fname
        self.prev = previous_stage

    def in_type(self):
        '''The type of image provided to the ImageFileWriter by the pipeline.
        '''
        return self.prev.out_type()

    def execute(self):
        '''Execute this pipeline stage--that is, write to file the itk Image
        provided by the input to this pipeline.'''
        from itk import ImageFileWriter  # pylint: disable=no-name-in-module

        writer = ImageFileWriter[self.in_type()].New()
        writer.SetFileName(self.fname)

        writer.SetInput(self.prev.execute())
        writer.Update()


class AnisoDiffStage(PipeStage):
    '''An itk PipeStage that implements
    CurvatureAnisotropicDiffusionImageFilter. Default values for parameters
    drawn from ITKExamples SegmentWithGeodesicActiveContourLevelSet.'''

    def __init__(self, previous_stage, timestep=0.01, iterations=50,
                 conductance=9.0):
        # pylint: disable=no-name-in-module
        from itk import CurvatureAnisotropicDiffusionImageFilter as templ

        params = {"SetTimeStep": timestep,
                  "SetNumberOfIterations": iterations,
                  "SetConductanceParameter": conductance}

        super(AnisoDiffStage, self).__init__(templ, previous_stage, params)


class GradMagStage(PipeStage):
    '''An itk PipeStage that implements GradientMagnitudeImageFilter'''

    def __init__(self, previous_stage):
                # pylint: disable=no-name-in-module
        from itk import GradientMagnitudeImageFilter as templ

        super(GradMagStage, self).__init__(templ, previous_stage, {})


class GradMagRecGaussStage(PipeStage):
    '''An itk PipeStage that implements
    GradientMagnitudeRecursiveGaussianImageFilter.'''

    def __init__(self, previous_stage, sigma):
        # pylint: disable=no-name-in-module
        from itk import GradientMagnitudeRecursiveGaussianImageFilter as templ

        params = {"SetSigma": sigma}

        super(GradMagRecGaussStage, self).__init__(templ, previous_stage,
                                                   params)


class SigmoidStage(PipeStage):
    '''An itk PipeStage that implements SigmoidImageFilter. Output min/max
    drawn from ITKExamples SegmentWithGeodesicActiveContourLevelSet.'''

    # pylint: disable=too-many-arguments
    def __init__(self, previous_stage, alpha, beta, out_max=1.0, out_min=0.0):
        # pylint: disable=no-name-in-module,no-member
        from itk import SigmoidImageFilter

        params = {"SetOutputMinimum": out_min,
                  "SetOutputMaximum": out_max,
                  "SetAlpha": alpha,
                  "SetBeta": beta}

        template = SigmoidImageFilter
        super(SigmoidStage, self).__init__(template, previous_stage, params)


class FastMarchingStage(PipeStage):
    '''An itk PipeStage that implements SigmoidImageFilter. It can be run as a
    pure distance calculator with the 'imageless' parameter set to true (an
    input pipe is still required to produce correct output size) or as a true
    image segmentation filter with 'imageless' set to false'''

    # pylint: disable=too-many-arguments
    def __init__(self, previous_stage, imageless, seeds, seed_value,
                 stopping_value=1000):
        # pylint: disable=no-name-in-module,no-member
        from itk import FastMarchingImageFilter

        params = {"SetStoppingValue": stopping_value}

        self.imageless = imageless
        if imageless:
            params["SetSpeedConstant"] = 1.0

        template = FastMarchingImageFilter
        super(FastMarchingStage, self).__init__(template, previous_stage,
                                                params)

        self.instance.SetTrialPoints(self.build_seeds(seeds, seed_value))

    def _bind_input(self):
        output = self.prev.execute()

        self.instance.SetOutputSize(output.GetBufferedRegion().GetSize())
        self.instance.SetOutputSpacing(output.GetSpacing())

        if not self.imageless:
            self.instance.SetInput(output)

    def build_seeds(self, seeds, seed_value):
        '''Construct an itk.VectorContainer of itk.LevelSetNode object from
        given input seeds.'''
        # pylint: disable=no-name-in-module,no-member
        from itk import LevelSetNode, VectorContainer, UI

        (px_type, dim) = extract_image_type(self.in_type())
        node_type = LevelSetNode[px_type, dim]

        seed_vect = VectorContainer[UI, node_type].New()
        seed_vect.Initialize()

        for i, seed in enumerate(seeds):
            node = node_type()
            node.SetValue(-seed_value)
            node.SetIndex(seed)
            seed_vect.InsertElement(i, node)

        return seed_vect


class LevelSetFilterStage(PipeStage):
    '''A base class for PipeStages wrapping ImageFilters that have two inputs:
    a usual input and a 'feature input'.'''

    def __init__(self, templ, previous_stage, feature_stage, params):
        self.prev_feature = feature_stage

        super(LevelSetFilterStage, self).__init__(templ, previous_stage,
                                                  params)

    def _bind_input(self):
        super(LevelSetFilterStage, self)._bind_input()
        self.instance.SetFeatureImage(self.prev_feature.execute())

    def _finished(self, instance):
        print instance.GetElapsedIterations()

    def _instantiate(self, template):
        # LevelSetImageFilters have an unusual 3-argument
        # templating, which is problematic for PipeStage's dynamic
        # instantiation protocol. Appropriate implementation here.

        img_type = self.in_type()
        feature_type = self.prev_feature.out_type()

        for avail_templ in self.template:
            if avail_templ[0] == img_type and avail_templ[1] == feature_type:
                return template[avail_templ].New()

        s = " ".join(["Could not instantiate", str(self.template), "because ",
                      "no valid template combination of", str(img_type), "and",
                      str(feature_type), "could be found. Possibilites were:",
                      str([t for t in self.template])])
        raise TypeError(s)


class ShapeDetectionStage(LevelSetFilterStage):
    '''An itk PipeStage that implements the ShapeDetectionLevelSetImageFilter
    '''

    def __init__(self, previous_stage, feature_stage, prop_curve_ratio):
        # pylint: disable=no-name-in-module,no-member
        from itk import ShapeDetectionLevelSetImageFilter as Shape

        (propagation, curvature) = (1.0, prop_curve_ratio)

        params = {"SetPropagationScaling": -propagation,
                  "SetCurvatureScaling": curvature,
                  # "SetMaximumRMSError": 2,
                  # "SetNumberOfIterations": 800,
                  "SetIsoSurfaceValue": -1,
                  }

        super(ShapeDetectionStage, self).__init__(Shape, previous_stage,
                                                  feature_stage, params)


class GeoContourLSetStage(LevelSetFilterStage):
    '''An itk PipeStage that implements a
    GeodesicActiveContourLevelSetImageFilter in the pipestage framework.'''

    def __init__(self, prev_stage, feature_stage, **kwargs):
        # pylint: disable=no-name-in-module,no-member
        from itk import GeodesicActiveContourLevelSetImageFilter as Geodesic

        params = {"SetPropagationScaling": kwargs['propagation_scaling'],
                  "SetNumberOfIterations": kwargs['iterations'],
                  "SetCurvatureScaling": kwargs.get('curvature_scaling', 1.0),
                  "SetMaximumRMSError": kwargs.get('max_rms', 0.02)}

        super(GeoContourLSetStage, self).__init__(Geodesic, prev_stage,
                                                  feature_stage, params)


class BinaryThreshStage(PipeStage):
    '''An itk PipeStage that implements the BinaryThresholdImageFilter.'''

    def __init__(self, previous_stage, threshold):
        # pylint: disable=no-name-in-module,no-member
        from itk import BinaryThresholdImageFilter as BinThresh
        from itk import NumericTraits

        params = {'SetLowerThreshold': threshold[0],
                  'SetUpperThreshold': threshold[1]}

        super(BinaryThreshStage, self).__init__(BinThresh, previous_stage,
                                                params)

        # the superclass configures the in_type and out_type stuff to work
        # so we have to wait until after the superclass constructor to access
        # the previous PipeStage in a nice way.
        px_type = extract_image_type(self.out_type())[0]

        # you'd think that this would lead to high-valued surround, but the
        # opposite turns out to be true. I flipped them and now it's working.
        # could investigate later...
        self.instance.SetOutsideValue(NumericTraits[px_type].max())
        self.instance.SetInsideValue(NumericTraits[px_type].min())

    def out_type(self):
        preferred = IMG_UC()

        avail_templ = [t[1] for t in self.template
                       if t[0] == self.in_type()]

        if preferred in avail_templ:
            return preferred
        elif len(avail_templ) > 1:
            print "Warning: automatically choosing", avail_templ[0], \
                  "as BinaryThreshStage output."
            return avail_templ[0]
        else:
            assert len(avail_templ) == 0
            return avail_templ[0]


class WatershedStage(PipeStage):
    '''An itk PipeStage that wraps WatershedImageFilter'''

    def __init__(self, previous_stage, level, threshold):
        # pylint: disable=no-name-in-module
        from itk import WatershedImageFilter

        params = {"SetLevel": level,
                  "SetThreshold": threshold}

        super(WatershedStage, self).__init__(WatershedImageFilter,
                                             previous_stage,
                                             params)

    def out_type(self):
        from itk import Image, UL

        # ugly hack -- could be determined on a  platform/compiletime basis
        return Image[UL, 3]
        # return extract_image_type(self.instance.GetOutput())

    def _instantiate(self, template):
        return template[self.in_type()].New()


class ConverterStage(PipeStage):
    '''An itk PipeStage that implements CastImageFilter to convert from the
    pipeline output type to the specified type. Dimensionality is determined
    dynamically based upon the input pipe stage.'''

    def __init__(self, previous_stage, type_out):
        # pylint: disable=no-name-in-module
        from itk import CastImageFilter

        self.type_out = type_out
        super(ConverterStage, self).__init__(CastImageFilter, previous_stage)

    def out_type(self):
        import itk

        assert self.out_type is not None

        if self.type_out is float:
            self.type_out = itk.F
        elif hasattr(itk, str(self.type_out)):
            self.type_out = getattr(itk, self.type_out)

        dim = extract_image_type(self.in_type())[1]

        return itk.Image[self.type_out, dim]
