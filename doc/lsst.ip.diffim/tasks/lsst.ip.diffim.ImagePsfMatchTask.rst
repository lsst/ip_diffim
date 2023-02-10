.. lsst-task-topic:: lsst.ip.diffim.ImagePsfMatchTask

##########################
ImagePsfMatchTask
##########################

Upon initialization, the kernel configuration is defined by self.config.kernel.active.
The task creates an lsst.afw.math.Warper from the subConfig self.config.kernel.active.warpingConfig.
A schema for the selection and measurement of candidate lsst.ip.diffim.KernelCandidates is
defined, and used to initize subTasks selectDetection (for candidate detection) and selectMeasurement
(for candidate measurement).

.. _lsst.ip.diffim.ImagePsfMatchTask-description:

Description
==================

Build a Psf-matching kernel using two input images, either as MaskedImages (in which case they need
to be astrometrically aligned) or Exposures (in which case astrometric alignment will happen by
default but may be turned off).  This requires a list of input Sources which may be provided
by the calling Task; if not, the Task will perform a coarse source detection
and selection for this purpose. Sources are vetted for signal-to-noise and masked pixels
(in both the template and science image), and substamps around each acceptable
source are extracted and used to create an instance of KernelCandidate.
Each KernelCandidate is then placed within a lsst.afw.math.SpatialCellSet, which is used by an ensemble of
lsst.afw.math.CandidateVisitor instances to build the Psf-matching kernel.   These visitors include, in
the order that they are called: BuildSingleKernelVisitor, KernelSumVisitor, BuildSpatialKernelVisitor,
and AssessSpatialKernelVisitor.

Sigma clipping of KernelCandidates is performed as follows:

- BuildSingleKernelVisitor, using the substamp diffim residuals from the per-source kernel fit,
    if PsfMatchConfig.singleKernelClipping is True
- KernelSumVisitor, using the mean and standard deviation of the kernel sum from all candidates,
    if PsfMatchConfig.kernelSumClipping is True
- AssessSpatialKernelVisitor, using the substamp diffim ressiduals from the spatial kernel fit,
    if PsfMatchConfig.spatialKernelClipping is True

The actual solving for the kernel (and differential background model) happens in
lsst.ip.diffim.PsfMatchTask._solve.  This involves a loop over the SpatialCellSet that first builds the
per-candidate matching kernel for the requested number of KernelCandidates per cell
(PsfMatchConfig.nStarPerCell).  The quality of this initial per-candidate difference image is examined,
using moments of the pixel residuals in the difference image normalized by the square root of the variance
(i.e. sigma); ideally this should follow a normal (0, 1) distribution,
but the rejection thresholds are set
by the config (PsfMatchConfig.candidateResidualMeanMax and PsfMatchConfig.candidateResidualStdMax).
All candidates that pass this initial build are then examined en masse to find the
mean/stdev of the kernel sums across all candidates.
Objects that are significantly above or below the mean,
typically due to variability or sources that are saturated in one image but not the other,
are also rejected.This threshold is defined by PsfMatchConfig.maxKsumSigma.
Finally, a spatial model is built using all currently-acceptable candidates,
and the spatial model used to derive a second set of (spatial) residuals
which are again used to reject bad candidates, using the same thresholds as above.

.. _lsst.ip.diffim.ImagePsfMatchTask-invoke:

Invoking the Task
==================

There is no run() method for this Task.  Instead there are 4 methods that
may be used to invoke the Psf-matching.  These are
`~lsst.ip.diffim.imagePsfMatch.ImagePsfMatchTask.matchMaskedImages`,
`~lsst.ip.diffim.imagePsfMatch.ImagePsfMatchTask.subtractMaskedImages`,
`~lsst.ip.diffim.imagePsfMatch.ImagePsfMatchTask.matchExposures`, and
`~lsst.ip.diffim.imagePsfMatch.ImagePsfMatchTask.subtractExposures`.

The methods that operate on lsst.afw.image.MaskedImage require that the images already be astrometrically
aligned, and are the same shape.  The methods that operate on lsst.afw.image.Exposure allow for the
input images to be misregistered and potentially be different sizes; by default a
lsst.afw.math.LanczosWarpingKernel is used to perform the astrometric alignment.  The methods
that "match" images return a Psf-matched image, while the methods that "subtract" images
return a Psf-matched and template subtracted image.

See each method's returned lsst.pipe.base.Struct for more details.

.. _lsst.ip.diffim.ImagePsfMatchTask--debug:

Debug Variables
==================

The ``pipetask`` command line interface supports a
flag ``--debug`` to import ref:  supports a flag `--debug`` to import debug.py from your PYTHONPATH; see :ref:`lsstDebug`
for more about debug.py files.

The available variables in ImagePsfMatchTask include:

.. code-block:: py

    import sys
    import lsstDebug
    def DebugInfo(name):
        di = lsstDebug.getInfo(name)
        if name == "lsst.ip.diffim.psfMatch":
            di.display = True                 # enable debug output
            di.maskTransparency = 80          # display mask transparency
            di.displayCandidates = True       # show all the candidates and residuals
            di.displayKernelBasis = False     # show kernel basis functions
            di.displayKernelMosaic = True     # show kernel realized across the image
            di.plotKernelSpatialModel = False # show coefficients of spatial model
            di.showBadCandidates = True       # show the bad candidates (red) along with good (green)
        elif name == "lsst.ip.diffim.imagePsfMatch":
            di.display = True                 # enable debug output
            di.maskTransparency = 30          # display mask transparency
            di.displayTemplate = True         # show full (remapped) template
            di.displaySciIm = True            # show science image to match to
            di.displaySpatialCells = True     # show spatial cells
            di.displayDiffIm = True           # show difference image
            di.showBadCandidates = True       # show the bad candidates (red) along with good (green)
        elif name == "lsst.ip.diffim.diaCatalogSourceSelector":
            di.display = False                # enable debug output
            di.maskTransparency = 30          # display mask transparency
            di.displayExposure = True         # show exposure with candidates indicated
            di.pauseAtEnd = False             # pause when done
        return di
    lsstDebug.Info = DebugInfo
    lsstDebug.frame = 1

Note that if you want addional logging info, you may add to your scripts:

.. code-block:: py

    import lsst.utils.logging as logUtils
    logUtils.trace_set_at("lsst.ip.diffim", 4)

.. _lsst.ip.diffim.ImagePsfMatchTask-examples:

Examples
==================

A complete example of using ImagePsfMatchTask

Create a subclass of ImagePsfMatchTask that allows us to either match exposures, or subtract exposures:

.. code-block:: none

    class MyImagePsfMatchTask(ImagePsfMatchTask):

        def __init__(self, args, kwargs):
            ImagePsfMatchTask.__init__(self, args, kwargs)

        def run(self, templateExp, scienceExp, mode):
            if mode == "matchExposures":
                return self.matchExposures(templateExp, scienceExp)
            elif mode == "subtractExposures":
                return self.subtractExposures(templateExp, scienceExp)

And allow the user the freedom to either run the script in default mode,
or point to their own images on disk.
Note that these images must be readable as an lsst.afw.image.Exposure.

We have enabled some minor display debugging in this script via the --debug option.  However, if you
have an lsstDebug debug.py in your PYTHONPATH you will get additional debugging displays.  The following
block checks for this script:

.. code-block:: py

    if args.debug:
        try:
            import debug
            # Since I am displaying 2 images here, set the starting frame number for the LSST debug LSST
            debug.lsstDebug.frame = 3
        except ImportError as e:
            print(e, file=sys.stderr)

Finally, we call a run method that we define below.
First set up a Config and modify some of the parameters.
E.g. use an "Alard-Lupton" sum-of-Gaussian basis,
fit for a differential background, and use low order spatial
variation in the kernel and background:

.. code-block:: py

    def run(args):
    #
    # Create the Config and use sum of gaussian basis
    #
    config = ImagePsfMatchTask.ConfigClass()
    config.kernel.name = "AL"
    config.kernel.active.fitForBackground = True
    config.kernel.active.spatialKernelOrder = 1
    config.kernel.active.spatialBgOrder = 0

Make sure the images (if any) that were sent to the script exist on disk and are readable.  If no images
are sent, make some fake data up for the sake of this example script (have a look at the code if you want
more details on generateFakeImages):

.. code-block:: py

    # Run the requested method of the Task
    if args.template is not None and args.science is not None:
        if not os.path.isfile(args.template):
            raise FileNotFoundError("Template image %s does not exist" % (args.template))
        if not os.path.isfile(args.science):
            raise FileNotFoundError("Science image %s does not exist" % (args.science))
        try:
            templateExp = afwImage.ExposureF(args.template)
        except Exception as e:
            raise RuntimeError("Cannot read template image %s" % (args.template))
        try:
            scienceExp = afwImage.ExposureF(args.science)
        except Exception as e:
            raise RuntimeError("Cannot read science image %s" % (args.science))
    else:
        templateExp, scienceExp = generateFakeImages()
        config.kernel.active.sizeCellX = 128
        config.kernel.active.sizeCellY = 128

Create and run the Task:

.. code-block:: py

    # Create the Task
    psfMatchTask = MyImagePsfMatchTask(config=config)
    # Run the Task
    result = psfMatchTask.run(templateExp, scienceExp, args.mode)

And finally provide some optional debugging displays:

.. code-block:: py

    if args.debug:
    # See if the LSST debug has incremented the frame number; if not start with frame 3
    try:
        frame = debug.lsstDebug.frame + 1
    except Exception:
        frame = 3
    afwDisplay.Display(frame=frame).mtv(result.matchedExposure,
                                        title="Example script: Matched Template Image")
    if "subtractedExposure" in result.getDict():
        afwDisplay.Display(frame=frame + 1).mtv(result.subtractedExposure,
                                                title="Example script: Subtracted Image")