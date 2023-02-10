.. lsst-task-topic:: lsst.ip.diffim.ModelPsfMatchTask

##########################
ModelPsfMatchTask
##########################

This Task differs from ``ImagePsfMatchTask`` in that it matches two Psf _models_ by realizing
them in an Exposure-sized SpatialCellSet and then inserting each Psf-image pair into ``KernelCandidates``.
Because none of the pairs of sources that are to be matched should be invalid, all sigma clipping is
turned off in ModelPsfMatchConfig.  And because there is no tracked _variance_ in the Psf images, the
debugging and logging QA info should be interpreted with caution.

One item of note is that the sizes of Psf models are fixed (e.g. its defined as a 21x21 matrix).  When the
Psf-matching kernel is being solved for, the Psf "image" is convolved with each kernel basis function,
leading to a loss of information around the borders.
This pixel loss will be problematic for the numerical
stability of the kernel solution if the size of the convolution kernel
(set by ModelPsfMatchConfig.kernelSize) is much bigger than: psfSize//2.
Thus the sizes of Psf-model matching kernels are typically smaller
than their image-matching counterparts.  If the size of the kernel is too small, the convolved stars will
look "boxy"; if the kernel is too large, the kernel solution will be "noisy".  This is a trade-off that
needs careful attention for a given dataset.

The primary use case for this Task is in matching an Exposure to a
constant-across-the-sky Psf model for the purposes of image coaddition.
It is important to note that in the code, the "template" Psf is the Psf
that the science image gets matched to.  In this sense the order of template and science image are
reversed, compared to ImagePsfMatchTask, which operates on the template image.

.. _lsst.ip.diffim.ModelPsfMatchTask-debug:

Debug variables
===============

The ``pipetask`` command line interface supports a
flag ``--debug`` to import ref:  supports a flag `--debug`` to import debug.py from your PYTHONPATH; see :ref:`lsstDebug`
for more about debug.py files.

The available variables in ModelsPsfMatchTask include:

.. code-block:: py

    import sys
    import lsstDebug
    def DebugInfo(name):
        di = lsstDebug.getInfo(name)
        if name == "lsst.ip.diffim.psfMatch":
            di.display = True                 # global
            di.maskTransparency = 80          # mask transparency
            di.displayCandidates = True       # show all the candidates and residuals
            di.displayKernelBasis = False     # show kernel basis functions
            di.displayKernelMosaic = True     # show kernel realized across the image
            di.plotKernelSpatialModel = False # show coefficients of spatial model
            di.showBadCandidates = True       # show the bad candidates (red) along with good (green)
        elif name == "lsst.ip.diffim.modelPsfMatch":
            di.display = True                 # global
            di.maskTransparency = 30          # mask transparency
            di.displaySpatialCells = True     # show spatial cells before the fit
        return di
    lsstDebug.Info = DebugInfo
    lsstDebug.frame = 1

Note that if you want additional logging info, you may add to your scripts:

.. code-block:: py

    import lsst.utils.logging as logUtils
    logUtils.trace_set_at("lsst.ip.diffim", 4)

.. _lsst.ip.diffim.ModelPsfMatchTask-examples:

Examples
=========
A complete example of using ModelPsfMatchTask

Create a subclass of ModelPsfMatchTask that accepts two exposures.
Note that the "template" exposure contains the Psf that will get matched to,
and the "science" exposure is the one that will be convolved:

.. code-block:: py

    class MyModelPsfMatchTask(ModelPsfMatchTask):
        def __init__(self, *args, **kwargs):
            ModelPsfMatchTask.__init__(self, *args, **kwargs)
        def run(self, templateExp, scienceExp):
            return ModelPsfMatchTask.run(self, scienceExp, templateExp.getPsf())

And allow the user the freedom to either run the script in default mode,
or point to their own images on disk. Note that these
images must be readable as an lsst.afw.image.Exposure:

.. code-block:: none

    if __name__ == "__main__":
        import argparse
        parser = argparse.ArgumentParser(description="Demonstrate the use of ModelPsfMatchTask")
        parser.add_argument("--debug", "-d", action="store_true", help="Load debug.py?", default=False)
        parser.add_argument("--template", "-t", help="Template Exposure to use", default=None)
        parser.add_argument("--science", "-s", help="Science Exposure to use", default=None)
        args = parser.parse_args()

We have enabled some minor display debugging in this script via the ``-–debug`` option.
However, if you have an :ref:`lsstDebug` debug.py in your PYTHONPATH you will get additional
debugging displays. The following block checks for this script:

.. code-block:: none

    if args.debug:
        try:
            import debug
            # Since I am displaying 2 images here, set the starting frame number for the LSST debug LSST
            debug.lsstDebug.frame = 3
        except ImportError as e:
            print(e, file=sys.stderr)

Finally, we call a run method that we define below.
First set up a Config and modify some of the parameters.
In particular we don't want to "grow" the sizes of the kernel or KernelCandidates,
since we are operating with fixed–size images (i.e. the size of the input Psf models).

.. code-block:: py

    def run(args):
        #
        # Create the Config and use sum of gaussian basis
        #
        config = ModelPsfMatchTask.ConfigClass()
        config.kernel.active.scaleByFwhm = False

Make sure the images (if any) that were sent to the script exist on disk and are readable.
If no images are sent, make some fake data up for the sake of this example script
(have a look at the code if you want more details on generateFakeData):

.. code-block:: none

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
        templateExp, scienceExp = generateFakeData()
        config.kernel.active.sizeCellX = 128
        config.kernel.active.sizeCellY = 128

.. code-block:: none

    if args.debug:
        afwDisplay.Display(frame=1).mtv(templateExp, title="Example script: Input Template")
        afwDisplay.Display(frame=2).mtv(scienceExp, title="Example script: Input Science Image")

Create and run the Task:

.. code-block:: none

    # Create the Task
    psfMatchTask = MyModelPsfMatchTask(config=config)
    # Run the Task
    result = psfMatchTask.run(templateExp, scienceExp)

And finally provide optional debugging display of the Psf-matched (via the Psf models) science image:

.. code-block:: none

    if args.debug:
        # See if the LSST debug has incremented the frame number; if not start with frame 3
        try:
            frame = debug.lsstDebug.frame + 1
        except Exception:
            frame = 3
        afwDisplay.Display(frame=frame).mtv(result.psfMatchedExposure,
                                            title="Example script: Matched Science Image")