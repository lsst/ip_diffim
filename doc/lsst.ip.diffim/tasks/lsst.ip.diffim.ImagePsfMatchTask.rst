.. lsst-task-topic:: lsst.ip.diffim.ImagePsfMatchTask

#################
ImagePsfMatchTask
#################

ImagePsfMatchTask creates a PSF-matching kernel for two images.

.. _lsst.ip.diffim.ImagePsfMatchTask-description:

Description
===========

Build a Psf-matching kernel using two input images, either as MaskedImages
(in which case they need to be astrometrically aligned) or Exposures (in which
case astrometric alignment will happen by default but may be turned off). This
requires a list of input Sources which may be provided by the calling Task; if
not, the Task will perform a coarse source detection and selection for this
purpose. Sources are vetted for signal-to-noise and masked pixels (in both the
template and science image), and substamps around each acceptable source are
extracted and used to create an instance of KernelCandidate. Each
KernelCandidate is then placed within a `~lsst.afw.math.SpatialCellSet`, which is
used by an ensemble of `~lsst.afw.math.CandidateVisitor` instances to build the
Psf-matching kernel. These visitors include, in the order that they are
called: ``BuildSingleKernelVisitor``, ``KernelSumVisitor``, ``BuildSpatialKernelVisitor``,
and ``AssessSpatialKernelVisitor``.

Upon initialization, the kernel configuration is defined by
``self.config.kernel.active``. The task creates an `~lsst.afw.math.Warper` from the
subConfig ``self.config.kernel.active.warpingConfig``. A schema for the selection
and measurement of candidate `~lsst.ip.diffim.KernelCandidates` is defined, and
used to initize subTasks selectDetection (for candidate detection) and
selectMeasurement(for candidate measurement).

Sigma clipping of KernelCandidates is performed as follows:

* BuildSingleKernelVisitor, using the substamp diffim residuals from the per-source kernel fit, if PsfMatchConfig.singleKernelClipping is True
* KernelSumVisitor, using the mean and standard deviation of the kernel sum from all candidates, if PsfMatchConfig.kernelSumClipping is True
* AssessSpatialKernelVisitor, using the substamp diffim ressiduals from the spatial kernel fit, if PsfMatchConfig.spatialKernelClipping is True

The actual solving for the kernel (and differential background model) happens
in ``lsst.ip.diffim.PsfMatchTask._solve``.  This involves a loop over the
SpatialCellSet that first builds the per-candidate matching kernel for the
requested number of KernelCandidates per cell
(``PsfMatchConfig.nStarPerCell``).  The quality of this initial per-candidate
difference image is examined, using moments of the pixel residuals in the
difference image normalized by the square root of the variance(i.e. sigma);
ideally this should follow a normal (0, 1) distribution, but the rejection
thresholds are set by the config (``PsfMatchConfig.candidateResidualMeanMax`` and
``PsfMatchConfig.candidateResidualStdMax``). All candidates that pass this initial
build are then examined en masse to find the mean/stdev of the kernel sums
across all candidates. Objects that are significantly above or below the mean,
typically due to variability or sources that are saturated in one image but
not the other, are also rejected.This threshold is defined by
``PsfMatchConfig.maxKsumSigma``. Finally, a spatial model is built using all
currently-acceptable candidates, and the spatial model used to derive a second
set of (spatial) residuals which are again used to reject bad candidates,
using the same thresholds as above.

.. _lsst.ip.diffim.ImagePsfMatchTask-invoke:

Invoking the Task
-----------------

There is no run() method for this Task.  Instead there are 4 methods that
may be used to invoke the Psf-matching.  These are
`~lsst.ip.diffim.imagePsfMatch.ImagePsfMatchTask.matchMaskedImages`,
`~lsst.ip.diffim.imagePsfMatch.ImagePsfMatchTask.subtractMaskedImages`,
`~lsst.ip.diffim.imagePsfMatch.ImagePsfMatchTask.matchExposures`, and
`~lsst.ip.diffim.imagePsfMatch.ImagePsfMatchTask.subtractExposures`.

The methods that operate on `~lsst.afw.image.MaskedImage` require that the images
already be astrometrically aligned, and are the same shape.  The methods that
operate on `~lsst.afw.image.Exposure` allow for the input images to be
misregistered and potentially be different sizes; by default a
`~lsst.afw.math.LanczosWarpingKernel` is used to perform the astrometric
alignment.  The methods that "match" images return a Psf-matched image, while
the methods that "subtract" images return a Psf-matched and template
subtracted image.

See each method's returned `~lsst.pipe.base.Struct` for more details.

.. _lsst.ip.diffim.ImagePsfMatchTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.ImagePsfMatchTask

.. _lsst.ip.diffim.ImagePsfMatchTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.ImagePsfMatchTask

.. _lsst.ip.diffim.ImagePsfMatchTask-config:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.ModelPsfMatchTask

.. _lsst.ip.diffim.ImagePsfMatchTask-debug:

Debugging
=========

The ``pipetask`` command line interface supports a ``--debug`` flag to import
``debug.py`` from your PYTHONPATH; see :ref:`lsstDebug` for more about ``debug.py``
files.
The available variables in ImagePsfMatchTask include:

display : `bool`
    Enable debug display output.
maskTransparency : `float`
    Transparency of mask planes in the output display.
displayCandidates : `bool`
    Show all the candidates and residuals.
displayKernelBasis : `bool`
    Show kernel basis functions.
displayKernelMosaic : `bool`
    Show kernel realized across the image.
plotKernelSpatialModel : `bool`
    Show coefficients of spatial model.
showBadCandidates : `bool`
    Show the bad candidates (red) along with good (green).
displayTemplate : `bool`
    Show full (remapped) template.
displaySciIm : `bool`
    Show science image to match to.
displaySpatialCells : `bool`
    Show spatial cells.
displayDiffIm : `bool`
    Show difference image.
