#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import numpy as np
import lsst.daf.base as dafBase
import lsst.pex.logging as pexLog
import lsst.pex.config as pexConfig
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.pipe.base as pipeBase
from lsst.meas.algorithms import SourceDetectionTask, SourceMeasurementTask, getBackground
from .makeKernelBasisList import makeKernelBasisList
from .psfMatch import PsfMatchTask, PsfMatchConfigDF, PsfMatchConfigAL
from . import utils as diUtils 
from . import diffimLib
from . import diffimTools
import lsst.afw.display.ds9 as ds9

sigma2fwhm = 2. * np.sqrt(2. * np.log(2.))

class ImagePsfMatchConfig(pexConfig.Config):
    """!Configuration for image-to-image Psf matching"""
    kernel = pexConfig.ConfigChoiceField(
        doc="kernel type",
        typemap=dict(
            AL=PsfMatchConfigAL,
            DF=PsfMatchConfigDF
        ),
        default="AL",
    )
    selectDetection = pexConfig.ConfigurableField(
        target=SourceDetectionTask,
        doc="Initial detections used to feed stars to kernel fitting",
    )
    selectMeasurement = pexConfig.ConfigurableField(
        target=SourceMeasurementTask,
        doc="Initial measurements used to feed stars to kernel fitting",
    )

    def setDefaults(self):
        # High sigma detections only
        self.selectDetection.reEstimateBackground = False
        self.selectDetection.thresholdValue = 10.0

        # Minimal set of measurments for star selection
        self.selectMeasurement.algorithms.names.clear()
        self.selectMeasurement.algorithms.names = ('flux.psf', 'flags.pixel', 'shape.sdss',
                                                   'flux.gaussian', 'skycoord')
        self.selectMeasurement.slots.modelFlux = None
        self.selectMeasurement.slots.apFlux = None

## \addtogroup LSST_task_documentation
## \{
## \page ImagePsfMatchTask
## \ref ImagePsfMatchTask_ "ImagePsfMatchTask"
## \copybrief ImagePsfMatchTask
## \}

class ImagePsfMatchTask(PsfMatchTask):
    """!
\anchor ImagePsfMatchTask_

\brief Psf-match two MaskedImages or Exposures using the sources in the images

\section ip_diffim_imagepsfmatch_Contents Contents

 - \ref ip_diffim_imagepsfmatch_Purpose
 - \ref ip_diffim_imagepsfmatch_Initialize
 - \ref ip_diffim_imagepsfmatch_IO
 - \ref ip_diffim_imagepsfmatch_Config
 - \ref ip_diffim_imagepsfmatch_Metadata
 - \ref ip_diffim_imagepsfmatch_Debug
 - \ref ip_diffim_imagepsfmatch_Example

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_imagepsfmatch_Purpose   Description

Build a Psf-matching kernel using two input images, either as MaskedImages (in which case they need
    to be astrometrically aligned) or Exposures (in which case astrometric alignment will happen by
    default but may be turned off).  This requires a list of input Sources which may be provided
by the calling Task; if not, the Task will perform a coarse source detection and selection for this purpose.
Sources are vetted for signal-to-noise and masked pixels (in both the template and science image), and 
substamps around each acceptable source are extracted and used to create an instance of KernelCandidate.  
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
(i.e. sigma); ideally this should follow a normal (0, 1) distribution, but the rejection thresholds are set 
by the config (PsfMatchConfig.candidateResidualMeanMax and PsfMatchConfig.candidateResidualStdMax).  
All candidates that pass this initial build are then examined en masse to find the 
mean/stdev of the kernel sums across all candidates.  Objects that are significantly above or below the mean, 
typically due to variability or sources that are saturated in one image but not the other, are also rejected.  
This threshold is defined by PsfMatchConfig.maxKsumSigma.  Finally, a spatial model is built using all
currently-acceptable candidates, and the spatial model used to derive a second set of (spatial) residuals
which are again used to reject bad candidates, using the same thresholds as above.

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_imagepsfmatch_Initialize    Task initialization

\copydoc \_\_init\_\_

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_imagepsfmatch_IO        Invoking the Task

There is no run() method for this Task.  Instead there are 4 methods that
may be used to invoke the Psf-matching.  These are
\link lsst.ip.diffim.imagePsfMatch.ImagePsfMatchTask.matchMaskedImages matchMaskedImages\endlink,
\link lsst.ip.diffim.imagePsfMatch.ImagePsfMatchTask.subtractMaskedImages subtractMaskedImages\endlink,
\link lsst.ip.diffim.imagePsfMatch.ImagePsfMatchTask.matchExposures matchExposures\endlink, and 
\link lsst.ip.diffim.imagePsfMatch.ImagePsfMatchTask.subtractExposures subtractExposures\endlink.

The methods that operate on lsst.afw.image.MaskedImage require that the images already be astrometrically
aligned, and are the same shape.  The methods that operate on lsst.afw.image.Exposure allow for the 
input images to be misregistered and potentially be different sizes; by default a 
lsst.afw.math.LanczosWarpingKernel is used to perform the astrometric alignment.  The methods 
that "match" images return a Psf-matched image, while the methods that "subtract" images 
return a Psf-matched and template subtracted image.

See each method's returned lsst.pipe.base.Struct for more details.

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_imagepsfmatch_Config       Configuration parameters

See \ref ImagePsfMatchConfig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_imagepsfmatch_Metadata   Quantities set in Metadata

See \ref ip_diffim_psfmatch_Metadata "PsfMatchTask"

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_imagepsfmatch_Debug     Debug variables

The \link lsst.pipe.base.cmdLineTask.CmdLineTask command line task\endlink interface supports a
flag \c -d/--debug to import \b debug.py from your \c PYTHONPATH.  The relevant contents of debug.py 
for this Task include:

\code{.py}
    import sys
    import lsstDebug
    def DebugInfo(name):
        di = lsstDebug.getInfo(name)   
        if name == "lsst.ip.diffim.psfMatch":
            di.display = True                 # enable debug output
            di.maskTransparency = 80          # ds9 mask transparency
            di.displayCandidates = True       # show all the candidates and residuals
            di.displayKernelBasis = False     # show kernel basis functions
            di.displayKernelMosaic = True     # show kernel realized across the image
            di.plotKernelSpatialModel = False # show coefficients of spatial model
            di.showBadCandidates = True       # show the bad candidates (red) along with good (green)
        elif name == "lsst.ip.diffim.imagePsfMatch":
            di.display = True                 # enable debug output
            di.maskTransparency = 30          # ds9 mask transparency
            di.displayTemplate = True         # show full (remapped) template
            di.displaySciIm = True            # show science image to match to
            di.displaySpatialCells = True     # show spatial cells
            di.displayDiffIm = True           # show difference image
            di.showBadCandidates = True       # show the bad candidates (red) along with good (green) 
        elif name == "lsst.ip.diffim.diaCatalogSourceSelector":
            di.display = False                # enable debug output
            di.maskTransparency = 30          # ds9 mask transparency
            di.displayExposure = True         # show exposure with candidates indicated
            di.pauseAtEnd = False             # pause when done
        return di
    lsstDebug.Info = DebugInfo
    lsstDebug.frame = 1      
\endcode

Note that if you want addional logging info, you may add to your scripts:
\code{.py}
import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)
\endcode

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_imagepsfmatch_Example   A complete example of using ImagePsfMatchTask

This code is imagePsfMatchTask.py in the examples directory, and can be run as \em e.g.
\code
examples/imagePsfMatchTask.py --debug 
examples/imagePsfMatchTask.py --debug --mode="matchExposures"
examples/imagePsfMatchTask.py --debug --template /path/to/templateExp.fits --science /path/to/scienceExp.fits
\endcode

\dontinclude imagePsfMatchTask.py
Create a subclass of ImagePsfMatchTask that allows us to either match exposures, or subtract exposures:
\skip MyImagePsfMatchTask
\until self.subtractExposures

And allow the user the freedom to either run the script in default mode, or point to their own images on disk.
Note that these images must be readable as an lsst.afw.image.Exposure:
\skip main
\until parse_args

We have enabled some minor display debugging in this script via the --debug option.  However, if you 
have an lsstDebug debug.py in your PYTHONPATH you will get additional debugging displays.  The following
block checks for this script:
\skip args.debug
\until sys.stderr

\dontinclude imagePsfMatchTask.py
Finally, we call a run method that we define below.  First set up a Config and modify some of the parameters.
E.g. use an "Alard-Lupton" sum-of-Gaussian basis, fit for a differential background, and use low order spatial
variation in the kernel and background:
\skip run(args)
\until spatialBgOrder

Make sure the images (if any) that were sent to the script exist on disk and are readable.  If no images
are sent, make some fake data up for the sake of this example script (have a look at the code if you want
more details on generateFakeImages):
\skip requested
\until sizeCellY

Create and run the Task:
\skip Create
\until args.mode

And finally provide some optional debugging displays:
\skip args.debug
\until result.subtractedExposure
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    """    
    ConfigClass = ImagePsfMatchConfig

    def __init__(self, *args, **kwargs):
        """!Create the ImagePsfMatchTask

        \param *args arguments to be passed to lsst.ip.diffim.PsfMatchTask.__init__
        \param **kwargs keyword arguments to be passed to lsst.ip.diffim.PsfMatchTask.__init__

        Upon initialization, the kernel configuration is defined by self.config.kernel.active.
        The task creates an lsst.afw.math.Warper from the subConfig self.config.kernel.active.warpingConfig.
        A schema for the selection and measurement of candidate lsst.ip.diffim.KernelCandidates is
        defined, and used to initize subTasks selectDetection (for candidate detection) and selectMeasurement
        (for candidate measurement).
        """
        PsfMatchTask.__init__(self, *args, **kwargs)
        self.kConfig = self.config.kernel.active
        self._warper = afwMath.Warper.fromConfig(self.kConfig.warpingConfig)
        self.selectSchema = afwTable.SourceTable.makeMinimalSchema()
        self.selectSchema.setVersion(0)
        self.selectAlgMetadata = dafBase.PropertyList()
        self.makeSubtask("selectDetection", schema=self.selectSchema)
        self.makeSubtask("selectMeasurement", schema=self.selectSchema, algMetadata=self.selectAlgMetadata)

    def getFwhmPix(self, psf):
        """!Return the FWHM in pixels of a Psf"""
        sigPix = psf.computeShape().getDeterminantRadius()
        return sigPix * sigma2fwhm

    @pipeBase.timeMethod
    def matchExposures(self, templateExposure, scienceExposure,
                       templateFwhmPix=None, scienceFwhmPix=None,
                       candidateList=None, doWarping=True, convolveTemplate=True):
        """!Warp and PSF-match an exposure to the reference

        Do the following, in order:
        - Warp templateExposure to match scienceExposure,
            if doWarping True and their WCSs do not already match
        - Determine a PSF matching kernel and differential background model
            that matches templateExposure to scienceExposure
        - Convolve templateExposure by PSF matching kernel

        @param templateExposure: Exposure to warp and PSF-match to the reference masked image
        @param scienceExposure: Exposure whose WCS and PSF are to be matched to
        @param templateFwhmPix: FWHM (in pixels) of the Psf in the template image (image to convolve)
        @param scienceFwhmPix: FWHM (in pixels) of the Psf in the science image
        @param candidateList: a list of footprints/maskedImages for kernel candidates; 
                              if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF
        @param doWarping: what to do if templateExposure's and scienceExposure's WCSs do not match:
            - if True then warp templateExposure to match scienceExposure
            - if False then raise an Exception
        @param convolveTemplate: convolve the template image or the science image
            - if True, templateExposure is warped if doWarping, templateExposure is convolved
            - if False, templateExposure is warped if doWarping, scienceExposure is convolved

        @return a pipeBase.Struct containing these fields:
        - matchedImage: the PSF-matched exposure =
            warped templateExposure convolved by psfMatchingKernel. This has:
            - the same parent bbox, Wcs and Calib as scienceExposure
            - the same filter as templateExposure
            - no Psf (because the PSF-matching process does not compute one)
        - psfMatchingKernel: the PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to solve for the PSF matching kernel

        Raise a RuntimeError if doWarping is False and templateExposure's and scienceExposure's
            WCSs do not match
        """
        if not self._validateWcs(templateExposure, scienceExposure):
            if doWarping:
                self.log.info("Astrometrically registering template to science image")
                templatePsf = templateExposure.getPsf()
                templateExposure = self._warper.warpExposure(scienceExposure.getWcs(),
                    templateExposure, destBBox=scienceExposure.getBBox())
                templateExposure.setPsf(templatePsf)
            else:
                pexLog.Trace(self.log.getName(), 1, "ERROR: Input images not registered")
                raise RuntimeError("Input images not registered")
        if templateFwhmPix is None:
            if not templateExposure.hasPsf():
                self.log.warn("No estimate of Psf FWHM for template image")
            else:
                templateFwhmPix = self.getFwhmPix(templateExposure.getPsf())

        if scienceFwhmPix is None:
            if not scienceExposure.hasPsf():
                self.log.warn("No estimate of Psf FWHM for science image")
            else:
                scienceFwhmPix = self.getFwhmPix(scienceExposure.getPsf())

        kernelSize = makeKernelBasisList(self.kConfig, templateFwhmPix, scienceFwhmPix)[0].getWidth()
        candidateList = self.makeCandidateList(templateExposure, scienceExposure, kernelSize, candidateList)

        if convolveTemplate:
            results = self.matchMaskedImages(
                templateExposure.getMaskedImage(), scienceExposure.getMaskedImage(), candidateList,
                templateFwhmPix=templateFwhmPix, scienceFwhmPix=scienceFwhmPix)
        else:
            results = self.matchMaskedImages(
                scienceExposure.getMaskedImage(), templateExposure.getMaskedImage(), candidateList,
                templateFwhmPix=scienceFwhmPix, scienceFwhmPix=templateFwhmPix)

        psfMatchedExposure = afwImage.makeExposure(results.matchedImage, scienceExposure.getWcs())
        psfMatchedExposure.setFilter(templateExposure.getFilter())
        psfMatchedExposure.setCalib(scienceExposure.getCalib())
        results.warpedExposure  = templateExposure
        results.matchedExposure = psfMatchedExposure
        return results

    @pipeBase.timeMethod
    def matchMaskedImages(self, templateMaskedImage, scienceMaskedImage, candidateList,
                          templateFwhmPix=None, scienceFwhmPix=None):
        """!PSF-match a MaskedImage (templateMaskedImage) to a reference MaskedImage (scienceMaskedImage)

        Do the following, in order:
        - Determine a PSF matching kernel and differential background model
            that matches templateMaskedImage to scienceMaskedImage
        - Convolve templateMaskedImage by the PSF matching kernel

        @param templateMaskedImage: masked image to PSF-match to the reference masked image;
            must be warped to match the reference masked image
        @param scienceMaskedImage: maskedImage whose PSF is to be matched to
        @param templateFwhmPix: FWHM (in pixels) of the Psf in the template image (image to convolve)
        @param scienceFwhmPix: FWHM (in pixels) of the Psf in the science image
        @param candidateList: a list of footprints/maskedImages for kernel candidates; 
                              if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF

        @return a pipeBase.Struct containing these fields:
        - psfMatchedMaskedImage: the PSF-matched masked image =
            templateMaskedImage convolved with psfMatchingKernel.
            This has the same xy0, dimensions and wcs as scienceMaskedImage.
        - psfMatchingKernel: the PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to solve for the PSF matching kernel

        Raise a RuntimeError if input images have different dimensions
        """

        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayTemplate = lsstDebug.Info(__name__).displayTemplate
        displaySciIm = lsstDebug.Info(__name__).displaySciIm
        displaySpatialCells = lsstDebug.Info(__name__).displaySpatialCells
        maskTransparency = lsstDebug.Info(__name__).maskTransparency
        if not maskTransparency:
            maskTransparency = 0
        ds9.setMaskTransparency(maskTransparency)

        if not candidateList:
            raise RuntimeError("Candidate list must be populated by makeCandidateList")

        if not self._validateSize(templateMaskedImage, scienceMaskedImage):
            pexLog.Trace(self.log.getName(), 1, "ERROR: Input images different size")
            raise RuntimeError("Input images different size")

        if display and displayTemplate:
            ds9.mtv(templateMaskedImage, frame=lsstDebug.frame, title="Image to convolve")
            lsstDebug.frame += 1

        if display and  displaySciIm:
            ds9.mtv(scienceMaskedImage, frame=lsstDebug.frame, title="Image to not convolve")
            lsstDebug.frame += 1

        kernelCellSet = self._buildCellSet(templateMaskedImage,
                                           scienceMaskedImage,
                                           candidateList)

        if display and displaySpatialCells:
            diUtils.showKernelSpatialCells(scienceMaskedImage, kernelCellSet,
                                           symb="o", ctype=ds9.CYAN, ctypeUnused=ds9.YELLOW, ctypeBad=ds9.RED,
                                           size=4, frame=lsstDebug.frame, title="Image to not convolve")
            lsstDebug.frame += 1

        if templateFwhmPix and scienceFwhmPix:
            self.log.info("Matching Psf FWHM %.2f -> %.2f pix" % (templateFwhmPix, scienceFwhmPix))

        if self.kConfig.useBicForKernelBasis:
            tmpKernelCellSet = self._buildCellSet(templateMaskedImage,
                                                  scienceMaskedImage,
                                                  candidateList)
            nbe = diffimTools.NbasisEvaluator(self.kConfig, templateFwhmPix, scienceFwhmPix)
            bicDegrees = nbe(tmpKernelCellSet, self.log)
            basisList = makeKernelBasisList(self.kConfig, templateFwhmPix, scienceFwhmPix,
                                            alardDegGauss=bicDegrees[0], metadata=self.metadata)
            del tmpKernelCellSet
        else:
            basisList = makeKernelBasisList(self.kConfig, templateFwhmPix, scienceFwhmPix,
                                            metadata=self.metadata)

        spatialSolution, psfMatchingKernel, backgroundModel = self._solve(kernelCellSet, basisList)




        psfMatchedMaskedImage = afwImage.MaskedImageF(templateMaskedImage.getBBox())
        doNormalize = False
        afwMath.convolve(psfMatchedMaskedImage, templateMaskedImage, psfMatchingKernel, doNormalize)
        return pipeBase.Struct(
            matchedImage=psfMatchedMaskedImage,
            psfMatchingKernel=psfMatchingKernel,
            backgroundModel=backgroundModel,
            kernelCellSet=kernelCellSet,
        )

    @pipeBase.timeMethod
    def subtractExposures(self, templateExposure, scienceExposure,
                          templateFwhmPix=None, scienceFwhmPix=None,
                          candidateList=None, doWarping=True, convolveTemplate=True):
        """!Register, Psf-match and subtract two Exposures

        Do the following, in order:
        - Warp templateExposure to match scienceExposure, if their WCSs do not already match
        - Determine a PSF matching kernel and differential background model
            that matches templateExposure to scienceExposure
        - PSF-match templateExposure to scienceExposure
        - Compute subtracted exposure (see return values for equation).

        @param templateExposure: exposure to PSF-match to scienceExposure
        @param scienceExposure: reference Exposure
        @param templateFwhmPix: FWHM (in pixels) of the Psf in the template image (image to convolve)
        @param scienceFwhmPix: FWHM (in pixels) of the Psf in the science image
        @param candidateList: a list of footprints/maskedImages for kernel candidates;
                              if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF
        @param doWarping: what to do if templateExposure's and scienceExposure's WCSs do not match:
            - if True then warp templateExposure to match scienceExposure
            - if False then raise an Exception
        @param convolveTemplate: convolve the template image or the science image
            - if True, templateExposure is warped if doWarping, templateExposure is convolved
            - if False, templateExposure is warped if doWarping, scienceExposure is convolved

        @return a pipeBase.Struct containing these fields:
        - subtractedExposure: subtracted Exposure = scienceExposure - (matchedImage + backgroundModel)
        - matchedImage: templateExposure after warping to match templateExposure (if doWarping true),
            and convolving with psfMatchingKernel
        - psfMatchingKernel: PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to determine PSF matching kernel
        """
        results = self.matchExposures(
            templateExposure=templateExposure,
            scienceExposure=scienceExposure,
            templateFwhmPix=templateFwhmPix,
            scienceFwhmPix=scienceFwhmPix,
            candidateList=candidateList,
            doWarping=doWarping,
            convolveTemplate=convolveTemplate
        )

        subtractedExposure = afwImage.ExposureF(scienceExposure, True)
        if convolveTemplate:
            subtractedMaskedImage  = subtractedExposure.getMaskedImage()
            subtractedMaskedImage -= results.matchedExposure.getMaskedImage()
            subtractedMaskedImage -= results.backgroundModel
        else:
            subtractedExposure.setMaskedImage(results.warpedExposure.getMaskedImage())
            subtractedMaskedImage  = subtractedExposure.getMaskedImage()
            subtractedMaskedImage -= results.matchedExposure.getMaskedImage()
            subtractedMaskedImage -= results.backgroundModel

            # Preserve polarity of differences
            subtractedMaskedImage *= -1

            # Place back on native photometric scale
            subtractedMaskedImage /= results.psfMatchingKernel.computeImage(
                afwImage.ImageD(results.psfMatchingKernel.getDimensions()), False)

        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayDiffIm = lsstDebug.Info(__name__).displayDiffIm
        maskTransparency = lsstDebug.Info(__name__).maskTransparency
        if not maskTransparency:
            maskTransparency = 0
        ds9.setMaskTransparency(maskTransparency)
        if display and displayDiffIm:
            ds9.mtv(templateExposure, frame=lsstDebug.frame, title="Template")
            lsstDebug.frame += 1
            ds9.mtv(results.matchedExposure, frame=lsstDebug.frame, title="Matched template")
            lsstDebug.frame += 1
            ds9.mtv(scienceExposure, frame=lsstDebug.frame, title="Science Image")
            lsstDebug.frame += 1
            ds9.mtv(subtractedExposure, frame=lsstDebug.frame, title="Difference Image")
            lsstDebug.frame += 1

        results.subtractedExposure = subtractedExposure
        return results

    @pipeBase.timeMethod
    def subtractMaskedImages(self, templateMaskedImage, scienceMaskedImage, candidateList,
            templateFwhmPix=None, scienceFwhmPix=None):
        """!Psf-match and subtract two MaskedImages

        Do the following, in order:
        - PSF-match templateMaskedImage to scienceMaskedImage
        - Determine the differential background
        - Return the difference: scienceMaskedImage -
            ((warped templateMaskedImage convolved with psfMatchingKernel) + backgroundModel)

        @param templateMaskedImage: MaskedImage to PSF-match to scienceMaskedImage
        @param scienceMaskedImage: reference MaskedImage
        @param templateFwhmPix: FWHM (in pixels) of the Psf in the template image (image to convolve)
        @param scienceFwhmPix: FWHM (in pixels) of the Psf in the science image
        @param candidateList: a list of footprints/maskedImages for kernel candidates;
                              if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF

        @return a pipeBase.Struct containing these fields:
        - subtractedMaskedImage = scienceMaskedImage - (matchedImage + backgroundModel)
        - matchedImage: templateMaskedImage convolved with psfMatchingKernel
        - psfMatchingKernel: PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to determine PSF matching kernel
        """
        if not candidateList:
            raise RuntimeError("Candidate list must be populated by makeCandidateList")

        results = self.matchMaskedImages(
            templateMaskedImage=templateMaskedImage,
            scienceMaskedImage=scienceMaskedImage,
            candidateList=candidateList,
            templateFwhmPix=templateFwhmPix,
            scienceFwhmPix=scienceFwhmPix,
            )

        subtractedMaskedImage  = afwImage.MaskedImageF(scienceMaskedImage, True)
        subtractedMaskedImage -= results.matchedImage
        subtractedMaskedImage -= results.backgroundModel
        results.subtractedMaskedImage = subtractedMaskedImage

        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayDiffIm = lsstDebug.Info(__name__).displayDiffIm
        maskTransparency = lsstDebug.Info(__name__).maskTransparency
        if not maskTransparency:
            maskTransparency = 0
        ds9.setMaskTransparency(maskTransparency)
        if display and displayDiffIm:
            ds9.mtv(subtractedMaskedImage, frame=lsstDebug.frame)
            lsstDebug.frame += 1

        return results

    def getSelectSources(self, exposure, sigma=None, doSmooth=True, idFactory=None):
        """!Get sources to use for Psf-matching

        This method runs detection and measurement on an exposure.
        The returned set of sources will be used as candidates for
        Psf-matching.

        @param exposure: Exposure on which to run detection/measurement
        @param sigma: Detection threshold
        @param doSmooth: Whether or not to smooth the Exposure with Psf before detection
        @param idFactory: Factory for the generation of Source ids

        @return source catalog containing candidates for the Psf-matching
        """

        if idFactory:
            table = afwTable.SourceTable.make(self.selectSchema, idFactory)
        else:
            table = afwTable.SourceTable.make(self.selectSchema)
        mi = exposure.getMaskedImage()

        imArr = mi.getImage().getArray()
        maskArr = mi.getMask().getArray()
        miArr = np.ma.masked_array(imArr, mask=maskArr)
        try:
            bkgd = getBackground(mi, self.kConfig.afwBackgroundConfig).getImageF()
        except:
            self.log.warn("Failed to get background model.  Falling back to median background estimation")
            bkgd = np.ma.extras.median(miArr)


        #Take off background for detection
        mi -= bkgd
        try:
            table.setMetadata(self.selectAlgMetadata) 
            detRet = self.selectDetection.makeSourceCatalog(
                table=table,
                exposure=exposure,
                sigma=sigma,
                doSmooth=doSmooth
                )
            selectSources = detRet.sources
            self.selectMeasurement.measure(exposure, selectSources)
        finally:
            # Put back on the background in case it is needed down stream
            mi += bkgd
            del bkgd
        return selectSources

    def makeCandidateList(self, templateExposure, scienceExposure, kernelSize, candidateList=None):
        """!Make a list of acceptable KernelCandidates

        Accept or generate a list of candidate sources for
        Psf-matching, and examine the Mask planes in both of the
        images for indications of bad pixels

        @param templateExposure: Exposure that will be convolved
        @param scienceExposure: Exposure that will be matched-to
        @param kernelSize: Dimensions of the Psf-matching Kernel, used to grow detection footprints
        @param candidateList: List of Sources to examine

        @return a list of dicts having a "source" and "footprint"
        field for the Sources deemed to be appropriate for Psf
        matching
        """
        if candidateList is None:
            candidateList = self.getSelectSources(scienceExposure)

        listTypes = set(type(x) for x in candidateList)
        if (not len(listTypes) == 1) or (type(listTypes.pop()) != type(afwTable.SourceRecord)):
            raise RuntimeError("Can only make candidate list from set of SourceRecords.  Got %s instead." \
                                   % (type(candidateList[0])))
        candidateList = diffimTools.sourceToFootprintList(candidateList,
                                                          templateExposure, scienceExposure,
                                                          kernelSize,
                                                          self.kConfig.detectionConfig,
                                                          self.log)
        if len(candidateList) == 0:
            raise RuntimeError("Cannot find any objects suitable for KernelCandidacy")

        return candidateList

    def _adaptCellSize(self, candidateList):
        """! NOT IMPLEMENTED YET"""
        nCand = len(candidateList)
        return self.kConfig.sizeCellX, self.kConfig.sizeCellY

    def _buildCellSet(self, templateMaskedImage, scienceMaskedImage, candidateList):
        """!Build a SpatialCellSet for use with the solve method

        @param templateMaskedImage: MaskedImage to PSF-matched to scienceMaskedImage
        @param scienceMaskedImage: reference MaskedImage
        @param candidateList: a list of footprints/maskedImages for kernel candidates;
                              if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF

        @return kernelCellSet: a SpatialCellSet for use with self._solve
        """
        if not candidateList:
            raise RuntimeError("Candidate list must be populated by makeCandidateList")

        sizeCellX, sizeCellY = self._adaptCellSize(candidateList)

        # Object to store the KernelCandidates for spatial modeling
        kernelCellSet = afwMath.SpatialCellSet(templateMaskedImage.getBBox(),
                                               sizeCellX, sizeCellY)

        policy = pexConfig.makePolicy(self.kConfig)
        # Place candidates within the spatial grid
        for cand in candidateList:
            bbox = cand['footprint'].getBBox()

            tmi  = afwImage.MaskedImageF(templateMaskedImage, bbox)
            smi  = afwImage.MaskedImageF(scienceMaskedImage, bbox)
            cand = diffimLib.makeKernelCandidate(cand['source'], tmi, smi, policy)

            self.log.logdebug("Candidate %d at %f, %f" % (cand.getId(), cand.getXCenter(), cand.getYCenter()))
            kernelCellSet.insertCandidate(cand)

        return kernelCellSet

    def _validateSize(self, templateMaskedImage, scienceMaskedImage):
        """!Return True if two image-like objects are the same size
        """
        return templateMaskedImage.getDimensions() == scienceMaskedImage.getDimensions()

    def _validateWcs(self, templateExposure, scienceExposure):
        """!Return True if the WCS of the two Exposures have the same origin and extent
        """
        templateWcs    = templateExposure.getWcs() 
        scienceWcs     = scienceExposure.getWcs()
        templateBBox   = templateExposure.getBBox()
        scienceBBox    = scienceExposure.getBBox()

        # LLC
        templateOrigin = templateWcs.pixelToSky(afwGeom.Point2D(templateBBox.getBegin()))
        scienceOrigin  = scienceWcs.pixelToSky(afwGeom.Point2D(scienceBBox.getBegin()))

        # URC
        templateLimit  = templateWcs.pixelToSky(afwGeom.Point2D(templateBBox.getEnd()))
        scienceLimit   = scienceWcs.pixelToSky(afwGeom.Point2D(scienceBBox.getEnd()))

        self.log.info("Template Wcs : %f,%f -> %f,%f" %
                      (templateOrigin[0], templateOrigin[1],
                       templateLimit[0], templateLimit[1]))
        self.log.info("Science Wcs : %f,%f -> %f,%f" %
                      (scienceOrigin[0], scienceOrigin[1],
                       scienceLimit[0], scienceLimit[1]))

        templateBBox = afwGeom.Box2D(templateOrigin.getPosition(), templateLimit.getPosition())
        scienceBBox  = afwGeom.Box2D(scienceOrigin.getPosition(), scienceLimit.getPosition())
        if not (templateBBox.overlaps(scienceBBox)):
            raise RuntimeError("Input images do not overlap at all")

        if ( (templateOrigin.getPosition() != scienceOrigin.getPosition()) or
             (templateLimit.getPosition()  != scienceLimit.getPosition())  or
             (templateExposure.getDimensions() != scienceExposure.getDimensions())):
            return False
        return True
