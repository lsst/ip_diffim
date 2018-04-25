# LSST Data Management System
# Copyright 2008-2016 LSST Corporation.
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

from . import diffimLib
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.log as log
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from .makeKernelBasisList import makeKernelBasisList
from .psfMatch import PsfMatchTask, PsfMatchConfigAL
from . import utils as dituils
import lsst.afw.display.ds9 as ds9

__all__ = ("ModelPsfMatchTask", "ModelPsfMatchConfig")

sigma2fwhm = 2. * np.sqrt(2. * np.log(2.))


def nextOddInteger(x):
    nextInt = int(np.ceil(x))
    return nextInt + 1 if nextInt % 2 == 0 else nextInt


class ModelPsfMatchConfig(pexConfig.Config):
    """!Configuration for model-to-model Psf matching"""

    kernel = pexConfig.ConfigChoiceField(
        doc="kernel type",
        typemap=dict(
            AL=PsfMatchConfigAL,
        ),
        default="AL",
    )
    doAutoPadPsf = pexConfig.Field(
        dtype=bool,
        doc=("If too small, automatically pad the science Psf? "
             "Pad to smallest dimensions appropriate for the matching kernel dimensions, "
             "as specified by autoPadPsfTo. If false, pad by the padPsfBy config."),
        default=True,
    )
    autoPadPsfTo = pexConfig.RangeField(
        dtype=float,
        doc=("Minimum Science Psf dimensions as a fraction of matching kernel dimensions. "
             "If the dimensions of the Psf to be matched are less than the "
             "matching kernel dimensions * autoPadPsfTo, pad Science Psf to this size. "
             "Ignored if doAutoPadPsf=False."),
        default=1.4,
        min=1.0,
        max=2.0
    )
    padPsfBy = pexConfig.Field(
        dtype=int,
        doc="Pixels (even) to pad Science Psf by before matching. Ignored if doAutoPadPsf=True",
        default=0,
    )

    def setDefaults(self):
        # No sigma clipping
        self.kernel.active.singleKernelClipping = False
        self.kernel.active.kernelSumClipping = False
        self.kernel.active.spatialKernelClipping = False
        self.kernel.active.checkConditionNumber = False

        # Variance is ill defined
        self.kernel.active.constantVarianceWeighting = True

        # Do not change specified kernel size
        self.kernel.active.scaleByFwhm = False


## @addtogroup LSST_task_documentation
## @{
## @page ModelPsfMatchTask
## @ref ModelPsfMatchTask_ "ModelPsfMatchTask"
## @copybrief ModelPsfMatchTask
## @}


class ModelPsfMatchTask(PsfMatchTask):
    """!
@anchor ModelPsfMatchTask_

@brief Matching of two model Psfs, and application of the Psf-matching kernel to an input Exposure

@section ip_diffim_modelpsfmatch_Contents Contents

 - @ref ip_diffim_modelpsfmatch_Purpose
 - @ref ip_diffim_modelpsfmatch_Initialize
 - @ref ip_diffim_modelpsfmatch_IO
 - @ref ip_diffim_modelpsfmatch_Config
 - @ref ip_diffim_modelpsfmatch_Metadata
 - @ref ip_diffim_modelpsfmatch_Debug
 - @ref ip_diffim_modelpsfmatch_Example

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

@section ip_diffim_modelpsfmatch_Purpose   Description

This Task differs from ImagePsfMatchTask in that it matches two Psf _models_, by realizing
them in an Exposure-sized SpatialCellSet and then inserting each Psf-image pair into KernelCandidates.
Because none of the pairs of sources that are to be matched should be invalid, all sigma clipping is
turned off in ModelPsfMatchConfig.  And because there is no tracked _variance_ in the Psf images, the
debugging and logging QA info should be interpreted with caution.

One item of note is that the sizes of Psf models are fixed (e.g. its defined as a 21x21 matrix).  When the
Psf-matching kernel is being solved for, the Psf "image" is convolved with each kernel basis function,
leading to a loss of information around the borders.  This pixel loss will be problematic for the numerical
stability of the kernel solution if the size of the convolution kernel (set by ModelPsfMatchConfig.kernelSize)
is much bigger than: psfSize//2.  Thus the sizes of Psf-model matching kernels are typically smaller
than their image-matching counterparts.  If the size of the kernel is too small, the convolved stars will
look "boxy"; if the kernel is too large, the kernel solution will be "noisy".  This is a trade-off that
needs careful attention for a given dataset.

The primary use case for this Task is in matching an Exposure to a constant-across-the-sky Psf model for the
purposes of image coaddition.  It is important to note that in the code, the "template" Psf is the Psf
that the science image gets matched to.  In this sense the order of template and science image are
reversed, compared to ImagePsfMatchTask, which operates on the template image.

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

@section ip_diffim_modelpsfmatch_Initialize    Task initialization

@copydoc \_\_init\_\_

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

@section ip_diffim_modelpsfmatch_IO        Invoking the Task

@copydoc run

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

@section ip_diffim_modelpsfmatch_Config       Configuration parameters

See @ref ModelPsfMatchConfig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

@section ip_diffim_modelpsfmatch_Metadata   Quantities set in Metadata

See @ref ip_diffim_psfmatch_Metadata "PsfMatchTask"

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

@section ip_diffim_modelpsfmatch_Debug     Debug variables

The @link lsst.pipe.base.cmdLineTask.CmdLineTask command line task@endlink interface supports a
flag @c -d/--debug to import @b debug.py from your @c PYTHONPATH.  The relevant contents of debug.py
for this Task include:

@code{.py}
    import sys
    import lsstDebug
    def DebugInfo(name):
        di = lsstDebug.getInfo(name)
        if name == "lsst.ip.diffim.psfMatch":
            di.display = True                 # global
            di.maskTransparency = 80          # ds9 mask transparency
            di.displayCandidates = True       # show all the candidates and residuals
            di.displayKernelBasis = False     # show kernel basis functions
            di.displayKernelMosaic = True     # show kernel realized across the image
            di.plotKernelSpatialModel = False # show coefficients of spatial model
            di.showBadCandidates = True       # show the bad candidates (red) along with good (green)
        elif name == "lsst.ip.diffim.modelPsfMatch":
            di.display = True                 # global
            di.maskTransparency = 30          # ds9 mask transparency
            di.displaySpatialCells = True     # show spatial cells before the fit
        return di
    lsstDebug.Info = DebugInfo
    lsstDebug.frame = 1
@endcode

Note that if you want addional logging info, you may add to your scripts:
@code{.py}
import lsst.log.utils as logUtils
logUtils.traceSetAt("ip.diffim", 4)
@endcode

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

@section ip_diffim_modelpsfmatch_Example   A complete example of using ModelPsfMatchTask

This code is modelPsfMatchTask.py in the examples directory, and can be run as @em e.g.
@code
examples/modelPsfMatchTask.py
examples/modelPsfMatchTask.py --debug
examples/modelPsfMatchTask.py --debug --template /path/to/templateExp.fits --science /path/to/scienceExp.fits
@endcode

@dontinclude modelPsfMatchTask.py
Create a subclass of ModelPsfMatchTask that accepts two exposures.  Note that the "template" exposure
contains the Psf that will get matched to, and the "science" exposure is the one that will be convolved:
@skip MyModelPsfMatchTask
@until return

And allow the user the freedom to either run the script in default mode, or point to their own images on disk.
Note that these images must be readable as an lsst.afw.image.Exposure:
@skip main
@until parse_args

We have enabled some minor display debugging in this script via the --debug option.  However, if you
have an lsstDebug debug.py in your PYTHONPATH you will get additional debugging displays.  The following
block checks for this script:
@skip args.debug
@until sys.stderr

@dontinclude modelPsfMatchTask.py
Finally, we call a run method that we define below.  First set up a Config and modify some of the parameters.
In particular we don't want to "grow" the sizes of the kernel or KernelCandidates, since we are operating with
fixed--size images (i.e. the size of the input Psf models).
@skip run(args)
@until False

Make sure the images (if any) that were sent to the script exist on disk and are readable.  If no images
are sent, make some fake data up for the sake of this example script (have a look at the code if you want
more details on generateFakeData):
@skip requested
@until sizeCellY

Display the two images if --debug:
@skip args.debug
@until Science

Create and run the Task:
@skip Create
@until result

And finally provide optional debugging display of the Psf-matched (via the Psf models) science image:
@skip args.debug
@until result.psfMatchedExposure

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    """
    ConfigClass = ModelPsfMatchConfig

    def __init__(self, *args, **kwargs):
        """!Create a ModelPsfMatchTask

        @param *args arguments to be passed to lsst.ip.diffim.PsfMatchTask.__init__
        @param **kwargs keyword arguments to be passed to lsst.ip.diffim.PsfMatchTask.__init__

        Upon initialization, the kernel configuration is defined by self.config.kernel.active.  This Task
        does have a run() method, which is the default way to call the Task.
        """
        PsfMatchTask.__init__(self, *args, **kwargs)
        self.kConfig = self.config.kernel.active

    @pipeBase.timeMethod
    def run(self, exposure, referencePsfModel, kernelSum=1.0):
        """!Psf-match an exposure to a model Psf

        @param exposure: Exposure to Psf-match to the reference Psf model;
            it must return a valid PSF model via exposure.getPsf()
        @param referencePsfModel: The Psf model to match to (an lsst.afw.detection.Psf)
        @param kernelSum: A multipicative factor to apply to the kernel sum (default=1.0)

        @return
        - psfMatchedExposure: the Psf-matched Exposure.  This has the same parent bbox, Wcs, Calib and
            Filter as the input Exposure but no Psf.  In theory the Psf should equal referencePsfModel but
            the match is likely not exact.
        - psfMatchingKernel: the spatially varying Psf-matching kernel
        - kernelCellSet: SpatialCellSet used to solve for the Psf-matching kernel
        - referencePsfModel: Validated and/or modified reference model used

        Raise a RuntimeError if the Exposure does not contain a Psf model
        """
        if not exposure.hasPsf():
            raise RuntimeError("exposure does not contain a Psf model")

        maskedImage = exposure.getMaskedImage()

        self.log.info("compute Psf-matching kernel")
        result = self._buildCellSet(exposure, referencePsfModel)
        kernelCellSet = result.kernelCellSet
        referencePsfModel = result.referencePsfModel
        fwhmScience = exposure.getPsf().computeShape().getDeterminantRadius() * sigma2fwhm
        fwhmModel = referencePsfModel.computeShape().getDeterminantRadius() * sigma2fwhm

        basisList = makeKernelBasisList(self.kConfig, fwhmScience, fwhmModel, metadata=self.metadata)
        spatialSolution, psfMatchingKernel, backgroundModel = self._solve(kernelCellSet, basisList)

        if psfMatchingKernel.isSpatiallyVarying():
            sParameters = np.array(psfMatchingKernel.getSpatialParameters())
            sParameters[0][0] = kernelSum
            psfMatchingKernel.setSpatialParameters(sParameters)
        else:
            kParameters = np.array(psfMatchingKernel.getKernelParameters())
            kParameters[0] = kernelSum
            psfMatchingKernel.setKernelParameters(kParameters)

        self.log.info("Psf-match science exposure to reference")
        psfMatchedExposure = afwImage.ExposureF(exposure.getBBox(), exposure.getWcs())
        psfMatchedExposure.setFilter(exposure.getFilter())
        psfMatchedExposure.setCalib(exposure.getCalib())
        psfMatchedExposure.getInfo().setVisitInfo(exposure.getInfo().getVisitInfo())
        psfMatchedExposure.setPsf(referencePsfModel)
        psfMatchedMaskedImage = psfMatchedExposure.getMaskedImage()

        # Normalize the psf-matching kernel while convolving since its magnitude is meaningless
        # when PSF-matching one model to another.
        doNormalize = True
        afwMath.convolve(psfMatchedMaskedImage, maskedImage, psfMatchingKernel, doNormalize)

        self.log.info("done")
        return pipeBase.Struct(psfMatchedExposure=psfMatchedExposure,
                               psfMatchingKernel=psfMatchingKernel,
                               kernelCellSet=kernelCellSet,
                               metadata=self.metadata,
                               )

    def _diagnostic(self, kernelCellSet, spatialSolution, spatialKernel, spatialBg):
        """!Print diagnostic information on spatial kernel and background fit

        The debugging diagnostics are not really useful here, since the images we are matching have
        no variance.  Thus override the _diagnostic method to generate no logging information"""
        return

    def _buildCellSet(self, exposure, referencePsfModel):
        """!Build a SpatialCellSet for use with the solve method

        @param exposure: The science exposure that will be convolved; must contain a Psf
        @param referencePsfModel: Psf model to match to

        @return
            -kernelCellSet: a SpatialCellSet to be used by self._solve
            -referencePsfModel: Validated and/or modified reference model used to populate the SpatialCellSet

        If the reference Psf model and science Psf model have different dimensions,
        adjust the referencePsfModel (the model to which the exposure PSF will be matched)
        to match that of the science Psf. If the science Psf dimensions vary across the image,
        as is common with a WarpedPsf, either pad or clip (depending on config.padPsf)
        the dimensions to be constant.
        """
        sizeCellX = self.kConfig.sizeCellX
        sizeCellY = self.kConfig.sizeCellY

        scienceBBox = exposure.getBBox()
        # Extend for proper spatial matching kernel all the way to edge, especially for narrow strips
        scienceBBox.grow(afwGeom.Extent2I(sizeCellX, sizeCellY))

        sciencePsfModel = exposure.getPsf()

        dimenR = referencePsfModel.getLocalKernel().getDimensions()
        psfWidth, psfHeight = dimenR

        regionSizeX, regionSizeY = scienceBBox.getDimensions()
        scienceX0, scienceY0 = scienceBBox.getMin()

        kernelCellSet = afwMath.SpatialCellSet(afwGeom.Box2I(scienceBBox), sizeCellX, sizeCellY)

        nCellX = regionSizeX//sizeCellX
        nCellY = regionSizeY//sizeCellY

        if nCellX == 0 or nCellY == 0:
            raise ValueError("Exposure dimensions=%s and sizeCell=(%s, %s). Insufficient area to match" %
                             (scienceBBox.getDimensions(), sizeCellX, sizeCellY))

        # Survey the PSF dimensions of the Spatial Cell Set
        # to identify the minimum enclosed or maximum bounding square BBox.
        widthList = []
        heightList = []
        for row in range(nCellY):
            posY = sizeCellY*row + sizeCellY//2 + scienceY0
            for col in range(nCellX):
                posX = sizeCellX*col + sizeCellX//2 + scienceX0
                widthS, heightS = sciencePsfModel.computeBBox(afwGeom.Point2D(posX, posY)).getDimensions()
                widthList.append(widthS)
                heightList.append(heightS)

        psfSize = max(max(heightList), max(widthList))

        if self.config.doAutoPadPsf:
            minPsfSize = nextOddInteger(self.kConfig.kernelSize*self.config.autoPadPsfTo)
            paddingPix = max(0, minPsfSize - psfSize)
        else:
            if self.config.padPsfBy % 2 != 0:
                raise ValueError("Config padPsfBy (%i pixels) must be even number." %
                                 self.config.padPsfBy)
            paddingPix = self.config.padPsfBy

        if paddingPix > 0:
            self.log.info("Padding Science PSF from (%s, %s) to (%s, %s) pixels" %
                          (psfSize, psfSize, paddingPix + psfSize, paddingPix + psfSize))
            psfSize += paddingPix

        # Check that PSF is larger than the matching kernel
        maxKernelSize = psfSize - 1
        if maxKernelSize % 2 == 0:
            maxKernelSize -= 1
        if self.kConfig.kernelSize > maxKernelSize:
            message = """
                Kernel size (%d) too big to match Psfs of size %d.
                Please reconfigure by setting one of the following:
                1) kernel size to <= %d
                2) doAutoPadPsf=True
                3) padPsfBy to >= %s
                """ % (self.kConfig.kernelSize, psfSize,
                       maxKernelSize, self.kConfig.kernelSize - maxKernelSize)
            raise ValueError(message)

        dimenS = afwGeom.Extent2I(psfSize, psfSize)

        if (dimenR != dimenS):
            try:
                referencePsfModel = referencePsfModel.resized(psfSize, psfSize)
                self.log.info("Adjusted dimensions of reference PSF model from %s to %s" % (dimenR, dimenS))
            except Exception as e:
                self.log.warn("Zero padding or clipping the reference PSF model of type %s and dimensions %s"
                              " to the science Psf dimensions %s because: %s",
                              referencePsfModel.__class__.__name__, dimenR, dimenS, e)
            dimenR = dimenS

        policy = pexConfig.makePolicy(self.kConfig)
        for row in range(nCellY):
            # place at center of cell
            posY = sizeCellY * row + sizeCellY//2 + scienceY0

            for col in range(nCellX):
                # place at center of cell
                posX = sizeCellX * col + sizeCellX//2 + scienceX0

                log.log("TRACE4." + self.log.getName(), log.DEBUG,
                        "Creating Psf candidate at %.1f %.1f", posX, posY)

                # reference kernel image, at location of science subimage
                referenceMI = self._makePsfMaskedImage(referencePsfModel, posX, posY, dimensions=dimenR)

                # kernel image we are going to convolve
                scienceMI = self._makePsfMaskedImage(sciencePsfModel, posX, posY, dimensions=dimenR)

                # The image to convolve is the science image, to the reference Psf.
                kc = diffimLib.makeKernelCandidate(posX, posY, scienceMI, referenceMI, policy)
                kernelCellSet.insertCandidate(kc)

        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displaySpatialCells = lsstDebug.Info(__name__).displaySpatialCells
        maskTransparency = lsstDebug.Info(__name__).maskTransparency
        if not maskTransparency:
            maskTransparency = 0
        if display:
            ds9.setMaskTransparency(maskTransparency)
        if display and displaySpatialCells:
            dituils.showKernelSpatialCells(exposure.getMaskedImage(), kernelCellSet,
                                           symb="o", ctype=ds9.CYAN, ctypeUnused=ds9.YELLOW, ctypeBad=ds9.RED,
                                           size=4, frame=lsstDebug.frame, title="Image to be convolved")
            lsstDebug.frame += 1
        return pipeBase.Struct(kernelCellSet=kernelCellSet,
                               referencePsfModel=referencePsfModel,
                               )

    def _makePsfMaskedImage(self, psfModel, posX, posY, dimensions=None):
        """! Return a MaskedImage of the a PSF Model of specified dimensions
        """
        rawKernel = psfModel.computeKernelImage(afwGeom.Point2D(posX, posY)).convertF()
        if dimensions is None:
            dimensions = rawKernel.getDimensions()
        if rawKernel.getDimensions() == dimensions:
            kernelIm = rawKernel
        else:
            # make image of proper size
            kernelIm = afwImage.ImageF(dimensions)
            bboxToPlace = afwGeom.Box2I(afwGeom.Point2I((dimensions.getX() - rawKernel.getWidth())//2,
                                                        (dimensions.getY() - rawKernel.getHeight())//2),
                                        rawKernel.getDimensions())
            kernelIm.assign(rawKernel, bboxToPlace)

        kernelMask = afwImage.Mask(dimensions, 0x0)
        kernelVar = afwImage.ImageF(dimensions, 1.0)
        return afwImage.MaskedImageF(kernelIm, kernelMask, kernelVar)
