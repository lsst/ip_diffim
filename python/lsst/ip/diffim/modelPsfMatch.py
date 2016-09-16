from __future__ import absolute_import
from builtins import range
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
import numpy as num
from . import diffimLib
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.logging as pexLog
import lsst.pex.config as pexConfig
import lsst.meas.algorithms as measAlg
import lsst.pipe.base as pipeBase
from .makeKernelBasisList import makeKernelBasisList
from .psfMatch import PsfMatchTask, PsfMatchConfigAL
from . import utils as diUtils
import lsst.afw.display.ds9 as ds9

sigma2fwhm = 2. * num.sqrt(2. * num.log(2.))


class ModelPsfMatchConfig(pexConfig.Config):
    """!Configuration for model-to-model Psf matching"""

    kernel = pexConfig.ConfigChoiceField(
        doc="kernel type",
        typemap=dict(
            AL=PsfMatchConfigAL,
        ),
        default="AL",
    )

    def setDefaults(self):
        # No sigma clipping
        self.kernel.active.singleKernelClipping = False
        self.kernel.active.kernelSumClipping = False
        self.kernel.active.spatialKernelClipping = False
        self.kernel.active.checkConditionNumber = False

        # Variance is ill defined
        self.kernel.active.constantVarianceWeighting = True

        # Psfs are typically small; reduce the kernel size
        self.kernel.active.kernelSizeMin = 11
        self.kernel.active.kernelSize = 11

## \addtogroup LSST_task_documentation
## \{
## \page ModelPsfMatchTask
## \ref ModelPsfMatchTask_ "ModelPsfMatchTask"
## \copybrief ModelPsfMatchTask
## \}


class ModelPsfMatchTask(PsfMatchTask):
    """!
\anchor ModelPsfMatchTask_

\brief Matching of two model Psfs, and application of the Psf-matching kernel to an input Exposure

\section ip_diffim_modelpsfmatch_Contents Contents

 - \ref ip_diffim_modelpsfmatch_Purpose
 - \ref ip_diffim_modelpsfmatch_Initialize
 - \ref ip_diffim_modelpsfmatch_IO
 - \ref ip_diffim_modelpsfmatch_Config
 - \ref ip_diffim_modelpsfmatch_Metadata
 - \ref ip_diffim_modelpsfmatch_Debug
 - \ref ip_diffim_modelpsfmatch_Example

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_modelpsfmatch_Purpose   Description

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

\section ip_diffim_modelpsfmatch_Initialize    Task initialization

\copydoc \_\_init\_\_

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_modelpsfmatch_IO        Invoking the Task

\copydoc run

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_modelpsfmatch_Config       Configuration parameters

See \ref ModelPsfMatchConfig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_modelpsfmatch_Metadata   Quantities set in Metadata

See \ref ip_diffim_psfmatch_Metadata "PsfMatchTask"

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_modelpsfmatch_Debug     Debug variables

The \link lsst.pipe.base.cmdLineTask.CmdLineTask command line task\endlink interface supports a
flag \c -d/--debug to import \b debug.py from your \c PYTHONPATH.  The relevant contents of debug.py
for this Task include:

\code{.py}
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
\endcode

Note that if you want addional logging info, you may add to your scripts:
\code{.py}
import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)
\endcode

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_modelpsfmatch_Example   A complete example of using ModelPsfMatchTask

This code is modelPsfMatchTask.py in the examples directory, and can be run as \em e.g.
\code
examples/modelPsfMatchTask.py
examples/modelPsfMatchTask.py --debug
examples/modelPsfMatchTask.py --debug --template /path/to/templateExp.fits --science /path/to/scienceExp.fits
\endcode

\dontinclude modelPsfMatchTask.py
Create a subclass of ModelPsfMatchTask that accepts two exposures.  Note that the "template" exposure
contains the Psf that will get matched to, and the "science" exposure is the one that will be convolved:
\skip MyModelPsfMatchTask
@until return

And allow the user the freedom to either run the script in default mode, or point to their own images on disk.
Note that these images must be readable as an lsst.afw.image.Exposure:
\skip main
@until parse_args

We have enabled some minor display debugging in this script via the --debug option.  However, if you
have an lsstDebug debug.py in your PYTHONPATH you will get additional debugging displays.  The following
block checks for this script:
\skip args.debug
@until sys.stderr

\dontinclude modelPsfMatchTask.py
Finally, we call a run method that we define below.  First set up a Config and modify some of the parameters.
In particular we don't want to "grow" the sizes of the kernel or KernelCandidates, since we are operating with
fixed--size images (i.e. the size of the input Psf models).
\skip run(args)
@until False

Make sure the images (if any) that were sent to the script exist on disk and are readable.  If no images
are sent, make some fake data up for the sake of this example script (have a look at the code if you want
more details on generateFakeData):
\skip requested
@until sizeCellY

Display the two images if --debug:
\skip args.debug
@until Science

Create and run the Task:
\skip Create
@until result

And finally provide optional debugging display of the Psf-matched (via the Psf models) science image:
\skip args.debug
@until result.psfMatchedExposure

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    """
    ConfigClass = ModelPsfMatchConfig

    def __init__(self, *args, **kwargs):
        """!Create a ModelPsfMatchTask

        \param *args arguments to be passed to lsst.ip.diffim.PsfMatchTask.__init__
        \param **kwargs keyword arguments to be passed to lsst.ip.diffim.PsfMatchTask.__init__

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

        Raise a RuntimeError if the Exposure does not contain a Psf model
        """
        if not exposure.hasPsf():
            raise RuntimeError("exposure does not contain a Psf model")

        maskedImage = exposure.getMaskedImage()

        self.log.info("compute Psf-matching kernel")
        kernelCellSet = self._buildCellSet(exposure, referencePsfModel)
        width, height = referencePsfModel.getLocalKernel().getDimensions()
        psfAttr1 = measAlg.PsfAttributes(exposure.getPsf(), width//2, height//2)
        psfAttr2 = measAlg.PsfAttributes(referencePsfModel, width//2, height//2)
        s1 = psfAttr1.computeGaussianWidth(psfAttr1.ADAPTIVE_MOMENT) # gaussian sigma in pixels
        s2 = psfAttr2.computeGaussianWidth(psfAttr2.ADAPTIVE_MOMENT) # gaussian sigma in pixels
        fwhm1 = s1 * sigma2fwhm # science Psf
        fwhm2 = s2 * sigma2fwhm # template Psf

        basisList = makeKernelBasisList(self.kConfig, fwhm1, fwhm2, metadata=self.metadata)
        spatialSolution, psfMatchingKernel, backgroundModel = self._solve(kernelCellSet, basisList)

        if psfMatchingKernel.isSpatiallyVarying():
            sParameters = num.array(psfMatchingKernel.getSpatialParameters())
            sParameters[0][0] = kernelSum
            psfMatchingKernel.setSpatialParameters(sParameters)
        else:
            kParameters = num.array(psfMatchingKernel.getKernelParameters())
            kParameters[0] = kernelSum
            psfMatchingKernel.setKernelParameters(kParameters)

        self.log.info("Psf-match science exposure to reference")
        psfMatchedExposure = afwImage.ExposureF(exposure.getBBox(), exposure.getWcs())
        psfMatchedExposure.setFilter(exposure.getFilter())
        psfMatchedExposure.setCalib(exposure.getCalib())
        psfMatchedMaskedImage = psfMatchedExposure.getMaskedImage()

        # Normalize the psf-matching kernel while convolving since its magnitude is meaningless
        # when PSF-matching one model to another.
        doNormalize = True
        afwMath.convolve(psfMatchedMaskedImage, maskedImage, psfMatchingKernel, doNormalize)

        self.log.info("done")
        return pipeBase.Struct(psfMatchedExposure=psfMatchedExposure,
                               psfMatchingKernel=psfMatchingKernel,
                               kernelCellSet=kernelCellSet,
                               metadata=self.metadata)

    def _diagnostic(self, kernelCellSet, spatialSolution, spatialKernel, spatialBg):
        """!Print diagnostic information on spatial kernel and background fit

        The debugging diagnostics are not really useful here, since the images we are matching have
        no variance.  Thus override the _diagnostic method to generate no logging information"""
        return

    def _buildCellSet(self, exposure, referencePsfModel):
        """!Build a SpatialCellSet for use with the solve method

        @param exposure: The science exposure that will be convolved; must contain a Psf
        @param referencePsfModel: Psf model to match to

        @return kernelCellSet: a SpatialCellSet to be used by self._solve

        Raise a RuntimeError if the reference Psf model and science Psf model have different dimensions
        """
        scienceBBox = exposure.getBBox()
        sciencePsfModel = exposure.getPsf()
        # The Psf base class does not support getKernel() in general, as there are some Psf
        # classes for which this is not meaningful.
        # Many Psfs we use in practice are KernelPsfs, and this algorithm will work fine for them,
        # but someday it should probably be modified to support arbitrary Psfs.
        referencePsfModel = measAlg.KernelPsf.swigConvert(referencePsfModel)
        sciencePsfModel = measAlg.KernelPsf.swigConvert(sciencePsfModel)
        if referencePsfModel is None or sciencePsfModel is None:
            raise RuntimeError("ERROR: Psf matching is only implemented for KernelPsfs")
        if (referencePsfModel.getKernel().getDimensions() != sciencePsfModel.getKernel().getDimensions()):
            pexLog.Trace(self.log.getName(), 1,
                         "ERROR: Dimensions of reference Psf and science Psf different; exiting")
            raise RuntimeError("ERROR: Dimensions of reference Psf and science Psf different; exiting")

        psfWidth, psfHeight = referencePsfModel.getKernel().getDimensions()
        maxKernelSize = min(psfWidth, psfHeight) - 1
        if maxKernelSize % 2 == 0:
            maxKernelSize -= 1
        if self.kConfig.kernelSize > maxKernelSize:
            raise ValueError("Kernel size (%d) too big to match Psfs of size %d; reduce to at least %d" % (
                self.kConfig.kernelSize, psfWidth, maxKernelSize))

        # Infer spatial order of Psf model!
        #
        # Infer from the number of spatial parameters.
        # (O + 1) * (O + 2) / 2 = N
        # O^2 + 3 * O + 2 * (1 - N) = 0
        #
        # Roots are [-3 +/- sqrt(9 - 8 * (1 - N))] / 2
        #
        nParameters = sciencePsfModel.getKernel().getNSpatialParameters()
        root = num.sqrt(9 - 8 * (1 - nParameters))
        if (root != root // 1):            # We know its an integer solution
            pexLog.Trace(self.log.getName(), 3, "Problem inferring spatial order of image's Psf")
        else:
            order = (root - 3) / 2
            if (order != order // 1):
                pexLog.Trace(self.log.getName(), 3, "Problem inferring spatial order of image's Psf")
            else:
                pexLog.Trace(self.log.getName(), 2, "Spatial order of Psf = %d; matching kernel order = %d" %
                             (order, self.kConfig.spatialKernelOrder))

        regionSizeX, regionSizeY = scienceBBox.getDimensions()
        scienceX0, scienceY0 = scienceBBox.getMin()

        sizeCellX = self.kConfig.sizeCellX
        sizeCellY = self.kConfig.sizeCellY

        kernelCellSet = afwMath.SpatialCellSet(
            afwGeom.Box2I(afwGeom.Point2I(scienceX0, scienceY0),
                          afwGeom.Extent2I(regionSizeX, regionSizeY)),
            sizeCellX, sizeCellY
        )

        nCellX = regionSizeX // sizeCellX
        nCellY = regionSizeY // sizeCellY
        dimenR = referencePsfModel.getKernel().getDimensions()
        dimenS = sciencePsfModel.getKernel().getDimensions()

        policy = pexConfig.makePolicy(self.kConfig)
        for row in range(nCellY):
            # place at center of cell
            posY = sizeCellY * row + sizeCellY // 2 + scienceY0

            for col in range(nCellX):
                # place at center of cell
                posX = sizeCellX * col + sizeCellX // 2 + scienceX0

                pexLog.Trace(self.log.getName(), 5, "Creating Psf candidate at %.1f %.1f" % (posX, posY))

                # reference kernel image, at location of science subimage
                kernelImageR = referencePsfModel.computeImage(afwGeom.Point2D(posX, posY)).convertF()
                kernelMaskR = afwImage.MaskU(dimenR)
                kernelMaskR.set(0)
                kernelVarR = afwImage.ImageF(dimenR)
                kernelVarR.set(1.0)
                referenceMI = afwImage.MaskedImageF(kernelImageR, kernelMaskR, kernelVarR)

                # kernel image we are going to convolve
                kernelImageS = sciencePsfModel.computeImage(afwGeom.Point2D(posX, posY)).convertF()
                kernelMaskS = afwImage.MaskU(dimenS)
                kernelMaskS.set(0)
                kernelVarS = afwImage.ImageF(dimenS)
                kernelVarS.set(1.0)
                scienceMI = afwImage.MaskedImageF(kernelImageS, kernelMaskS, kernelVarS)

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
            diUtils.showKernelSpatialCells(exposure.getMaskedImage(), kernelCellSet,
                                           symb="o", ctype=ds9.CYAN, ctypeUnused=ds9.YELLOW, ctypeBad=ds9.RED,
                                           size=4, frame=lsstDebug.frame, title="Image to be convolved")
            lsstDebug.frame += 1
        return kernelCellSet
