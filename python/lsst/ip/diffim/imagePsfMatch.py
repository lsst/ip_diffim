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
import numpy as num
import lsst.pex.logging as pexLog
import lsst.pex.config as pexConfig
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.detection as afwDetect
import lsst.pipe.base as pipeBase
import lsst.meas.algorithms as measAlg
from .makeKernelBasisList import makeKernelBasisList
from .psfMatch import PsfMatch, PsfMatchConfigDF, PsfMatchConfigAL
from . import utils as diUtils 
from . import diffimLib
from . import diffimTools
import lsst.afw.display.ds9 as ds9

sigma2fwhm = 2. * num.sqrt(2. * num.log(2.))

class ImagePsfMatchConfig(pexConfig.Config):
    kernel = pexConfig.ConfigChoiceField(
        doc = "kernel type",
        typemap = dict(
            AL = PsfMatchConfigAL,
            DF = PsfMatchConfigDF
        ),
        default = "AL",
    )

class ImagePsfMatchTask(PsfMatch):
    """PSF-match images to reference images

    Fits the following model:
    image to not convolve = (image to convolve convolved with PSF matching kernel) + background model
    """
    ConfigClass = ImagePsfMatchConfig

    def __init__(self, *args, **kwargs):
        """Create a PsfMatchToImage
        
        @param config: see lsst.ip.diffim.PsfMatchConfig
        @param logName: name by which messages are logged
        """
        PsfMatch.__init__(self, *args, **kwargs)
        self.kconfig = self.config.kernel.active
        self._warper = afwMath.Warper.fromConfig(self.kconfig.warpingConfig)

    @pipeBase.timeMethod
    def run(self, imageToConvolve, imageToNotConvolve, mode, **kwargs):
        self.log.warn("run is deprecated; call the appropriate method directly")

        if mode == "matchExposures":
            return self.matchExposures(imageToConvolve, imageToNotConvolve, **kwargs)

        elif mode == "matchMaskedImages":
            return self.matchMaskedImages(imageToConvolve, imageToNotConvolve, **kwargs)

        elif mode == "subtractExposures":
            return self.subtractExposures(imageToConvolve, imageToNotConvolve, **kwargs)

        elif mode == "subtractMaskedImages":
            return self.subtractMaskedImages(imageToConvolve, imageToNotConvolve, **kwargs)
        else:
            raise ValueError("Invalid mode requested")
        
    @pipeBase.timeMethod
    def matchExposures(self, exposureToConvolve, exposureToNotConvolve,
                       psfFwhmPixTc = None, psfFwhmPixTnc = None,
                       candidateList = None, doWarping = True, swapImageToConvolve = False):
        """Warp and PSF-match an exposure to the reference

        Do the following, in order:
        - Warp exposureToConvolve to match exposureToNotConvolve,
            if doWarping True and their WCSs do not already match
        - Determine a PSF matching kernel and differential background model
            that matches exposureToConvolve to exposureToNotConvolve
        - Convolve exposureToConvolve by PSF matching kernel
        
        @param exposureToConvolve: Exposure to warp and PSF-match to the reference masked image
        @param exposureToNotConvolve: Exposure whose WCS and PSF are to be matched
        @param psfFwhmPixTc: FWHM (in pixels) of the Psf in the template image (image to convolve)
        @param psfFwhmPixTnc: FWHM (in pixels) of the Psf in the science image
        @param candidateList: a list of footprints/maskedImages for kernel candidates; if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF
        @param doWarping: what to do if exposureToConvolve's and exposureToNotConvolve's WCSs do not match:
            - if True then warp exposureToConvolve to match exposureToNotConvolve
            - if False then raise an Exception
        @param swapImageToConvolve: switch which image is used as the refernce, in case of e.g. deconvolution
            - if True, exposureToConvolve is warped if doWarping, exposureToNotConvolve is convolved
            - if False, exposureToConvolve is warped if doWarping, exposureToConvolve is convolved
        
        @return a pipeBase.Struct containing these fields:
        - matchedImage: the PSF-matched exposure =
            warped exposureToConvolve convolved by psfMatchingKernel. This has:
            - the same parent bbox, Wcs and Calib as exposureToNotConvolve
            - the same filter as exposureToConvolve
            - no Psf (because the PSF-matching process does not compute one)
        - psfMatchingKernel: the PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to solve for the PSF matching kernel

        @raise RuntimeError if doWarping is False and exposureToConvolve's and exposureToNotConvolve's
            WCSs do not match
        """
        if not self._validateWcs(exposureToConvolve, exposureToNotConvolve):
            if doWarping:
                pexLog.Trace(self.log.getName(), 1, "Astrometrically registering template to science image")
                exposureToConvolve = self._warper.warpExposure(exposureToNotConvolve.getWcs(), 
                    exposureToConvolve, destBBox = exposureToNotConvolve.getBBox(afwImage.PARENT))
            else:
                pexLog.Trace(self.log.getName(), 1, "ERROR: Input images not registered")
                raise RuntimeError, "Input images not registered"

        if psfFwhmPixTc is None:
            if not exposureToConvolve.hasPsf():
                pexLog.Trace(self.log.getName(), 1, "WARNING: no estimate of Psf FWHM for template image")
            else:
                psf = exposureToConvolve.getPsf()
                width, height = psf.getKernel().getDimensions()
                psfAttr = measAlg.PsfAttributes(psf, width//2, height//2)
                psfSigPixTc = psfAttr.computeGaussianWidth(psfAttr.ADAPTIVE_MOMENT)
                psfFwhmPixTc = psfSigPixTc * sigma2fwhm 
        if psfFwhmPixTnc is None:
            if not exposureToNotConvolve.hasPsf():
                pexLog.Trace(self.log.getName(), 1, "WARNING: no estimate of Psf FWHM for science image")
            else:
                psf = exposureToNotConvolve.getPsf()
                width, height = psf.getKernel().getDimensions()
                psfAttr = measAlg.PsfAttributes(psf, width//2, height//2)
                psfSigPixTnc = psfAttr.computeGaussianWidth(psfAttr.ADAPTIVE_MOMENT)
                psfFwhmPixTnc = psfSigPixTnc * sigma2fwhm 

        if candidateList != None:
            if type(candidateList[0]) == afwTable.SourceRecord:
                candidateList = diffimTools.sourceToFootprintList(candidateList, exposureToConvolve, exposureToNotConvolve, 
                                                                  self.kconfig.detectionConfig, self.log)
        if swapImageToConvolve:
            results = self.matchMaskedImages(
                exposureToNotConvolve.getMaskedImage(), exposureToConvolve.getMaskedImage(),
                psfFwhmPixTc = psfFwhmPixTnc, psfFwhmPixTnc = psfFwhmPixTc,
                candidateList = candidateList)
        else:
            results = self.matchMaskedImages(
                exposureToConvolve.getMaskedImage(), exposureToNotConvolve.getMaskedImage(),
                psfFwhmPixTc = psfFwhmPixTc, psfFwhmPixTnc = psfFwhmPixTnc,
                candidateList = candidateList)
        
        psfMatchedExposure = afwImage.makeExposure(results.matchedImage, exposureToNotConvolve.getWcs())
        psfMatchedExposure.setFilter(exposureToConvolve.getFilter())
        psfMatchedExposure.setCalib(exposureToNotConvolve.getCalib())
        results.warpedExposure  = exposureToConvolve
        results.matchedExposure = psfMatchedExposure
        return results

    @pipeBase.timeMethod
    def matchMaskedImages(self, maskedImageToConvolve, maskedImageToNotConvolve, 
                          psfFwhmPixTc = None, psfFwhmPixTnc = None, 
                          candidateList = None):
        """PSF-match a MaskedImage to a reference MaskedImage

        Do the following, in order:
        - Determine a PSF matching kernel and differential background model
            that matches maskedImageToConvolve to maskedImageToNotConvolve
        - Convolve maskedImageToConvolve by the PSF matching kernel
        
        @param maskedImageToConvolve: masked image to PSF-match to the reference masked image;
            must be warped to match the reference masked image
        @param maskedImageToNotConvolve: maskedImage whose PSF is to be matched
        @param psfFwhmPixTc: FWHM (in pixels) of the Psf in the template image (image to convolve)
        @param psfFwhmPixTnc: FWHM (in pixels) of the Psf in the science image
        @param candidateList: a list of footprints/maskedImages for kernel candidates; if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF
        
        @return a pipeBase.Struct containing these fields:
        - psfMatchedMaskedImage: the PSF-matched masked image =
            maskedImageToConvolve convolved with psfMatchingKernel.
            This has the same xy0, dimensions and wcs as maskedImageToNotConvolve.
        - psfMatchingKernel: the PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to solve for the PSF matching kernel
        
        @raise RuntimeError if input images have different dimensions
        """


        if not self._validateSize(maskedImageToConvolve, maskedImageToNotConvolve):
            pexLog.Trace(self.log.getName(), 1, "ERROR: Input images different size")
            raise RuntimeError, "Input images different size"
            
        self.log.log(pexLog.Log.INFO, "compute PSF-matching kernel")
        if psfFwhmPixTc and psfFwhmPixTnc:
            pexLog.Trace(self.log.getName(), 2, "Matching Psf FWHM %.2f -> %.2f pix" % (psfFwhmPixTc, psfFwhmPixTnc))




        if self.kconfig.useBicForKernelBasis:
            tmpKernelCellSet = self._buildCellSet(maskedImageToConvolve,
                                                  maskedImageToNotConvolve,
                                                  candidateList = candidateList)
            nbe = diffimTools.NbasisEvaluator(self.kconfig, psfFwhmPixTc, psfFwhmPixTnc)
            bicDegrees = nbe(tmpKernelCellSet, self.log)
            basisList = makeKernelBasisList(self.kconfig, psfFwhmPixTc, psfFwhmPixTnc, bicDegrees[0])
            del tmpKernelCellSet
        else:
            basisList = makeKernelBasisList(self.kconfig, psfFwhmPixTc, psfFwhmPixTnc)


        kernelCellSet = self._buildCellSet(maskedImageToConvolve,
                                           maskedImageToNotConvolve,
                                           candidateList = candidateList)

        spatialSolution, psfMatchingKernel, backgroundModel = self._solve(kernelCellSet, basisList)
        conditionNum = spatialSolution.getConditionNumber(eval("diffimLib.KernelSolution.%s" % (self.kconfig.conditionNumberType)))        
        self.metadata.set("spatialConditionNum", conditionNum)

        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayTemplate = lsstDebug.Info(__name__).displayTemplate
        displaySciIm = lsstDebug.Info(__name__).displaySciIm
        displaySpatialCells = lsstDebug.Info(__name__).displaySpatialCells
        maskTransparency = lsstDebug.Info(__name__).maskTransparency   
        if not maskTransparency:
            maskTransparency = 0
        ds9.setMaskTransparency(maskTransparency)

        if display and displayTemplate:
            ds9.mtv(maskedImageToConvolve, frame=lsstDebug.frame, title="Image to convolve")
            lsstDebug.frame += 1

        if display and displaySpatialCells:
            diUtils.showKernelSpatialCells(maskedImageToNotConvolve, kernelCellSet, 
                                           symb="o", ctype=ds9.CYAN, ctypeUnused=ds9.YELLOW, ctypeBad=ds9.RED,
                                           size=4, frame=lsstDebug.frame)
            lsstDebug.frame += 1
        elif display and  displaySciIm:
            ds9.mtv(maskedImageToNotConvolve, frame=lsstDebug.frame, title="Image to not convolve")
            lsstDebug.frame += 1

        psfMatchedMaskedImage = afwImage.MaskedImageF(maskedImageToConvolve.getBBox(afwImage.PARENT))
        doNormalize = False
        afwMath.convolve(psfMatchedMaskedImage, maskedImageToConvolve, psfMatchingKernel, doNormalize)
        self.log.log(pexLog.Log.INFO, "done")
        return pipeBase.Struct(
            matchedImage = psfMatchedMaskedImage,
            psfMatchingKernel = psfMatchingKernel,
            backgroundModel = backgroundModel,
            kernelCellSet = kernelCellSet,
        )

    @pipeBase.timeMethod
    def subtractExposures(self, exposureToConvolve, exposureToNotConvolve,
                          psfFwhmPixTc = None, psfFwhmPixTnc = None,
                          candidateList = None, doWarping = True, swapImageToConvolve = False):
        """Subtract two Exposures
        
        Do the following, in order:
        - Warp exposureToConvolve to match exposureToNotConvolve, if their WCSs do not already match
        - Determine a PSF matching kernel and differential background model
            that matches exposureToConvolve to exposureToNotConvolve
        - PSF-match exposureToConvolve to exposureToNotConvolve
        - Compute subtracted exposure (see return values for equation).

        @param exposureToConvolve: exposure to PSF-matched to exposureToNotConvolve
        @param exposureToNotConvolve: reference Exposure
        @param psfFwhmPixTc: FWHM (in pixels) of the Psf in the template image (image to convolve)
        @param psfFwhmPixTnc: FWHM (in pixels) of the Psf in the science image
        @param candidateList: a list of footprints/maskedImages for kernel candidates; if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF
        @param doWarping: what to do if exposureToConvolve's and exposureToNotConvolve's WCSs do not match:
            - if True then warp exposureToConvolve to match exposureToNotConvolve
            - if False then raise an Exception
        @param swapImageToConvolve: switch which image is used as the refernce, in case of e.g. deconvolution
            - if True, exposureToConvolve is warped if doWarping, exposureToNotConvolve is convolved
            - if False, exposureToConvolve is warped if doWarping, exposureToConvolve is convolved
        
        @return a pipeBase.Struct containing these fields:
        - subtractedExposure: subtracted Exposure = exposureToNotConvolve - (matchedImage + backgroundModel)
        - matchedImage: exposureToConvolve after warping to match exposureToConvolve (if doWarping true),
            and convolving with psfMatchingKernel
        - psfMatchingKernel: PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to determine PSF matching kernel
        """     
        results = self.matchExposures(
            exposureToConvolve = exposureToConvolve,
            exposureToNotConvolve = exposureToNotConvolve,
            psfFwhmPixTc = psfFwhmPixTc,
            psfFwhmPixTnc = psfFwhmPixTnc,
            candidateList = candidateList,
            doWarping = doWarping,
            swapImageToConvolve = swapImageToConvolve
        )
        
        subtractedExposure = afwImage.ExposureF(exposureToNotConvolve, True)
        if swapImageToConvolve:
            subtractedExposure.setMaskedImage(results.warpedExposure.getMaskedImage())
            subtractedMaskedImage  = subtractedExposure.getMaskedImage()
            subtractedMaskedImage -= results.matchedExposure.getMaskedImage()
            subtractedMaskedImage -= results.backgroundModel

            # Preserve polaity of differences
            subtractedMaskedImage *= -1

            # Place back on native photometric scale
            subtractedMaskedImage /= results.psfMatchingKernel.computeImage(afwImage.ImageD(results.psfMatchingKernel.getDimensions()), False)

        else:
            subtractedMaskedImage  = subtractedExposure.getMaskedImage()
            subtractedMaskedImage -= results.matchedExposure.getMaskedImage()
            subtractedMaskedImage -= results.backgroundModel

        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayDiffIm = lsstDebug.Info(__name__).displayDiffIm
        maskTransparency = lsstDebug.Info(__name__).maskTransparency   
        if not maskTransparency:
            maskTransparency = 0
        ds9.setMaskTransparency(maskTransparency)
        if display and displayDiffIm:
            ds9.mtv(subtractedExposure, frame=lsstDebug.frame)
            lsstDebug.frame += 1

        results.subtractedExposure = subtractedExposure
        return results

    @pipeBase.timeMethod
    def subtractMaskedImages(self, maskedImageToConvolve, maskedImageToNotConvolve,
                             psfFwhmPixTc = None, psfFwhmPixTnc = None,
                             candidateList = None):
        """Subtract two MaskedImages
        
        Do the following, in order:
        - PSF-match maskedImageToConvolve to maskedImageToNotConvolve
        - Determine the differential background
        - Return the difference: maskedImageToNotConvolve -
            ((warped maskedImageToConvolve convolved with psfMatchingKernel) + backgroundModel)
        
        @param maskedImageToConvolve: MaskedImage to PSF-matched to maskedImageToNotConvolve
        @param maskedImageToNotConvolve: reference MaskedImage
        @param psfFwhmPixTc: FWHM (in pixels) of the Psf in the template image (image to convolve)
        @param psfFwhmPixTnc: FWHM (in pixels) of the Psf in the science image
        @param candidateList: a list of footprints/maskedImages for kernel candidates; if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF
        
        @return a pipeBase.Struct containing these fields:
        - subtractedMaskedImage = maskedImageToNotConvolve - (matchedImage + backgroundModel)
        - matchedImage: maskedImageToConvolve convolved with psfMatchingKernel
        - psfMatchingKernel: PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to determine PSF matching kernel
        """
        results = self.matchMaskedImages(
            maskedImageToConvolve = maskedImageToNotConvolve,
            maskedImageToNotConvolve = maskedImageToConvolve,
            psfFwhmPixTc = psfFwhmPixTnc,
            psfFwhmPixTnc = psfFwhmPixTc,
            candidateList = candidateList,
            )

        subtractedMaskedImage  = afwImage.MaskedImageF(maskedImageToNotConvolve, True)
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

    def _adaptCellSize(self, candidateList):
        """ NOT IMPLEMENTED YET"""
        nCand = len(candidateList)
        return self.kconfig.sizeCellX, self.kconfig.sizeCellY

    def _buildCellSet(self, maskedImageToConvolve, maskedImageToNotConvolve, candidateList = None):
        """Build a SpatialCellSet for use with the solve method

        @param maskedImageToConvolve: MaskedImage to PSF-matched to maskedImageToNotConvolve
        @param maskedImageToNotConvolve: reference MaskedImage
        @param candidateList: a list of footprints/maskedImages for kernel candidates; if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF
        
        @return kernelCellSet: a SpatialCellSet for use with self._solve
        """
        # Candidate source footprints to use for Psf matching
        if candidateList == None:
            self.log.log(pexLog.Log.INFO, "temporarily subtracting backgrounds for detection")
            mi1 = maskedImageToConvolve.Factory(maskedImageToConvolve, True)
            mi2 = maskedImageToNotConvolve.Factory(maskedImageToNotConvolve, True)
            tmp = diffimTools.backgroundSubtract(self.kconfig.afwBackgroundConfig, [mi1, mi2])

            detConfig = self.kconfig.detectionConfig
            kcDetect = diffimLib.KernelCandidateDetectionF(pexConfig.makePolicy(detConfig))
            kcDetect.apply(mi1, mi2)
            candidateList = kcDetect.getFootprints()

        sizeCellX, sizeCellY = self._adaptCellSize(candidateList)

        # Object to store the KernelCandidates for spatial modeling
        kernelCellSet = afwMath.SpatialCellSet(maskedImageToConvolve.getBBox(afwImage.PARENT),
                                               sizeCellX, sizeCellY)
        
            
        policy = pexConfig.makePolicy(self.kconfig)
        # Place candidates within the spatial grid
        for cand in candidateList:
            bbox = cand.getBBox()
            
            # Grab the centers in the parent's coordinate system
            xC   = 0.5 * ( bbox.getMinX() + bbox.getMaxX() )
            yC   = 0.5 * ( bbox.getMinY() + bbox.getMaxY() )
            
            tmi  = afwImage.MaskedImageF(maskedImageToConvolve, bbox, afwImage.PARENT)
            smi  = afwImage.MaskedImageF(maskedImageToNotConvolve, bbox, afwImage.PARENT)
            cand = diffimLib.makeKernelCandidate(xC, yC, tmi, smi, policy)
            
            pexLog.Trace(self.log.getName(), 5,
                         "Candidate %d at %f, %f" % (cand.getId(), cand.getXCenter(), cand.getYCenter()))
            kernelCellSet.insertCandidate(cand)

        return kernelCellSet

    def _validateSize(self, maskedImageToConvolve, maskedImageToNotConvolve):
        """Return True if two image-like objects are the same size
        """
        return maskedImageToConvolve.getDimensions() == maskedImageToNotConvolve.getDimensions()
    
    def _validateWcs(self, exposureToConvolve, exposureToNotConvolve):
        """Return True if two Exposures have the same WCS
        """
        templateWcs    = exposureToConvolve.getWcs() 
        scienceWcs     = exposureToNotConvolve.getWcs()
        templateBBox   = exposureToConvolve.getBBox(afwImage.PARENT)
        scienceBBox    = exposureToNotConvolve.getBBox(afwImage.PARENT)

        # LLC
        templateOrigin = templateWcs.pixelToSky(afwGeom.Point2D(templateBBox.getBegin()))
        scienceOrigin  = scienceWcs.pixelToSky(afwGeom.Point2D(scienceBBox.getBegin()))

        # URC
        templateLimit  = templateWcs.pixelToSky(afwGeom.Point2D(templateBBox.getEnd()))
        scienceLimit   = scienceWcs.pixelToSky(afwGeom.Point2D(scienceBBox.getEnd()))
        
        print ("Template limits : %f,%f -> %f,%f" %
                     (templateOrigin[0], templateOrigin[1],
                      templateLimit[0], templateLimit[1]))
        print ("Science limits : %f,%f -> %f,%f" %
                     (scienceOrigin[0], scienceOrigin[1],
                      scienceLimit[0], scienceLimit[1]))

        templateBBox = afwGeom.Box2D(templateOrigin.getPosition(), templateLimit.getPosition())
        scienceBBox  = afwGeom.Box2D(scienceOrigin.getPosition(), scienceLimit.getPosition())
        if not (templateBBox.overlaps(scienceBBox)):
            raise RuntimeError, "Input images do not overlap at all"
            

        if ( (templateOrigin.getPosition() != scienceOrigin.getPosition()) or \
             (templateLimit.getPosition()  != scienceLimit.getPosition())  or \
             (exposureToConvolve.getDimensions() != exposureToNotConvolve.getDimensions())):
            return False
        return True
