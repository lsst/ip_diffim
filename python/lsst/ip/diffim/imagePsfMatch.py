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
    def run(self, templateImage, scienceImage, mode, **kwargs):
        self.log.warn("run is deprecated; call the appropriate method directly")

        if mode == "matchExposures":
            return self.matchExposures(templateImage, scienceImage, **kwargs)

        elif mode == "matchMaskedImages":
            return self.matchMaskedImages(templateImage, scienceImage, **kwargs)

        elif mode == "subtractExposures":
            return self.subtractExposures(templateImage, scienceImage, **kwargs)

        elif mode == "subtractMaskedImages":
            return self.subtractMaskedImages(templateImage, scienceImage, **kwargs)
        else:
            raise ValueError("Invalid mode requested")
        
    @pipeBase.timeMethod
    def matchExposures(self, templateExposure, scienceExposure,
                       templateFwhmPix = None, scienceFwhmPix = None,
                       candidateList = None, doWarping = True, convolveTemplate = True):
        """Warp and PSF-match an exposure to the reference

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
        @param candidateList: a list of footprints/maskedImages for kernel candidates; if None then source detection is run.
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

        @raise RuntimeError if doWarping is False and templateExposure's and scienceExposure's
            WCSs do not match
        """
        if not self._validateWcs(templateExposure, scienceExposure):
            if doWarping:
                pexLog.Trace(self.log.getName(), 1, "Astrometrically registering template to science image")
                templateExposure = self._warper.warpExposure(scienceExposure.getWcs(), 
                    templateExposure, destBBox = scienceExposure.getBBox(afwImage.PARENT))
            else:
                pexLog.Trace(self.log.getName(), 1, "ERROR: Input images not registered")
                raise RuntimeError, "Input images not registered"

        if templateFwhmPix is None:
            if not templateExposure.hasPsf():
                pexLog.Trace(self.log.getName(), 1, "WARNING: no estimate of Psf FWHM for template image")
            else:
                psf = templateExposure.getPsf()
                width, height = psf.getKernel().getDimensions()
                psfAttr = measAlg.PsfAttributes(psf, width//2, height//2)
                templateSigPix = psfAttr.computeGaussianWidth(psfAttr.ADAPTIVE_MOMENT)
                templateFwhmPix = templateSigPix * sigma2fwhm 
        if scienceFwhmPix is None:
            if not scienceExposure.hasPsf():
                pexLog.Trace(self.log.getName(), 1, "WARNING: no estimate of Psf FWHM for science image")
            else:
                psf = scienceExposure.getPsf()
                width, height = psf.getKernel().getDimensions()
                psfAttr = measAlg.PsfAttributes(psf, width//2, height//2)
                scienceSigPix = psfAttr.computeGaussianWidth(psfAttr.ADAPTIVE_MOMENT)
                scienceFwhmPix = scienceSigPix * sigma2fwhm 

        if candidateList != None:
            if type(candidateList[0]) == afwTable.SourceRecord:
                candidateList = diffimTools.sourceToFootprintList(candidateList, templateExposure, scienceExposure, 
                                                                  self.kconfig.detectionConfig, self.log)
        if convolveTemplate:
            results = self.matchMaskedImages(
                templateExposure.getMaskedImage(), scienceExposure.getMaskedImage(),
                templateFwhmPix = templateFwhmPix, scienceFwhmPix = scienceFwhmPix,
                candidateList = candidateList)
        else:
            results = self.matchMaskedImages(
                scienceExposure.getMaskedImage(), templateExposure.getMaskedImage(),
                templateFwhmPix = scienceFwhmPix, scienceFwhmPix = templateFwhmPix,
                candidateList = candidateList)
        
        psfMatchedExposure = afwImage.makeExposure(results.matchedImage, scienceExposure.getWcs())
        psfMatchedExposure.setFilter(templateExposure.getFilter())
        psfMatchedExposure.setCalib(scienceExposure.getCalib())
        results.warpedExposure  = templateExposure
        results.matchedExposure = psfMatchedExposure
        return results

    @pipeBase.timeMethod
    def matchMaskedImages(self, templateMaskedImage, scienceMaskedImage, 
                          templateFwhmPix = None, scienceFwhmPix = None, 
                          candidateList = None):
        """PSF-match a MaskedImage (templateMaskedImage) to a reference MaskedImage (scienceMaskedImage)

        Do the following, in order:
        - Determine a PSF matching kernel and differential background model
            that matches templateMaskedImage to scienceMaskedImage
        - Convolve templateMaskedImage by the PSF matching kernel
        
        @param templateMaskedImage: masked image to PSF-match to the reference masked image;
            must be warped to match the reference masked image
        @param scienceMaskedImage: maskedImage whose PSF is to be matched to
        @param templateFwhmPix: FWHM (in pixels) of the Psf in the template image (image to convolve)
        @param scienceFwhmPix: FWHM (in pixels) of the Psf in the science image
        @param candidateList: a list of footprints/maskedImages for kernel candidates; if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF
        
        @return a pipeBase.Struct containing these fields:
        - psfMatchedMaskedImage: the PSF-matched masked image =
            templateMaskedImage convolved with psfMatchingKernel.
            This has the same xy0, dimensions and wcs as scienceMaskedImage.
        - psfMatchingKernel: the PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to solve for the PSF matching kernel
        
        @raise RuntimeError if input images have different dimensions
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


        if not self._validateSize(templateMaskedImage, scienceMaskedImage):
            pexLog.Trace(self.log.getName(), 1, "ERROR: Input images different size")
            raise RuntimeError, "Input images different size"

        if display and displayTemplate:
            ds9.mtv(templateMaskedImage, frame=lsstDebug.frame, title="Image to convolve")
            lsstDebug.frame += 1

        if display and  displaySciIm:
            ds9.mtv(scienceMaskedImage, frame=lsstDebug.frame, title="Image to not convolve")
            lsstDebug.frame += 1

        kernelCellSet = self._buildCellSet(templateMaskedImage,
                                           scienceMaskedImage,
                                           candidateList = candidateList)

        if display and displaySpatialCells:
            diUtils.showKernelSpatialCells(scienceMaskedImage, kernelCellSet, 
                                           symb="o", ctype=ds9.CYAN, ctypeUnused=ds9.YELLOW, ctypeBad=ds9.RED,
                                           size=4, frame=lsstDebug.frame, title="Image to not convolve")
            lsstDebug.frame += 1

        if templateFwhmPix and scienceFwhmPix:
            pexLog.Trace(self.log.getName(), 2, "Matching Psf FWHM %.2f -> %.2f pix" % (templateFwhmPix, scienceFwhmPix))

        if self.kconfig.useBicForKernelBasis:
            tmpKernelCellSet = self._buildCellSet(templateMaskedImage,
                                                  scienceMaskedImage,
                                                  candidateList = candidateList)
            nbe = diffimTools.NbasisEvaluator(self.kconfig, templateFwhmPix, scienceFwhmPix)
            bicDegrees = nbe(tmpKernelCellSet, self.log)
            basisList = makeKernelBasisList(self.kconfig, templateFwhmPix, scienceFwhmPix, bicDegrees[0])
            del tmpKernelCellSet
        else:
            basisList = makeKernelBasisList(self.kconfig, templateFwhmPix, scienceFwhmPix)

        spatialSolution, psfMatchingKernel, backgroundModel = self._solve(kernelCellSet, basisList)




        psfMatchedMaskedImage = afwImage.MaskedImageF(templateMaskedImage.getBBox(afwImage.PARENT))
        doNormalize = False
        afwMath.convolve(psfMatchedMaskedImage, templateMaskedImage, psfMatchingKernel, doNormalize)
        self.log.log(pexLog.Log.INFO, "done")
        return pipeBase.Struct(
            matchedImage = psfMatchedMaskedImage,
            psfMatchingKernel = psfMatchingKernel,
            backgroundModel = backgroundModel,
            kernelCellSet = kernelCellSet,
        )

    @pipeBase.timeMethod
    def subtractExposures(self, templateExposure, scienceExposure,
                          templateFwhmPix = None, scienceFwhmPix = None,
                          candidateList = None, doWarping = True, convolveTemplate = True):
        """Subtract two Exposures
        
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
        @param candidateList: a list of footprints/maskedImages for kernel candidates; if None then source detection is run.
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
            templateExposure = templateExposure,
            scienceExposure = scienceExposure,
            templateFwhmPix = templateFwhmPix,
            scienceFwhmPix = scienceFwhmPix,
            candidateList = candidateList,
            doWarping = doWarping,
            convolveTemplate = convolveTemplate
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
            subtractedMaskedImage /= results.psfMatchingKernel.computeImage(afwImage.ImageD(results.psfMatchingKernel.getDimensions()), False)

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
    def subtractMaskedImages(self, templateMaskedImage, scienceMaskedImage,
                             templateFwhmPix = None, scienceFwhmPix = None,
                             candidateList = None):
        """Subtract two MaskedImages
        
        Do the following, in order:
        - PSF-match templateMaskedImage to scienceMaskedImage
        - Determine the differential background
        - Return the difference: scienceMaskedImage -
            ((warped templateMaskedImage convolved with psfMatchingKernel) + backgroundModel)
        
        @param templateMaskedImage: MaskedImage to PSF-match to scienceMaskedImage
        @param scienceMaskedImage: reference MaskedImage
        @param templateFwhmPix: FWHM (in pixels) of the Psf in the template image (image to convolve)
        @param scienceFwhmPix: FWHM (in pixels) of the Psf in the science image
        @param candidateList: a list of footprints/maskedImages for kernel candidates; if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF
        
        @return a pipeBase.Struct containing these fields:
        - subtractedMaskedImage = scienceMaskedImage - (matchedImage + backgroundModel)
        - matchedImage: templateMaskedImage convolved with psfMatchingKernel
        - psfMatchingKernel: PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to determine PSF matching kernel
        """
        results = self.matchMaskedImages(
            templateMaskedImage = templateMaskedImage,
            scienceMaskedImage = scienceMaskedImage,
            templateFwhmPix = templateFwhmPix,
            scienceFwhmPix = scienceFwhmPix,
            candidateList = candidateList,
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

    def _adaptCellSize(self, candidateList):
        """ NOT IMPLEMENTED YET"""
        nCand = len(candidateList)
        return self.kconfig.sizeCellX, self.kconfig.sizeCellY

    def _buildCellSet(self, templateMaskedImage, scienceMaskedImage, candidateList = None):
        """Build a SpatialCellSet for use with the solve method

        @param templateMaskedImage: MaskedImage to PSF-matched to scienceMaskedImage
        @param scienceMaskedImage: reference MaskedImage
        @param candidateList: a list of footprints/maskedImages for kernel candidates; if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF
        
        @return kernelCellSet: a SpatialCellSet for use with self._solve
        """
        # Candidate source footprints to use for Psf matching
        if candidateList == None:
            self.log.log(pexLog.Log.INFO, "temporarily subtracting backgrounds for detection")
            mi1 = templateMaskedImage.Factory(templateMaskedImage, True)
            mi2 = scienceMaskedImage.Factory(scienceMaskedImage, True)
            tmp = diffimTools.backgroundSubtract(self.kconfig.afwBackgroundConfig, [mi1, mi2])

            detConfig = self.kconfig.detectionConfig
            kcDetect = diffimLib.KernelCandidateDetectionF(pexConfig.makePolicy(detConfig))
            kcDetect.apply(mi1, mi2)
            candidateList = kcDetect.getFootprints()

        sizeCellX, sizeCellY = self._adaptCellSize(candidateList)

        # Object to store the KernelCandidates for spatial modeling
        kernelCellSet = afwMath.SpatialCellSet(templateMaskedImage.getBBox(afwImage.PARENT),
                                               sizeCellX, sizeCellY)
        
            
        policy = pexConfig.makePolicy(self.kconfig)
        # Place candidates within the spatial grid
        for cand in candidateList:
            bbox = cand.getBBox()
            
            # Grab the centers in the parent's coordinate system
            xC   = 0.5 * ( bbox.getMinX() + bbox.getMaxX() )
            yC   = 0.5 * ( bbox.getMinY() + bbox.getMaxY() )
            
            tmi  = afwImage.MaskedImageF(templateMaskedImage, bbox, afwImage.PARENT)
            smi  = afwImage.MaskedImageF(scienceMaskedImage, bbox, afwImage.PARENT)
            cand = diffimLib.makeKernelCandidate(xC, yC, tmi, smi, policy)
            
            pexLog.Trace(self.log.getName(), 5,
                         "Candidate %d at %f, %f" % (cand.getId(), cand.getXCenter(), cand.getYCenter()))
            kernelCellSet.insertCandidate(cand)

        return kernelCellSet

    def _validateSize(self, templateMaskedImage, scienceMaskedImage):
        """Return True if two image-like objects are the same size
        """
        return templateMaskedImage.getDimensions() == scienceMaskedImage.getDimensions()
    
    def _validateWcs(self, templateExposure, scienceExposure):
        """Return True if the WCS of the two Exposures have the same origin and extent
        """
        templateWcs    = templateExposure.getWcs() 
        scienceWcs     = scienceExposure.getWcs()
        templateBBox   = templateExposure.getBBox(afwImage.PARENT)
        scienceBBox    = scienceExposure.getBBox(afwImage.PARENT)

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
            raise RuntimeError, "Input images do not overlap at all"

        if ( (templateOrigin.getPosition() != scienceOrigin.getPosition()) or \
             (templateLimit.getPosition()  != scienceLimit.getPosition())  or \
             (templateExposure.getDimensions() != scienceExposure.getDimensions())):
            return False
        return True
