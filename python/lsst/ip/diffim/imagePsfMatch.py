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
import diffimLib 
import lsst.pex.logging as pexLog
import lsst.pex.config as pexConfig
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.pipe.base as pipeBase
import diffimTools
from .makeKernelBasisList import makeKernelBasisList
from .psfMatch import PsfMatch, PsfMatchConfigDF, PsfMatchConfigAL

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
                       footprints = None, doWarping = True):
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
        @param footprints: a list of footprints of sources; if None then source detection is run
        @param doWarping: what to do if exposureToConvolve's and exposureToNotConvolve's WCSs do not match:
            - if True then warp exposureToConvolve to match exposureToNotConvolve
            - if False then raise an Exception
        
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


        results = self.matchMaskedImages(
            exposureToConvolve.getMaskedImage(), exposureToNotConvolve.getMaskedImage(),
            psfFwhmPixTc = psfFwhmPixTc, psfFwhmPixTnc = psfFwhmPixTnc,
            footprints = footprints)
        
        psfMatchedExposure = afwImage.makeExposure(results.matchedImage, exposureToNotConvolve.getWcs())
        psfMatchedExposure.setFilter(exposureToConvolve.getFilter())
        psfMatchedExposure.setCalib(exposureToNotConvolve.getCalib())
        results.matchedImage = psfMatchedExposure
        return results

    @pipeBase.timeMethod
    def matchMaskedImages(self, maskedImageToConvolve, maskedImageToNotConvolve, 
                          psfFwhmPixTc = None, psfFwhmPixTnc = None, 
                          footprints = None):
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
        @param footprints: a list of footprints of sources; if None then source detection is run
        
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
        kernelCellSet = self._buildCellSet(maskedImageToConvolve,
                                           maskedImageToNotConvolve,
                                           footprints = footprints)
        basisList = makeKernelBasisList(self.kconfig, psfFwhmPixTc, psfFwhmPixTnc)
        spatialSolution, psfMatchingKernel, backgroundModel = self._solve(kernelCellSet, basisList)
        conditionNum = spatialSolution.getConditionNumber(eval("diffimLib.KernelSolution.%s" % (self.kconfig.conditionNumberType)))        
        self.metadata.set("spatialConditionNum", conditionNum)
    
        self.log.log(pexLog.Log.INFO, "PSF-match science MaskedImage to reference")
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
                          footprints = None, doWarping = True):
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
        @param footprints: a list of footprints of sources; if None then source detection is run
        @param doWarping: what to do if exposureToConvolve's and exposureToNotConvolve's WCSs do not match:
            - if True then warp exposureToConvolve to match exposureToNotConvolve
            - if False then raise an Exception
        
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
            footprints = footprints,
            doWarping = doWarping,
        )

        subtractedExposure  = afwImage.ExposureF(exposureToNotConvolve, True)
        smi  = subtractedExposure.getMaskedImage()
        smi -= results.matchedImage.getMaskedImage()
        smi -= results.backgroundModel
        results.subtractedImage = subtractedExposure
        return results

    @pipeBase.timeMethod
    def subtractMaskedImages(self, maskedImageToConvolve, maskedImageToNotConvolve,
                             psfFwhmPixTc = None, psfFwhmPixTnc = None,
                             footprints = None):
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
        @param footprints: a list of footprints of sources; if None then source detection is run
        
        @return a pipeBase.Struct containing these fields:
        - subtractedMaskedImage = maskedImageToNotConvolve - (matchedImage + backgroundModel)
        - matchedImage: maskedImageToConvolve convolved with psfMatchingKernel
        - psfMatchingKernel: PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to determine PSF matching kernel
        """
        results = self.matchMaskedImages(
            maskedImageToConvolve = maskedImageToConvolve,
            maskedImageToNotConvolve = maskedImageToNotConvolve,
            psfFwhmPixTc = psfFwhmPixTc,
            psfFwhmPixTnc = psfFwhmPixTnc,
            footprints = footprints,
        )
        
        subtractedMaskedImage  = afwImage.MaskedImageF(maskedImageToNotConvolve, True)
        subtractedMaskedImage -= results.matchedImage
        subtractedMaskedImage -= results.backgroundModel
        results.subtractedMaskedImage = subtractedMaskedImage
        return results

    def _adaptCellSize(self, footprints):
        """ NOT IMPLEMENTED YET"""
        nFp = len(footprints)
        return self.kconfig.sizeCellX, self.kconfig.sizeCellY

    def _buildCellSet(self, maskedImageToConvolve, maskedImageToNotConvolve, footprints = None):
        """Build a SpatialCellSet for use with the solve method

        @param maskedImageToConvolve: MaskedImage to PSF-matched to maskedImageToNotConvolve
        @param maskedImageToNotConvolve: reference MaskedImage
        @param footprints: a list of footprints of sources; if None then source detection is run
        
        @return kernelCellSet: a SpatialCellSet for use with self._solve
        """
        
        # Candidate source footprints to use for Psf matching
        if footprints == None:
            self.log.log(pexLog.Log.INFO, "temporarily subtracting backgrounds for detection")
            mi1 = maskedImageToConvolve.Factory(maskedImageToConvolve, True)
            mi2 = maskedImageToNotConvolve.Factory(maskedImageToNotConvolve, True)
            tmp = diffimTools.backgroundSubtract(self.kconfig.afwBackgroundConfig, [mi1, mi2])
            
            detConfig = self.kconfig.detectionConfig
            kcDetect = diffimLib.KernelCandidateDetectionF(pexConfig.makePolicy(detConfig))
            kcDetect.apply(mi1, mi2)
            footprints = kcDetect.getFootprints()

        sizeCellX, sizeCellY = self._adaptCellSize(footprints)

        # Object to store the KernelCandidates for spatial modeling
        kernelCellSet = afwMath.SpatialCellSet(maskedImageToConvolve.getBBox(afwImage.PARENT),
                                               sizeCellX, sizeCellY)
        
            
        policy = pexConfig.makePolicy(self.kconfig)
        # Place candidate footprints within the spatial grid
        for fp in footprints:
            bbox = fp.getBBox()
            
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
