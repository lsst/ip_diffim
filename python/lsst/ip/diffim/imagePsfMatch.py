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

class ImagePsfMatchConfig(PsfMatchConfig):
    """The default config is basically designed for Image Psf matching"""
    def __init__(self):
        PsfMatchConfig.__init__(self)

class SnapPsfMatchConfigDF(PsfMatchConfigDF):
    """Version of Psf Matching optimized for snap subtraction"""
    def __init__(self):
        PsfMatchConfigDF.__init__(self)
        
        # No spatial variation in model
        self.spatialKernelOrder = 0
        
        # Don't fit for differential background
        self.fitForBackground = False

        # Small kernel size
        self.kernelSize = 7

        # With zero spatial order don't worry about spatial clipping
        self.spatialKernelClipping = False

        # No regularization
        self.useRegularization = False
            
        
class SnapPsfMatchConfigAL(PsfMatchConfigAL):
    """Version of Psf Matching optimized for snap subtraction"""
    def __init__(self):
        PsfMatchConfigAL.__init__(self)
        
        # No spatial variation in model
        self.spatialKernelOrder = 0
        
        # Don't fit for differential background
        self.fitForBackground = False

        # Small kernel size
        self.kernelSize = 7

        # With zero spatial order don't worry about spatial clipping
        self.spatialKernelClipping = False

        # Simple basis set
        self.kernelBasisSet = "alard-lupton"
        self.alardNGauss = 2
        self.alardDegGauss = (4, 2)
        self.alardSigGauss = (1.0, 2.5)
            
        

class ImagePsfMatch(PsfMatch):
    """PSF-match images to reference images

    Fits the following model:
    image to not convolve = (image to convolve convolved with PSF matching kernel) + background model
    """
    def __init__(self, policy, logName="lsst.ip.diffim.ImagePsfMatch"):
        """Create a PsfMatchToImage
        
        @param policy: see lsst/ip/diffim/policy/PsfMatchingDictionary.paf
        @param logName: name by which messages are logged
        """
        PsfMatch.__init__(self, policy, logName)
        self._warper = afwMath.Warper.fromPolicy(policy.getPolicy("warpingPolicy"))

    def _validateSize(self, maskedImageToConvolve, maskedImageToNotConvolve):
        """Return True if two image-like objects are the same size
        """
        return maskedImageToConvolve.getDimensions() == maskedImageToNotConvolve.getDimensions()
    
    def _validateWcs(self, exposureToConvolve, exposureToNotConvolve):
        """Return True if two Exposures have the same WCS
        """
        templateWcs    = exposureToConvolve.getWcs() 
        scienceWcs     = exposureToNotConvolve.getWcs()
        
        # LLC
        templateOrigin = templateWcs.pixelToSky(0, 0)
        scienceOrigin  = scienceWcs.pixelToSky(0, 0)
        # URC
        templateLimit  = templateWcs.pixelToSky(exposureToConvolve.getWidth(),
                                                exposureToConvolve.getHeight())
        scienceLimit   = scienceWcs.pixelToSky(exposureToNotConvolve.getWidth(),
                                               exposureToNotConvolve.getHeight())
        
        pexLog.Trace(self._log.getName(), 2,
                     "Template limits : %f,%f -> %f,%f" %
                     (templateOrigin[0], templateOrigin[1],
                      templateLimit[0], templateLimit[1]))
        pexLog.Trace(self._log.getName(), 2,
                     "Science limits : %f,%f -> %f,%f" %
                     (scienceOrigin[0], scienceOrigin[1],
                      scienceLimit[0], scienceLimit[1]))

        if ( (templateOrigin.getPosition() != scienceOrigin.getPosition()) or \
             (templateLimit.getPosition()  != scienceLimit.getPosition())  or \
             (exposureToConvolve.getDimensions() != exposureToNotConvolve.getDimensions())):
            return False
        return True
        
    def matchExposures(self, exposureToConvolve, exposureToNotConvolve,
                       footprints = None, doWarping = True):
        """Warp and PSF-match an exposure to the reference

        Do the following, in order:
        - Warp exposureToConvolve to match exposureToNotConvolve,
            if doWarping True and their WCSs do not already match
        - Determine a PSF matching kernel and differential background model
            that matches exposureToConvolve to exposureToNotConvolve
        - Convolve exposureToConvolve by PSF matching kernel
        
        @param exposure: Exposure to warp and PSF-match to the reference masked image
        @param referenceMaskedImage: maskedImage whose PSF is to be matched
        @param footprints: a list of footprints of sources; if None then source detection is run
        @param doWarping: what to do if exposureToConvolve's and exposureToNotConvolve's WCSs do not match:
            - if True then warp exposureToConvolve to match exposureToNotConvolve
            - if False then raise an Exception
        
        @return
        - psfMatchedExposure: the PSF-matched exposure =
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
                pexLog.Trace(self._log.getName(), 1, "Astrometrically registering template to science image")
                exposureToConvolve = self._warper.warpExposure(exposureToNotConvolve.getWcs(), 
                    exposureToConvolve, destBBox = exposureToNotConvolve.getBBox(afwImage.PARENT))
            else:
                pexLog.Trace(self._log.getName(), 1, "ERROR: Input images not registered")
                raise RuntimeError, "Input images not registered"
                
        psfMatchedMaskedImage, psfMatchingKernel, backgroundModel, kernelCellSet = self.matchMaskedImages(
            exposureToConvolve.getMaskedImage(), exposureToNotConvolve.getMaskedImage(),
            footprints = footprints)
        
        psfMatchedExposure = afwImage.makeExposure(psfMatchedMaskedImage, exposureToNotConvolve.getWcs())
        psfMatchedExposure.setFilter(exposureToConvolve.getFilter())
        psfMatchedExposure.setCalib(exposureToNotConvolve.getCalib())

        return (psfMatchedExposure, psfMatchingKernel, backgroundModel, kernelCellSet)
    
    def matchMaskedImages(self, maskedImageToConvolve, maskedImageToNotConvolve, footprints = None):
        """PSF-match a MaskedImage to a reference MaskedImage

        Do the following, in order:
        - Determine a PSF matching kernel and differential background model
            that matches maskedImageToConvolve to maskedImageToNotConvolve
        - Convolve maskedImageToConvolve by the PSF matching kernel
        
        @param maskedImageToConvolve: masked image to PSF-match to the reference masked image;
            must be warped to match the reference masked image
        @param maskedImageToNotConvolve: maskedImage whose PSF is to be matched
        @param footprints: a list of footprints of sources; if None then source detection is run
        
        @return
        - psfMatchedMaskedImage: the PSF-matched masked image =
            maskedImageToConvolve convolved with psfMatchingKernel.
            This has the same xy0, dimensions and wcs as maskedImageToNotConvolve.
        - psfMatchingKernel: the PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to solve for the PSF matching kernel
        
        @raise RuntimeError if input images have different dimensions
        """
        if not self._validateSize(maskedImageToConvolve, maskedImageToNotConvolve):
            pexLog.Trace(self._log.getName(), 1, "ERROR: Input images different size")
            raise RuntimeError, "Input images different size"
            
        self._log.log(pexLog.Log.INFO, "compute PSF-matching kernel")
        kernelCellSet = self._buildCellSet(maskedImageToConvolve,
                                           maskedImageToNotConvolve,
                                           footprints = footprints)
        psfMatchingKernel, backgroundModel = self._solve(kernelCellSet)
    
        self._log.log(pexLog.Log.INFO, "PSF-match science MaskedImage to reference")
        psfMatchedMaskedImage = afwImage.MaskedImageF(maskedImageToConvolve.getBBox(afwImage.PARENT))
        doNormalize = False
        afwMath.convolve(psfMatchedMaskedImage, maskedImageToConvolve, psfMatchingKernel, doNormalize)
        self._log.log(pexLog.Log.INFO, "done")
        return (psfMatchedMaskedImage, psfMatchingKernel, backgroundModel, kernelCellSet)

    def subtractExposures(self, exposureToConvolve, exposureToNotConvolve,
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
        
        @return
        - subtractedExposure: subtracted Exposure = exposureToNotConvolve -
            ((warped exposureToConvolve convolved with psfMatchingKernel) + backgroundModel)
        - psfMatchingKernel: PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to determine PSF matching kernel
        """
        results = self.matchExposures(exposureToConvolve, exposureToNotConvolve,
                                      footprints = footprints,
                                      doWarping = doWarping)

        psfMatchedExposure, psfMatchingKernel, backgroundModel, kernelCellSet = results
        subtractedExposure  = afwImage.ExposureF(exposureToNotConvolve, True)
        smi  = subtractedExposure.getMaskedImage()
        smi -= psfMatchedExposure.getMaskedImage()
        smi -= backgroundModel
        return (subtractedExposure, psfMatchingKernel, backgroundModel, kernelCellSet)

    def subtractMaskedImages(self, maskedImageToConvolve, maskedImageToNotConvolve,
                             footprints = None):
        """Subtract two MaskedImages
        
        Do the following, in order:
        - PSF-match maskedImageToConvolve to maskedImageToNotConvolve
        - Determine the differential background
        - Return the difference: maskedImageToNotConvolve -
            ((warped maskedImageToConvolve convolved with psfMatchingKernel) + backgroundModel)
        
        @param maskedImageToConvolve: MaskedImage to PSF-matched to maskedImageToNotConvolve
        @param maskedImageToNotConvolve: reference MaskedImage
        @param footprints: a list of footprints of sources; if None then source detection is run
        
        @return
        - subtractedMaskedImage = maskedImageToNotConvolve -
            ((maskedImageToConvolve convolved with psfMatchingKernel) + backgroundModel)
        - psfMatchingKernel: PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to determine PSF matching kernel
        """

        results = self.matchMaskedImages(maskedImageToConvolve,
                                         maskedImageToNotConvolve,
                                         footprints = footprints)
        
        psfMatchedMaskedImage, psfMatchingKernel, backgroundModel, kernelCellSet = results 
        subtractedMaskedImage  = afwImage.MaskedImageF(maskedImageToNotConvolve, True)
        subtractedMaskedImage -= psfMatchedMaskedImage
        subtractedMaskedImage -= backgroundModel
        return (subtractedMaskedImage, psfMatchingKernel, backgroundModel, kernelCellSet)

    def _buildCellSet(self, maskedImageToConvolve, maskedImageToNotConvolve, footprints = None):
        """Build a SpatialCellSet for use with the solve method

        @param maskedImageToConvolve: MaskedImage to PSF-matched to maskedImageToNotConvolve
        @param maskedImageToNotConvolve: reference MaskedImage
        @param footprints: a list of footprints of sources; if None then source detection is run
        
        @return kernelCellSet: a SpatialCellSet for use with self._solve
        """
        # Object to store the KernelCandidates for spatial modeling
        kernelCellSet = afwMath.SpatialCellSet(maskedImageToConvolve.getBBox(afwImage.PARENT),
                                               self._policy.getInt("sizeCellX"),
                                               self._policy.getInt("sizeCellY"))
        
        # Candidate source footprints to use for Psf matching
        if footprints == None:
            # If you need to fit for background in ip_diffim, we need
            # to subtract it off before running detection
            if self._policy.get("fitForBackground"):
                self._log.log(pexLog.Log.INFO, "temporarily subtracting backgrounds for detection")
                bkgds = diffimTools.backgroundSubtract(self._policy.getPolicy("afwBackgroundPolicy"),
                                                       [maskedImageToConvolve, maskedImageToNotConvolve])
            
            kcDetect = diffimLib.KernelCandidateDetectionF(self._policy.getPolicy("detectionPolicy"))
            kcDetect.apply(maskedImageToConvolve, maskedImageToNotConvolve)
            footprints = kcDetect.getFootprints()

            if self._policy.get("fitForBackground"):
                maskedImageToConvolve += bkgds[0]
                maskedImageToNotConvolve += bkgds[1]

        # Place candidate footprints within the spatial grid
        for fp in footprints:
            bbox = fp.getBBox()
            
            # Grab the centers in the parent's coordinate system
            xC   = 0.5 * ( bbox.getMinX() + bbox.getMaxX() )
            yC   = 0.5 * ( bbox.getMinY() + bbox.getMaxY() )
            
            tmi  = afwImage.MaskedImageF(maskedImageToConvolve, bbox, afwImage.PARENT)
            smi  = afwImage.MaskedImageF(maskedImageToNotConvolve, bbox, afwImage.PARENT)
            cand = diffimLib.makeKernelCandidate(xC, yC, tmi, smi, self._policy)
            
            pexLog.Trace(self._log.getName(), 5,
                         "Candidate %d at %f, %f" % (cand.getId(), cand.getXCenter(), cand.getYCenter()))
            kernelCellSet.insertCandidate(cand)

        return kernelCellSet
        

