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
import sys
import diffimLib
import lsst.pex.logging as pexLog
import lsst.pex.exceptions as pexExcept
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

from .warpTemplateExposure import warpTemplateExposure
import diffimTools

__all__ = ["PsfMatch", "ImagePsfMatch", "ModelPsfMatch"]

class PsfMatch(object):
    """Base class for PSF matching
    """
    def __init__(self, policy, logName="lsst.ip.diffim.PsfMatch"):
        """Create a PsfMatchToImage
        
        @param policy: see lsst/ip/diffim/policy/PsfMatchingDictionary.paf
        @param logName: name by which messages are logged
        """
        self._policy = policy
        self._log = pexLog.Log(pexLog.Log.getDefaultLog(), logName)

    def solve(self, kernelCellSet, returnOnExcept = False):
        # Create the Psf matching kernel
        try:
            kb = diffimLib.fitSpatialKernelFromCandidates(kernelCellSet, self._policy)
        except pexExcept.LsstCppException, e:
            pexLog.Trace(self._log.getName(), 1, "ERROR: Unable to calculate psf matching kernel")
            pexLog.Trace(self._log.getName(), 2, e.args[0].what())
    
            if returnOnExcept:
                return None, None
            else:
                raise
        else:
            spatialKernel = kb.first
            spatialBg     = kb.second
    
        # What is the status of the processing?
        nGood = 0
        for cell in kernelCellSet.getCellList():
            for cand in cell.begin(True):
                cand = diffimLib.cast_KernelCandidateF(cand)
                if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                    nGood += 1
        if nGood == 0:
            pexLog.Trace(self._log.getName(), 1, "WARNING")
        pexLog.Trace(self._log.getName(), 1, "Used %d kernels for spatial fit" % (nGood))
    
        return spatialKernel, spatialBg
            
class ImagePsfMatch(PsfMatch):
    """PSF-match images to reference images
    """
    def __init__(self, policy, logName="lsst.ip.diffim.ImagePsfMatch"):
        """Create a PsfMatchToImage
        
        @param policy: see lsst/ip/diffim/policy/PsfMatchingDictionary.paf
        @param logName: name by which messages are logged
        """
        PsfMatch.__init__(self, policy, logName)

    def validateSize(self, maskedImageToConvolve, maskedImageToNotConvolve):
        # Make sure they are the same size
        return maskedImageToConvolve.getDimensions() == maskedImageToNotConvolve.getDimensions()
    
    def validateWcs(self, exposureToConvolve, exposureToNotConvolve):
        # Make sure they end up the same dimensions on the sky
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
             (exposureToConvolve.getHeight() != exposureToNotConvolve.getHeight()) or \
             (exposureToConvolve.getWidth()  != exposureToNotConvolve.getWidth()) ):
            return False
        return True

    def subtractExposures(self, exposureToConvolve, exposureToNotConvolve,
                          footprints = None, doWarping = True):

        results = self.matchExposures(exposureToConvolve, exposureToNotConvolve,
                                      footprints = footprints,
                                      doWarping = doWarping)

        psfMatchedExposure, psfMatchingKernel, backgroundModel, kernelCellSet = results
        subtractedExposure  = afwImage.ExposureF(exposureToNotConvolve, True)
        smi  = subtractedExposure.getMaskedImage()
        smi -= psfMatchedExposure.getMaskedImage()
        smi -= backgroundModel
        return (subtractedExposure, psfMatchingKernel, backgroundModel, kernelCellSet)
        
    def matchExposures(self, exposureToConvolve, exposureToNotConvolve,
                       footprints = None, doWarping = True):
        """Match an exposure to the reference
        
        @param exposure: Exposure to PSF-match to the reference masked image;
            must be warped to match the reference masked image
        @param referenceMaskedImage: maskedImage whose PSF is to be matched
        
        @returns
        - psfMatchedExposure: the PSF-matched exposure.
            This has the same xy0, dimensions and wcs as exposure but no psf.
            There is no psf because the PSF-matching process does not compute one.
        - psfMatchingKernel: the PSF matching kernel
        """

        if not self.validateWcs(exposureToConvolve, exposureToNotConvolve):
            if doWarping:
                pexLog.Trace(self._log.getName(), 1,
                             "Astrometrically registering template to science image")
                exposureToConvolve = warpTemplateExposure(exposureToConvolve,
                                                          exposureToNotConvolve,
                                                          self._policy.getPolicy("warpingPolicy"))
            else:
                pexLog.Trace(self._log.getName(), 1, "ERROR: Input images not registered")
                raise RuntimeError, "Input images not registered"
                
        psfMatchedMaskedImage, psfMatchingKernel, backgroundModel, kernelCellSet = self.matchMaskedImages(
            exposureToConvolve.getMaskedImage(), exposureToNotConvolve.getMaskedImage(),
            footprints = footprints
            )
        
        psfMatchedExposure = afwImage.makeExposure(psfMatchedMaskedImage, exposureToNotConvolve.getWcs())
        return (psfMatchedExposure, psfMatchingKernel, backgroundModel, kernelCellSet)

    def subtractMaskedImages(self, maskedImageToConvolve, maskedImageToNotConvolve,
                             footprints = None):

        results = self.matchMaskedImages(maskedImageToConvolve,
                                         maskedImageToNotConvolve,
                                         footprints = footprints)
        
        psfMatchedMaskedImage, psfMatchingKernel, backgroundModel, kernelCellSet = results 
        subtractedMaskedImage  = afwImage.MaskedImageF(maskedImageToNotConvolve, True)
        subtractedMaskedImage -= psfMatchedMaskedImage
        subtractedMaskedImage -= backgroundModel
        return (subtractedMaskedImage, psfMatchingKernel, backgroundModel, kernelCellSet)
    
    def matchMaskedImages(self, maskedImageToConvolve, maskedImageToNotConvolve, footprints = None):
        """PSF-match a Masked Image to a reference MaskedImage
        
        @param maskedImage: masked image to PSF-match to the reference masked image;
            must be warped to match the reference masked image
        @param referenceMaskedImage: maskedImage whose PSF is to be matched
        
        @returns
        - psfMatchedMaskedImage: the PSF-matched masked image
        - psfMatchingKernel: the PSF matching kernel
        """
        if not self.validateSize(maskedImageToConvolve, maskedImageToNotConvolve):
            pexLog.Trace(self._log.getName(), 1, "ERROR: Input images different size")
            raise RuntimeError, "Input images different size"
            
        self._log.log(pexLog.Log.INFO, "compute PSF-matching kernel")
        kernelCellSet = self._buildCellSet(maskedImageToConvolve,
                                           maskedImageToNotConvolve,
                                           footprints = footprints)
        psfMatchingKernel, backgroundModel = self.solve(kernelCellSet)
    
        self._log.log(pexLog.Log.INFO, "PSF-match science MaskedImage to reference")
        psfMatchedMaskedImage = afwImage.MaskedImageF(maskedImageToConvolve.getBBox(afwImage.PARENT))
        doNormalize = False
        afwMath.convolve(psfMatchedMaskedImage, maskedImageToConvolve, psfMatchingKernel, doNormalize)
        self._log.log(pexLog.Log.INFO, "done")
        return (psfMatchedMaskedImage, psfMatchingKernel, backgroundModel, kernelCellSet)


    def _buildCellSet(self, maskedImageToConvolve, maskedImageToNotConvolve, footprints = None):
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
        
        

class ModelPsfMatch(PsfMatch):
    """PSF-match PSF models to reference PSF models
    """
    def __init__(self, policy, logName="lsst.ip.diffim.ModelPsfMatch"):
        """Create a PsfMatchToModel
        
        @param policy: see lsst/ip/diffim/policy/PsfMatchingDictionary.paf
        @param logName: name by which messages are logged
        """
        PsfMatch.__init__(self, policy, logName)
    
    def matchExposure(self, exposure, referencePsfModel):
        """PSF-match an exposure to a model.
        
        @param exposure: Exposure to PSF-match to the reference masked image;
            must contain a PSF model and must be warped to match the reference masked image
        @param referencePsfModel: PSF model to match (an afwDetection.Psf)
        
        @returns
        - psfMatchedExposure: the PSF-matched exposure.
            This has the same xy0, dimensions and wcs as exposure but no psf.
            In theory the psf should equal referencePsfModel but the match is likely not exact.
        - psfMatchingKernel: the PSF matching kernel
        """
        if not exposure.hasPsf():
            raise RuntimeError("exposure does not contain a PSF model")

        maskedImage = exposure.getMaskedImage()

        self._log.log(pexLog.Log.INFO, "compute PSF-matching kernel")
        kernelCellSet = self._buildCellSet(referencePsfModel,
                                           exposure.getBBox(afwImage.PARENT),
                                           exposure.getPsf())
        psfMatchingKernel, backgroundModel = self.solve(kernelCellSet)
        
        self._log.log(pexLog.Log.INFO, "PSF-match science exposure to reference")
        psfMatchedExposure = afwImage.ExposureF(exposure.getBBox(afwImage.PARENT), exposure.getWcs())
        psfMatchedMaskedImage = psfMatchedExposure.getMaskedImage()

        # Normalize the psf-matching kernel while convolving since its magnitude is meaningless
        # when PSF-matching one model to another.
        doNormalize = True
        afwMath.convolve(psfMatchedMaskedImage, maskedImage, psfMatchingKernel, doNormalize)
        
        self._log.log(pexLog.Log.INFO, "done")
        return (psfMatchedExposure, psfMatchingKernel, kernelCellSet)

    def _buildCellSet(self, referencePsfModel, scienceBBox, sciencePsfModel, mergePolicy = False):
        
        if (referencePsfModel.getKernel().getDimensions() != sciencePsfModel.getKernel().getDimensions()):
            pexLog.Trace(self._log.getName(), 1,
                         "ERROR: Dimensions of reference Psf and science Psf different; exiting")
            raise RuntimeError, "ERROR: Dimensions of reference Psf and science Psf different; exiting"
    
        kernelWidth, kernelHeight = referencePsfModel.getKernel().getDimensions()
        maxPsfMatchingKernelSize = 1 + (min(kernelWidth - 1, kernelHeight - 1) // 2)
        if maxPsfMatchingKernelSize % 2 == 0:
            maxPsfMatchingKernelSize -= 1
        if self._policy.get('kernelSize') > maxPsfMatchingKernelSize:
            pexLog.Trace(self._log.getName(), 1,
                         "WARNING: Resizing matching kernel to size %d x %d" % (maxPsfMatchingKernelSize,
                                                                                maxPsfMatchingKernelSize))
            self._policy.set('kernelSize', maxPsfMatchingKernelSize)
    
        # Chanes to policy particular for matchPsfModels
        if mergePolicy:
            policyFile = pexPolicy.DefaultPolicyFile("ip_diffim", "MatchPsfModels.paf", "policy")
            matchPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)
            matchPolicy.mergeDefaults(self._policy.getDictionary())
            self._policy = matchPolicy
        
        regionSizeX, regionSizeY = scienceBBox.getDimensions()
        scienceX0,   scienceY0   = scienceBBox.getMin()
    
        sizeCellX = self._policy.get("sizeCellX")
        sizeCellY = self._policy.get("sizeCellY")
    
        kernelCellSet = afwMath.SpatialCellSet(
            afwGeom.Box2I(afwGeom.Point2I(scienceX0, scienceY0),
                          afwGeom.Extent2I(regionSizeX, regionSizeY)),
            sizeCellX, sizeCellY
            )
    
        nCellX    = regionSizeX // sizeCellX
        nCellY    = regionSizeY // sizeCellY
        dimenR    = referencePsfModel.getKernel().getDimensions()
        dimenS    = sciencePsfModel.getKernel().getDimensions()
    
        for row in range(nCellY):
            # place at center of cell
            posY = sizeCellY * row + sizeCellY // 2 + scienceY0
            
            for col in range(nCellX):
                # place at center of cell
                posX = sizeCellX * col + sizeCellX // 2 + scienceX0
    
                pexLog.Trace(self._log.getName(), 5, "Creating Psf candidate at %.1f %.1f" % (posX, posY))
    
                # reference kernel image, at location of science subimage
                kernelImageR = referencePsfModel.computeImage(afwGeom.Point2D(posX, posY), True).convertF()
                sum = afwMath.makeStatistics(kernelImageR, afwMath.SUM).getValue(afwMath.SUM)
                kernelImageR /= sum
                kernelMaskR   = afwImage.MaskU(dimenR)
                kernelMaskR.set(0)
                kernelVarR    = afwImage.ImageF(kernelImageR, True)
                referenceMI   = afwImage.MaskedImageF(kernelImageR, kernelMaskR, kernelVarR)
     
                kernelImageS = sciencePsfModel.computeImage(afwGeom.Point2D(posX, posY), True).convertF()
                sum = afwMath.makeStatistics(kernelImageS, afwMath.SUM).getValue(afwMath.SUM)
                kernelImageS /= sum
                kernelMaskS   = afwImage.MaskU(dimenS)
                kernelMaskS.set(0)
                kernelVarS    = afwImage.ImageF(kernelImageS, True)
                scienceMI     = afwImage.MaskedImageF(kernelImageS, kernelMaskS, kernelVarS)
    
                # The image to convolve is the science image, to the reference Psf.
                kc = diffimLib.makeKernelCandidate(posX, posY, scienceMI, referenceMI, self._policy)
                kernelCellSet.insertCandidate(kc)

        return kernelCellSet
            
            
    
