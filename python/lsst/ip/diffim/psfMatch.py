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
import lsst.pex.logging as pexLog
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim

__all__ = ["PsfMatchToImage", "PsfMatchToModel"]

class PsfMatchToImage(object):
    """PSF-match images to reference images
    """
    def __init__(self, policy, logName="ip.diffim.psfMatchToImage"):
        """Create a PsfMatchToImage
        
        @param policy: see lsst/ip/diffim/policy/PsfMatchingDictionary.paf
        @param logName: name by which messages are logged
        """
        self._log = pexLog.Log(pexLog.Log.getDefaultLog(), "ip.diffim.psfMatchToImage")
        self._policy = policy
    
    def matchExposure(self, exposure, referenceMaskedImage):
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
        psfMatchedMaskedImage, psfMatchingKernel = self.matchMaskedImage(
            exposure.getMaskedImage(), referenceMaskedImage)
        psfMatchedExposure = afwImage.makeExposure(psfMatchedMaskedImage, exposure.getWcs())
        return (psfMatchedExposure, psfMatchingKernel)
    
    def matchMaskedImage(self, maskedImage, referenceMaskedImage):
        """PSF-match a Masked Image to a reference MaskedImage
        
        @param maskedImage: masked image to PSF-match to the reference masked image;
            must be warped to match the reference masked image
        @param referenceMaskedImage: maskedImage whose PSF is to be matched
        
        @returns
        - psfMatchedMaskedImage: the PSF-matched masked image
        - psfMatchingKernel: the PSF matching kernel
        """
        self._log.log(pexLog.Log.INFO, "compute PSF-matching kernel")
        psfMatchingKernel, backgroundModel, kernelCellSet = ipDiffim.psfMatchImageToImage(
            maskedImage, referenceMaskedImage, self._policy)
    
        self._log.log(pexLog.Log.INFO, "PSF-match science MaskedImage to reference")
        psfMatchedMaskedImage = afwImage.MaskedImageF(maskedImage.getBBox(afwImage.PARENT))
        doNormalize = False
        afwMath.convolve(psfMatchedMaskedImage, maskedImage, psfMatchingKernel, doNormalize)
        self._log.log(pexLog.Log.INFO, "done")
        return (psfMatchedMaskedImage, psfMatchingKernel)


class PsfMatchToModel(object):
    """PSF-match PSF models to reference PSF models
    """
    def __init__(self, policy, logName="ip.diffim.psfMatchToImage"):
        """Create a PsfMatchToModel
        
        @param policy: see lsst/ip/diffim/policy/PsfMatchingDictionary.paf
        @param logName: name by which messages are logged
        """
        self._log = pexLog.Log(pexLog.Log.getDefaultLog(), "ip.diffim.psfMatchToImage")
        self._policy = policy
    
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
        psfMatchingKernel, backgroundModel, kernelCellSet = ipDiffim.psfMatchModelToModel(
            referencePsfModel, exposure.getBBox(afwImage.PARENT), exposure.getPsf(), self._policy)
        
        self._log.log(pexLog.Log.INFO, "PSF-match science exposure to reference")
        psfMatchedExposure = afwImage.ExposureF(exposure.getBBox(afwImage.PARENT), exposure.getWcs())
        psfMatchedMaskedImage = psfMatchedExposure.getMaskedImage()

        # Normalize the psf-matching kernel while convolving since its magnitude is meaningless
        # when PSF-matching one model to another.
        doNormalize = True
        afwMath.convolve(psfMatchedMaskedImage, maskedImage, psfMatchingKernel, doNormalize)
        
        self._log.log(pexLog.Log.INFO, "done")
        return (psfMatchedExposure, psfMatchingKernel)
