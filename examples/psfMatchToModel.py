#!/usr/bin/env python

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

"""PSF-match a masked image to model PSF
"""
import os
import sys
import math

import lsst.daf.base as dafBase
import lsst.daf.persistence as dafPersist
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.ip.diffim as ipDiffim
import lsst.coadd.utils as coaddUtils

PsfTypeName = "pcaPsf"

PackageName = "ip_diffim"
PolicyDictName = "PsfMatchingDictionary.paf"

FWHMPerSigma = 2 * math.sqrt(2 * math.log(2))

# set BBox to an afwGeom.BoxI to process a subframe of each input exposure
# set to an empty bbox to process the whole image
BBox = afwImage.BBox(afwImage.PointI(0, 0), 1024, 1024)

def psfFromBoost(psfPath, typeName):
    """Unpersist a PSF from a boost file.
    
    @param psfPath: path to boost file
    @param typeName: type of PSF, e.g. "dgPsf", "pcaPsf"
    """
    loc = dafPersist.LogicalLocation(psfPath)
    storageList = dafPersist.StorageList()
    additionalData = dafBase.PropertySet()
    persistence = dafPersist.Persistence.getPersistence(pexPolicy.Policy())
    storageList.append(persistence.getRetrieveStorage("BoostStorage", loc))
    psfptr = persistence.unsafeRetrieve(typeName, storageList, additionalData)
    psf = afwDet.Psf.swigConvert(psfptr)
    return psf

def psfMatchToModel(exposurePath, psfPath, desFwhm, policy):
    exposure = afwImage.ExposureF(exposurePath, 0, BBox)
    maskedImage = exposure.getMaskedImage()

    exposurePsf = psfFromBoost(psfPath, PsfTypeName)
    exposurePsfKernel = exposurePsf.getKernel()

    ctrXPos = maskedImage.indexToPosition((maskedImage.getWidth()  - 1) / 2, afwImage.X)
    ctrYPos = maskedImage.indexToPosition((maskedImage.getHeight() - 1) / 2, afwImage.Y)
    print "Dimensions=%d, %d; xy0=%d, %d, ctrPos=%0.1f, %0.1f" % (
        maskedImage.getWidth(), maskedImage.getHeight(),
        maskedImage.getX0(), maskedImage.getY0(),
        ctrXPos, ctrYPos)
    
    kernelWidth = exposurePsfKernel.getWidth()
    kernelHeight = exposurePsfKernel.getHeight()
    print "Create double Gaussian PSF model with core fwhm %0.1f and size %dx%d" % \
        (desFwhm, kernelWidth, kernelHeight)
    coreSigma = desFwhm / FWHMPerSigma
    modelPsf = afwDet.createPsf("DoubleGaussian", kernelWidth, kernelHeight,
        coreSigma, coreSigma * 2.5, 0.1)
    kernelImage = afwImage.ImageD(kernelWidth, kernelHeight)
    
    # make sure PSF-matching kernel is small enough to work with the exposure's PSF kernel
    # it must be less than half as big
    maxPsfMatchingKernelSize = 1 + (min(kernelWidth - 1, kernelHeight - 1) // 2)
    if maxPsfMatchingKernelSize%2 == 0:
        maxPsfMatchingKernelSize -= 1
        
    if maxPsfMatchingKernelSize < 5:
        raise RuntimeError("PSF is too small to do a useful job")
    desKernelSize = policy.get("kernelSize")
    if desKernelSize > maxPsfMatchingKernelSize:
        print "Reducing size of PSF matching kernel from %s to %s" % \
            (desKernelSize, maxPsfMatchingKernelSize)
        policy.set("kernelSize", maxPsfMatchingKernelSize)

    kernelSum = modelPsf.getKernel().computeImage(kernelImage, False, ctrXPos, ctrYPos)
    print "Model PSF kernel sum=%s" % (kernelSum,)
    kernelImage.writeFits("modelPsf.fits")

    kernelSum = exposurePsfKernel.computeImage(kernelImage, False, ctrXPos, ctrYPos)
    print "Exposure PSF kernel sum=%s at center of exposure" % (kernelSum,)
    kernelImage.writeFits("exposurePsfKernel.fits")

    print "Compute PSF-matching kernel"
    psfMatchingKernel, backgroundModel, kernelCellSet = ipDiffim.psfMatchModelToModel(
        modelPsf, coaddUtils.bboxFromImage(exposure), exposurePsf, policy)

    psfMatchingKernelImage = afwImage.ImageD(
        psfMatchingKernel.getWidth(), psfMatchingKernel.getHeight())
    kernelSum = psfMatchingKernel.computeImage(psfMatchingKernelImage, False, ctrXPos, ctrYPos)
    print "PSF matching kernel sum = %s at center of exposure" % (kernelSum,)
    kernelImage.writeFits("psfMatchingKernel.fits")
    
    print "Convolve exposure with PSF matching kernel"
    psfMatchedMaskedImage = afwImage.MaskedImageF(maskedImage.getDimensions())
    # Normalize the psf-matching kernel while convolving since its magnitude is meaningless
    # when PSF-matching one model to another.
    doNormalize = True
    afwMath.convolve(psfMatchedMaskedImage, maskedImage, psfMatchingKernel, doNormalize)
    
    psfMatchedExposure = afwImage.makeExposure(psfMatchedMaskedImage, exposure.getWcs())
    psfMatchedExposure.setXY0(exposure.getXY0())
    psfMatchedExposure.writeFits("psfMatched.fits")


if __name__ == "__main__":
    pexLog.Trace.setVerbosity('lsst.coadd', 5)
    pexLog.Trace.setVerbosity('lsst.ip.diffim', 5)

    helpStr = """Usage: warpPsfMatchToModelAndCoadd.py exposurePath psfPath desFwhm [policy]

where:
- exposurePath: path to desired exposure
- psfPath: path to boost-persisted psf for exposure
- desFwhm: FWHM of desired PSF (pixels);
    the desired PSF model is a double Gaussian:
    - the core has the user-specified FWHM
    - the wings have a Gaussian with 2.5 the core's and amplitude 0.1 the core
- policy: path to a policy file

The policy dictionary is: policy/%s
""" % (PolicyDictName,)
    if len(sys.argv) not in (3, 4):
        print helpStr
        sys.exit(0)
    
    exposurePath = sys.argv[1]
    psfPath = sys.argv[2]
    desFwhm = float(sys.argv[3])

    if len(sys.argv) > 4:
        policyPath = sys.argv[4]
        policy = pexPolicy.Policy(policyPath)
    else:
        policy = pexPolicy.Policy()

    policyFile = pexPolicy.DefaultPolicyFile(PackageName, PolicyDictName, "policy")
    defPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)
    policy.mergeDefaults(defPolicy.getDictionary())

    psfMatchToModel(exposurePath, psfPath, desFwhm, policy)    
