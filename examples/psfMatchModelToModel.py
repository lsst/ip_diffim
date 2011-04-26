import lsst.ip.diffim                    as ipDiffim

import sys, re, os
import numpy as num

import lsst.daf.base                     as dafBase
import lsst.daf.persistence              as dafPersist
import lsst.afw.detection                as afwDet
import lsst.afw.math                     as afwMath
import lsst.afw.geom                     as afwGeom
import lsst.afw.image                    as afwImage
import lsst.pex.policy                   as pexPolicy
import lsst.pex.logging                  as pexLog

import lsst.afw.display.ds9              as ds9

def pcapsf_read_boost(fn):
    print '# Reading', fn, 
    loc = dafPersist.LogicalLocation(fn)
    storageList = dafPersist.StorageList()
    additionalData = dafBase.PropertySet()
    persistence = dafPersist.Persistence.getPersistence(pexPolicy.Policy())
    storageList.append(persistence.getRetrieveStorage("BoostStorage", loc))
    psfptr = persistence.unsafeRetrieve("pcaPsf", storageList, additionalData)
    psf = afwDet.Psf.swigConvert(psfptr)
    return psf

if __name__ == '__main__':
    import lsst.ip.diffim.diffimTools as diffimTools
    
    pexLog.Trace_setVerbosity("lsst.ip.diffim", 5)

    calexpPath = sys.argv[1]
    boostPath  = sys.argv[2]
    
    psf    = pcapsf_read_boost(boostPath)
    calexp = afwImage.ExposureF(calexpPath) 

    # match to this
    sigGauss = 5
    gaussPsf = afwDet.createPsf("DoubleGaussian", psf.getKernel().getWidth(),
                                psf.getKernel().getHeight(), sigGauss)

    policy = ipDiffim.makeDefaultPolicy()

    # Some policy specializations for model matching

    # We don't want to do any sigma clipping since there are no
    # "outliers"
    policy.set("singleKernelClipping", False)
    policy.set("kernelSumClipping", False)
    policy.set("spatialKernelClipping", False)

    # If doing any deconvolvution, this frequently gets triggered
    policy.set("checkConditionNumber", False)

    # The Kernel size is limited by the Psf size
    kernelSize = psf.getKernel().getWidth() // 2
    if (kernelSize // 2) != 1:
        kernelSize -= 1
    policy.set("kernelSize", kernelSize)

    # We have no background
    policy.set("fitForBackground", False)

    # Use AL basis set
    policy.set("kernelBasisSet", "alard-lupton")
    policy.set("usePcaForSpatialKernel", False)

    # We need to make sure that not too many bases are being used,
    # compared to the number of non-EDGE pixels in "Psf image
    # convolved with AL basis" image.
    nPixelsUnmasked = (psf.getKernel().getWidth() - kernelSize // 2) * \
                      (psf.getKernel().getHeight() - kernelSize // 2)
    basisList = ipDiffim.makeKernelBasisList(policy)

    if nPixelsUnmasked < (kernelSize * kernelSize):
        pexLog.Trace("lsst.ip.diffim.psfMatchModelToModel", 2,
                     "Warning : fewer unmasked pixels in image (%d) than in kernel (%d)" \
                     % (nPixelsUnmasked,
                        len(basisList)))
        
    if nPixelsUnmasked < len(basisList):
        pexLog.Trace("lsst.ip.diffim.psfMatchModelToModel", 2,
                     "Warning : fewer unmasked pixels in image (%d) than basis components (%d)" \
                     % (nPixelsUnmasked,
                        len(basisList)))
            
    # Make sure that there are not too many bases for the number of
    # unmasked pixels in "Psf convolved with the Kernel" image;
    # consider reducing basis set size
    #policy.set("alardNGauss", 2)
    #policy.set("alardSigGauss", 0.5)
    #policy.add("alardSigGauss", 3.)
    #policy.set("alardDegGauss", 3)
    #policy.add("alardDegGauss", 2)

    # Spatial order of Psf.
    #
    # Infer from the number of spatial parameters.
    # (O + 1) * (O + 2) / 2 = N
    # O^2 + 3 * O + 2 * (1 - N) = 0
    #
    # Roots are [-3 +/- sqrt(9 - 8 * (1 - N))] / 2
    #
    nParameters = psf.getKernel().getNSpatialParameters()
    root        = num.sqrt(9 - 8 * (1 - nParameters))
    assert(root == root // 1)   # We know its an integer solution
    order       = (root - 3) / 2
    assert(order == order // 1) # Ditto
    policy.set("spatialKernelOrder", int(order))


    imageBBox = afwGeom.Box2I(afwGeom.Point2I(calexp.getX0(), calexp.getY0()),
                              afwGeom.Extent2I(calexp.getWidth(), calexp.getHeight()))
    sk, sb, kcs = ipDiffim.psfMatchModelToModel(gaussPsf, imageBBox, psf, policy)

    cim = afwImage.MaskedImageF(calexp.getMaskedImage().getDimensions())
    afwMath.convolve(cim, calexp.getMaskedImage(), sk, False)
    cexp = afwImage.ExposureF(cim, calexp.getWcs())
    cexp.writeFits(re.sub('.fits', '_remap.fits', calexpPath))
    
    if 0:
        diffimTools.displaySpatialKernelQuality(kcs, sk, sb, frame = 1)
    else:
        ds9.setMaskPlaneVisibility("DETECTED", False)
        ds9.mtv(calexp, frame = 1)
        diffimTools.displaySpatialKernelMosaic(psf.getKernel(), calexp.getWidth(), calexp.getHeight(),
                                               frame = 2)
        diffimTools.displaySpatialKernelMosaic(gaussPsf.getKernel(), calexp.getWidth(), calexp.getHeight(),
                                               frame = 3)
        diffimTools.displayKernelMosaic(kcs, frame = 4)
        diffimTools.displayCandidateMosaic(kcs, frame = 5)
        diffimTools.displaySpatialKernelMosaic(sk, calexp.getWidth(), calexp.getHeight(), frame = 6)
        ds9.mtv(cim, frame = 7)
    

# foreach i ( 85501910 85563260 85597903 85656377 85661696 )
# python examples/psfMatchModelToModel.py ~/LSST/PT1/psfMatch/pt1prod_im0133/update/calexp/v"$i"-fr/R13/S11.fits ~/LSST/PT1/psfMatch/pt1prod_im0133/update/psf/v"$i"-fr/R13/S11.boost
# end
