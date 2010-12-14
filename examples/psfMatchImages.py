import sys, os
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as pexLog

verbosity = 3
pexLog.Trace_setVerbosity('lsst.ip.diffim', verbosity)

imageToConvolve     = afwImage.MaskedImageF(sys.argv[1])
imageToNotConvolve  = afwImage.MaskedImageF(sys.argv[2])
outputImage         = sys.argv[3]

policy              = ipDiffim.createDefaultPolicy()
policy.set("detThreshold", 5.0)
policy.set("kernelBasisSet", "alard-lupton")
policy.set("usePcaForSpatialKernel", True)
policy.set("spatialKernelOrder", 1)

spatialKernel, spatialBg, kernelCellSet = ipDiffim.psfMatchImageToImage(imageToConvolve,
                                                                        imageToNotConvolve,
                                                                        policy)
 
cMi = afwImage.MaskedImageF(imageToConvolve.getDimensions())
afwMath.convolve(cMi, imageToConvolve, spatialKernel, False)
cMi.writeFits(outputImage)

