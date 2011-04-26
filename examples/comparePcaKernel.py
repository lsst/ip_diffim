#!/usr/bin/env python
import os
import sys
import eups
import lsst.afw.geom as afwGeom
import lsst.afw.image.imageLib as afwImage
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as pexLogging
import lsst.afw.display.ds9 as ds9

display = True

verbosity = 5
pexLogging.Trace_setVerbosity("lsst.ip.diffim", verbosity)

defDataDir   = eups.productDir("afwdata") 
imageProcDir = eups.productDir("ip_diffim")

if len(sys.argv) == 1:
    defTemplatePath = os.path.join(defDataDir, "CFHT", "D4", "cal-53535-i-797722_2_tmpl")
    defSciencePath  = os.path.join(defDataDir, "CFHT", "D4", "cal-53535-i-797722_2")
    templateMaskedImage = afwImage.MaskedImageF(defTemplatePath)
    scienceMaskedImage  = afwImage.MaskedImageF(defSciencePath)
    bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(512, 512))
    templateMaskedImage = afwImage.MaskedImageF(templateMaskedImage, bbox, afwImage.LOCAL)
    scienceMaskedImage  = afwImage.MaskedImageF(scienceMaskedImage, bbox, afwImage.LOCAL)
    
elif len(sys.argv) == 3:
    defTemplatePath = sys.argv[1]
    defSciencePath  = sys.argv[2]
    templateMaskedImage = afwImage.MaskedImageF(defTemplatePath)
    scienceMaskedImage  = afwImage.MaskedImageF(defSciencePath)
else:
    sys.exit(1)
    
policy = ipDiffim.makeDefaultPolicy()


# same for all kernels
policy.set("singleKernelClipping", True)
policy.set("kernelSumClipping", True)
policy.set("spatialKernelClipping", False)
policy.set("spatialKernelOrder", 0)
policy.set("spatialBgOrder", 0)
policy.set("usePcaForSpatialKernel", True)
policy.set("fitForBackground", True)


kcDetect = ipDiffim.KernelCandidateDetectionF(policy.getPolicy("detectionPolicy"))
kcDetect.apply(templateMaskedImage, scienceMaskedImage)
footprints = kcDetect.getFootprints()

# specific to delta function
policy.set("kernelBasisSet", "delta-function")
policy.set("useRegularization", False)
spatialKernel1, spatialBg1, kernelCellSet1 = ipDiffim.psfMatchImageToImage(templateMaskedImage,
                                                                           scienceMaskedImage,
                                                                           policy,
                                                                           footprints)

# alard lupton
policy.set("kernelBasisSet", "alard-lupton")
policy.set("useRegularization", False)
spatialKernel2, spatialBg2, kernelCellSet2 = ipDiffim.psfMatchImageToImage(templateMaskedImage,
                                                                           scienceMaskedImage,
                                                                           policy,
                                                                           footprints)

# regularized delta function
policy.set("kernelBasisSet", "delta-function")
policy.set("useRegularization", True)
spatialKernel3, spatialBg3, kernelCellSet3 = ipDiffim.psfMatchImageToImage(templateMaskedImage,
                                                                           scienceMaskedImage,
                                                                           policy,
                                                                           footprints)

basisList1 = spatialKernel1.getKernelList()
basisList2 = spatialKernel2.getKernelList()
basisList3 = spatialKernel3.getKernelList()

frame = 1
for idx in range(min(5, len(basisList1))):
    kernel = basisList1[idx]
    im     = afwImage.ImageD(spatialKernel1.getDimensions())
    ksum   = kernel.computeImage(im, False)    
    ds9.mtv(im, frame=frame)
    frame += 1

for idx in range(min(5, len(basisList2))):
    kernel = basisList2[idx]
    im     = afwImage.ImageD(spatialKernel2.getDimensions())
    ksum   = kernel.computeImage(im, False)    
    ds9.mtv(im, frame=frame)
    frame += 1


for idx in range(min(5, len(basisList3))):
    kernel = basisList3[idx]
    im     = afwImage.ImageD(spatialKernel3.getDimensions())
    ksum   = kernel.computeImage(im, False)    
    ds9.mtv(im, frame=frame)
    frame += 1

