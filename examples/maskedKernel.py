import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools

import lsst.afw.display.ds9 as ds9
ds9.ds9Cmd("mask transparency 50")

import lsst.pex.logging as pexLogging
verbosity = 5
pexLogging.Trace_setVerbosity('lsst.ip.diffim', verbosity)

policy1  = ipDiffim.makeDefaultPolicy()
policy1.set("kernelBasisSet", "delta-function")
policy1.set("useRegularization", False)
kList1 = ipDiffim.makeKernelBasisList(policy1)
soln1  = ipDiffim.MaskedKernelSolutionF(kList1, False)

policy2  = ipDiffim.makeDefaultPolicy()
policy2.set("kernelBasisSet", "delta-function")
policy2.set("useRegularization", True)
kList2 = ipDiffim.makeKernelBasisList(policy2)
soln2  = ipDiffim.MaskedKernelSolutionF(kList2, False)

policy3  = ipDiffim.makeDefaultPolicy()
policy3.set("kernelBasisSet", "alard-lupton")
kList3 = ipDiffim.makeKernelBasisList(policy3)
soln3  = ipDiffim.MaskedKernelSolutionF(kList3, False)

#scienceImage  = afwImage.ExposureF('s2L007173_0100g4TANSIPwInv.fits')
#templateImage = afwImage.ExposureF('oneTemplate100006_0072g4.fits')
#xCenter       = 713
#yCenter       = 641

scienceImage  = afwImage.ExposureF('s2L006417_0516g4TANSIPwInv.fits')
templateImage = afwImage.ExposureF('s2L100006_05000501g4.fits')
xCenter       = 572
yCenter       = 106


##### TO DO; DEAL WITH EDGE NANS; REALLY NEED TO SPREAD MASK...
stampSize     = 200
maskSize      = 10
candBBox      = afwGeom.Box2I(afwGeom.Point2I(xCenter - stampSize//2, yCenter - stampSize//2),
                              afwGeom.Extent2I(stampSize, stampSize))
maskBBox      = afwGeom.Box2I(afwGeom.Point2I(xCenter - maskSize//2, yCenter - maskSize//2),
                              afwGeom.Extent2I(maskSize, maskSize))

diffimTools.backgroundSubtract(policy1.getPolicy("afwBackgroundPolicy"),
                               [scienceImage.getMaskedImage(),])
templateImage = ipDiffim.warpTemplateExposure(templateImage,
                                              scienceImage,
                                              policy1.getPolicy("warpingPolicy"))

templateImage = templateImage.getMaskedImage()
scienceImage = scienceImage.getMaskedImage()

# Create deep copy since this will be modified
pixelMask  = afwImage.MaskU(scienceImage.getMask(), True)
# Add in template mask
pixelMask |= templateImage.getMask()
# And mask out the variability!
# just hack it for now since its late
for y in range(yCenter - maskSize//2, yCenter + maskSize//2):
    for x in range(xCenter - maskSize//2, xCenter + maskSize//2):
        maskVal  = pixelMask.get(x,y)
        maskVal |= afwImage.MaskU_getPlaneBitMask("BAD")
        pixelMask.set(x, y, maskVal)

tsi = afwImage.MaskedImageF(templateImage, candBBox, afwImage.LOCAL)
ssi = afwImage.MaskedImageF(scienceImage, candBBox, afwImage.LOCAL)

ds9.mtv(tsi, frame = 1)
ds9.mtv(ssi, frame = 2)

for soln in (soln1, soln2, soln3): # soln2, soln3):
    soln.buildSingleMask(tsi.getImage(),
                         ssi.getImage(),
                         ssi.getVariance(),
                         maskBBox)
    soln.solve()

    #soln.build(tsi.getImage(),
    #           ssi.getImage(),
    #           ssi.getVariance(),
    #           afwImage.MaskU(pixelMask, candBBox))
    #soln.solve()
    
k1 = soln1.getKernel()
k2 = soln2.getKernel()
k3 = soln3.getKernel()

d1 = ipDiffim.convolveAndSubtract(tsi, ssi, k1, 0)
d2 = ipDiffim.convolveAndSubtract(tsi, ssi, k2, 0)
d3 = ipDiffim.convolveAndSubtract(tsi, ssi, k3, 0)
                                  

ds9.mtv(tsi, frame = 1)
ds9.mtv(ssi, frame = 2)
ds9.mtv(soln1.makeKernelImage(), frame = 3)
ds9.mtv(d1, frame = 4)
#
ds9.mtv(tsi, frame = 5)
ds9.mtv(ssi, frame = 6)
ds9.mtv(soln2.makeKernelImage(), frame = 7)
ds9.mtv(d2, frame = 8)
#
ds9.mtv(tsi, frame = 9)
ds9.mtv(ssi, frame = 10)
ds9.mtv(soln3.makeKernelImage(), frame = 11)
ds9.mtv(d3, frame = 12)
