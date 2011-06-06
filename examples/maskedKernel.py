import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
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

#scienceExposure  = afwImage.ExposureF('s2L007173_0100g4TANSIPwInv.fits')
#templateExposure = afwImage.ExposureF('oneTemplate100006_0072g4.fits')
#xCenter       = 713
#yCenter       = 641

scienceExposure  = afwImage.ExposureF('s2L006417_0516g4TANSIPwInv.fits')
templateExposure = afwImage.ExposureF('s2L100006_05000501g4.fits')
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
                               [scienceExposure.getMaskedImage(),])

warper = afwMath.Warper.fromPolicy(policy1.getPolicy("warpingPolicy"))
templateExposure = warper.warpExposure(scienceExposure.getWcs(), templateExposure,
                destBBox = scienceExposure.getBBox(afwImage.PARENT))

templateExposure = templateExposure.getMaskedImage()
scienceExposure = scienceExposure.getMaskedImage()

# Create deep copy since this will be modified
pixelMask  = afwImage.MaskU(scienceExposure.getMask(), True)
# Add in template mask
pixelMask |= templateExposure.getMask()
# And mask out the variability!
afwDet.setMaskFromFootprint(pixelMask, afwDet.Footprint(maskBBox), afwImage.MaskU_getPlaneBitMask("BAD"))

# Grab subimages
tsi = afwImage.MaskedImageF(templateExposure, candBBox, afwImage.PARENT)
ssi = afwImage.MaskedImageF(scienceExposure, candBBox, afwImage.PARENT)
msi = afwImage.MaskU(pixelMask, candBBox, afwImage.PARENT)

#import pdb; pdb.set_trace()
bbox   = ssi.getBBox(afwImage.PARENT)
fp     = afwDet.Footprint(bbox)
isbad  = afwImage.MaskU_getPlaneBitMask("BAD")
ds9.mtv(afwImage.ImageU(msi.getArray()), frame=0)

mask   = afwDet.footprintAndMask(fp, msi, 0x2)
m2 = msi.Factory(msi.getBBox(afwImage.PARENT))
afwDet.setMaskFromFootprint(m2, mask, isbad)
ds9.mtv(afwImage.ImageU(m2.getArray()), frame=1)

print 'cand', candBBox
print 'mask', maskBBox

ds9.mtv(tsi, frame = 1)
ds9.mtv(ssi, frame = 2)

for soln in (soln1, soln2, soln3): # soln2, soln3):
    soln.build(tsi.getImage(),
               ssi.getImage(),
               ssi.getVariance(),
               msi)
    soln.solve()
    
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
