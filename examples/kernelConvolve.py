import lsst.afw.image as afwImage
import lsst.afw.math  as afwMath
import sys

tRoot  = sys.argv[1]
iRoot  = sys.argv[2]
kernel = sys.argv[3]

tMi   = afwImage.MaskedImageF(tRoot)
iMi   = afwImage.MaskedImageF(iRoot)
kImg  = afwImage.ImageD(kernel)
k     = afwMath.FixedKernel(kImg)

cMi   = afwImage.MaskedImageF(tMi.getDimensions())
afwMath.convolve(cMi, tMi, k, False)

goodBBox = k.shrinkBBox(cMi.getBBox(afwImage.LOCAL))
cMi2 = afwImage.MaskedImageF(cMi, goodBBox)
cMi2.writeFits('conv')
