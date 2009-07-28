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

bbox      = afwImage.BBox(afwImage.PointI(k.getCtrX(),
                                          k.getCtrY()) ,
                          afwImage.PointI(cMi.getWidth() - (k.getWidth() - k.getCtrX()),
                                          cMi.getHeight() - (k.getHeight() - k.getCtrY())))
cMi2 = afwImage.MaskedImageF(cMi, bbox)
cMi2.writeFits('conv')
