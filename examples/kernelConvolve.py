import eups
import lsst.afw.image as afwImage
import lsst.afw.math  as afwMath
import numpy as num
import sys, os

id    = sys.argv[1]
tRoot = 'tFoot_c_%s' % (id)
iRoot = 'iFoot_c_%s' % (id)
kRoot = 'kernel_c_%s.fits' % (id)

tMi   = afwImage.MaskedImageF(tRoot)
iMi   = afwImage.MaskedImageF(iRoot)
kImg  = afwImage.ImageD(kRoot)
k     = afwMath.FixedKernel(kImg)

cMi   = afwImage.MaskedImageF(tMi.getDimensions())
afwMath.convolve(cMi, tMi, k, False)
cMi.writeFits('c1')

iMi  -= cMi
iMi.writeFits('d1')

gaussFunction = afwMath.GaussianFunction2D(2, 3)
gaussKernel   = afwMath.AnalyticKernel(19, 19, gaussFunction)
kImageIn      = afwImage.ImageD(19, 19)
gaussKernel.computeImage(kImageIn, False)
afwMath.convolve(cMi, tMi, gaussKernel, False)
cMi.writeFits('c2')
