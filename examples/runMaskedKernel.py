import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.pex.logging as pexLog
import sys, os, re

import lsst.afw.display.ds9 as ds9
ds9.ds9Cmd("mask transparency 50")

import lsst.pex.logging as pexLogging
verbosity = 5
pexLogging.Trace_setVerbosity('lsst.ip.diffim', verbosity)

infile = sys.argv[1]

# Diffim policy
policy  = ipDiffim.makeDefaultPolicy()
#policy.set("kernelBasisSet", "delta-function")
#policy.set("useRegularization", False)
policy.set("kernelBasisSet", "alard-lupton")
policy.set("fitForBackground", True)
kList = ipDiffim.makeKernelBasisList(policy)
soln  = ipDiffim.MaskedKernelSolutionF(kList, policy.get("fitForBackground"))

# Size of footprint, and mask around candidate variability
stampSize     = 200
maskSize      = 10

for line in open(infile).readlines():
    if line.startswith('#'):
        fields = line.split()
        ra     = float(fields[1])
        decl   = float(fields[2])
        continue
    
    fields = line.split()
    sciencePath      = os.path.basename(fields[0])
    templatePath     = os.path.basename(fields[1])
    templateExposure = afwImage.ExposureF(templatePath)
    scienceExposure  = afwImage.ExposureF(sciencePath)
    
    #diffimTools.backgroundSubtract(policy.getPolicy("afwBackgroundPolicy"),
    #                               [scienceExposure.getMaskedImage(),])
    templateExposure = ipDiffim.warpTemplateExposure(templateExposure,
                                                  scienceExposure,
                                                  policy.getPolicy("warpingPolicy"))

    xCenter, yCenter = map(int, scienceExposure.getWcs().skyToPixel(ra, decl))

    if ((xCenter < 0) or (yCenter < 0) or
        (xCenter > scienceExposure.getWidth()) or (yCenter > scienceExposure.getHeight())):
        pexLog.Trace("lsst.ip.diffim.runMaskedKernel", 1,
                     "%s : Candidate off image (%.3f, %.3f)" % (sciencePath, xCenter, yCenter))
        continue
    pexLog.Trace("lsst.ip.diffim.runMaskedKernel", 1,
                 "%s : Candidate at %.3f, %.3f" % (sciencePath, xCenter, yCenter))

    
    minX = max(xCenter - stampSize//2, 0)
    minY = max(yCenter - stampSize//2, 0)
    maxX = min(xCenter + stampSize//2, scienceExposure.getWidth())
    maxY = min(yCenter + stampSize//2, scienceExposure.getHeight())
    candBBox = afwGeom.Box2I(afwGeom.Point2I(minX, minY),
                             afwGeom.Point2I(maxX, maxY))
    
    minX = max(xCenter - maskSize//2, 0)
    minY = max(yCenter - maskSize//2, 0)
    maxX = min(xCenter + maskSize//2, scienceExposure.getWidth())
    maxY = min(yCenter + maskSize//2, scienceExposure.getHeight())

    templateImage = templateExposure.getMaskedImage()
    scienceImage  = scienceExposure.getMaskedImage()
    ds9.mtv(afwImage.MaskedImageF(templateImage, candBBox, afwImage.LOCAL), frame = 1)
    ds9.mtv(afwImage.MaskedImageF(scienceImage, candBBox, afwImage.LOCAL), frame = 2)

    # Create deep copy since this will be modified
    pixelMask  = afwImage.MaskU(scienceImage.getMask(), True)
    # Add in template mask
    pixelMask |= templateImage.getMask()
    # And mask out the variability!
    for y in range(minY, maxY + 1):
        for x in range(minX, maxX + 1):
            maskVal  = pixelMask.get(x,y)
            maskVal |= afwImage.MaskU_getPlaneBitMask("BAD")
            pixelMask.set(x, y, maskVal)

    # Sub image for speed
    templateSubimage = afwImage.MaskedImageF(templateImage, candBBox, afwImage.LOCAL)
    scienceSubimage  = afwImage.MaskedImageF(scienceImage, candBBox, afwImage.LOCAL)
    maskSubmask      = afwImage.MaskU(pixelMask, candBBox, afwImage.LOCAL)
            
    # Spread it here
    #fooMi  = afwImage.MaskedImageF(scienceSubimage.getImage(),
    #                               maskSubmask,
    #                               scienceSubimage.getVariance())
    #foocMi = afwImage.MaskedImageF(fooMi.getDimensions())
    #afwMath.convolve(foocMi, fooMi, kList[0], False)
    #finalMask = foocMi.getMask()
    finalMask = maskSubmask

    soln.build(templateSubimage.getImage(), scienceSubimage.getImage(),
               scienceSubimage.getVariance(), finalMask)
    try:
        soln.solve()
    except Exception, e:
        pexLog.Trace("lsst.ip.diffim.runMaskedKernel", 1,
                     "%s : Exception caught solving kernel" % (sciencePath))
        continue
    
    kernel = soln.getKernel()
    bg     = soln.getBackground()
    diffim = ipDiffim.convolveAndSubtract(templateSubimage, scienceSubimage,
                                          kernel, bg)
    diffim.writeFits(re.sub(".fits.gz", "_diff.fits", sciencePath))
    ds9.mtv(soln.makeKernelImage(), frame = 3)
    ds9.mtv(diffim, frame = 4)    
    
