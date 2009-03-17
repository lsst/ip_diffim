#!/usr/bin/env python
import os
import optparse
import eups
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim

DefScienceImageName = "medswarp1lanczos2.fits"
DefTemplateImageName = "med_img.fits"
DefBorderWidth = 5
DefVerbosity = 0

def getArg(ind, defValue):
    if ind < len(args):
        return args[ind]
    return defValue

def printBBox(descr, bbox):
    print "%s: [%s, %s], [%s, %s]" % (descr, bbox.getX0(), bbox.getY0(), bbox.getX1(), bbox.getY1())

def computeBBox(scienceImagePath, templateImagePath, borderWidth):
    """Compute and print template bounding box
    
    Inputs:
    - scienceImagePath: path to science image with Wcs info in header
    - templateImagePath: path to template image with Wcs info in header
    - borderWidth: size of border for template bounding box (pixels)
    """
    sciHdr = dafBase.PropertySet()
    sciIm = afwImage.ImageF(scienceImagePath, 0, sciHdr)
    sciWcs = afwImage.Wcs(sciHdr)
    sciDim = sciIm.getDimensions()
    sciBBox = afwImage.BBox(afwImage.PointI(0, 0), sciDim[0], sciDim[1])
    printBBox("Science  BBox", sciBBox)
    
    tmpHdrData = dafBase.PropertySet()
    tmpIm = afwImage.ImageF(templateImagePath, 0, tmpHdrData)
    tmpWcs = afwImage.Wcs(tmpHdrData)
    tmpDim = tmpIm.getDimensions()
    
    borderWidth = 5
    tmpBBox = ipDiffim.computeTemplateBbox(sciWcs, sciDim, tmpWcs, tmpDim, borderWidth)
    printBBox("Template BBox", tmpBBox)

if __name__ == "__main__":
    dataDir = eups.productDir("afwdata")
    if not dataDir:
        defScienceImagePath = None
        defTemplateImagePath = None
    else:
        defScienceImagePath = os.path.join(dataDir, DefScienceImageName)
        defTemplateImagePath = os.path.join(dataDir, DefTemplateImageName)

    usage = """usage: %%prog [options] [scienceImage templateImage

    Computes bounding box on a template image corresponding to a science image
    
    Warning: the two FITS images must contain WCS information in their header

    Note:
    - image arguments are paths to Image fits files; they must include WCS information
    - default scienceImage = %s
    - default templateImage = %s
    """ % (defScienceImagePath, defTemplateImagePath)
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-b", "--border",
                      type=int, default=DefBorderWidth,
                      help="width of border to be added to template exposure; default=%s" % (DefBorderWidth,))
    parser.add_option("-v", "--verbosity",
                      type=int, default=DefVerbosity,
                      help="verbosity of diagnostic trace messages; 1 for just warnings, more for more" + \
                      " information; default=%s" % (DefVerbosity,))
    
    (opt, args) = parser.parse_args()

    scienceImagePath = getArg(0, defScienceImagePath)
    templateImagePath = getArg(1, defTemplateImagePath)
    if None in (scienceImagePath, templateImagePath):
        print "Must setup afwdata or else provide image paths"
    borderWidth = opt.border

    print "Science  Image =", scienceImagePath
    print "Template Image =", templateImagePath
    print "Border Width =", borderWidth

    if opt.verbosity > 0:
        print "Verbosity =", opt.verbosity
        lsst.pex.logging.Trace_setVerbosity("lsst.ip.diffim", opt.verbosity)

    computeBBox(scienceImagePath, templateImagePath, borderWidth)