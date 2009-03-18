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

def computeBBox(scienceImagePath, templateImagePath, borderWidth, templateSubImagePath=None):
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
    
    tmpHdrData = afwImage.readMetadata(templateImagePath)
    tmpDim = afwImage.PointI(tmpHdrData.getInt("NAXIS1"), tmpHdrData.getInt("NAXIS2"))
    tmpWcs = afwImage.Wcs(tmpHdrData)
    
    borderWidth = 5
    tmpBBox = ipDiffim.computeTemplateBbox(sciWcs, sciDim, tmpWcs, tmpDim, borderWidth)
    printBBox("Template BBox", tmpBBox)

    print "Reading sub template exposure"
    tmpSubHdr = dafBase.PropertySet()
    tmpSubIm = afwImage.ImageF(templateImagePath, 0, tmpSubHdr, tmpBBox)
    
    if templateSubImagePath:
        print "Saving templateSubImage to:", templateSubImagePath
        tmpSubIm.writeFits(templateSubImagePath, tmpSubHdr)

if __name__ == "__main__":
    dataDir = eups.productDir("afwdata")
    if not dataDir:
        defScienceImagePath = None
        defTemplateImagePath = None
    else:
        defScienceImagePath = os.path.join(dataDir, DefScienceImageName)
        defTemplateImagePath = os.path.join(dataDir, DefTemplateImageName)

    usage = """usage: %%prog [options] scienceImage templateImage [templateSubImage]

    Computes bounding box on a template image corresponding to a science image.
    If templateSubImage is specified, saves the corresponding template subimage as a FITS file.
    
    Warning: the two input FITS images must contain WCS information in their header

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
    templateSubImagePath = getArg(2, None)
    if None in (scienceImagePath, templateImagePath):
        print "Must setup afwdata or else provide image paths"
    borderWidth = opt.border

    print "Science  Image =", scienceImagePath
    print "Template Image =", templateImagePath
    print "Border Width =", borderWidth

    if opt.verbosity > 0:
        print "Verbosity =", opt.verbosity
        lsst.pex.logging.Trace_setVerbosity("lsst.ip.diffim", opt.verbosity)

    computeBBox(scienceImagePath, templateImagePath, borderWidth, templateSubImagePath)
        
    