#!/usr/bin/env python
import os
import optparse
import eups
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim

DefScienceExposureName = "medswarp1lanczos2.fits"
DefTemplateExposureName = "med_img.fits"
DefBorderWidth = 5
DefVerbosity = 0

def getArg(ind, defValue):
    if ind < len(args):
        return args[ind]
    return defValue

def printBBox(descr, bbox):
    print "%s: [%s, %s], [%s, %s]" % (descr, bbox.getX0(), bbox.getY0(), bbox.getX1(), bbox.getY1())

def computeBBox(scienceExposurePath, templateExposurePath, borderWidth, templateSubExposurePath=None):
    """Compute and print template bounding box
    
    Inputs:
    - scienceExposurePath: path to science image with Wcs info in header
    - templateExposurePath: path to template image with Wcs info in header
    - borderWidth: size of border for template bounding box (pixels)
    - templateSubExposure: path for output template SubExposure; not written if None
    
    All exposure paths must omit the final "_img.fits"
    """
    sciExp = afwImage.ExposureF(scienceExposurePath)
    if not sciExp.hasWcs():
        raise RuntimeError("Science exposure %r has no Wcs" % (scienceExposurePath,))
    sciWcs = sciExp.getWcs()
    sciDim = sciExp.getMaskedImage().getDimensions()
    sciBBox = afwImage.BBox(afwImage.PointI(0, 0), sciDim[0], sciDim[1])
    printBBox("Science BBox", sciBBox)
    
    # instead of reading in the template Exposure (which may be huge), just use the header
    tmpHdrData = afwImage.readMetadata(templateExposurePath + "_img.fits")
    tmpDim = afwImage.PointI(tmpHdrData.getInt("NAXIS1"), tmpHdrData.getInt("NAXIS2"))
    tmpWcs = afwImage.Wcs(tmpHdrData)
    tmpFullBBox = afwImage.BBox(afwImage.PointI(0, 0), tmpDim[0], tmpDim[1])
    printBBox("Template Full BBox", tmpFullBBox)
    
    borderWidth = 5
    tmpBBox = ipDiffim.computeTemplateBbox(sciWcs, sciDim, tmpWcs, tmpDim, borderWidth)
    printBBox("Template BBox", tmpBBox)

    if templateSubExposurePath:
        print "Reading sub template exposure"
        tmpSubExp = afwImage.ExposureF(templateExposurePath, 0, tmpBBox)

        print "Saving templateSubExposure to:", templateSubExposurePath
        tmpSubExp.writeFits(templateSubExposurePath)

if __name__ == "__main__":
    dataDir = eups.productDir("afwdata")
    if not dataDir:
        defScienceExposurePath = None
        defTemplateExposurePath = None
    else:
        defScienceExposurePath = os.path.join(dataDir, DefScienceExposureName)
        defTemplateExposurePath = os.path.join(dataDir, DefTemplateExposureName)

    usage = """usage: %%prog [options] scienceExposure templateExposure [templateSubExposure]

    Computes bounding box on a template exposure corresponding to a science exposure.
    If templateSubExposure is specified, saves the template sub-exposure as a FITS file.
    
    Warning: the two input FITS images must contain WCS information in their header

    Note:
    - image arguments are paths to Exposure fits files; they must not include the final "_img.fits"
      and must include WCS information
    - default scienceExposure = %s
    - default templateExposure = %s
    """ % (defScienceExposurePath, defTemplateExposurePath)
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-b", "--border",
                      type=int, default=DefBorderWidth,
                      help="width of border to be added to template exposure; default=%s" % (DefBorderWidth,))
    parser.add_option("-v", "--verbosity",
                      type=int, default=DefVerbosity,
                      help="verbosity of diagnostic trace messages; 1 for just warnings, more for more" + \
                      " information; default=%s" % (DefVerbosity,))
    
    (opt, args) = parser.parse_args()

    scienceExposurePath = getArg(0, defScienceExposurePath)
    templateExposurePath = getArg(1, defTemplateExposurePath)
    templateSubExposurePath = getArg(2, None)
    if None in (scienceExposurePath, templateExposurePath):
        print "Must setup afwdata or else provide image paths"
    borderWidth = opt.border

    print "Science  Exposure =", scienceExposurePath
    print "Template Exposure =", templateExposurePath
    print "Border Width =", borderWidth

    if opt.verbosity > 0:
        print "Verbosity =", opt.verbosity
        lsst.pex.logging.Trace_setVerbosity("lsst.ip.diffim", opt.verbosity)

    computeBBox(scienceExposurePath, templateExposurePath, borderWidth, templateSubExposurePath)
        
    