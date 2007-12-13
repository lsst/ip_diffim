#!/usr/bin/env python
import os
import sys
import optparse

import eups
import lsst.mwi.data as mwiData
import lsst.fw.Core.fwLib as fw
import lsst.imageproc.imageprocLib
import lsst.mwi.utils

def main():
    defDataDir = os.environ.get("FWDATA_DIR", "")
    
    defSciencePath = os.path.join(defDataDir, "871034p_1_MI")
    defTemplatePath = os.path.join(defDataDir, "871034p_1_MI")
    defOutputPath = "remappedImage"
    defKernelType = "lanczos"
    defKernelSize = 5
    defVerbosity = 5 # change to 0 once this all works to hide all messages
    
    usage = """usage: %%prog [options] [originalImage [remapImage [outputImage]]]
    Note:
    - image arguments are paths to MaskedImage fits files;
      they must NOT include the final _img.fits|_var.fits|_msk.fits
    - output image is the original image remapped to match the remap image's WCS and size
    - default originalImage = %s
    - default remapImage = %s
    - default outputImage = %s
    """ % (defSciencePath, defTemplatePath, defOutputPath)
    
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--kernelSize",
                      type=int, default=defKernelSize,
                      help="kernel size (cols=rows); should be even; default=%s" % (defKernelSize,))
    parser.add_option("", "--kernelType",
                      type=str, default=defKernelType,
                      help="kernel type; default=%s" % (defKernelType,))
    parser.add_option("-v", "--verbosity",
                      type=int, default=defVerbosity,
                      help="verbosity of diagnostic trace messages; 1 for just warnings, more for more" + \
                      " information; default=%s" % (defVerbosity,))
    
    (opt, args) = parser.parse_args()
    
    def getArg(ind, defValue):
        if ind < len(args):
            return args[ind]
        return defValue
    
    originalPath = getArg(0, defTemplatePath)
    remapPath = getArg(1, defSciencePath)
    outputPath = getArg(2, defOutputPath)
    print "Remapping masked image  ", originalPath
    print "to match wcs and size of", remapPath
    
    originalExposure = fw.ExposureD()
    originalExposure.readFits(originalPath)
    
    remapExposure  = fw.ExposureD()
    remapExposure.readFits(remapPath)
    
    if opt.verbosity > 0:
        print "Verbosity =", opt.verbosity
        lsst.mwi.utils.Trace_setVerbosity("lsst.imageproc", opt.verbosity)
    
    numEdgePixels = lsst.imageproc.imageprocLib.wcsMatch(
        remapExposure, originalExposure, opt.kernelType, opt.kernelSize, opt.kernelSize)
    print "Remapped masked image has %s edge pixels" % (numEdgePixels)
    
    print "Writing remapped masked image to %s" % (outputPath,)
    remapExposure.getMaskedImage().writeFits(outputPath)

if __name__ == "__main__":
    memId0 = mwiData.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if mwiData.Citizen_census(0, memId0) != 0:
        print mwiData.Citizen_census(0, memId0), "Objects leaked:"
        print mwiData.Citizen_census(mwiData.cout, memId0)
