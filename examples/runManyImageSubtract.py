#!/usr/bin/env python
"""Subtract multiple of pairs of images as specified in a file.
"""
from __future__ import with_statement

import os
import sys
import optparse

import eups
import lsst.mwi.data as mwiData
import lsst.fw.Core.fwLib as fw
import lsst.mwi.utils
import lsst.imageproc

def main():
    imageProcDir = eups.productDir("imageproc", "setup")
    if imageProcDir == None:
        print "Error: imageproc not setup"
        sys.exit(1)
    moduleDir = os.path.dirname(os.path.abspath(__file__))

    defPolicyPath = os.path.join(imageProcDir, "pipeline", "ImageSubtractStageDictionary.paf")
    defFileList = os.path.join(moduleDir, "fileList.txt")
    defVerbosity = 0
    
    usage = """usage: %%prog [options] [fileList]

Notes:
- fileList is a list of image files to subtract; default=%r
  Each line consists of:
  sciencePath templatePath [differencePath]
  blank lines and lines beginning with # are ignored
- image paths must NOT include the final _img.fits
- environment variables (such as $FW_DATA) are expanded in each of the three paths
- the result is science image - template image
- the template image is convolved, the science image is not
- default difference image is <scienceName>_diff
  where <scienceName> is the name portion of sciencePath
- default --policy=%s
    """ % (defFileList, defPolicyPath)
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-p", "--policy", default=defPolicyPath, help="policy file")
    parser.add_option("-t", "--trial", action="store_true", default=False,
        help="trial run: show what images would be subtracted, but don't subtract")
    parser.add_option("-v", "--verbosity", type=int, default=defVerbosity,
        help="verbosity of diagnostic trace messages; 1 for just warnings, more for more information")
    (options, args) = parser.parse_args()
    
    def getArg(ind, defValue):
        if ind < len(args):
            return args[ind]
        return defValue

    fileListPath = os.path.abspath(getArg(0, defFileList))
    print "File list:", fileListPath

    if options.verbosity > 0:
        print "Verbosity =", options.verbosity
        lsst.mwi.utils.Trace_setVerbosity("lsst.imageproc", options.verbosity)

    policyPath = options.policy
    policy = lsst.mwi.policy.Policy.createPolicy(policyPath)

    with file(fileListPath, "rU") as fileList:    
        for lineNum, dataStr in enumerate(fileList):
            dataStr = dataStr.strip()
            if not dataStr:
                continue
            if dataStr.startswith("#"):
                continue
            dataList = dataStr.split()
            if len(dataList) < 2 or len(dataList) > 3:
                print "Cannot parse line %s: %r" % (lineNum, dataStr)
            sciencePath, templatePath = dataList[0:2]
            if len(dataList) > 2:
                differencePath = dataList[2]
            else:
                # use default = <science_image_name>_diff
                differencePath = "%s_diff" % (os.path.basename(sciencePath),)
            sciencePath = os.path.abspath(os.path.expandvars(sciencePath))
            templatePath = os.path.abspath(os.path.expandvars(templatePath))
            differencePath = os.path.abspath(os.path.expandvars(differencePath))
            print "Compute %r = \n  %r - %r" % (differencePath, sciencePath, templatePath)
            
            templateMaskedImage = fw.MaskedImageD()
            templateMaskedImage.readFits(templatePath)
            
            scienceMaskedImage  = fw.MaskedImageD()
            scienceMaskedImage.readFits(sciencePath)

            if not options.trial:
                differenceImage, psfMatchKernelPtr, backgroundFunctionPtr = lsst.imageproc.imageSubtract(
                    imageToConvolve = templateMaskedImage,
                    imageToNotConvolve = scienceMaskedImage,
                    policy = policy,
                )
                differenceImage.writeFits(differencePath)

if __name__ == "__main__":
    memId0 = mwiData.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if mwiData.Citizen_census(0, memId0) != 0:
        print mwiData.Citizen_census(0, memId0), "Objects leaked:"
        print mwiData.Citizen_census(mwiData.cout, memId0)
