#!/usr/bin/env python
"""Subtract one pair of images.
"""
import os
import sys
import optparse

import eups
import lsst.daf.base as dafBase
import lsst.afw.Core.afwLib as afw
import lsst.ip.diffim
import lsst.pex.logging

def main():
    defDataDir = os.environ.get("FWDATA_DIR", "")
    imageProcDir = eups.productDir("ip_diffim", "setup")
    if imageProcDir == None:
        print "Error: ip_diffim not setup"
        sys.exit(1)

    defSciencePath = os.path.join(defDataDir, "CFHT", "D4", "cal-53535-i-797722_1")
    defTemplatePath = os.path.join(defDataDir, "CFHT", "D4", "cal-53535-i-797722_1_tmpl")
    defPolicyPath = os.path.join(imageProcDir, "pipeline", "ImageSubtractStageDictionary.paf")
    defOutputPath = "diffImage"
    defVerbosity = 0
    
    usage = """usage: %%prog [options] [scienceImage [templateImage [outputImage]]]]

Notes:
- image arguments are paths to MaskedImage fits files
- image arguments must NOT include the final _img.fits
- the result is science image - template image
- the template image is convolved, the science image is not
- default scienceMaskedImage=%s
- default templateMaskedImage=%s
- default outputImage=%s 
- default --policy=%s
""" % (defSciencePath, defTemplatePath, defOutputPath, defPolicyPath)
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-p", "--policy", default=defPolicyPath, help="policy file")
    parser.add_option("-d", "--debugIO", action="store_true", default=False,
        help="write diagnostic intermediate files")
    parser.add_option("-s", "--switchConvolve", action="store_true", default=False,
        help="switch which image is convolved; still detect on template")
    parser.add_option("-v", "--verbosity", type=int, default=defVerbosity,
        help="verbosity of diagnostic trace messages; 1 for just warnings, more for more information")
    (options, args) = parser.parse_args()
    
    def getArg(ind, defValue):
        if ind < len(args):
            return args[ind]
        return defValue
    
    sciencePath = getArg(0, defSciencePath)
    templatePath = getArg(1, defTemplatePath)
    outputPath = getArg(2, defOutputPath)
    policyPath = options.policy
    
    print "Science image: ", sciencePath
    print "Template image:", templatePath
    print "Output image:  ", outputPath
    print "Policy file:   ", policyPath
    
    templateMaskedImage = afw.MaskedImageF()
    templateMaskedImage.readFits(templatePath)
    
    scienceMaskedImage  = afw.MaskedImageF()
    scienceMaskedImage.readFits(sciencePath)
    
    policy = lsst.pex.policy.Policy.createPolicy(policyPath)
    if options.debugIO:
        policy.set("debugIO", True)
    if options.switchConvolve:
        policy.set("switchConvolve", True)
    
    if options.verbosity > 0:
        print "Verbosity =", options.verbosity
        lsst.pex.logging.Trace_setVerbosity("lsst.ip.diffim", options.verbosity)
    
    # compute difference image
    differenceImage, psfMatchKernelPtr, backgroundFunctionPtr = lsst.ip.diffim.imageSubtract(
        imageToConvolve = templateMaskedImage,
        imageToNotConvolve = scienceMaskedImage,
        policy = policy,
    )
    differenceImage.writeFits(outputPath)

def run():
    from lsst.pex.logging import Log
    Log.getDefaultLog()                 # leaks a DataProperty
    
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), "Objects leaked:"
        print dafBase.Citizen_census(dafBase.cout, memId0)

if __name__ == "__main__":
    run()
