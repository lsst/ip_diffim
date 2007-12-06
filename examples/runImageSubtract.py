#!/usr/bin/env python
import os
import sys
import optparse

import lsst.mwi.data as mwiData
import lsst.fw.Core.fwLib as fw
import lsst.imageproc
import lsst.mwi.utils

def runImageSubtract():
    defDataDir = os.environ.get("FWDATA_DIR", "")
    try:
        imageProcDir = os.environ["IMAGEPROC_DIR"]
    except KeyError:
        print "Error: imageproc not setup"
        sys.exit(1)
    pipelineDir = os.path.join(imageProcDir, "pipeline", "examples", "imageSubtractPipeline")

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
    
    templateMaskedImage = fw.MaskedImageD()
    templateMaskedImage.readFits(templatePath)
    
    scienceMaskedImage  = fw.MaskedImageD()
    scienceMaskedImage.readFits(sciencePath)
    
    policy = lsst.mwi.policy.Policy.createPolicy(policyPath)
    if options.debugIO:
        policy.set("debugIO", True)
    
    if options.verbosity > 0:
        print "Verbosity =", options.verbosity
        lsst.mwi.utils.Trace_setVerbosity("lsst.imageproc", options.verbosity)
    
    # compute difference image
    differenceImage, psfMatchKernelPtr, backgroundFunctionPtr = lsst.imageproc.imageSubtract(
        templateMaskedImage, scienceMaskedImage, policy)
    differenceImage.writeFits(outputPath)

if __name__ == "__main__":
    runImageSubtract()

    # check for memory leaks
    memId0 = 0
    if mwiData.Citizen_census(0, memId0) != 0:
        print mwiData.Citizen_census(0, memId0), "Objects leaked:"
        print mwiData.Citizen_census(mwiData.cout, memId0)
