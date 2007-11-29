#!/usr/bin/env python
import os
import sys
import optparse

import lsst.mwi.data as mwiData
import lsst.fw.Core.fwLib as fw
import lsst.imageproc
import lsst.mwi.utils

def runImageSubtract():
    defInDir = os.environ.get("FWDATA_DIR", "")
    moduleDir = os.path.split(__file__)[0]
    
    defSciencePath = os.path.join(defInDir, "871034p_1_MI")
    defTemplatePath = os.path.join(defInDir, "871034p_1_MI")
    defPolicyPath = os.path.join(moduleDir, "imageSubtractPipeline", "policy", "imageSubtract_policy.paf")
    defOutputPath = "diffImage"
    defVerbosity = 5 # change to 0 once this all works to hide all messages
    
    usage = """usage: %%prog [options] [scienceImage [templateImage [policyFile [outputImage]]]]
    Note:
    - image arguments are paths to MaskedImage fits files
    - image arguments must NOT include the final _img.fits
    - the result is science image - template image
    - the template image is convolved, the science image is not
    - default templateMaskedImage = %s
    - default scienceMaskedImage = %s
    - default policy = %s
    - default outputImage = %s 
    """ % (defSciencePath, defTemplatePath, defPolicyPath, defOutputPath)
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-v", "--verbosity",
                      type=int, default=defVerbosity,
                      help="verbosity of diagnostic trace messages; 1 for just warnings, more for more information")
    parser.add_option("-d", "--debugIO",
                      action="store_true", default=False,
                      help="write diagnostic intermediate files")
    (options, args) = parser.parse_args()
    
    def getArg(ind, defValue):
        if ind < len(args):
            return args[ind]
        return defValue
    
    templatePath = getArg(0, defTemplatePath)
    sciencePath = getArg(1, defSciencePath)
    policyPath = getArg(2, defPolicyPath)
    outputPath = getArg(3, defOutputPath)
    
    print "Template image:", templatePath
    print "Science image: ", sciencePath
    print "Policy file:   ", policyPath
    print "Output image:  ", outputPath
    
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
