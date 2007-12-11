#!/usr/bin/env python
"""Run the image subtraction pipeline two subtract two images
"""
from __future__ import with_statement

import os
import sys
import optparse
import shutil
import subprocess
import socket

import lsst.mwi.data as mwiData
import lsst.mwi.utils
import startPipeline

def main():
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
    
    usage = """usage: %%prog [options] [scienceImage [templateImage [outputImage]]]

Notes:
- image arguments are paths to MaskedImage fits files
- image arguments must NOT include the final _img.fits
- the result is science image - template image
- the template image is convolved, the science image is not
- default scienceMaskedImage=%s
- default templateMaskedImage=%s
- default outputImage=%s 
- default --policy=%s""" % (defSciencePath, defTemplatePath, defOutputPath, defPolicyPath)
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-p", "--policy", default=defPolicyPath, help="policy file")
    parser.add_option("-v", "--verbosity", type=int, default=defVerbosity,
        help="verbosity of diagnostic trace messages; default=%s" % (defVerbosity,))
    (options, args) = parser.parse_args()
    
    def getArg(ind, defValue):
        if ind < len(args):
            return args[ind]
        return defValue
    
    sciencePath = os.path.abspath(getArg(0, defSciencePath))
    templatePath = os.path.abspath(getArg(1, defTemplatePath))
    outputPath = os.path.abspath(getArg(2, defOutputPath))
    policyPath = options.policy
    
    print "Science image: ", sciencePath
    print "Template image:", templatePath
    print "Output image:  ", outputPath
    print "Policy file:   ", policyPath
    
    def copyTemplatedConfigFile(templateName, templateDict):
        """Read a templated configuration file, fill it in and write it out.
        templateName is a path relative to pipelineDir
        templateDict contains the items to substitute in the template file
        """
        #print "copyTemplatedConfigFile(%r, %r)" % (templateName, templateDict)
        templateBaseName, templateExt = os.path.splitext(templateName)
        if not templateBaseName.endswith("_template"):
            raise RuntimeError("Invalid templateName %r; must end with _template" % (templateName,))
        inPath = os.path.join(pipelineDir, templateName)
        with file(inPath, "rU") as templateFile:
            templateText = templateFile.read()
            destText = templateText % templateDict
        outName = templateBaseName[:-len("_template")] + templateExt
        outPath = os.path.join(pipelineDir, outName)
        with file(outPath, "w") as destFile:
            destFile.write(destText)
    
    # write configuration files, filling in templates as required
    copyTemplatedConfigFile(
        "nodelist_template.scr",
        dict(
            ipaddress = socket.gethostname(),
        ),
    )
    copyTemplatedConfigFile(
        "input_policy_template.paf",
        dict(
            scienceExposurePath = sciencePath,
            templateExposurePath = templatePath,
        ),
    )
    copyTemplatedConfigFile(
        "output_policy_template.paf",
        dict(
            differenceExposurePath = outputPath,
        ),
    )
    copyTemplatedConfigFile(
        "pipeline_policy_template.paf",
        dict(
            imageSubtractPolicyPath = policyPath,
        ),
    )
    
    if options.verbosity > 0:
        print "Verbosity =", options.verbosity
        lsst.mwi.utils.Trace_setVerbosity("lsst.imageproc", options.verbosity)
    lsst.mwi.utils.Trace_setVerbosity("dps", 3)
    
    nodeList = os.path.join(pipelineDir, "nodelist.scr")
    startPipeline.startPipeline(nodeList, "pipeline_policy.paf", "imageSubtractId")
    # or if you prefer to use a private copy of run.sh...
    #subprocess.call(os.path.join(pipelineDir, "run.sh"), cwd=pipelineDir)


if __name__ == "__main__":
    main()
    # check for memory leaks
    memId0 = 0
    if mwiData.Citizen_census(0, memId0) != 0:
        print mwiData.Citizen_census(0, memId0), "Objects leaked:"
        print mwiData.Citizen_census(mwiData.cout, memId0)
