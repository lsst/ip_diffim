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

defInDir = os.environ.get("FWDATA_DIR", "")
try:
    imageProcDir = os.environ["IMAGEPROC_DIR"]
except KeyError:
    print "Error: imageproc not setup"
    sys.exit(1)
pipelineDir = os.path.join(imageProcDir, "pipeline", "examples", "imageSubtractPipeline")

defSciencePath = os.path.join(defInDir, "871034p_1_MI")
defTemplatePath = os.path.join(defInDir, "871034p_1_MI")
defPolicyPath = os.path.join(imageProcDir, "pipeline", "ImageSubtractStageDictionary.paf")
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
                  help="verbosity of diagnostic trace messages; 9 for just warnings, less for less information")
(options, args) = parser.parse_args()

def getArg(ind, defValue):
    if ind < len(args):
        return args[ind]
    return defValue

templatePath = os.path.abspath(getArg(0, defTemplatePath))
sciencePath = os.path.abspath(getArg(1, defSciencePath))
policyPath = getArg(2, defPolicyPath)
outputPath = os.path.abspath(getArg(3, defOutputPath))

print "Template image:", templatePath
print "Science image: ", sciencePath
print "Policy file:   ", policyPath
print "Output image:  ", outputPath

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
    os.path.join("policy", "input_policy_template.paf"),
    dict(
        scienceExposurePath = sciencePath,
        templateExposurePath = templatePath,
    ),
)
copyTemplatedConfigFile(
    os.path.join("policy", "output_policy_template.paf"),
    dict(
        differenceExposurePath = outputPath,
    ),
)
copyTemplatedConfigFile(
    os.path.join("policy", "pipeline_policy_template.json"),
    dict(
        imageSubtractPolicyPath = policyPath,
    ),
)

if options.verbosity > 0:
    print "Verbosity =", options.verbosity
    lsst.mwi.utils.Trace_setVerbosity("lsst.imageproc", options.verbosity)
    lsst.mwi.utils.Trace_setVerbosity("dps", options.verbosity)

nodeList = os.path.join(pipelineDir, "nodelist.scr")
startPipeline.startPipeline(nodeList)
# or if you prefer to use a private copy of run.sh...
#subprocess.call(os.path.join(pipelineDir, "run.sh"), cwd=pipelineDir)

# check for memory leaks
memId0 = 0
if mwiData.Citizen_census(0, memId0) != 0:
    print mwiData.Citizen_census(0, memId0), "Objects leaked:"
    print mwiData.Citizen_census(mwiData.cout, memId0)
