#!/usr/bin/env python
"""Run the image copy pipeline to copy one image
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

defInDir = os.environ.get("FWDATA_DIR", "")
moduleDir = os.path.split(__file__)[0]

defInputPath = os.path.join(defInDir, "871034p_1_MI")
defOutputPath = "imageCopy"

ConfigRelPath = os.path.join(moduleDir, "imageCopyPipeline")

usage = """usage: %%prog [inputImage [outputImage]]
Note:
- image arguments are paths to MaskedImage fits files
- image arguments must NOT include the final _img.fits
- default inputMaskedImage = %s
- default outputImage = %s 
""" % (defInputPath, defOutputPath)

parser = optparse.OptionParser(usage)
(options, args) = parser.parse_args()

def getArg(ind, defValue):
    if ind < len(args):
        return args[ind]
    return defValue

inputImagePath = os.path.abspath(getArg(0, defInputPath))
outputImagePath = os.path.abspath(getArg(3, defOutputPath))

def copyTemplatedConfigFile(templateName, templateDict):
    """Read a templated configuration file, fill it in and write it out.
    templateName is a path relative to ConfigRelPath
    templateDict contains the items to substitute in the template file
    """
    #print "copyTemplatedConfigFile(%r, %r)" % (templateName, templateDict)
    templateBaseName, templateExt = os.path.splitext(templateName)
    if not templateBaseName.endswith("_template"):
        raise RuntimeError("Invalid templateName %r; must end with _template" % (templateName,))
    inPath = os.path.join(ConfigRelPath, templateName)
    with file(inPath, "rU") as templateFile:
        templateText = templateFile.read()
        destText = templateText % templateDict
    outName = templateBaseName[:-len("_template")] + templateExt
    outPath = os.path.join(ConfigRelPath, outName)
    with file(outPath, "w") as destFile:
        destFile.write(destText)

# write configuration files, filling in inputs as required
copyTemplatedConfigFile(
    "nodelist_template.scr",
    dict(
        ipaddress = socket.gethostname(),
    ),
)
copyTemplatedConfigFile(
    os.path.join("policy", "input_policy_template.paf"),
    dict(
        inputImagePath = inputImagePath,
    ),
)
copyTemplatedConfigFile(
    os.path.join("policy", "output_policy_template.paf"),
    dict(
        outputImagePath = outputImagePath,
    ),
)

configAbsPath = os.path.abspath(ConfigRelPath)
subprocess.call(os.path.join(configAbsPath, "run.sh"), cwd=configAbsPath)

# check for memory leaks
memId0 = 0
if mwiData.Citizen_census(0, memId0) != 0:
    print mwiData.Citizen_census(0, memId0), "Objects leaked:"
    print mwiData.Citizen_census(mwiData.cout, memId0)
