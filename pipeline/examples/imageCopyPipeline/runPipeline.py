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

def main():
    defInDir = os.environ.get("FWDATA_DIR", "")
    try:
        imageProcDir = os.environ["IMAGEPROC_DIR"]
    except KeyError:
        print "Error: imageproc not setup"
        sys.exit(1)
    pipelineDir = os.path.dirname(os.path.abspath(__file__))
    sys.path += [os.path.dirname(pipelineDir)]
    import startPipeline
    
    defInputPath = os.path.join(defInDir, "871034p_1_MI")
    defOutputPath = "imageCopy"
    
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
    
    # write configuration files, filling in inputs as required
    copyTemplatedConfigFile(
        "nodelist_template.scr",
        dict(
            ipaddress = socket.gethostname(),
        ),
    )
    copyTemplatedConfigFile(
        "input_policy_template.paf",
        dict(
            inputImagePath = inputImagePath,
        ),
    )
    copyTemplatedConfigFile(
        "output_policy_template.paf",
        dict(
            outputImagePath = outputImagePath,
        ),
    )

    lsst.mwi.utils.Trace_setVerbosity("dps", 3)
    
    nodeList = os.path.join(pipelineDir, "nodelist.scr")
    startPipeline.startPipeline(nodeList, "pipeline_policy.paf", "copyPipelineId")

if __name__ == "__main__":
    memId0 = mwiData.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if mwiData.Citizen_census(0, memId0) != 0:
        print mwiData.Citizen_census(0, memId0), "Objects leaked:"
        print mwiData.Citizen_census(mwiData.cout, memId0)
