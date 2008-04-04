#!/usr/bin/env python
"""Run the image subtraction pipeline to perform image subtraction on a series of images.
"""
from __future__ import with_statement

import os
import sys
import optparse
import subprocess
import socket

import eups
import lsst.daf.base as dafBase
import lsst.pex.logging
import lsst.pex.harness.startPipeline

def main():
    imageProcDir = eups.productDir("ip_diffim", "setup")
    if imageProcDir == None:
        print "Error: ip_diffim not setup"
        sys.exit(1)
    pipelineDir = os.path.dirname(os.path.abspath(__file__))

    defPolicyPath = os.path.join(imageProcDir, "pipeline", "ImageSubtractStageDictionary.paf")
    defVerbosity = 0
    
    usage = """usage: %%prog [options]

Notes:
- default --policy=%s""" % (defPolicyPath,)
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-p", "--policy", default=defPolicyPath, help="image subtraction policy file")
    parser.add_option("-v", "--verbosity", type=int, default=defVerbosity,
        help="verbosity of diagnostic trace messages; default=%s" % (defVerbosity,))
    (options, args) = parser.parse_args()
    
    subtractionPolicyPath = options.policy

    print "Image Subtraction Policy file:", subtractionPolicyPath
    
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
        "pipeline_policy_template.paf",
        dict(
            imageSubtractPolicyPath = subtractionPolicyPath,
        ),
    )
    
    if options.verbosity > 0:
        print "Verbosity =", options.verbosity
        lsst.pex.logging.Trace_setVerbosity("lsst.ip.diffim", options.verbosity)
    lsst.pex.logging.Trace_setVerbosity("pex.harness", 3)
    
    print """Starting the pipeline.
Once you see a message like:
  Python Slice handleEvents rank :  0  - waiting on receive...
then run feedManySubtractPipeline.py from a new process
to feed images to the image subtraction pipeline.

Control-C the pipeline when it is done (or you have had enough).
"""
    nodeList = os.path.join(pipelineDir, "nodelist.scr")
    lsst.pex.harness.startPipeline.startPipeline(nodeList, "pipeline_policy.paf", "imageManySubtractId")

if __name__ == "__main__":
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), "Objects leaked:"
        print dafBase.Citizen_census(dafBase.cout, memId0)
