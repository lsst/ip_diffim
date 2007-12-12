#!/usr/bin/env python
"""Run the image subtraction and detection pipeline on a series of images
"""
from __future__ import with_statement

import os
import sys
import optparse
import subprocess
import socket

import lsst.mwi.data as mwiData
import lsst.mwi.utils

def main():
    try:
        imageProcDir = os.environ["IMAGEPROC_DIR"]
    except KeyError:
        print "Error: imageproc not setup"
        sys.exit(1)
    try:
        detectionDir = os.environ["DETECTION_DIR"]
    except KeyError:
        print "Error: detection not setup"
        sys.exit(1)
    pipelineDir = os.path.dirname(os.path.abspath(__file__))
    sys.path += [os.path.dirname(pipelineDir)]
    import startPipeline

    defDetectionPolicyPath = os.path.join(detectionDir, "pipeline", "DetectionStagePolicy.paf")
    defVerbosity = 0
    
    usage = """usage: %%prog [options]

Notes:
- default --detpolicy=%s""" % (defDetectionPolicyPath,)
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-d", "--detpolicy", default=defDetectionPolicyPath, help="detection policy file")
    parser.add_option("-v", "--verbosity", type=int, default=defVerbosity,
        help="verbosity of diagnostic trace messages; default=%s" % (defVerbosity,))
    (options, args) = parser.parse_args()
    
    detectionPolicyPath = options.detpolicy

    print "Detection Policy file:", detectionPolicyPath
    
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
            detectionPolicyPath = detectionPolicyPath,
        ),
    )
    
    if options.verbosity > 0:
        print "Verbosity =", options.verbosity
        lsst.mwi.utils.Trace_setVerbosity("lsst.imageproc", options.verbosity)
    lsst.mwi.utils.Trace_setVerbosity("dps", 3)
    
    print """Starting the pipeline.
Once you see a message like:
  Python Slice handleEvents rank :  0  - waiting on receive...
then from a new shell run one of:
    eventgenerator.py < dataFile
or (for one sample event):
    triggervisitevent.py
to feed images to the image subtraction pipeline.

Control-C the pipeline when it is done (or you have had enough).
"""
    nodeList = os.path.join(pipelineDir, "nodelist.scr")
    startPipeline.startPipeline(nodeList, "pipeline_policy.paf", "RUN0001")

if __name__ == "__main__":
    main()
    # check for memory leaks
    memId0 = 0
    if mwiData.Citizen_census(0, memId0) != 0:
        print mwiData.Citizen_census(0, memId0), "Objects leaked:"
        print mwiData.Citizen_census(mwiData.cout, memId0)
