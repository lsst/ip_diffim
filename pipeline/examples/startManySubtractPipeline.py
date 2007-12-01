#!/usr/bin/env python
"""Run the image subtraction pipeline to perform image subtraction on a series of images.
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
    moduleDir = os.path.split(__file__)[0]
    configRelPath = os.path.join(moduleDir, "imageManySubtractPipeline")
    
    defPolicyPath = os.path.join("policy", "imageSubtract_policy.paf")
    defVerbosity = 5 # change to 0 once this all works to hide all messages
    
    usage = """usage: %%prog [options] [policyFile]
    Note:
    - default policy = %s
    """ % (defPolicyPath,)
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-v", "--verbosity",
                      type=int, default=defVerbosity,
                      help="verbosity of diagnostic trace messages; 9 for just warnings, less for less information")
    (options, args) = parser.parse_args()
    
    def getArg(ind, defValue):
        if ind < len(args):
            return args[ind]
        return defValue
    
    policyPath = getArg(1, defPolicyPath)

    print "Policy file:   ", policyPath
    
    def copyTemplatedConfigFile(templateName, templateDict):
        """Read a templated configuration file, fill it in and write it out.
        templateName is a path relative to configRelPath
        templateDict contains the items to substitute in the template file
        """
        #print "copyTemplatedConfigFile(%r, %r)" % (templateName, templateDict)
        templateBaseName, templateExt = os.path.splitext(templateName)
        if not templateBaseName.endswith("_template"):
            raise RuntimeError("Invalid templateName %r; must end with _template" % (templateName,))
        inPath = os.path.join(configRelPath, templateName)
        with file(inPath, "rU") as templateFile:
            templateText = templateFile.read()
            destText = templateText % templateDict
        outName = templateBaseName[:-len("_template")] + templateExt
        outPath = os.path.join(configRelPath, outName)
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
        os.path.join("policy", "pipeline_policy_template.json"),
        dict(
            imageSubtractPolicyPath = policyPath,
        ),
    )
    
    if options.verbosity > 0:
        print "Verbosity =", options.verbosity
        lsst.mwi.utils.Trace_setVerbosity("lsst.imageproc", options.verbosity)
    
    print """Starting the pipeline.
Once it is running feed it images to subtract by running feedManySubtractPipeline.py
from a new process.
Type control-c to kill the pipeline when you are finished.
"""
    configAbsPath = os.path.abspath(configRelPath)
    subprocess.call(os.path.join(configAbsPath, "run.sh"), cwd=configAbsPath)

if __name__ == "__main__":
    main()
    # check for memory leaks
    memId0 = 0
    if mwiData.Citizen_census(0, memId0) != 0:
        print mwiData.Citizen_census(0, memId0), "Objects leaked:"
        print mwiData.Citizen_census(mwiData.cout, memId0)
