#!/usr/bin/env python
"""Feed the image subtraction pipeline with a series of images.

To do:
- My attempts to create an eventTransmitter fails with "Connection refused"
  (but this is printed to stdout; there's no obvious way to get that info in my program!)
- Try creating an event by supplying a policy containing:
  topicName = "triggerAssociationEvent"
  useLocalSockets = True
"""
from __future__ import with_statement

import os
import sys
import optparse
import socket
import time

import lsst.mwi.data as mwiData
import lsst.mwi.policy
import lsst.mwi.utils
import lsst.events

def sendEvent(templatePath, sciencePath, differencePath, eventTransmitter, trialRun=False):
    rootProperty = mwiData.SupportFactory.createPropertyNode("root");

    rootProperty.addProperty(mwiData.DataProperty("visitId", 1)) # this may be required
    rootProperty.addProperty(mwiData.DataProperty("sciencePath", sciencePath))
    rootProperty.addProperty(mwiData.DataProperty("templatePath", templatePath))
    rootProperty.addProperty(mwiData.DataProperty("differencePath", differencePath))

    if not trialRun:
        eventTransmitter.publish("imageSubtractEventType", rootProperty)


def main():
    moduleDir = os.path.split(__file__)[0]
    defFileList = os.path.join(moduleDir, "fileList.txt")
    
    defVerbosity = 5 # change to 0 once this all works to hide all messages
    
    usage = """usage: %%prog [options] [fileList]
    Note:
    - fileList is a list of image files to subtract; default=%r
      Each line consists of:
      sciencePath templatePath [differencePath]
      blank lines and lines beginning with # are ignored
    - image paths must NOT include the final _img.fits
    - environment variables (such as $FW_DATA) are expanded in each of the three paths
    - the result is science image - template image
    - the template image is convolved, the science image is not
    - default difference image is <scienceName>_diff
      where <scienceName> is the name portion of sciencePath
    """ % (defFileList,)
    
    parser = optparse.OptionParser(usage)
    parser.add_option("-t", "--trial",
                      action="store_true", default=False,
                      help="trial run: show what would be done but don't run the pipeline")
    (options, args) = parser.parse_args()
    
    def getArg(ind, defValue):
        if ind < len(args):
            return args[ind]
        return defValue

    policyRelPath = os.path.join(moduleDir, "imageManySubtractPipeline", "policy")
    eventPolicy = lsst.mwi.policy.Policy(os.path.join(policyRelPath, "event_policy.paf"))
    eventTransmitter = lsst.events.EventTransmitter(eventPolicy)

    fileListPath = os.path.abspath(getArg(0, defFileList))
    print "File list:", fileListPath

    with file(fileListPath, "rU") as fileList:    
        for lineNum, dataStr in enumerate(fileList):
            dataStr = dataStr.strip()
            if not dataStr:
                continue
            if dataStr.startswith("#"):
                continue
            dataList = dataStr.split()
            if len(dataList) < 2 or len(dataList) > 3:
                print "Cannot parse line %s: %r" % (lineNum, dataStr)
            sciencePath, templatePath = dataList[0:2]
            if len(dataList) > 2:
                differencePath = dataList[2]
            else:
                differencePath = "%s_diff" % (os.path.split(sciencePath)[1],)
            sciencePath = os.path.abspath(os.path.expandvars(sciencePath))
            templatePath = os.path.abspath(os.path.expandvars(templatePath))
            differencePath = os.path.abspath(os.path.expandvars(differencePath))
            print "Computing %r = \n  %r - %r" % (differencePath, sciencePath, templatePath)
            sendEvent(templatePath, sciencePath, differencePath, eventTransmitter, options.trial)

if __name__ == "__main__":
    main()
    # check for memory leaks
    memId0 = 0
    if mwiData.Citizen_census(0, memId0) != 0:
        print mwiData.Citizen_census(0, memId0), "Objects leaked:"
        print mwiData.Citizen_census(mwiData.cout, memId0)
