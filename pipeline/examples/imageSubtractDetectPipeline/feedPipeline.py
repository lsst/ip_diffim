#!/usr/bin/env python
"""Feed the image subtraction pipeline with a series of images.

To do:
- Once pipelines can receive events using a local socket modify this example to work that way.
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

EventHost = "lsst8.ncsa.uiuc.edu"

def sendEvent(templatePath, sciencePath, differencePath, eventTransmitter):
    rootProperty = mwiData.SupportFactory.createPropertyNode("root");

    rootProperty.addProperty(mwiData.DataProperty("visitId", 1)) # this may be required
    rootProperty.addProperty(mwiData.DataProperty("sciencePath", sciencePath))
    rootProperty.addProperty(mwiData.DataProperty("templatePath", templatePath))
    rootProperty.addProperty(mwiData.DataProperty("differencePath", differencePath))

    eventTransmitter.publish("imageSubtractEventType", rootProperty)


def main():
    try:
        imageProcDir = os.environ["IMAGEPROC_DIR"]
    except KeyError:
        print "Error: imageproc not setup"
        sys.exit(1)
    pipelineDir = os.path.join(imageProcDir, "pipeline", "examples", "imageManySubtractPipeline")
    defFileList = os.path.join(pipelineDir, "fileList.txt")
    
    defVerbosity = 0
    
    usage = """usage: %%prog [options] [fileList]

Notes:
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
    parser.add_option("-t", "--trial", action="store_true", default=False,
        help="trial run: show what images would be subtracted, but don't run the pipeline")
    (options, args) = parser.parse_args()
    
    def getArg(ind, defValue):
        if ind < len(args):
            return args[ind]
        return defValue

    triggerEventTransmitter = lsst.events.EventTransmitter(EventHost, "triggerImageSubtraction")

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
                # use default = <science_image_name>_diff
                differencePath = "%s_diff" % (os.path.basename(sciencePath),)
            sciencePath = os.path.abspath(os.path.expandvars(sciencePath))
            templatePath = os.path.abspath(os.path.expandvars(templatePath))
            differencePath = os.path.abspath(os.path.expandvars(differencePath))
            print "Compute %r = \n  %r - %r" % (differencePath, sciencePath, templatePath)
            if not options.trial:
                sendEvent(templatePath, sciencePath, differencePath, triggerEventTransmitter)
#    if not options.trial:
#        print "Sending shutdown event"
# the technique has changed and I don't know how to do it...

if __name__ == "__main__":
    main()
    # check for memory leaks
    memId0 = 0
    if mwiData.Citizen_census(0, memId0) != 0:
        print mwiData.Citizen_census(0, memId0), "Objects leaked:"
        print mwiData.Citizen_census(mwiData.cout, memId0)
