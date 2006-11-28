#!/usr/bin/env python
"""
test1.fic
                                                                                
Description
    First test of basic python module access within WCS pipeline
                                                                                
"""
from lsst.apps.fw.Image import CCDImage
import os
import sys

print "Test access of FIC with a single target: lbcb.20050812.072401_2.fic.\nStart test1.fic....."
try:
    testData = os.environ['LSSTProto'] + \
                          '/tests/WCS/data/lbcb.20050812.072401_2.fic'
except:
    testData = "./lbcb.20050812.072401_2.fic"

try:
    im=CCDImage(testData)
except:
    print 'Test1.fic: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'Test1.fic FAILED.'
    sys.exit(1)
print "End test1.fic."
