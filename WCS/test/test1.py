#!/usr/bin/env python
"""
test1
                                                                                
Description
    First test of basic python module access within WCS pipeline
                                                                                
"""
from lsst.fw.Image import CCDImage
import os
import sys

print "Start test1..... access single image file"
try:
    testData = os.environ['LSSTProto'] + \
                          '/SampleData/data/lbcb.20050812.072401_2.fits'
except:
    testData = os.path.abspath("./data/lbcb.20050812.072401_2.fits")

if  not os.path.exists(testData):
    print('Test1 FAILED; unable to find data file: %s.\n' %(testData))
    sys.exit(1)

print "Image file: %s" %(testData)
try:
    im=CCDImage(testData,os.path.abspath("./conf/CCDImage.conf"))
except:
    print 'Test1: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'Test1 FAILED.'
    sys.exit(1)
print "End test1."
