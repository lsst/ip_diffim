#!/usr/bin/env python
"""
test4
                                                                                
Description
    Fourth test of basic python module access within WCS pipeline
                                                                                
"""
from lsst.fw.Image import CCDImage
import RO.DS9
import os
import sys

print "Start test4...access of single image file then do source extraction."

# Find the image, either in system location or local directory
try:
    testData = os.environ['LSSTProto'] +\
                          '/SampleData/data/lbcb.20050812.072401_2.fits'
except:
    testData = "./lbcb.20050812.072401_2.fits"

if  not os.path.exists(testData):
    print('Test4 FAILED; unable to find data file: %s.\n' %(testData))
    sys.exit(1)
                                                                                
print "Image file: %s" %(testData)


# Initialization of a single Image fits file or single HDU
try:
    ccdim=CCDImage(testData,os.path.abspath("./conf/CCDImage.conf"))
except:
    print 'test4: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test4 FAILED.'
    sys.exit(1) 

# Using sextractor, extract sources; then build Images's Source Collection
try:
    sc=ccdim.ExtractSources()
except:
    print 'test4: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test4 FAILED.'
    sys.exit(1)

# Use ds9 to display the Source Collection
ds9=RO.DS9.DS9Win()
ccdim.Display(ds9)
sc.DisplaySources(ds9,flagMask=255, starGalCut=0.8)

print "End test4"
