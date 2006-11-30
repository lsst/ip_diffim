#!/usr/bin/env python
"""
test4.fic
                                                                                
Description
    Fourth test of basic python module access within WCS pipeline
                                                                                
"""
from lsst.fw.Image import CCDImage
import RO.DS9
import os
import sys

print "Test access of single image file: lbcb.20050812.072401_2 and its source extraction.\nStart test4.fic...."

# Find the image, either in system location or local directory
try:
    testData = os.environ['LSSTProto'] +\
                          '/SampleData/data/lbcb.20050812.072401_2.fic'
except:
    testData = "./lbcb.20050812.072401_2.fic"

# Initialization of a single Image fits file or single HDU
try:
    ccdim=CCDImage(testData)
except:
    print 'test4.fic: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test4.fic FAILED.'
    sys.exit(1) 

# Using sextractor, extract sources; then build Images's Source Collection
try:
    sc=ccdim.ExtractSources()
except:
    print 'test4.fic: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test4.fic FAILED.'
    sys.exit(1)

# Use ds9 to display the Source Collection
#ds9=RO.DS9.DS9Win()
#ccdim.Display(ds9)
#sc.DisplaySources(ds9,flagMask=255, starGalCut=0.8)

print "End test4.fic"
