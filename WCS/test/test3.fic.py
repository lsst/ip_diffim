#!/usr/bin/env python
"""
test3.fic
                                                                                
Description
    Third test of basic python module access within WCS pipeline
                                                                                
"""
from lsst.fw.Image import CCDImage
import RO.DS9
import sys
import os

print "Test source collection of an image at pole: polaris_20sec_2.fic\nStart test3.fic...."

# Find the image, either in system location or local directory
try:
    testData = os.environ['LSSTProto'] + '/SampleData/data/polaris_20sec_2.fic'
except:
    testData = "./polaris_20sec_2.fic"
                                                                                
# Initialization of a single Image fits file or single HDU
try:
    ccdim=CCDImage(testData)
except:
    print 'test3.fic: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test3.fic FAILED.'
    sys.exit(1)

# Using sextractor, extract sources; then build Images's Source Collection
try:
    sc=ccdim.ExtractSources()
except:
    print 'test3.fic: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test3.fic FAILED.'
    sys.exit(1)

# Use ds9 to display the Source Collection
#ds9=RO.DS9.DS9Win()
#ccdim.Display(ds9)
#sc.DisplaySources(ds9,flagMask=255, starGalCut=0.8)

print "End test3.fic"
