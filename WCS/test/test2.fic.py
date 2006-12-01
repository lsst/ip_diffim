#!/usr/bin/env python
"""
test2.fic
                                                                                
Description
    Second test of basic python module access within WCS pipeline
                                                                                
"""

import sys
import os
from lsst.fw.Image import Image
from lsst.fw.Image import CCDImage
import RO.DS9

print "Test access of target 7 of FIC 642538p.fic and collection of its sources.\nStart test2.fic...."

# Find the image, either in system location or local directory
try:
    testData = os.environ['LSSTProto'] + '/SampleData/data/642538p.fic'
except:
    testData = "./642538p.fic"
                                                                                
# Initialization of an MEF or single Image fits file
try:
    immef=Image(testData)
except:
    print 'test2.fic: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test2.fic FAILED.'
    sys.exit(1)

# Initialization of a single Image fits file or single HDU 
try:
    ccd7=CCDImage(immef.GetCCDHDU(7))
except:
    print 'test2.fic: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test2.fic FAILED.'
    sys.exit(1)

# Using sextractor, extract sources; then build Images's Source Collection
try:
    sc7=ccd7.ExtractSources()
except:
    print 'test2.fic: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test2.fic FAILED.'
    sys.exit(1)

# Use ds9 to display the Source Collection
ds9=RO.DS9.DS9Win()
ccd7.Display(ds9)
sc7.DisplaySources(ds9,flagMask=255, starGalCut=0.8)

print "End test2.fic"
