#!/usr/bin/env python
"""
test2
                                                                                
Description
    Second test of basic python module access within WCS pipeline
                                                                                
"""

import sys
import os
from lsst.apps.fw.Image import Image
from lsst.apps.fw.Image import CCDImage
import RO.DS9

print "Start test2....access HDU 7 of MEF image then collect its sources"

# Find the image, either in system location or local directory
try:
    testData = os.environ['LSSTProto'] + '/tests/WCS/data/642538p.fits'
except:
    testData = os.path.abspath("./642538p.fits")

if  not os.path.exists(testData):
    print('Test2 FAILED; unable to find data file: %s.\n' %(testData))
    sys.exit(1)
                                                                                
print "Image file: %s" %(testData)
                                                                                
# Initilization of an MEF or single Image fits file
try:
    immef=Image(testData)
except:
    print 'test2: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test2 FAILED.'
    sys.exit(1)

# Initialization of a single Image fits file or single HDU 
try:
    ccd7=CCDImage(immef.GetCCDHDU(7),os.path.abspath("./conf/CCDImage.conf"))
except:
    print 'test2: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test2 FAILED.'
    sys.exit(1)

# Using sextractor, extract sources; then build Images's Source Collection
try:
    sc7=ccd7.ExtractSources()
except:
    print 'test2: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test2 FAILED.'
    sys.exit(1)

# Use ds9 to display the Source Collection
#ds9=RO.DS9.DS9Win()
#ccd7.Display(ds9)
#sc7.DisplaySources(ds9,flagMask=255, starGalCut=0.8)

print "End test2"

