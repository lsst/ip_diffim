#!/usr/bin/env python
"""
Fifth test of basic python module access within WCS pipeline
                                                                                
"""

import os
import sys
from lsst.fw.Image import CCDImage
from lsst.fw.Collection import *
from lsst.imageproc.WCS import *
import RO.DS9

print "Start test5...star matching on image taken at the pole"

# Find the image, either in system location or local directory
#                          '/SampleData/data/polaris_20sec_2.fits'
try:
    testData = os.environ['LSSTProto'] +\
                          '/SampleData/data/lbcb.20050812.072401_2.fits'
except:
    testData = "./lbcb.20050812.072401_2.fits"

if  not os.path.exists(testData):
    print('Test5 FAILED; unable to find data file: %s.\n' %(testData))
    sys.exit(1)
                                                                                
print "Image file: %s" %(testData)

# Initialization of a single Image fits file or single HDU
try:
    ccdim=CCDImage(testData,os.path.abspath("./conf/CCDImage.conf"))
except:
    print 'test5: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test5 FAILED.'
    sys.exit(1)

# Using sextractor, extract sources; then build Images's Source Collection
try:
    sc=ccdim.ExtractSources()
except:
    print 'test5: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test5 FAILED.'
    sys.exit(1)

# Use ds9 to display the Source Collection
ds9=RO.DS9.DS9Win()
ccdim.Display(ds9)
sc.DisplaySources(ds9,flagMask=255, starGalCut=0.8)

# Acquire a fiducial stars for GSC catalog in region overlapping source image
try:
    scat=StarCollection(ccdim.GetSkyRegion(), 'gsc', nStarMax=500)
except:
    print 'Test5: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
    sys.exit(1)
                                                                                
# Initialize the match class with the source and fiducial star collections
try:
	match=StarSourceMatchCollection(scat, sc, ccdim.GetMetaData(),policyFile=os.path.abspath("./conf/StarSourceMatchCollection.conf"))
except:
    print 'Test5: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
    sys.exit(1)
                                                                                
# Diagnostics to display the matched stars
#match.printDefaults()
#match.printInput()

# Match the source and fiducial stars; results remain within match instance
try:
    match.StarMatch()
    print "End test5"
except:
    print 'Test5: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
    print "Found no matched stars"

