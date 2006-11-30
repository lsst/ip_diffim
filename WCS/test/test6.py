#!/usr/bin/env python
"""
test6
                                                                                
Description
    Sixth test of basic python module access within WCS pipeline
                                                                                
"""

import sys
import os
from lsst.fw.Image import CCDImage
from lsst.fw.Collection import *
from lsst.imageproc.WCS import *
import RO.DS9

print "Test star match and WCS build on single image: lbcb.20050812.072401_2.fits\nStart test6...."

# Acquire test image from system or local directory
try:
    testData = os.environ['LSSTProto'] +\
                          '/SampleData/data/lbcb.20050812.072401_2.fits'
except:
    testData = "./lbcb.20050812.072401_2.fits"
                                                                                
if  not os.path.exists(testData):
    print('Test6 FAILED; unable to find data file: %s.\n' %(testData))
    sys.exit(1)

# Acquire the CCD Image and relevant HDU info
ccdim=CCDImage(testData,os.path.abspath("./conf/CCDImage.conf"))

# Extract sources from CCD
sc=ccdim.ExtractSources()

# Display if desired
#ds9=RO.DS9.DS9Win()
#ccdim.Display(ds9)
#sc.DisplaySources(ds9,flagMask=255, starGalCut=0.8)

# Acquire Fiducial Stars from gsc catalog
try:
    scat=StarCollection(ccdim.GetSkyRegion(), 'gsc', nStarMax=500)
except:
    print 'Test6: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
    sys.exit(1)

# Initialize the match class with the source and fiducial star collections
match=StarSourceMatchCollection(scat, sc, ccdim.GetMetaData(),policyFile=os.path.abspath("./conf/StarSourceMatchCollection.conf"))

# Perform the star match and acquire the transformation parameters
# match the source and fiducial stars; results remain within match
try:
    match.StarMatch()
    try:
        # Create the WCS for the CCD
        wcs=WCS(match)
    except:
        print 'Test6: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
        print ("Test6: Failure: no WCS constructed\n")
        sys.exit(1)

except:
    print 'Test6: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
    print ("Test6: Failure: no WCS constructed\n")
    sys.exit(1)

print "End Test6"
