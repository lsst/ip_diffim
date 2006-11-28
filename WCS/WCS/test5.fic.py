#!/usr/bin/env python
"""
Fifth test of basic python module access within WCS pipeline
                                                                                
"""

import os
import sys
from lsst.apps.fw.Image import CCDImage
from lsst.apps.fw.Collection import *
from lsst.apps.imageproc.WCS import *
import RO.DS9

print "Test star matching on image taken at the pole: polaris_20sec_2.fic\nStart test5.fic...."

# Find the image, either in system location or local directory
try:
    testData = os.environ['LSSTProto'] +\
                          '/tests/WCS/data/polaris_20sec_2.fic'
#                          '/tests/WCS/data/lbcb.20050812.072401_2.fic'
except:
    testData = "./lbcb.20050812.072401_2.fic"

# Initialization of a single target FIC Image file
try:
    ccdim=CCDImage(testData)
except:
    print 'test5.fic: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test5.fic FAILED.'
    sys.exit(1)

# Using sextractor, extract sources; then build Images's Source Collection
try:
    sc=ccdim.ExtractSources()
except:
    print 'test5.fic: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    print 'test5.fic FAILED.'
    sys.exit(1)

# Use ds9 to display the Source Collection
#ds9=RO.DS9.DS9Win()
#ccdim.Display(ds9)
#sc.DisplaySources(ds9,flagMask=255, starGalCut=0.8)

# Acquire a fiducial stars for GSC catalog in region overlapping source image
try:
    scat=StarCollection(ccdim.GetSkyRegion(), 'gsc', nStarMax=500)
except:
    print 'Test5.fic: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
    sys.exit(1)
                                                                                
# Initialize the match class with the source and fiducial star collections
try:
	match=StarSourceMatchCollection(scat, sc, ccdim.GetMetaData())
except:
    print 'Test5.fic: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
    sys.exit(1)
                                                                                
# Diagnostics to display the matched stars
#match.printDefaults()
#match.printInput()

# Match the source and fiducial stars; results remain within match instance
try:
    match.StarMatch()
    print "End test5.fic"
except:
    print 'Test5.fic: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
    print "Found no matched stars"
