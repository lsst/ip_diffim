#!/usr/bin/env python

import sys
import os
from lsst.fw.Image import MosaicImage
from lsst.fw.Collection import *
from lsst.imageproc.WCS import *
import RO.DS9

print "Start test8.fic: \nStar match and WCS build on specified CCD (from execute line) from FIC"

if( len(sys.argv) < 2) :
    print("Syntax:    ./test8.fic.py <ccd#> \n")
    sys.exit(1)

# Acquire test image from system or local directory
try:
    testData = os.environ['LSSTProto'] +\
                          '/SampleData/data/642538p.fic'
except:
    testData = "./642538p.fic"

# Acquire mosaic configuration file from system or local directory
#
#    N O T E: both the mosaic and ccd conf files need to be tuned to the
#             detector specifics
#
#             The mosaic conf file defines the CCD conf file to be used.
try:
    mosaicConfFile = os.environ['LSST_POLICY_DIR'] +\
                          '/CFHT12K_Mosaic.conf'
    ccdConfFile = os.environ['LSST_POLICY_DIR'] +\
                          '/CFHT12K_CCD.conf'
except:
    mosaicConfFile = "./CFHT12K_Mosaic.conf"
    ccdConfFile = "./CFHT12K_CCD.conf"
                                                                                

# Acquire the Mosaic Image (and all relevant HDU info)
immef=MosaicImage(testData,mosaicConfFile,ccdConfFile)

# Select CCD to process
try:
    ccd=immef.GetCCDImage(int(sys.argv[1]))
except:
    print 'Test8.fic: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
    sys.exit(1)

# Extract sources from that CCD
sc=ccd.ExtractSources()

# Extract fiducial stars from ascii download of USNO-B named "txtcat"
#                           and located in $WCS_CAT/txtcat
try:
    scat=StarCollection(ccd.GetSkyRegion(),\
                                                    'txtcat', nStarMax=5000)
except:
    print 'Test8.fic: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
    sys.exit(1)

try:
    scat.SortByMag(truncLen=300)
except:
    print 'Test8.fic: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
    sys.exit(1)

# Display if desired
#ds9=RO.DS9.DS9Win()
#ccd.Display(ds9)
#sc.DisplaySources(ds9,flagMask=0, starGalCut=0)

# Initialize the match class with the source and fiducial star collections
match=StarSourceMatchCollection(scat, sc, ccd.GetMetaData())

# Perform the star match and acquire the transformation parameters
# match the source and fiducial stars; results remain within match instance
try:
    match.StarMatch()
    # Create the WCS for the CCD
    try:
        wcs=WCS(match)
    except:
        print 'Test8.fic: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
        print ("Test8.fic: Failure: no WCS constructed\n")
        sys.exit(1)
except:
    print 'Test8.fic: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
    print ("Test8.fic: Failure: no stars matched, no WCS constructed\n")
    sys.exit(1)

print "End test8.fic"
