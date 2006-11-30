#!/usr/bin/env python

import sys
import os
from lsst.fw.Image import MosaicImage
from lsst.fw.Collection import *
from lsst.imageproc.WCS import *
import RO.DS9

print "Star match and WCS build on all CCD's from MEF image.\nStart test9...."

# Acquire test image from system or local directory
try:
    testData = os.environ['LSSTProto'] +\
                          '/SampleData/data/642538p.fits'
except:
    testData = "./642538p.fits"

# Acquire mosaic configuration file from system or local directory
# 
#    N O T E: both the mosaic and ccd conf files need to be tuned to the
#	      detector specifics
# 
#	      The mosaic conf file defines the CCD conf file to be used.
try:
    mosaicConfFile = os.environ['LSST_POLICY_DIR'] +\
                          '/CFHT12K_Mosaic.conf'
    ccdConfFile = os.environ['LSST_POLICY_DIR'] +\
                          '/CFHT12K_CCD.conf'
except:
    mosaicConfFile = "./CFHT12K_Mosaic.conf"
    ccdConfFile = "./CFHT12K_CCD.conf"


# Acquire the Mosaic Image (and relevant HDU info)
immef=MosaicImage(testData,mosaicConfFile,ccdConfFile)


# Loop over all HDU in MEF file
for i in range(immef.NumCCDS()):
    # Select CCD to process
    try:
        ccd=immef.GetCCDImage(i)
    except:
        print 'Test9: ',sys.exc_type,"\n Cause: ",sys.exc_value,"\n"
        sys.exit(1)


    # At the moment, all ccd's have the SkyRegion of the full mosaic
    # so just get the star catalog once
    if i==0:
        # Extract fiducial stars from ascii download of USNO-B named "txtcat"
        #                           and located in $WCS_CAT/txtcat
        try:
            scat=StarCollection( \
                                   ccd.GetSkyRegion(),'txtcat', nStarMax=5000)
        except:
            print 'Test9: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
            sys.exit(1)

        try:
            scat.SortByMag(truncLen=300)
        except:
            print 'Test9: ',sys.exc_type,"\nCause: ",sys.exc_value,"\n"
            sys.exit(1)

    # Extract sources from that CCD
    sc=ccd.ExtractSources()

    # Display if desired
    #ds9=RO.DS9.DS9Win()
    #ccd.Display(ds9)
    #sc.DisplaySources(ds9,flagMask=0, starGalCut=0)

    # Initialize the match class with the source and fiducial star collections
    match=StarSourceMatchCollection(scat, sc, ccd.GetMetaData())

    # Perform the star match and acquire the transformation parameters
    try:
        match.StarMatch()
        # Create the WCS for the CCD
        try:
            wcs=WCS(match)
        except:
            print ("Test9: failed during WCS construction for HDU:%d\n" % (i))
    except:
        print ("Test9: No stars matched, so no WCS constructed for HDU:%d\n"% (i))

print "End test9"



