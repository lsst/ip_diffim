#!/usr/bin/env python
"""
doWCS
                                                                                
Description
     Determine and then update the WCS for a fits image
                                                                                
Input
    FitsFile    full pathname to the FITS image to be updated;
                Default: none

Interactive Output
    match       StarSourceMatchCollection python object or None if no match
    wcs         WCS python object

Shell Output
    none

Transient Debris
    lots.....from sextractor invocation
                                                                                
"""
import sys
from lsst.fw.Image.Image import Image
from lsst.fw.Image.CCDImage import CCDImage
from lsst.fw.Collection import StarCollection
import RO.DS9
from lsst.imageproc.WCS import *

if (len(sys.argv) < 3):
    print 'doWCS syntax:  ./doWCS.py <single FITS image file> <CCDImage configuration>'
    sys.exit(1)

# Select CCD to process
try:
    ccdim=CCDImage(sys.argv[1], sys.argv[2])
except:
    print 'doWCS: ' ,sys.exc_type," Cause: ",sys.exc_value
    print 'doWCS: Failed to acquire HDU: %s' % (sys.argv[1])
    sys.exit(1)

# Extract sources from that CCD
sc=ccdim.ExtractSources()

# Display if desired
ds9=RO.DS9.DS9Win()
ccdim.Display(ds9)
sc.DisplaySources(ds9,flagMask=0, starGalCut=0)


try:
    scat=StarCollection(ccdim.GetSkyRegion(), 'gsc', nStarMax=500)
except:
    print 'doWCS: ',sys.exc_type," Cause: ",sys.exc_value
    print 'doWCS: failed to acquire GSC fiducial stars\n'
    sys.exit(1)

try:
    scat.SortByMag()
except:
    print 'doWCS: ',sys.exc_type," Cause: ",sys.exc_value
    print 'doWCS: failed to sort by magnitude'
    sys.exit(1)

# Initialize the match class with the source and fiducial star collections
match=StarSourceMatchCollection(scat, sc, ccdim.GetMetaData())


# Perform the star match and acquire the transformation parameters
try: 
    match.StarMatch()
    # Create the WCS for the CCD
    try:
        wcs=WCS(match)
    except:
        print 'doWCS: Failed during construction of the WCS'
        sys.exit(1)
except:
    print 'doWCS: ',sys.exc_type," Cause: ",sys.exc_value
    print 'doWCS: No stars matched \n'
    sys.exit(1)

