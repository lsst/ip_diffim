#!/usr/bin/env python
"""
doMosaicWCS
                                                                                
Description
     Determine the WCS for an MEF image
                                                                                
Input
    FitsFile    full pathname to the MEF image to be updated
                Format: string; Default: none
    CCDId       Id of CCD within Mosaic Image to be processed
                Format: integer Default: none; 
Interactive Output
    match       StarSourceMatchCollection python object  for the CCD
                or None if no match
    wcs         WCS python object for the CCD
                                                                                
Shell Output
    none

Transient Debris
    lots... from sextractor invocation

"""
import sys
import RO.DS9
from lsst.fw.Image.Image import Image
from lsst.fw.Image.CCDImage import CCDImage
from lsst.fw.Image.MosaicImage import MosaicImage
from lsst.fw.Collection.StarCollection import StarCollection
from lsst.imageproc.WCS.StarSourceMatchCollection import StarSourceMatchCollection
from lsst.imageproc.WCS.WCS import WCS

if (len(sys.argv) < 5):
    print 'doMosaicWCS syntax:  ./doMosaicWCS.py <MEF file> <CCD #> <Mosaic Config> <CCD Config>'
    sys.exit(1)
                                                                                
# Select MEF to process
try:
    mosim=MosaicImage(sys.argv[1],sys.argv[3],sys.argv[4])
except:
    print 'doMosaicWCS: ',sys.exc_type," Cause: ",sys.exc_value
    print 'doMosaicWCS: Failed to acquire MEF: %s' % (sys.argv[1])
    sys.exit(1)
                                                                                
# Extract specific HDU from that MEF
try:
    #ccdim=mosim.GetCCDImage(int(sys.argv[2]),RegionFuzz=0.5)
    ccdim=mosim.GetCCDImage(int(sys.argv[2]))
except:
    print 'doMosaicWCS: ',sys.exc_type," Cause: ",sys.exc_value
    print 'doMosaicWCS: failed to acquire specified HDU\n'
    sys.exit(1)

# Extract sources from that HDU
sc=ccdim.ExtractSources()

# Display if desired
ds9=RO.DS9.DS9Win()
ccdim.Display(ds9)
sc.DisplaySources(ds9,flagMask=0, starGalCut=0)

try:
    scat=StarCollection(ccdim.GetSkyRegion(), \
    'ub1', nStarMax=1500)
except:
    print 'doMosaicWCS: ',sys.exc_type," Cause: ",sys.exc_value
    print 'doMosaicWCS: failed to acquire fiducial stars\n'
    sys.exit(1)

try:
    scat.SortByMag(truncLen=200)
except:
    print 'doMosaicWCS: ',sys.exc_type," Cause: ",sys.exc_value
    print 'doMosaicWCS: failed to sort by magnitude\n'
    sys.exit(1)

# Initialize the match class with the source and fiducial star collections
match=StarSourceMatchCollection(scat, sc, ccdim.GetMetaData())

try:
    match.StarMatch()
    wcs=WCS(match)
except:
    print ("doMosaicWCS: No stars matched. Bailing\n")
