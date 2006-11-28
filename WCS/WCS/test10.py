#!/usr/bin/env python

import sys
import os
from lsst.apps.imageproc.WCS import *

print "Start test10: \nConstruct and use WCS from FITS keywords"

wcsKW = {}
wcsKW['CTYPE1']  = 'RA---TAN'            
wcsKW['CTYPE2']  = 'DEC--TAN'            
wcsKW['CRVAL1']  =        209.845842793  
wcsKW['CRVAL2']  =         39.343830375  
wcsKW['CRPIX1']  =            6338.2370  
wcsKW['CRPIX2']  =            4069.5796  
wcsKW['LIM1_1']  =                  1.0 
wcsKW['LIM2_2']  =                  1.0 
wcsKW['WAT0_001']= 'system=image'       
wcsKW['WAT1_001']= 'wtype=tan axtype=ra' 
wcsKW['WAT2_001']= 'wtype=tan axtype=dec' 
wcsKW['CD1_1']   =          0.000056218
wcsKW['CD1_2']   =         -0.000000660
wcsKW['CD2_1']   =         -0.000000659
wcsKW['CD2_2']   =         -0.000056833

wcs=WCS(kwWCS=wcsKW)
(ra, dec) = wcs.Pix2WCS(6338.237, 4069.5796)

print ra, dec
print "End test10"

