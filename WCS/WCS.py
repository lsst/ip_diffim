"""
WCS

Creates a WCS object from a StarSourceMatchCollection
"""
__all__ = ["WCS"]

from lsst.apps.imageproc.WCS import StarSourceMatchCollection
from lsst.apps.support.PyWCSlib import PyWCSlib
import numarray 
import os
import re
import string
import sys
import tempfile
import logging
import pyfits
from RO.StringUtil import *
from lsst.apps.fw.Policy import Policy


#*********************************************************************
class WCS:
    def __init__ (self, match=None, kwWCS=None, policyFile=None, **kw):
        """
        Input
            match           StarSourceMatchCollection; Default: none
            policyFile      Filename for class Policy specification;
                            Format: string    
                            Default: 
                                if none provided, use "WCS.conf"
            kw              Keyword dictionary which is used, along with 
                            Policy, to manage processing; Format: dictionary
        Return
            none
        """
        if ( not policyFile ):
            conf_file = "WCS.conf"
        else:
            conf_file = policyFile
        self.policyFile = policyFile
        self.config = Policy(self.policyFile).conf
        # activeKW = kwMerge(kw, policy)
        #

        if ( match ):
            self._BuildWCSFromMatch(match)
        elif (kwWCS):
            self.wcsHeader = kwWCS
            self.worldCoor = PyWCSlib.wcsinit(2, self.wcsHeader)
        else:
            raise EnvironmentError, 'Insufficient args supplied to WCS constructor'
        
    def _BuildWCSFromMatch(self, match):
#
# generate astrom input file
#
        try:
            astromInp = self.config['astromInpFile']
        except:
            astromInp = './TmpAstromConfig'
        try:
            astromFITS = self.config['astromFITSFile']
        except:
            astromFITS = './TmpAstromFITS'
        
        match.WriteMatchesAstrom(astromInp)
        #
        # run astrom
        #
        astromCmd = 'astrom %s fits=%s summary=astrom.info log=astrom.log' % (astromInp, astromFITS)
        # print astromCmd

        os.system('rm -f ' + astromFITS + '*.fits')
        try:
            os.system(astromCmd)
        except:
            raise EnvironmentError, 'Failed to run ASTROM'
        #
        # slurp up the wcs file produced by astrom, and create a struct WorldCoor
        #
        wcsFITS = astromFITS + '01.fits'
        self.wcsIm = pyfits.open(wcsFITS)
        self.wcsHeader = self.wcsIm[0].header
        self.worldCoor = PyWCSlib.wcsinit(2, self.wcsHeader)


    def SetFITSWCS(self, ccdim):
        """
        SetFITSWCS updates FITS header in given image to incorporate WCS.
           Might be good idea to delete ALL wcs-related fits keywords first.
                                                                                
        Input
            ccdim           CCDImage ; Default: none
        Return
            none
        """
        destHeader = ccdim.GetHeader()
        ignoreKW = ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2']
        for kwVal in self.wcsHeader.items():
            kw = kwVal[0]
            val = kwVal[1]
            if ignoreKW.count(kw) == 0:
                destHeader.update(kw, val)
                print 'update', kw, val

        ccdim.hdus.flush()
    
    def WCS2pix(self, ra, dec):
        """
        WCS2pix generates and returns pixel x, y coordinates from ra, dec

        Input
            ra          RA in source image; 
                        Format: Float; Units: decimal degrees
            dec         DEC in source image; Units: decimal degrees
                        Format: Float; Units: decimal degrees
        Return
            pixcoord    (x,y) pixel coordinate in source image; 
                        Format: list; Units: pixels
        """
        radec = numarray.array((ra, dec), 'Float64')
        (pixcoord, phi, theta, imgcrd, stat) = PyWCSlib.wcss2p(self.worldCoor, radec)
        return pixcoord


    def Pix2WCS(self, x, y):
        """
        Pix2WCS generates and returns ra, dec coordinates from pixel x,y

        Input
            x            x pixel coordinate in source image; 
                         Format: float; Units: pixels
            y            y pixel coordinate in source image; 
                         Format: Float; Units: pixels
        Return
            worldcoordinate
                         (ra,dec) world coordinate;
                         Format: list; Units: decimal degrees
        """
        xy = numarray.array((x,y), 'Float64')
        (worldcoord, imgcrd, phi, theta, stat) = PyWCSlib.wcsp2s(self.worldCoor, xy)
        return worldcoord

