#!/usr/bin/env python
"""
StarSourceMatchCollection

    The StarSourceMatchCollection class uses a star matching algorithm 
    to determine the best match between an image's sources and a set of 
    fiducial stars. The resultant array of matches may be used to 
    construct a WCS for the source image
"""
__all__ = ["StarSourceMatchCollection"]

import numarray 
import os
import re
import string
import sys
import tempfile
import logging
import fileinput
from RO.StringUtil import *
from lsst.apps.fw.Policy import Policy


#*************************************************
class StarSourceMatchCollection(object):
    #------------------------------------------------------------------------
    def __init__ (self,
                stars,
                sources,
                metaData,
                logFile="StarMatch.log",
                policyFile=None,
                **kws):
        """
        Initialization of StarSourceMatchCollection

        Input
            stars           StarCollection containing stars from catalog
            sources         ImageSourceCollection containing detected sources 
                            in image; Default: none
            logFile         filename for logging; Default is "Log.out"
                            Filename format: ascii string
            policyFile      filename of Policy rules; 
                            1st default: "$LSSTProto/StarSourceMatchCollection.conf"
                            2cd default: "./StarSourceMatchCollection.conf"
                            Filename format: ascii string
        Return
            none
        """

        if ( not policyFile ):
            conf_file = "StarSourceMatchCollection.conf"
        else:
            conf_file = policyFile

        self.policy = Policy(conf_file, kws)
        self.logFile = logFile
        self.stars = stars
        self.sources = sources
        self.metaData = metaData
        
    #------------------------------------------------------------------------
    
    #--------------------------------------------------------------------
    def StarMatch(self):
        """
        Run the StarLink program FINDOFF. Creates two output files which are
        subsequently merged by  a later user call to
        StarSourceMatchCollection.GetMatches()

        Input  
            None

        Return 
            None

        Input Files
            starFindoffIn       file named by Policy, created locally from
                                StarCollection, then read by StarLink:FINDOFF 
            sourceFindoffIn     file named by Policy, created locally from
                                ImageSourceCollection, then read by 
                                StarLink:FINDOFF 

        Output Files
            starFindoffOut      file named by Policy, created by 
                                StarLink:FINDOFF execution, contains matched
                                stars from StarCollection
            sourceFindoffOut    file named by Policy, created by 
                                StarLink:FINDOFF execution, contains matched
                                stars from ImageSourceCollection

        """

        self.starFindoffInp = self.policy.Get('starFindoffInp')
        self.sourceFindoffInp = self.policy.Get('sourceFindoffInp')
        self.starFindoffOut = self.policy.Get('starFindoffOut')
        self.sourceFindoffOut = self.policy.Get('sourceFindoffOut')
        self.maxSources = self.policy.Get('maxSources')
        measError = self.policy.Get('measError')
        maxOffset = self.policy.Get('maxOffset')
        minMatch = self.policy.Get('minMatch')
        minComplete = self.policy.Get('minComplete')

        # delete any FINDOFF files that exist

        for fileName in (self.starFindoffInp, self.sourceFindoffInp, self.starFindoffOut, self.sourceFindoffOut):
            if (os.access(fileName,os.R_OK|os.W_OK)):
                os.unlink(fileName)

        #--------------------------------------------------------------------
        # write xi/eta/mag/id from StarCollection to FINDOFF input file
        #
        self.stars.WriteStarsXiEta(self.starFindoffInp)

        #--------------------------------------------------------------------
        # write source xi/eta/mag/id from ImageSourceCollection to FINDOFF input file

        self.sources.WriteSourcesXiEta(self.sourceFindoffInp, takeTop=self.maxSources)

        # At this point, we're ready to match stars. Run FINDOFF

        findoffCmd = 'findoff inlist=\'"%s,%s"\' ndfnames=false outlist=\'"%s,%s"\' error=%d maxdisp=%d minmatch=%d complete=%f fast=TRUE failsafe=TRUE' % \
                     (self.starFindoffInp, self.sourceFindoffInp, self.starFindoffOut, self.sourceFindoffOut, measError, maxOffset, minMatch, minComplete )

        try:
            os.system(findoffCmd)
        except:
            raise EnvironmentError, 'Failed to run FINDOFF'

 
    #-----------------------------------------------------------------------
    def GetMatches(self):
        """
        Return the matched set of fiducial/source stars.  Reads the output
        files generated by StarSourceMatchCollection.StarMatch

        Input
            none
        Return
            nstars  count of matched stars
            ra      numarray of ra of fiducial stars; 
                    Format: float; Units: decimal degrees
            dec     numarray of dec of fiducial stars; 
                    Format: float; Units: decimal degrees
            xf2     numarray of x pixel location of source stars
            yf2     numarray of y pixel location of source stars

        Input Files:
            starFindoffOut      file named by Policy, created by 
                                StarLink:FINDOFF execution, contains matched
                                stars from StarCollection
            sourceFindoffOut    file named by Policy, created by 
                                StarLink:FINDOFF execution, contains matched
                                stars from ImageSourceCollection
        """

        starsMatched = self._readarray(self.starFindoffOut)
        sourcesMatched = self._readarray(self.sourceFindoffOut)

        self.nMatch = starsMatched.getshape()[0]

        self.idStars = starsMatched[:,4].astype('Int32')
        self.idSources = sourcesMatched[:,4].astype('Int32')
        
        # Note the unsymmetric way that sorces and stars are accessed - ugly!
        return(self.nMatch,self.stars.raArray[self.idStars],self.stars.decArray[self.idStars],self.sources.sourceArray[self.idSources,0],self.sources.sourceArray[self.idSources,1])

    #-----------------------------------------------------------------------

    def WriteMatches(self, fileName):
        """
        Write (RA,DEC) world coordinates and (x,y)  pixel coordinates for
        matched stars.

        Input
            filename    name of output file; Format: string; Default: none

        Return
            none

        Output File    
            World Coords and Pixel Coords for all matched stars
            File Format: '%f %f %f %f\n' % (ra dec x y)
        """
        (n, ra, dec, x, y) = self.GetMatches()

        f=os.open(fileName,os.O_CREAT|os.O_RDWR)
        for i in numarray.arange(n):
            os.write(f,'%f %f %f %f\n' % (x[i],y[i],ra[i],dec[i]))

        os.close(f)
                     
    #-----------------------------------------------------------------------
        

    def WriteMatchesAstrom(self, fileName):
        """
        Write (RA,DEC) world coordinates and proper motions for
        matched stars in a form suitable for StarLink Astrom program.

        Input
            filename    name of output file; Format: string; Default: none

        Return
            none

        Output File    
            ASTR file format
        """
        (n, ra, dec, x, y) = self.GetMatches()
        raPM = self.stars.raPMArray[self.idStars]
        decPM = self.stars.decPMArray[self.idStars]

        # must translate : into blanks for astrom
        xlate = string.maketrans(':', ' ')

        # must delete file if it exists - else may leave garbage at end

        if (os.access(fileName,os.R_OK|os.W_OK)):
            os.unlink(fileName)
            
        f=os.open(fileName,os.O_CREAT|os.O_RDWR)

        # print out file header before match data

        os.write(f,'ASTR\n')
        fieldRa = self.metaData.GetKW('ra')
        fieldDec = self.metaData.GetKW('dec')
        epoch = self.metaData.GetKW('epoch')
        os.write(f,'%s %s J2000.0 %s\n' % (fieldRa.translate(xlate), fieldDec.translate(xlate), epoch))

        # now the match data
        for i in numarray.arange(self.nMatch):
            raStr= dmsStrFromDeg(ra[i]/15.0).translate(xlate)
            decStr = dmsStrFromDeg(dec[i]).translate(xlate)
            os.write(f,'%s %s %f %f J2000.0\n%f %f\n' % (raStr,decStr,raPM[i],decPM[i],x[i],y[i]))

        os.close(f)
                                         
     #-----------------------------------------------------------------------
    def WriteMatchesIRAF(self, fileName):
        """
        Write (RA,DEC) world coordinates and proper motions for
        matched stars in a form suitable for input to iraf ccmap

        Input
            filename    name of output file; Format: string; Default: none

        Return
            none

        Output File    
            iraf file format
        """
        (n, ra, dec, x, y) = self.GetMatches()

        # must delete file if it exists - else may leave garbage at end

        if (os.access(fileName,os.R_OK|os.W_OK)):
            os.unlink(fileName)
            
        f=os.open(fileName,os.O_CREAT|os.O_RDWR)


        # write the match data
        for i in numarray.arange(self.nMatch):
            os.write(f,'%f %f %f %f\n' % (x[i], y[i], ra[i], dec[i]))

        os.close(f)
                                         
    #-----------------------------------------------------------------------
        
    def _readarray(self, fileName):
        """
        Input
            filename    name of file written by StarLink program FINDOFF;
                        Format: string; Default: none
        Return
            data        2-D numarray containing data extracted from input file.
        """
        fileinp = fileinput.input(fileName)
        inplist = []
        for l in fileinp:
            if l.find('#') == -1:
                inplist = inplist + [l.split()]

        fileinp.close()

        nrow = len(inplist)
        ncol = len(inplist[0])
        dat = numarray.zeros([nrow, ncol],'Float32')
        for i in range(nrow):
            for j in range(ncol):
                dat[i,j] = float(inplist[i][j])

        return dat
