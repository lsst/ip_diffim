#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#
from __future__ import division

import os
from optparse import OptionParser
import lsst.afw.image                    as afwImage
import lsst.afw.geom                     as afwGeom
import lsst.afw.detection                as afwDet
import lsst.pex.policy                   as pexPolicy
import lsst.pex.logging                  as pexLog 
import lsst.daf.persistence              as dafPersist
import lsst.daf.base                     as dafBase
import numpy as num
import lsst.afw.display.ds9 as ds9
import lsst.pex.config as pexConfig

scaling = 5

def parseOptions():
    """Parse the command line options."""
    parser = OptionParser(
            usage="""%prog cdDiffSources crDiffExposure
            
Read in sources and test for junk""")
    options, args = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")
    return options, args

def readSourceSet(boostFile):
    loc = dafPersist.LogicalLocation(boostFile)
    storageList = dafPersist.StorageList()
    additionalData = dafBase.PropertySet()
    persistence = dafPersist.Persistence.getPersistence(pexPolicy.Policy())
    storageList.append(persistence.getRetrieveStorage("BoostStorage", loc))
    psvptr = persistence.unsafeRetrieve("PersistableSourceVector", storageList, additionalData)
    psv = afwDet.PersistableSourceVector.swigConvert(psvptr)
    return psv.getSources()

class DiaSourceAnalystConfig(pexConfig.Config):
    srcBadMaskPlanes = pexConfig.ListField(
        dtype = str,
        doc = """Mask planes that lead to an invalid detection.
                 Options: EDGE SAT BAD CR INTRP
                 E.g. : EDGE SAT BAD allows CR-masked and interpolated pixels""",
        default = ("EDGE", "SAT", "BAD")
    )
    fBadPixels = pexConfig.Field(
        dtype = float,
        doc = "Fraction of bad pixels allowed in footprint",
        default = 0.1
    )
    fluxPolarityRatio = pexConfig.Field(
        dtype = float,
        doc = "Minimum fraction of flux in correct-polarity pixels",
        default = 0.75
    )
    nPolarityRatio = pexConfig.Field(
        dtype = float,
        doc = "Minimum fraction of correct-polarity pixels in unmasked subset",
        default = 0.7
    )
    nMaskedRatio = pexConfig.Field(
        dtype = float,
        doc = "Minimum fraction of correct-polarity unmasked to masked pixels",
        default = 0.6,
    )
    nGoodRatio =  pexConfig.Field(
        dtype = float,
        doc = "Minimum fraction of correct-polarity unmasked to all pixels",
        default = 0.5
    )


class DiaSourceAnalyst(object):
    def __init__(self, config):
        self.config = config
            
        self.bitMask = 0
        srcBadMaskPlanes = self.config.srcBadMaskPlanes
        for maskPlane in srcBadMaskPlanes:
            self.bitMask |= afwImage.MaskU_getPlaneBitMask(maskPlane) 
        
    def countDetected(self, mask):
        idxP = num.where(mask & afwImage.MaskU_getPlaneBitMask("DETECTED"))
        idxN = num.where(mask & afwImage.MaskU_getPlaneBitMask("DETECTED_NEGATIVE"))
        return len(idxP[0]), len(idxN[0])
        
    def countMasked(self, mask):
        idxM = num.where(mask & self.bitMask)
        return len(idxM[0])
        
    def countPolarity(self, mask, pixels):
        unmasked = ((mask & self.bitMask) == 0)
        idxP  = num.where( (pixels >= 0) & unmasked)
        idxN  = num.where( (pixels < 0)  & unmasked)
        fluxP = num.sum(pixels[idxP])
        fluxN = num.sum(pixels[idxN])
        #import pdb; pdb.set_trace()
        return len(idxP[0]), len(idxN[0]), fluxP, fluxN
        
    def testSource(self, source, subMi):
        imArr, maArr, varArr   = subMi.getArrays()    
        flux                   = source.getApFlux()
        
        nPixels                = subMi.getWidth() * subMi.getHeight()
        nPos, nNeg, fPos, fNeg = self.countPolarity(maArr, imArr)
        nDetPos, nDetNeg       = self.countDetected(maArr)
        nMasked                = self.countMasked(maArr)
        assert (nPixels == (nMasked + nPos + nNeg))
        
        # 1) Too many pixels in the detection are masked  
        fMasked    = (nMasked / nPixels)
        fMaskedTol = self.config.fBadPixels
        if fMasked > fMaskedTol:
            pexLog.Trace("lsst.ip.diffim.DiaSourceAnalysis", 1,
                "Candidate %d : BAD fBadPixels %.2f > %.2f" % (source.getId(), fMasked, fMaskedTol))
            return False
            
        if flux > 0:
            # positive-going source
            fluxRatio  = fPos / (fPos + abs(fNeg))
            ngoodRatio = nPos / nPixels
            maskRatio  = nPos / (nPos + nMasked)
            npolRatio  = nPos / (nPos + nNeg)
        else:
            # negative-going source
            fluxRatio  = abs(fNeg)/ (fPos + abs(fNeg))
            ngoodRatio = nNeg / nPixels
            maskRatio  = nNeg / (nNeg + nMasked)
            npolRatio  = nNeg / (nNeg + nPos)
        
        # 2) Not enough flux in unmasked correct-polarity pixels
        fluxRatioTolerance = self.config.fluxPolarityRatio
        if fluxRatio < fluxRatioTolerance:
            pexLog.Trace("lsst.ip.diffim.DiaSourceAnalysis", 1,
                "Candidate %d : BAD flux polarity %.2f < %.2f (pos=%.2f neg=%.2f)" % (source.getId(), 
                fluxRatio, fluxRatioTolerance, fPos, fNeg))
            return False
        
        # 3) Not enough unmasked pixels of correct polarity
        polarityTolerance = self.config.nPolarityRatio
        if npolRatio < polarityTolerance:
            pexLog.Trace("lsst.ip.diffim.DiaSourceAnalysis", 1,
                "Candidate %d : BAD polarity count %.2f < %.2f (pos=%d neg=%d)" % (source.getId(), 
                npolRatio, polarityTolerance, nPos, nNeg))
            return False
            
        # 4) Too many masked vs. correct polarity pixels
        maskedTolerance = self.config.nMaskedRatio
        if maskRatio < maskedTolerance:
            pexLog.Trace("lsst.ip.diffim.DiaSourceAnalysis", 1,
                "Candidate %d : BAD unmasked count %.2f < %.2f (pos=%d neg=%d mask=%d)" % (source.getId(), 
                maskRatio, maskedTolerance, nPos, nNeg, nMasked))
            return False

        # 5) Too few unmasked, correct polarity pixels
        ngoodTolerance = self.config.nGoodRatio
        if ngoodRatio < ngoodTolerance:
            pexLog.Trace("lsst.ip.diffim.DiaSourceAnalysis", 1,
                "Candidate %d : BAD good pixel count %.2f < %.2f (pos=%d neg=%d tot=%d)" % (source.getId(), 
                ngoodRatio, ngoodTolerance, nPos, nNeg, nPixels))
            return False
                    
        pexLog.Trace("lsst.ip.diffim.DiaSourceAnalysis", 1,
            "Candidate %d : OK flux=%.2f nPos=%d nNeg=%d nTot=%d nDetPos=%d nDetNeg=%d fPos=%.2f fNeg=%2f" % (source.getId(),
            flux, nPos, nNeg, nPixels, nDetPos, nDetNeg, fPos, fNeg))
        return True                                        
        

def main():
    """Main program"""
    options, args = parseOptions()
    (crDiffSourceFile, crDiffExposureFile) = args

    crDiffSources = readSourceSet(crDiffSourceFile)
    crDiffExposure = afwImage.ExposureF(crDiffExposureFile)
    #import pdb; pdb.set_trace()
    
    analyst = DiaSourceAnalyst()
    
    expX0 = crDiffExposure.getX0()
    expY0 = crDiffExposure.getY0()
    expX1 = expX0 + crDiffExposure.getWidth() - 1
    expY1 = expY0 + crDiffExposure.getHeight() - 1
    
    for i in range(crDiffSources.size()):
        crDiffSource = crDiffSources[i]

        # This segfaults; revisit once the stack stabilizes
        #footprint    = crDiffSource.getFootprint()
        #bbox         = footprint.getBBox()
        
        xAstrom      = crDiffSource.getXAstrom()
        yAstrom      = crDiffSource.getYAstrom()
        Ixx          = max(1.0, crDiffSource.getIxx())
        Iyy          = max(1.0, crDiffSource.getIyy())
        x0           = max(expX0, int(xAstrom - scaling * Ixx))
        x1           = min(expX1, int(xAstrom + scaling * Ixx))
        y0           = max(expY0, int(yAstrom - scaling * Iyy))
        y1           = min(expY1, int(yAstrom + scaling * Iyy))
        bbox         = afwGeom.Box2I(afwGeom.Point2I(x0, y0),
                                     afwGeom.Point2I(x1, y1))
        subExp       = afwImage.ExposureF(crDiffExposure, bbox)
        subMi        = subExp.getMaskedImage()
        imArr, maArr, varArr = subMi.getArrays()

        if analyst.testSource(crDiffSource, subMi):
            ds9.mtv(subExp, frame=1)
            raw_input('Next: ')
if __name__ == "__main__":
    main()
