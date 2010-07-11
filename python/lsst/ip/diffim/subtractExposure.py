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

from .subtractMaskedImage import subtractMaskedImage
from .warpTemplateExposure import warpTemplateExposure

import lsst.pex.logging as pexLog
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9

def subtractExposure(exposureToConvolve, exposureToNotConvolve, policy, display=False):
    # Make sure they end up the same dimensions on the sky
    templateWcs    = exposureToConvolve.getWcs() 
    scienceWcs     = exposureToNotConvolve.getWcs()

    # LLC
    templateOrigin = templateWcs.xyToRaDec(0, 0)
    scienceOrigin  = scienceWcs.xyToRaDec(0, 0)
    # URC
    templateLimit  = templateWcs.xyToRaDec(exposureToConvolve.getWidth(),
                                           exposureToConvolve.getHeight())
    scienceLimit   = scienceWcs.xyToRaDec(exposureToNotConvolve.getWidth(),
                                          exposureToNotConvolve.getHeight())

    pexLog.Trace("lsst.ip.diffim.subtractExposure", 1,
                 "Template limits : %f,%f -> %f,%f" %
                 (templateOrigin[0], templateOrigin[1],
                  templateLimit[0], templateLimit[1]))
    pexLog.Trace("lsst.ip.diffim.subtractExposure", 1,
                 "Science limits : %f,%f -> %f,%f" %
                 (scienceOrigin[0], scienceOrigin[1],
                  scienceLimit[0], scienceLimit[1]))
                  

    # Remap if they do not line up in size on the sky, and size in pixels
    if ( (templateOrigin != scienceOrigin) or \
         (templateLimit  != scienceLimit)  or \
         (exposureToConvolve.getHeight() != exposureToNotConvolve.getHeight()) or \
         (exposureToConvolve.getWidth()  != exposureToNotConvolve.getWidth()) ):
        pexLog.Trace("lsst.ip.diffim.subtractExposure", 1,
                     "Astrometrically registering template to science image")
        exposureToConvolve = warpTemplateExposure(exposureToConvolve, exposureToNotConvolve, policy)
    
    maskedImageToConvolve = exposureToConvolve.getMaskedImage()
    maskedImageToNotConvolve = exposureToNotConvolve.getMaskedImage()

    if display:
        ds9.mtv(maskedImageToConvolve, frame=0)
        ds9.mtv(maskedImageToNotConvolve, frame=1)

    # Subtract their MaskedImages
    try:
        result = subtractMaskedImage(maskedImageToConvolve,
                                     maskedImageToNotConvolve,
                                     policy)
    except:
        pexLog.Trace("lsst.ip.diffim.subtractExposure", 1,
                     "ERROR: Unable to calculate psf matching kernel")
        raise
    else:
        differenceMaskedImage, spatialKernel, spatialBg, kernelCellSet = result
        
    if display:
        ds9.mtv(differenceMaskedImage, frame=2)
        
    
    # Generate an exposure from the results
    differenceExposure = afwImage.ExposureF(differenceMaskedImage, scienceWcs)


    return differenceExposure, spatialKernel, spatialBg, kernelCellSet




