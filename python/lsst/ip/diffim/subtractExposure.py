from subtractMaskedImage import subtractMaskedImage
from warpTemplateExposure import warpTemplateExposure

import lsst.pex.logging as pexLog
import lsst.pex.exceptions as pexExcept
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9

def subtractExposure(templateExposure, scienceExposure, policy, display=False):
    # Make sure they end up the same dimensions on the sky
    templateWcs    = templateExposure.getWcs() 
    scienceWcs     = scienceExposure.getWcs()

    # LLC
    templateOrigin = templateWcs.xyToRaDec(0,0)
    scienceOrigin  = scienceWcs.xyToRaDec(0,0)
    # URC
    templateLimit  = templateWcs.xyToRaDec(templateExposure.getWidth(),
                                           templateExposure.getHeight())
    scienceLimit   = scienceWcs.xyToRaDec(scienceExposure.getWidth(),
                                          scienceExposure.getHeight())

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
         (templateExposure.getHeight() != scienceExposure.getHeight()) or \
         (templateExposure.getWidth()  != scienceExposure.getWidth()) ):
        pexLog.Trace("lsst.ip.diffim.subtractExposure", 1,
                     "Astrometrically registering template to science image")
        templateExposure = warpTemplateExposure(templateExposure, scienceExposure, policy)
    
    templateMaskedImage = templateExposure.getMaskedImage()
    scienceMaskedImage  = scienceExposure.getMaskedImage()

    if display:
        ds9.mtv(templateMaskedImage, frame=0)
        ds9.mtv(scienceMaskedImage, frame=1)

    # Subtract their MaskedImages
    try:
        result = subtractMaskedImage(templateMaskedImage,
                                     scienceMaskedImage,
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




