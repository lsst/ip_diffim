from subtractMaskedImage import subtractMaskedImage
from warpTemplateExposure import warpTemplateExposure
import lsst.afw.image as afwImage

def subtractExposure(templateExposure, scienceExposure, policy):
    # Make sure they end up the same dimensions on the sky
    templateWcs    = templateExposure.getWcs() 
    scienceWcs     = scienceExposure.getWcs()

    # LLC
    templateOrigin = templateWcs.xyToRaDec(0,0)
    scienceOrigin  = scienceWcs.xyToRaDec(0,0)
    # URC
    templateLimit  = templateWcs.xyToRaDec(templateMaskedImage.getHeight(),
                                           templateMaskedImage.getWidth())
    scienceLimit   = scienceWcs.xyToRaDec(scienceMaskedImage.getHeight(),
                                          scienceMaskedImage.getWidth())

    # Remap if they do not line up in size on the sky, and size in pixels
    if ( (templateOrigin != scienceOrigin) or \
         (templateLimit  != scienceLimit)  or \
         (templateMaskedImage.getHeight() != scienceMaskedImage.getHeight()) or \
         (templateMaskedImage.getWidth()  != scienceMaskedImage.getWidth()) ):
        templateExposure = warpTemplateExposure(templateExposure, scienceExposure, policy)
    
    templateMaskedImage = templateExposure.getMaskedImage()
    scienceMaskedImage  = scienceExposure.getMaskedImage()

    # Subtract their MaskedImages
    differenceMaskedImage, spatialKernel, spatialBg = \
                           subtractMaskedImage(templateMaskedImage,
                                               scienceMaskedImage,
                                               policy)

    # Generate an exposure from the results
    differenceExposure = afwImage.ExposureF(differenceMaskedImage, scienceWcs)

    return differenceExposure, spatialKernel, spatialBg




