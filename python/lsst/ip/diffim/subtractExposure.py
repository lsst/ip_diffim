from subtractMaskedImage import subtractMaskedImage
import lsst.afw.image as afwImage

def subtractExposure(templateExposure, scienceExposure, policy, log, invert=False):
    # Make sure they end up the same dimensions on the sky
    templateWcs = templateExposure.getWcs() 
    scienceWcs = scienceExposure.getWcs()

    templateMaskedImage = templateExposure.getMaskedImage()
    scienceMaskedImage = scienceExposure.getMaskedImage()

    templateOrigin = templateWcs.xyToRaDec(0,0)
    scienceOrigin = scienceWcs.xyToRaDec(0,0)

    # Within some tolerance; do we have sky distance methods?
    #assert(templateOrigin[0] == scienceOrigin[0])
    #assert(templateOrigin[1] == scienceOrigin[1])
    assert(templateOrigin == scienceOrigin)

    templateLimit = templateWcs.xyToRaDec(templateMaskedImage.getHeight(),
            templateMaskedImage.getWidth())
    scienceLimit = scienceWcs.xyToRaDec(scienceMaskedImage.getHeight(),
            scienceMaskedImage.getWidth())

    # Within some tolerance; do we have sky distance methods?
    #assert(templateLimit[0]  == scienceLimit[0])
    #assert(templateLimit[1]  == scienceLimit[1])
    assert(templateLimit == scienceLimit)

    # Subtract their MaskedImages
    differenceMaskedImage, spatialKernel, backgroundModel, sdqaList = \
                           subtractMaskedImage(templateMaskedImage,
                                               scienceMaskedImage,
                                               policy,
                                               log, invert=invert)

    # Note : we assume that the Template is warped to the science image's WCS
    #      : meaning that the scienceWcs is the correct one to store in 
    #      : the diffim
    differenceExposure = afwImage.ExposureF(differenceMaskedImage, scienceWcs)

    return differenceExposure, spatialKernel, backgroundModel, sdqaList




