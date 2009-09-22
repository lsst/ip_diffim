import lsst.afw.math as afwMath

def warpTemplateExposure(templateExposure, scienceExposure, policy):
    # The destination Wcs is in scienceExposure
    # Create the warping Kernel according to policy
    warpingKernelName = policy.getString("warpingKernelName")
    warpingKernel     = afwMath.makeWarpingKernel(warpingKernelName)

    # create a blank exposure to hold the remaped template exposure
    remapedTemplateExposure = templateExposure.Factory(
                scienceExposure.getWidth(), 
                scienceExposure.getHeight(),
                scienceExposure.getWcs())
    scienceMaskedImage = scienceExposure.getMaskedImage()
    remapedMaskedImage = remapedTemplateExposure.getMaskedImage()
    remapedMaskedImage.setXY0(scienceMaskedImage.getXY0())

    # warp the template exposure
    afwMath.warpExposure(remapedTemplateExposure, 
                         templateExposure, 
                         warpingKernel)
        
    return remapedTemplateExposure
    
