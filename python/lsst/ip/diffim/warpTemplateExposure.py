import time
import lsst.afw.math as afwMath
import lsst.pex.logging as pexLog

def warpTemplateExposure(templateExposure, scienceExposure, policy):
    t0 = time.time()

    # The destination Wcs is in scienceExposure
    # Create the warping Kernel according to policy
    warpingKernelName = policy.getString("warpingKernelName")
    warpingKernel     = afwMath.makeWarpingKernel(warpingKernelName)

    # create a blank exposure to hold the remapped template exposure
    remappedTemplateExposure = templateExposure.Factory(
        scienceExposure.getWidth(), 
        scienceExposure.getHeight(),
        scienceExposure.getWcs())
    scienceMaskedImage = scienceExposure.getMaskedImage()
    remappedMaskedImage = remappedTemplateExposure.getMaskedImage()
    remappedMaskedImage.setXY0(scienceMaskedImage.getXY0())

    # warp the template exposure
    afwMath.warpExposure(remappedTemplateExposure, 
                         templateExposure, 
                         warpingKernel)
        
    t1 = time.time()
    pexLog.Trace("lsst.ip.diffim.warpTemplateExposure", 1,
                 "Total time for warping : %.2f s" % (t1-t0))

    return remappedTemplateExposure
    
