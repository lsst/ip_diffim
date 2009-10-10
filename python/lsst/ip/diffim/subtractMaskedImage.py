import time

# all the c++ level classes and routines
import diffimLib

# all the other diffim routines
from createPsfMatchingKernel import createPsfMatchingKernel

# all the other LSST packages
import lsst.afw.image as afwImage
import lsst.pex.logging as pexLog
import lsst.sdqa as sdqa

def subtractMaskedImage(templateMaskedImage, 
                        scienceMaskedImage, 
                        policy, 
                        footprints=None):
    
    # Make sure they are the EXACT same dimensions in pixels
    # This is non-negotiable
    assert (templateMaskedImage.getDimensions() == \
            scienceMaskedImage.getDimensions())

    # We also assume that at this stage, they are aligned at the pixel level
    # Assign to the coordinate system of the science image
    templateMaskedImage.setXY0(scienceMaskedImage.getXY0())

    
    spatialKernel, spatialBg, kernelCellSet = createPsfMatchingKernel(templateMaskedImage,
                                                                      scienceMaskedImage,
                                                                      policy,
                                                                      footprints)
    
    # no need to subtract a background in subtraction as we'll do so in a moment
    if policy.exists("backgroundPolicy"):
        background = 0.                  
    else:
        background = spatialBg

    t0 = time.time()
    differenceMaskedImage = diffimLib.convolveAndSubtract(templateMaskedImage,
                                                          scienceMaskedImage,
                                                          spatialKernel,
                                                          background)
    t1 = time.time()
    pexLog.Trace("lsst.ip.diffim.subtractMaskedImage", 1,
                 "Total time for final convolve and subtract : %.2f s" % (t1-t0))

    #
    # Maybe subtract a background model from the difference image
    #
    if policy.exists("backgroundPolicy"):
        algorithm = policy.get("backgroundPolicy.algorithm")
        binsize   = policy.get("backgroundPolicy.binsize")

        if algorithm == "NATURAL_SPLINE":
            bctrl = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE)
        else:
            raise RuntimeError, "Unknown backgroundPolicy.algorithm: %s" % (algorithm)

        bctrl.setNxSample(int(differenceMaskedImage.getWidth()//binsize) + 1)
        bctrl.setNySample(int(differenceMaskedImage.getHeight()//binsize) + 1)
        
        image   = differenceMaskedImage.getImage() 
        backobj = afwMath.makeBackground(image, bctrl)
        image  -= backobj.getImageF()
        del image; del backobj

    #
    # Place holder for Sqda on diffim; diffim stats and kernel sum
    #

    return differenceMaskedImage, spatialKernel, spatialBg, kernelCellSet
