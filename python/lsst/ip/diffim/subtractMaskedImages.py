import time

# all the c++ level classes and routines
import diffimLib

# all the other diffim routines
from .makePsfMatchingKernel import makePsfMatchingKernel

# all the other LSST packages
import lsst.afw.math as afwMath
import lsst.pex.logging as pexLog

def subtractMaskedImages(maskedImageToConvolve, 
                         maskedImageToNotConvolve, 
                         policy, 
                         footprints=None):
    
    # Make sure they are the EXACT same dimensions in pixels
    # This is non-negotiable
    assert (maskedImageToConvolve.getDimensions() == \
            maskedImageToNotConvolve.getDimensions())

    # We also assume that at this stage, they are aligned at the pixel level
    # Assign to the coordinate system of the science image
    maskedImageToConvolve.setXY0(maskedImageToNotConvolve.getXY0())


    try:
        result = makePsfMatchingKernel(maskedImageToConvolve,
                                       maskedImageToNotConvolve,
                                       policy,
                                       footprints)
    except:
        pexLog.Trace("lsst.ip.diffim.subtractMaskedImage", 1,
                     "ERROR: Unable to calculate psf matching kernel")
        raise
    else:
        spatialKernel, spatialBg, kernelCellSet = result
        
    # no need to subtract a background in subtraction as we'll do so in a moment
    if policy.get("useAfwBackground"):
        background = 0.                  
    else:
        background = spatialBg

    t0 = time.time()
    differenceMaskedImage = diffimLib.convolveAndSubtract(maskedImageToConvolve,
                                                          maskedImageToNotConvolve,
                                                          spatialKernel,
                                                          background)
    t1 = time.time()
    pexLog.Trace("lsst.ip.diffim.subtractMaskedImage", 1,
                 "Total time for final convolve and subtract : %.2f s" % (t1-t0))

    #
    # Instead subtract afw's background model from the difference image
    #
    if policy.get("useAfwBackground"):
        algorithm   = policy.get("backgroundPolicy.algorithm")
        binsize     = policy.get("backgroundPolicy.binsize")
        undersample = policy.get("backgroundPolicy.undersample")
        bctrl       = afwMath.BackgroundControl(algorithm)
        bctrl.setNxSample(differenceMaskedImage.getWidth()//binsize + 1)
        bctrl.setNySample(differenceMaskedImage.getHeight()//binsize + 1)
        bctrl.setUndersampleStyle(undersample)

        image   = differenceMaskedImage.getImage() 
        backobj = afwMath.makeBackground(image, bctrl)
        image  -= backobj.getImageF()
        del image
        del backobj

    #
    # Place holder for Sqda on diffim; diffim stats and kernel sum
    #

    return differenceMaskedImage, spatialKernel, spatialBg, kernelCellSet
