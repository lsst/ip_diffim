# all the c++ level classes and routines
import diffimLib

# all the other diffim routines
from createKernelFunctor import createKernelFunctor

# all the other LSST packages
import lsst.afw.image.imageLib as afwImage
import lsst.afw.math.mathLib as afwMath
import lsst.pex.logging as pexLog
import lsst.pex.exceptions as pexExcept

# Most general routine
def createPsfMatchingKernel(maskedImageToConvolve,
                            maskedImageToNotConvolve,
                            policy,
                            footprints=None):

    
    # Object to store the KernelCandidates for spatial modeling
    kernelCellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(maskedImageToConvolve.getX0(),
                                                                         maskedImageToConvolve.getY0()),
                                                         maskedImageToConvolve.getWidth(),
                                                         maskedImageToConvolve.getHeight()),
                                           policy.getInt("sizeCellX"),
                                           policy.getInt("sizeCellY"))
    
    # Candidate source footprints to use for Psf matching
    if footprints == None:
        footprints = diffimLib.getCollectionOfFootprintsForPsfMatching(maskedImageToConvolve,
                                                                       maskedImageToNotConvolve,
                                                                       policy)

    # Place candidate footprints within the spatial grid
    for fp in footprints:
        bbox = fp.getBBox()
        xC   = 0.5 * ( bbox.getX0() + bbox.getX1() )
        yC   = 0.5 * ( bbox.getY0() + bbox.getY1() )
        tmi  = afwImage.MaskedImageF(maskedImageToConvolve,  bbox)
        smi  = afwImage.MaskedImageF(maskedImageToNotConvolve, bbox)
        
        cand = diffimLib.makeKernelCandidate(xC, yC, tmi, smi)
        kernelCellSet.insertCandidate(cand)

    # Object to perform the Psf matching on a source-by-source basis
    kFunctor = createKernelFunctor(policy)

    # Create the Psf matching kernel
    try:
        KB = diffimLib.fitSpatialKernelFromCandidates(kFunctor, kernelCellSet, policy)
    except pexExcept.LsstCppException, e:
        pexLog.Trace("lsst.ip.diffim.createPsfMatchingKernel", 1,
                     "ERROR: Unable to calculate psf matching kernel")
        pexLog.Trace("lsst.ip.diffim.createPsfMatchingKernel", 2,
                     e.args[0].what())
        raise
    else:
        spatialKernel = KB.first
        spatialBg     = KB.second

    # What is the status of the processing?
    nGood = 0
    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(True):
            cand = diffimLib.cast_KernelCandidateF(cand)
            if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                nGood += 1
    if nGood == 0:
        pexLog.Trace("lsst.ip.diffim.createPsfMatchingKernel", 1, "WARNING")
    pexLog.Trace("lsst.ip.diffim.createPsfMatchingKernel", 1,
                 "Used %d kernels for spatial fit" % (nGood))
                 
    return spatialKernel, spatialBg, kernelCellSet


# Matching an image to a Gaussian model
# This requires sending footprints since we can only accept stars
# In this case, maskedImageToConvolve is played by the science image
def createPsfMatchingKernelToGaussian(maskedImageToConvolve,
                                      sigGauss,
                                      policy,
                                      footprints):

    # Idealized function you will match all the stars to
    gaussFunc   = afwMath.GaussianFunction2D(sigGauss, sigGauss)
    
    # Object to store the KernelCandidates for spatial modeling
    kernelCellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(maskedImageToConvolve.getX0(),
                                                                         maskedImageToConvolve.getY0()),
                                                         maskedImageToConvolve.getWidth(),
                                                         maskedImageToConvolve.getHeight()),
                                           policy.getInt("sizeCellX"),
                                           policy.getInt("sizeCellY"))

    # Place candidate footprints within the spatial grid
    for fp in footprints:
        bbox = fp.getBBox()
        xC   = 0.5 * ( bbox.getX0() + bbox.getX1() )
        yC   = 0.5 * ( bbox.getY0() + bbox.getY1() )
        tmi  = afwImage.MaskedImageF(maskedImageToConvolve,  bbox)

        # A bit of a hack, but to turn this function into an image I
        # need to go through kernel
        gaussKernel = afwMath.AnalyticKernel(tmi.getWidth(), tmi.getHeight(), gaussFunc)
        kImage      = afwImage.ImageD(tmi.getWidth(), tmi.getHeight())

        # Do we need to normalize the kernel here?  No; want the peak
        # value the same in all realizations of it
        gaussKernel.computeImage(kImage, False)

        # Create a masked image that your input image will be matched to
        # Empty mask and zero variance
        smi  = afwImage.MaskedImageF(kImage)
        
        cand = diffimLib.makeKernelCandidate(xC, yC, tmi, smi)
        kernelCellSet.insertCandidate(cand)

    # Object to perform the Psf matching on a source-by-source basis
    kFunctor = createKernelFunctor(policy)

    # Create the Psf matching kernel
    try:
        KB = diffimLib.fitSpatialKernelFromCandidates(kFunctor, kernelCellSet, policy)
    except pexExcept.LsstCppException, e:
        pexLog.Trace("lsst.ip.diffim.createPsfMatchingKernel", 1,
                     "ERROR: Unable to calculate psf matching kernel")
        pexLog.Trace("lsst.ip.diffim.createPsfMatchingKernel", 2,
                     e.args[0].what())
        raise
    else:
        spatialKernel = KB.first
        spatialBg     = KB.second

    # What is the status of the processing?
    nGood = 0
    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(True):
            cand = diffimLib.cast_KernelCandidateF(cand)
            if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                nGood += 1
    if nGood == 0:
        pexLog.Trace("lsst.ip.diffim.createPsfMatchingKernel", 1, "WARNING")
    pexLog.Trace("lsst.ip.diffim.createPsfMatchingKernel", 1,
                 "Used %d kernels for spatial fit" % (nGood))
                 
    return spatialKernel, spatialBg, kernelCellSet


# Specialized routines where I tweak the policy based on what you want done
def createMeanPsfMatchingKernel(maskedImageToConvolve,
                                maskedImageToNotConvolve,
                                policy):

    policy.set("spatialKernelOrder", 0)
    policy.set("singleKernelClipping", True)
    policy.set("kernelSumClipping", True)
    policy.set("spatialKernelClipping", False)

    return createPsfMatchingKernel(maskedImageToConvolve, maskedImageToNotConvolve, policy)
