# all the c++ level classes and routines
import diffimLib

# all the other LSST packages
import lsst.afw.image.imageLib as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.pex.policy as pexPolicy
import lsst.pex.logging as pexLog
import lsst.pex.exceptions as pexExcept
import lsst.meas.algorithms as measAlgorithms

display = True
if display:
    import lsst.afw.display.ds9 as ds9
    import diffimTools


# Most general routine
def psfMatchImageToImage(maskedImageToConvolve,
                         maskedImageToNotConvolve,
                         policy,
                         footprints = None,
                         returnOnExcept = False):

    # Reminder about background subtraction
    median1 = afwMath.makeStatistics(maskedImageToConvolve, afwMath.MEDIAN).getValue(afwMath.MEDIAN)
    median2 = afwMath.makeStatistics(maskedImageToNotConvolve, afwMath.MEDIAN).getValue(afwMath.MEDIAN)
    doBg    = policy.get("fitForBackground")
    if (not doBg) and abs(median1 - median2) > 10.0:
        pexLog.Trace("lsst.ip.diffim.psfMatchImageToImage", 1,
                     "WARNING: fitForBackground = False, but background1 = %.1f, background2 = %.1f" %
                     median1, median2)
        

    
    # Object to store the KernelCandidates for spatial modeling
    kernelCellSet = afwMath.SpatialCellSet(maskedImageToConvolve.getBBox(afwImage.PARENT),
                                           policy.getInt("sizeCellX"),
                                           policy.getInt("sizeCellY"))
    
    # Candidate source footprints to use for Psf matching
    if footprints == None:
        kcDetect = diffimLib.KernelCandidateDetectionF(policy.getPolicy("detectionPolicy"))
        kcDetect.apply(maskedImageToConvolve, maskedImageToNotConvolve)
        footprints = kcDetect.getFootprints()

    # Place candidate footprints within the spatial grid
    for fp in footprints:
        bbox = fp.getBBox()

        # Grab the centers in the parent's coordinate system
        xC   = 0.5 * ( bbox.getMinX() + bbox.getMaxX() )
        yC   = 0.5 * ( bbox.getMinY() + bbox.getMaxY() )

        tmi  = afwImage.MaskedImageF(maskedImageToConvolve, bbox, afwImage.PARENT)
        smi  = afwImage.MaskedImageF(maskedImageToNotConvolve, bbox, afwImage.PARENT)

        cand = diffimLib.makeKernelCandidate(xC, yC, tmi, smi, policy)

        pexLog.Trace("lsst.ip.diffim.psfMatchImageToImage", 7,
                     "Candidate %d at %f, %f" % (cand.getId(), cand.getXCenter(), cand.getYCenter()))
        
        kernelCellSet.insertCandidate(cand)


    # Create the Psf matching kernel
    try:
        kb = diffimLib.fitSpatialKernelFromCandidates(kernelCellSet, policy)
    except pexExcept.LsstCppException, e:
        pexLog.Trace("lsst.ip.diffim.psfMatchImageToImage", 1,
                     "ERROR: Unable to calculate psf matching kernel")
        pexLog.Trace("lsst.ip.diffim.psfMatchImageToImage", 2,
                     e.args[0].what())

        if returnOnExcept:
            return None, None, kernelCellSet
        else:
            raise
    else:
        spatialKernel = kb.first
        spatialBg     = kb.second

    # What is the status of the processing?
    nGood = 0
    for cell in kernelCellSet.getCellList():
        for cand in cell.begin(True):
            cand = diffimLib.cast_KernelCandidateF(cand)
            if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                nGood += 1
    if nGood == 0:
        pexLog.Trace("lsst.ip.diffim.psfMatchImageToImage", 1, "WARNING")
    pexLog.Trace("lsst.ip.diffim.psfMatchImageToImage", 1,
                 "Used %d kernels for spatial fit" % (nGood))

    return spatialKernel, spatialBg, kernelCellSet


# Requires images to have been pre-registered.  Once this happens,
# both the reference and the science image will have same XY0
def psfMatchModelToModel(referencePsfModel,
                         scienceBBox, sciencePsfModel,
                         policy, mergePolicy = False):

    if (referencePsfModel.getKernel().getDimensions() != sciencePsfModel.getKernel().getDimensions()):
        pexLog.Trace("lsst.ip.diffim.psfMatchModelToModel", 1,
                     "ERROR: Dimensions of reference Psf and science Psf different; exiting")
        raise RuntimeError, "ERROR: Dimensions of reference Psf and science Psf different; exiting"

    kernelWidth, kernelHeight = referencePsfModel.getKernel().getDimensions()
    maxPsfMatchingKernelSize = 1 + (min(kernelWidth - 1, kernelHeight - 1) // 2)
    if maxPsfMatchingKernelSize % 2 == 0:
        maxPsfMatchingKernelSize -= 1
    if policy.get('kernelSize') > maxPsfMatchingKernelSize:
        pexLog.Trace("lsst.ip.diffim.psfMatchModelToModel", 1,
                     "WARNING: Resizing matching kernel to size %d x %d" % (maxPsfMatchingKernelSize,
                                                                            maxPsfMatchingKernelSize))
        policy.set('kernelSize', maxPsfMatchingKernelSize)

    # Chanes to policy particular for matchPsfModels
    if mergePolicy:
        policyFile = pexPolicy.DefaultPolicyFile("ip_diffim", "MatchPsfModels.paf", "policy")
        matchPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)
        matchPolicy.mergeDefaults(policy.getDictionary())
        policy = matchPolicy
    
    regionSizeX, regionSizeY = scienceBBox.getDimensions()
    scienceX0,   scienceY0   = scienceBBox.getMin()

    sizeCellX = policy.get("sizeCellX")
    sizeCellY = policy.get("sizeCellY")

    kernelCellSet = afwMath.SpatialCellSet(afwGeom.Box2I(afwGeom.Point2I(scienceX0, scienceY0),
                                                         afwGeom.Extent2I(regionSizeX, regionSizeY)),
                                           sizeCellX, sizeCellY)

    nCellX    = regionSizeX // sizeCellX
    nCellY    = regionSizeY // sizeCellY

    dimenR    = referencePsfModel.getKernel().getDimensions()
    dimenS    = sciencePsfModel.getKernel().getDimensions()

    for row in range(nCellY):
        # place at center of cell
        posY = sizeCellY * row + sizeCellY // 2 + scienceY0
        
        for col in range(nCellX):
            # place at center of cell
            posX = sizeCellX * col + sizeCellX // 2 + scienceX0

            pexLog.Trace("lsst.ip.diffim.psfMatchModelToModel", 5,
                         "Creating Psf candidate at %.1f %.1f" % (posX, posY))

            # reference kernel image, at location of science subimage
            kernelImageR = referencePsfModel.computeImage(afwGeom.Point2D(posX, posY), True).convertF()
            sum = afwMath.makeStatistics(kernelImageR, afwMath.SUM).getValue(afwMath.SUM)
            kernelImageR /= sum
            kernelMaskR  = afwImage.MaskU(dimenR)
            kernelMaskR.set(0)
            kernelVarR   = afwImage.ImageF(dimenR)
            kernelVarR.set(0.01) # Total flux = 1, so this is order of magnitude
            referenceMI = afwImage.MaskedImageF(kernelImageR, kernelMaskR, kernelVarR)

            kernelImageS = sciencePsfModel.computeImage(afwGeom.Point2D(posX, posY), True).convertF()
            sum = afwMath.makeStatistics(kernelImageS, afwMath.SUM).getValue(afwMath.SUM)
            kernelImageS /= sum
            kernelMaskS  = afwImage.MaskU(dimenS)
            kernelMaskS.set(0)
            kernelVarS   = afwImage.ImageF(dimenS)
            kernelVarS.set(0.01) 
            scienceMI = afwImage.MaskedImageF(kernelImageS, kernelMaskS, kernelVarS)

            # The image to convolve is the science image, to the reference Psf.
            kc = diffimLib.makeKernelCandidate(posX, posY, scienceMI, referenceMI, policy)
            kernelCellSet.insertCandidate(kc)

    # Create the Psf matching kernel
    try:
        kb = diffimLib.fitSpatialKernelFromCandidates(kernelCellSet, policy)
    except pexExcept.LsstCppException, e:
        pexLog.Trace("lsst.ip.diffim.psfMatchModelToModel", 1,
                     "ERROR: Unable to calculate psf matching kernel")
        pexLog.Trace("lsst.ip.diffim.psfMatchModelToModel", 2,
                     e.args[0].what())
        raise
    else:
        spatialKernel = kb.first
        spatialBg     = kb.second    

    return spatialKernel, spatialBg, kernelCellSet

