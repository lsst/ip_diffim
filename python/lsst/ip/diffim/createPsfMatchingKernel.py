import lsst.ip.diffim as ipDiffim
import lsst.afw.image.imageLib as afwImage
import lsst.afw.math.mathLib as afwMath

# Most general routine
def createPsfMatchingKernel(templateMaskedImage,
                            scienceMaskedImage,
                            policy):

    # Object to store the KernelCandidates for spatial modeling
    kernelCellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(templateMaskedImage.getX0(),
                                                                         templateMaskedImage.getY0()),
                                                         templateMaskedImage.getWidth(),
                                                         templateMaskedImage.getHeight()),
                                           policy.getInt("sizeCellX"),
                                           policy.getInt("sizeCellY"))
    
    # Object to perform the Psf matching on a source-by-source basis
    kFunctor = ipDiffim.createKernelFunctor(policy)

    # Candidate source footprints to use for Psf matching
    footprints = ipDiffim.getCollectionOfFootprintsForPsfMatching(templateMaskedImage,
                                                                  scienceMaskedImage,
                                                                  policy)

    # Place candidate footprints within the spatial grid
    for fp in footprints:
        bbox = fp.getBBox()
        xC   = 0.5 * ( bbox.getX0() + bbox.getX1() )
        yC   = 0.5 * ( bbox.getY0() + bbox.getY1() )
        tmi  = afwImage.MaskedImageF(templateMaskedImage,  bbox)
        smi  = afwImage.MaskedImageF(scienceMaskedImage, bbox)
        
        cand = ipDiffim.makeKernelCandidate(xC, yC, tmi, smi)
        kernelCellSet.insertCandidate(cand)

    # Create the Psf matching kernel
    spatialKernel, spatialBg = ipDiffim.fitSpatialKernelFromCandidates(kFunctor, kernelCellSet, policy)

    return spatialKernel, spatialBg


# Specialized routines where I tweak the policy based on what you want done
def createMeanPsfMatchingKernel(templateMaskedImage,
                                scienceMaskedImage,
                                policy):

    policy.set("spatialKernelOrder", 0)
    policy.set("singleKernelClipping", True)
    policy.set("kernelSumClipping", True)
    policy.set("spatialKernelClipping", False)

    return createPsfMatchingKernel(templateMaskedImage, scienceMaskedImage, policy)
