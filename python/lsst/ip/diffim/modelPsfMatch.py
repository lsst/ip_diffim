
class ModelPsfMatchConfig(PsfMatchConfig):

class ModelPsfMatch(PsfMatch):
    """PSF-match PSF models to reference PSF models
    """
    def __init__(self, policy, logName="lsst.ip.diffim.ModelPsfMatch", mergePolicy = False):
        """Create a PsfMatchToModel
        
        @param policy: see lsst/ip/diffim/policy/PsfMatchingDictionary.paf
        @param logName: name by which messages are logged
        @param mergePolicy: merge ip_diffim/policy/MatchPsfModels.paf into existing policy;
            this is a temporary hack
        """
        PsfMatch.__init__(self, policy, logName)

        # Changes to policy particular to matchPsfModels
        if mergePolicy:
            policyFile = pexPolicy.DefaultPolicyFile("ip_diffim", "MatchPsfModels.paf", "policy")
            matchPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)
            matchPolicy.mergeDefaults(self._policy.getDictionary())
            self._policy = matchPolicy
    
    def matchExposure(self, exposure, referencePsfModel, kernelSum=1.0):
        """PSF-match an exposure to a PSF model.
        
        @param exposure: Exposure to PSF-match to the reference masked image;
            must contain a PSF model and must be warped to match the reference masked image
        @param referencePsfModel: PSF model to match (an afwDetection.Psf)
        @param kernelSum: A multipicative factor reflecting the difference in 
            zeropoints between the images; kernelSum = zpt(science) / zpt(ref)
        
        @return
        - psfMatchedExposure: the PSF-matched exposure.
            This has the same parent bbox, Wcs, Calib and Filter as exposure but no psf.
            In theory the psf should equal referencePsfModel but the match is likely not exact.
        - psfMatchingKernel: the PSF matching kernel
        - kernelCellSet: SpatialCellSet used to solve for the PSF matching kernel
        
        @raise RuntimeError if exposure does not contain a PSF model"
        """
        if not exposure.hasPsf():
            raise RuntimeError("exposure does not contain a PSF model")

        maskedImage = exposure.getMaskedImage()

        self._log.log(pexLog.Log.INFO, "compute PSF-matching kernel")
        kernelCellSet = self._buildCellSet(referencePsfModel,
                                           exposure.getBBox(afwImage.PARENT),
                                           exposure.getPsf(), kernelSum)
        psfMatchingKernel, backgroundModel = self._solve(kernelCellSet)
        
        self._log.log(pexLog.Log.INFO, "PSF-match science exposure to reference")
        psfMatchedExposure = afwImage.ExposureF(exposure.getBBox(afwImage.PARENT), exposure.getWcs())
        psfMatchedExposure.setFilter(exposure.getFilter())
        psfMatchedExposure.setCalib(exposure.getCalib())
        psfMatchedMaskedImage = psfMatchedExposure.getMaskedImage()

        # Normalize the psf-matching kernel while convolving since its magnitude is meaningless
        # when PSF-matching one model to another.
        doNormalize = True
        afwMath.convolve(psfMatchedMaskedImage, maskedImage, psfMatchingKernel, doNormalize)
        
        self._log.log(pexLog.Log.INFO, "done")
        return (psfMatchedExposure, psfMatchingKernel, kernelCellSet)

    def _buildCellSet(self, referencePsfModel, scienceBBox, sciencePsfModel, kernelSum = 1.0):
        """Build a SpatialCellSet for use with the solve method

        @param referencePsfModel: PSF model to match (an afwDetection.Psf)
        @param scienceBBox: parent bounding box on science image
        @param sciencePsfModel: PSF model for science image
        
        @return kernelCellSet: a SpatialCellSet for use with self._solve
        
        @raise RuntimeError if reference PSF model and science PSF model have different dimensions
        """
        if (referencePsfModel.getKernel().getDimensions() != sciencePsfModel.getKernel().getDimensions()):
            pexLog.Trace(self._log.getName(), 1,
                         "ERROR: Dimensions of reference Psf and science Psf different; exiting")
            raise RuntimeError, "ERROR: Dimensions of reference Psf and science Psf different; exiting"
    
        kernelWidth, kernelHeight = referencePsfModel.getKernel().getDimensions()
        maxPsfMatchingKernelSize = 1 + (min(kernelWidth - 1, kernelHeight - 1) // 2)
        if maxPsfMatchingKernelSize % 2 == 0:
            maxPsfMatchingKernelSize -= 1
        if self._policy.get('kernelSize') > maxPsfMatchingKernelSize:
            pexLog.Trace(self._log.getName(), 1,
                         "WARNING: Resizing matching kernel to size %d x %d" % (maxPsfMatchingKernelSize,
                                                                                maxPsfMatchingKernelSize))
            self._policy.set('kernelSize', maxPsfMatchingKernelSize)
        
        regionSizeX, regionSizeY = scienceBBox.getDimensions()
        scienceX0,   scienceY0   = scienceBBox.getMin()
    
        sizeCellX = self._policy.get("sizeCellX")
        sizeCellY = self._policy.get("sizeCellY")
    
        kernelCellSet = afwMath.SpatialCellSet(
            afwGeom.Box2I(afwGeom.Point2I(scienceX0, scienceY0),
                          afwGeom.Extent2I(regionSizeX, regionSizeY)),
            sizeCellX, sizeCellY
            )
    
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
    
                pexLog.Trace(self._log.getName(), 5, "Creating Psf candidate at %.1f %.1f" % (posX, posY))
    
                # reference kernel image, at location of science subimage
                kernelImageR = referencePsfModel.computeImage(afwGeom.Point2D(posX, posY), True).convertF()
                imsum = afwMath.makeStatistics(kernelImageR, afwMath.SUM).getValue(afwMath.SUM)
                kernelImageR /= imsum         # image sums to 1.0 
                kernelImageR /= kernelSum     # image sums to 1/kernelSum
                kernelMaskR   = afwImage.MaskU(dimenR)
                kernelMaskR.set(0)
                kernelVarR    = afwImage.ImageF(dimenR)
                kernelVarR.set(1.0)
                referenceMI   = afwImage.MaskedImageF(kernelImageR, kernelMaskR, kernelVarR)
     
                kernelImageS = sciencePsfModel.computeImage(afwGeom.Point2D(posX, posY), True).convertF()
                imsum = afwMath.makeStatistics(kernelImageS, afwMath.SUM).getValue(afwMath.SUM)
                kernelImageS /= imsum
                kernelMaskS   = afwImage.MaskU(dimenS)
                kernelMaskS.set(0)
                kernelVarS    = afwImage.ImageF(dimenS)
                kernelVarS.set(1.0)
                scienceMI     = afwImage.MaskedImageF(kernelImageS, kernelMaskS, kernelVarS)

                #referenceMI.writeFits('ref_%d_%d.fits' % (row, col))
                #scienceMI.writeFits('sci_%d_%d.fits' % (row, col))
 
                # The image to convolve is the science image, to the reference Psf.
                kc = diffimLib.makeKernelCandidate(posX, posY, scienceMI, referenceMI, self._policy)
                kernelCellSet.insertCandidate(kc)

        return kernelCellSet