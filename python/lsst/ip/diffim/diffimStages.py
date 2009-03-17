from lsst.pex.harness.Stage import Stage
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import diffimLib as ipDiffim

class DiffimStage(Stage):
    def process(self):
        self.activeClipboard = self.inputQueue.getNextDataset()
        
        scienceExposureKey = self._policy.get('scienceExposureKey')
        templateExposureKey = self._policy.get('templateExposureKey')
        
        scienceExposure = self.activeClipboard.get(scienceExposureKey)
        templateExposure = self.activeClipboard.get(templateExposureKey)
       
        diffimPolicy = self._policy.get('diffimPolicy')
        # step 1
        remapedTemplateExposure = warpTemplateExposure(templateExposure,
                scienceExposure, 
                diffimPolicy)
        
        # step 2
        products = subtractExposure(remapedTemplateExposure, 
                scienceExposure, 
                diffimPolicy)

        if products == None:
            raise RuntimeException("DiffimStage.subtractExposure failed")

        differenceExposure, spatialKernel, backgroundModel, sdqaSet = products

        exposureKey = self._policy.getString('differenceExposureKey')
        self.activeClipboard.put(exposureKey, differenceExposure)

        clipboardKey = self._policy.getString('sdqaRatingSetKey')
        self.activeClipboard.put(sdqaKey, sdqaSet)

        self.outputQueue.addDataset(self.activeClipboard)

def warpTemplateExposure(templateExposure, scienceExposure, policy):
    # Create the warping Kernel according to policy
    warpingKernelOrder = policy.getInt("warpingKernelOrder")
    warpingKernel = afwMath.LanczosWarpingKernel(warpingKernelOrder)

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
    
def subtractExposure(templateExposure, scienceExposure, policy):
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
                    policy)

    # Note : we assume that the Template is warped to the science image's WCS
    #      : meaning that the scienceWcs is the correct one to store in 
    #      : the diffim
    differenceExposure = afwImage.ExposureF(differenceMaskedImage, scienceWcs)

    return differenceExposure, spatialKernel, backgroundModel, sdqaList



def subtractMaskedImage(templateMaskedImage, 
        scienceMaskedImage, 
        policy, 
        fpList=None):
    # Make sure they are the EXACT same dimensions in pixels
    # This is non-negotiable
    assert (templateMaskedImage.getDimensions() == \
            scienceMaskedImage.getDimensions())
    
    kCols = policy.get('kernelCols')
    kRows = policy.get('kernelRows')

    kBasisList = ipDiffim.generateDeltaFunctionKernelSet(kCols, kRows)
    kFunctor   = ipDiffim.PsfMatchingFunctorF(kBasisList)

    if fpList == None:
        # Need to find own footprints
        fpList = ipDiffim.getCollectionOfFootprintsForPsfMatching(
                templateMaskedImage,
                scienceMaskedImage,
                policy)

    # Set up grid for spatial model
    spatialCells = ipDiffim.createSpatialModelKernelCells(templateMaskedImage,
            scienceMaskedImage,
            fpList,
            kFunctor,
            policy)

    # Set up fitting loop 
    maxSpatialIterations = policy.getInt('maxSpatialIterations')
    rejectKernels = policy.getBool('spatialKernelRejection')
    nRejected = -1
    nIter =  0
    
    # And fit spatial kernel model
    if policy.get('spatialKernelModel') == 'pca':
        # Fit spatial variation of principal components

        minPrincipalComponents = policy.getInt('minPrincipalComponents')
        maxPrincipalComponents = policy.getInt('maxPrincipalComponents')
        fracEigenVal = policy.getDouble('fracEigenVal')
        
        while (nRejected != 0) and (nIter < maxSpatialIterations):
            # Run the PCA
            mKernel, eKernelVector, eVal, eCoeff = \
                    ipDiffim.spatialModelKernelPca(spatialCells, policy)

            # Make the decision on how many components to use
            eFrac = numpy.cumsum(eVal)
            eFrac /= eFrac[-1]
            nEval = len(numpy.where(eFrac < fracEigenVal)[0])
            nEval = min(nEval, maxPrincipalComponents)
            nEval = max(nEval, minPrincipalComponents)

            # do spatial fit here by Principal Component
            sKernel, bgFunction = ipDiffim.spatialModelByPca(spatialCells,
                    mKernel,
                    eKernelVector,
                    eCoeff,
                    nEval,
                    policy)

            # Evaluate quality of spatial fit
            nRejected, sdqaList = ipDiffim.evaluateModelByPca(spatialCells,
                    bgFunction, 
                    sKernel,
                    policy, 
                    reject=rejectKernels)
               
            nIter += 1

    elif policy.get('spatialKernelModel') == 'pixel':
        # Fit function to each pixel

        while (nRejected != 0) and (nIter < maxSpatialIterations):
            # do spatial fit here pixel by pixel
            sKernel, bgFunction = ipDiffim.spatialModelByPixel(spatialCells,
                    kBasisList,
                    policy)

            # and check quality
            nRejected, sdqaList  = ipDiffim.evaluateModelByPixel(spatialCells,
                    bgFunction, 
                    sKernel, 
                    policy, 
                    reject=rejectKernels)

            nIter += 1
        
    else:
        # All that is supported
        # Throw exception!
        pass

    differenceMaskedImage = ipDiffim.convolveAndSubtract(templateMaskedImage,
            scienceMaskedImage,
            sKernel,
            bgFunction)

    #
    # Lets do some more Sdqa here
    #
    imStats = ipDiffim.ImageStatisticsF()
    imStats.apply(differenceMaskedImage)
    residualsRating = sdqa.SdqaRating("ip_diffim.residuals",
            imStats.getMean(), 
            imStats.getRms(), 
            sdqa.SdqaRating.AMP) 
    sdqaList.append(residualsRating)

    # Check kernel sum in the corners
    kSums  = []
    kImage = afwImage.ImageD(sKernel.getDimensions())
    for nRow in [0, templateMaskedImage.getHeight()]:
        for nCol in [0, templateMaskedImage.getWidth()]:
            kSums.append( sKernel.computeImage(kImage, False, nCol, nRow) )
    kSumArray = numpy.array(kSums)
    Trace('lsst.ip.diffim.subtractMaskedImage', 3, 
            'Final Kernel Sum from Image Corners : %0.3f (%0.3f)' % 
            (kSumArray.mean(), kSumArray.std()))
    kernelSumRating =  sdqa.SdqaRating("ip_diffim.kernelSum",
            kSumArray.mean(),
            kSumArray.std(), 
            sdqa.SdqaRating.AMP)  
    sdqaList.append(kernelSumRating)

    # What kind of metadata do we add here to MaskedImage?
    
    return differenceMaskedImage, sKernel, bgFunction, sdqaList

