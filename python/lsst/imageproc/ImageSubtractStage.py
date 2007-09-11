#! /usr/bin/env python

from Stage import Stage
import lsst.fw.Core.fwLib as fw
import lsst.imageproc.Core.imageprocLib as imageproc

class ImageSubtractStage(Stage):
    #------------------------------------------------------------------------
    def preprocess(self):
        self.activeClipboard = self.inputQueue.getNextDataset()

        #keys = self.activeClipboard.getKeys()
        keys = ['Policy', ]

        # get policy
        for key in keys:
            value = self.activeClipboard.get(key)
            print 'Python ImageSubtractStage preprocess(): stageId %i key %s value %s' % (self.stageId, key, value)

        # from policy generate kernelBasisSet
        # e.g. imageproc.generateDeltaFunctionKernelSet()

        #value = self.activeClipboard.put('BasisSet', kernelBasisSet)
        self.activeClipboard['BasisSet'] = kernelBasisSet

    #------------------------------------------------------------------------
    def processA(self):
        print 'Python ImageSubtractStage process : _rank %i stageId %d' % (self._rank, self.stageId)
        self.activeClipboard = self.inputQueue.getNextDataset()

        templateImage = self.activeClipboard['TemplateImage']
        inputImage    = self.activeClipboard['InputImage']
        basisSet      = self.activeClipboard['BasisSet']
        policy        = self.activeClipboard['Policy']

        # generate kernelPtr, kernelFunctionPtr, and backgroundFunctionPtr from Policy
        
        imageproc.computePsfMatchingKernelForMaskedImage(templateImage, inputImage, basisSet,
                                                         kernelPtr, kernelFunctionPtr, backgroundFunctionPtr, policy)

        differenceImage = inputImage - fw.convolve(templateImage, kernelPtr)

        # Ptr or not Ptr, that is the question.  The answer is Not
        
        self.activeClipboard['DifferenceImage'] = differenceImage
        self.activeClipboard['OutputKernel'] = kernel
        self.activeClipboard['BackgroundModel'] = backgroundFunction

    #------------------------------------------------------------------------
    def processB(self):
        print 'Python ImageSubtractStage process : _rank %i stageId %d' % (self._rank, self.stageId)
        self.activeClipboard = self.inputQueue.getNextDataset()

        ###########
        #
        # Get objects from clipboard
        #
        templateImage = self.activeClipboard['TemplateImage']
        inputImage    = self.activeClipboard['InputImage']
        basisSet      = self.activeClipboard['BasisSet']
        policy        = self.activeClipboard['Policy']

        ###########
        #
        # Get directives from policy
        #
        edgeMaskBit = policy.get('edgeMaskBit')
        convolveThreshold = policy.get('convolveThreshold')
        kernelRows = policy.get('kernelRows')
        kernelCols = policy.get('kernelCols')
        kernelSpatialOrder = policy.get('kernelSpatialOrder')
        backgroundSpatialOrder = policy.get('backgroundSpatialOrder')

        ###########
        #
        # Generate objects from policy directives
        #
        kernelBasisVec = imageproc.generateDeltaFunctionKernelSet(kernelRows, kernelCols)
        kernelPtr = fw.LinearCombinationKernel()
        kernelFunctionPtr = fw.PolynomialFunction2(kernelSpatialOrder)
        backgroundFunctionPtr = fw.PolynomialFunction2(backgroundSpatialOrder)

        ###########
        #
        # Get good footprints
        #
        if self.activeClipboard.has_key('FootprintList'):
            footprintList = self.activeClipboard['FootprintList']
        else:
            footprintList = imageproc.getCollectionOfFootprintsForPsfMatching(templateImage, intputImage, policy)

        ###########
        #
        # Calculate all individual Kernels
        #
        diffImContainerList = []
        nFootprint = 0
        for iFootprint in footprintList:
            footprintBBox = iFootprint.getBBox()

            imageToConvolveStampPtr = imageToConvolve.getSubImage(footprintBBox)
            imageToNotConvolveStampPtr = imageToNotConvolve.getSubImage(footprintBBox)

            if DEBUG_IO:
                imageToConvolveStampPtr.writeFits('csFits_%d' % (nFootprint))
                imageToNotConvolveStampPtr.writeFits('ncsFits_%d' % (nFootprint))

            kernelCoeffs, background = imageproc.computePsfMatchingKernelForPostageStamp(imageToConvolveStampPtr,
                                                                                         imageToNotConvolveStampPtr,
                                                                                         basisSet)
            footprintKernelPtr = fw.LinearCombinationKernel(basisSet, kernelCoeffs)

            diffImFootprintContainer = imageproc.DiffImContainer()
            diffImFootprintContainer.id = nFootprint
            diffImFootprintContainer.isGood = True
            diffImFootprintContainer.diffImFootprintPtr = iFootprint
            diffImFootprintContainer.diffImKernelPtr = footprintKernelPtr
            diffImFootprintContainer.background = background

            # renormalize the coordinates between -1 and 1?
            # for now, no
            center = footprintBBox.center()
            diffImFootprintContainer.colcNorm = center[0]
            diffImFootprintContainer.rowcNorm = center[1]

            # calculate the residual of the subtracted image here
            convolvedImageStamp = fw.convolve(imageToConvolveStampPtr,
                                              imageToNotConvolveStampPtr,
                                              footprintKernelPtr,
                                              convolveThreshold,
                                              edgeMaskBit)
            differenceImageStamp  = imageToNotConvolveStampPtr - convolvedImageStamp
            differenceImageStamp -= background
            nGoodPixels, meanOfResiduals, varianceOfResiduals = imageproc.calculateMaskedImageResiduals(differenceImageStamp)

            diffImFootprintContainer.footprintResidualMean = meanOfResiduals
            diffImFootprintContainer.footprintResidualVariance = varianceOfResiduals

            if abs(meanOfResiduals) > maximumFootprintResidualMean:
                diffImFootprintContainer.isGood = false
            if varianceOfResiduals > maximumFootprintResidualVariance:
                diffImFootprintContainer.isGood = false

            diffImContainerList.append(diffImFootprintContainer)

            if DEBUG_IO:
                 kImage = footprintKernelPtr->computeNewImage()
                 kImage.writeFits('kFits_%d.fits' % (nFootprint))
                 convolvedImageStamp.writeFits('d1Fits_%d' (nFootprint))

            nFootprint += 1

        #
        # Calculate all individual Kernels
        #
        ###########

        # Test if you do PCA from Policy
        # For now just do it
        kernelOutBasisList = imageproc.computePcaKernelBasis(diffImContainerList)

        # Compute spatial variation of the kernel
        spatialKernelPtr = computeSpatiallyVaryingPsfMatchingKernel(diffImContainerList, kernelOutBasisList, kernelFunctionPtr)

        ###########
        #
        # Compute spatial variation of the background
        #
        backgrounds = []
        variances = []
        position1 = []
        position2 = []

        nGood = 0
        backgroundSum = 0.0
        for iDiffImContainer in diffImContainerList:
            if iDiffImContainer.isGood == False:
                continue

            backgrounds.append(iDiffImContainer.background)
            variances.append(iDiffImContainer.footprintResidualVariance) # This is not entirely correct
            position1.append(iDiffImContainer.colcNorm)
            position2.append(iDiffImContainer.rowcNorm)

            backgroundSum += iDiffImContainer.background
            nGood         += 1

        if nGood == 0:
            # Throw execption; oops!
            return

        # Set up object to fit
        nSigmaSquared = 1.0
        backgroundFunction = fw.MinimizerFunctionBase2(backgrounds, variances, position1, position2, nSigmaSquared, backgroundFunctionPtr)

        # Set up initial guesses; 0th order term is average background, spatial terms zero
        nParameters = backgroundFunctionPtr.getNParameters()
        parameters = num.zeros((nParameters,))
        parameters[0] = backgroundSum / nGood

        # Set up error matrix
        errors = num.zeros((nParameters, 2))
        
        # Minimize!
        outParameters = backgroundFunction.minimize(inParameters, errors)

        # Set actual function parameters
        backgroundFunctionPtr.setParameters(outParameters)

        # Debugging info
        for i in range(len(outParameters)):
            utils.Trace(4, 'Fit Background Parameter %d : %f (%f,%f)' % (i, parameters[i], errors[i][0], errors[i][1]))

        #
        # Compute spatial variation of the background
        #
        ###########

        ###########
        #
        # Create final difference image
        #
        differenceImage  = inputImage - fw.convolve(templateImage, spatialKernelPtr)
        differenceImage -= backgroundFunctionPtr

        # NOTE - Do I do object detection here or not?  Suspect not

        ###########
        #
        # Post results to clipboard
        #
        self.activeClipboard['DifferenceImage'] = differenceImage
        self.activeClipboard['OutputKernel'] = spatialKernelPtr
        self.activeClipboard['BackgroundModel'] = backgroundFunctionPtr
        
    #------------------------------------------------------------------------
    def postprocess(self):
        print 'Python ImageSubtractStage postprocess : stageId %d' % self.stageId
        self.outputQueue.addDataset(self.activeClipboard)

        # potentially delete 'BasisSet'
