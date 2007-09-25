import lsst.dps.Stage 
import lsst.fw.Core.fwLib as fw
import lsst.imageproc.imageprocLib as imageproc

class ImageSubtractStage(lsst.dps.Stage):
    #------------------------------------------------------------------------
    def preprocess(self):
        self.activeClipboard = self.inputQueue.getNextDataset()
        keys = self.activeClipboard.getKeys()

        # get policy
        for key in keys:
            value = self.activeClipboard.get(key)
            print 'Python ImageSubtractStage preprocess(): stageId %i key %s value %s' % (self.stageId, key, value)

        # Possible courses of action:

        # Get policy
        # policy = self.activeClipboard.get('Policy')

        # Use policy to create a set of basis functions to be used in all Slices
        # kernelBasisVec = fw.vectorKernelPtrD()
        # imageproc.generateDeltaFunctionKernelSet_D(kernelRows, kernelCols, kernelBasisVec)
        # self.activeClipboard.put('BasisSet', kernelBasisVec)

    #------------------------------------------------------------------------
    def process(self):
        print 'Python ImageSubtractStage process : _rank %i stageId %d' % (self._rank, self.stageId)
        self.activeClipboard = self.inputQueue.getNextDataset()

        ###########
        #
        # Get objects from clipboard
        #
        templateMaskedImage = self.activeClipboard.get('TemplateImage')
        scienceMaskedImage  = self.activeClipboard.get('InputImage')
        basisSet            = self.activeClipboard.get('BasisSet')
        policy              = self.activeClipboard.get('Policy')

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
        computeInC = policy.get('computeInC')

        ###########
        #
        # Generate objects from policy directives
        #

        # create basis vectors
        kernelBasisVec = fw.vectorKernelPtrD()
        imageproc.generateDeltaFunctionKernelSet_D(kernelRows, kernelCols, kernelBasisVec)
        
        # create output kernel pointer
        kernelPtr = imageproc.LinearCombinationKernelPtrTypeD(fw.LinearCombinationKernelD())
        
        # and its function for spatial variation
        kernelFunctionPtr = fw.Function2PtrTypeD(fw.PolynomialFunction2D(kernelSpatialOrder))
        
        # and background function
        backgroundFunctionPtr = fw.Function2PtrTypeD(fw.PolynomialFunction2D(backgroundSpatialOrder))
        
        ###########
        #
        # Get good footprints
        #
        footprintList = detection.FootprintContainerT()
        
        if self.activeClipboard.has_key('FootprintList'):
            footprintList = self.activeClipboard.get('FootprintList')
        else:
            footprintList = imageproc.getCollectionOfFootprintsForPsfMatching_F(templateMaskedImage, scienceMaskedImage, footprintList, policy)

        if computeInC:
            imageproc.computePsfMatchingKernelForMaskedImage_FDD(templateMaskedImage,
                                                                 scienceMaskedImage,
                                                                 kernelBasisVec,
                                                                 footprintList,
                                                                 kernelPtr,
                                                                 kernelFunctionPtr,
                                                                 backgroundFunctionPtr,
                                                                 policy)
        else:
            
            ###########
            #
            # Calculate all individual Kernels
            #
            
            diffImContainerList = imageproc.vectorDiffImContainer_D()
            nFootprint = 0
            for iFootprint in footprintList:
                footprintBBox = iFootprint.getBBox()
        
                templateMaskedImageStampPtr = templateMaskedImage.getSubImage(footprintBBox)
                scienceMaskedImageStampPtr  = scienceMaskedImage.getSubImage(footprintBBox)
        
                if DEBUG_IO:
                    templateMaskedImageStampPtr.writeFits('tFits_%d' % (nFootprint))
                    scienceMaskedImageStampPtr.writeFits('sFits_%d' % (nFootprint))
        
        
                kernelCoeffs = fw.vectorD()
        
                # background is a single number; SWIG returns it here.
                background = imageproc.computePsfMatchingKernelForPostageStamp_FD(templateMaskedImageStampPtr.get(),
                                                                                  scienceMaskedImageStampPtr.get(),
                                                                                  kernelBasisVec,
                                                                                  kernelCoeffs,
                                                                                  policy)
                # Best kernel for this footprint
                footprintKernelPtr = fw.LinearCombinationKernelD(kernelBasisVec, kernelCoeffs)
        
                # Structure holding information about this footprint and its fit to a kernel
                diffImFootprintContainer = imageproc.DiffImContainer_D()
                diffImFootprintContainer.id = nFootprint
                diffImFootprintContainer.isGood = True
                diffImFootprintContainer.diffImFootprintPtr = iFootprint
                diffImFootprintContainer.diffImKernelPtr = footprintKernelPtr
                diffImFootprintContainer.background = background
        
                # renormalize the coordinates between -1 and 1?
                # for now, no
                # NOTE - check and make sure these are the correct coords
                center = footprintBBox.center()
                diffImFootprintContainer.colcNorm = center.y()
                diffImFootprintContainer.rowcNorm = center.x()
        
                # calculate the residual of the subtracted image here
                convolvedImageStamp = fw.convolveD(templateMaskedImageStampPtr,
                                                   footprintKernelPtr,
                                                   convolveThreshold,
                                                   edgeMaskBit)
                differenceImageStamp  = scienceMaskedImageStampPtr - convolvedImageStamp
                differenceImageStamp -= background
                
                nGoodPixels = 0
                meanOfResiduals = 0.0
                varianceOfResiduals = 0.0
                imageproc.calculateMaskedImageResiduals_F(differenceImageStamp, nGoodPixels, meanOfResiduals, varianceOfResiduals)
        
                diffImFootprintContainer.footprintResidualMean = meanOfResiduals
                diffImFootprintContainer.footprintResidualVariance = varianceOfResiduals
                
                if abs(meanOfResiduals) > maximumFootprintResidualMean:
                    diffImFootprintContainer.isGood = False
                if varianceOfResiduals > maximumFootprintResidualVariance:
                    diffImFootprintContainer.isGood = False
                    
                diffImContainerList.append(diffImFootprintContainer)
                
                if DEBUG_IO:
                    kImage = footprintKernelPtr.computeNewImage()
                    kImage.writeFits('kFits_%d.fits' % (nFootprint))
                    convolvedImageStamp.writeFits('dFits_%d' (nFootprint))
        
                nFootprint += 1
        
            #
            # Calculate all individual Kernels
            #
            ###########
            #
            # Calculate the basis kernels and their spatial variation
            #
        
            # Test if you do PCA from Policy
            # For now just do it
            kernelOutBasisList = fw.vectorKernelPtrD()
            imageproc.computePcaKernelBasis_D(diffImContainerList, kernelOutBasisList, policy)
            
            # Compute spatial variation of the kernel
            computeSpatiallyVaryingPsfMatchingKernel_DD(diffImContainerList,
                                                        kernelOutBasisList,
                                                        kernelPtr,
                                                        kernelFunctionPtr,
                                                        policy)
            
            #
            # Calculate the basis kernels and their spatial variation
            #
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
                sys.exit(1)
        
            backgrounds = num.array(backgrounds)
            variances   = num.array(variances)
            position1   = num.array(position1)
            position2   = num.array(position2)
        
            # Set up object to fit
            nSigmaSquared = 1.0
            backgroundFunction = fw.MinimizerFunctionBase2_D(backgrounds, variances, position1, position2, nSigmaSquared, backgroundFunctionPtr)
        
            # Set up initial guesses; 0th order term is average background, spatial terms zero
            nParameters = backgroundFunctionPtr.getNParameters()
            parameters = num.zeros((nParameters,))
            parameters[0] = backgroundSum / nGood
        
            # Step size
            stepsize  = num.ones((nParameters,))
            stepsize *= 0.1
            
            # Set up error matrix
            errors = num.zeros((nParameters, 2))
                
            # Minimize!
            fw.minimize_D(backgroundFunction, parameters, stepsize, errors)
        
            # Set actual function parameters
            backgroundFunctionPtr.setParameters(parameters)
        
            # Debugging info
            for i in range(len(parameters)):
                utils.Trace(4, 'Fit Background Parameter %d : %f (%f,%f)' % (i, parameters[i], errors[i][0], errors[i][1]))
                
            #
            # Compute spatial variation of the background
            #
            ###########
    
        ###########
        #
        # Create final difference image
        #
        if type(spatialKernelPtr) == type(fw.LinearCombinationKernel()):
            differenceImage  = inputImage - fw.convolveLinearD(templateImage, spatialKernelPtr, edgeMaskBit)
        else:
            differenceImage  = inputImage - fw.convolveD(templateImage, spatialKernelPtr, convolveThreshold, edgeMaskBit)
                
        #differenceImage -= backgroundFunctionPtr


        # NOTE - Do I do object detection here or not?  Suspect not

        ###########
        #
        # Post results to clipboard
        #
        self.activeClipboard.put('DifferenceImage', differenceImage)
        self.activeClipboard.put('OutputKernel', spatialKernelPtr)
        self.activeClipboard.put('BackgroundModel', backgroundFunctionPtr)
        
    #------------------------------------------------------------------------
    def postprocess(self):
        print 'Python ImageSubtractStage postprocess : stageId %d' % self.stageId
        self.outputQueue.addDataset(self.activeClipboard)

        # Possible courses of action:

        # Delete anything posted to the Clipboard in preprocess(), e.g. BasisSet
