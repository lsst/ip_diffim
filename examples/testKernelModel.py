import sys, os, optparse
import numpy
import ip_diffim_plot
import eups

# python
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.daf.base as dafBase
import lsst.ip.diffim as ipDiffim

from lsst.pex.logging import Log
from lsst.pex.logging import Trace
from lsst.pex.policy import Policy

DOKRIGING = False

# Temporary functions until we formalize this in the build system somewhere

### NOTE THESE ARE ALSO IN lsst.afw.image.testUtils
### Which to use?
def vectorToImage(inputVector, nCols, nRows):
    assert len(inputVector) == nCols * nRows
    outputImage = afwImage.ImageF(nCols, nRows)
    nVec = 0
    for nCol in range(nCols):
        for nRow in range(nRows):
            outputImage.setVal(nCol, nRow, inputVector[nVec])
            nVec += 1
    return outputImage

def imageToVector(inputImage):
    nCols = inputImage.getCols()
    nRows = inputImage.getRows()
    outputVector = numpy.zeros(nCols * nRows)
    nVec = 0    
    for nCol in range(nCols):
        for nRow in range(nRows):
            outputVector[nVec] = inputImage.getVal(nCol, nRow)
            nVec += 1
    return outputVector

def imageToMatrix(inputImage):
    nCols = inputImage.getCols()
    nRows = inputImage.getRows()
    outputMatrix = numpy.zeros((nCols,nRows))
    for nCol in range(nCols):
        for nRow in range(nRows):
            outputMatrix[nCol,nRow] = inputImage.getVal(nCol, nRow)
    return outputMatrix

def rejectKernelOutliers(difiList, policy):
    # Compare kernel sums; things that are bad (e.g. variable stars, moving objects, cosmic rays)
    # will have a kernel sum far from the mean; reject these.
    
    maxOutlierIterations = policy.get('lsst.ip.diffim.rejectKernelOutliers').get('maxOutlierIterations')
    maxOutlierSigma      = policy.get('lsst.ip.diffim.rejectKernelOutliers').get('maxOutlierSigma')

    # So are we using pexLog.Trace or pexLog.Log
    Trace('lsst.ip.diffim.rejectKernelOutliers', 3,
          'Rejecting kernels with deviant kSums');
    
    for nIter in xrange(maxOutlierIterations):
        goodDifiList = ipDiffim.getGoodFootprints(difiList)
        nFootprint   = len(goodDifiList)

        if nFootprint == 0:
            raise pex_ex.LsstOutOfRange('No good kernels found')

        ksumVector    = numpy.zeros(nFootprint)
        for i in range(nFootprint):
            ksumVector[i] = goodDifiList[i].getSingleKernelPtr().computeNewImage(False)[1]

        ksumMean = ksumVector.mean()
        ksumStd  = ksumVector.std()

        # reject kernels with aberrent statistics
        numRejected = 0
        for i in range(nFootprint):
            if numpy.fabs( (ksumVector[i]-ksumMean)/ksumStd ) > maxOutlierSigma:
                goodDifiList[i].setIsGood(False)
                numRejected += 1

                Trace('lsst.ip.diffim.rejectKernelOutliers', 5,
                      '# Kernel %d (kSum=%.3f) REJECTED due to bad kernel sum (mean=%.3f, std=%.3f)' %
                      (goodDifiList[i].getID(), ksumVector[i], ksumMean, ksumStd)
                      )
                
        Trace('lsst.ip.diffim.rejectKernelOutliers', 3,
              'Kernel Sum Iteration %d, rejected %d kernels : Kernel Sum = %0.3f (%0.3f)' %
              (nIter, numRejected, ksumMean, ksumStd)
              )

        if numRejected == 0:
            break

    if nIter == (maxOutlierIterations-1) and numRejected != 0:
        Trace('lsst.ip.diffim.rejectKernelOutliers', 1,
              'Detection of kernels with deviant kSums reached its limit of %d iterations'
              (maxOutlierIterations)
              )


    # final loop to report values
    goodDifiList = ipDiffim.getGoodFootprints(difiList)
    nFootprint   = len(goodDifiList)
    ksumVector   = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        ksumVector[i] = goodDifiList[i].getSingleKernelPtr().computeNewImage(False)[1]
    ksumMean = ksumVector.mean()
    ksumStd  = ksumVector.std()
    Trace('lsst.ip.diffim.rejectKernelOutliers', 5,
          'Kernel Sum : %0.3f +/- %0.3f, %d kernels' % (ksumMean, ksumStd, nFootprint));

def fitSpatialFunction(spatialFunction, values, variances, col, row, policy):
    nSigmaSq = policy.get('lsst.afw.math.minimize.nSigmaSq')
    stepsize = policy.get('lsst.afw.math.minimize.stepsize')

    # initialize fit parameters
    nParameters   = spatialFunction.getNParameters()
    parameters    = numpy.zeros(nParameters)           # start with no spatial variation
    parameters[0] = 1.0                                # except for the constant term
    stepsize      = stepsize * numpy.ones(nParameters)
    
    spatialFit    = afwMath.minimize(spatialFunction,
                                     parameters,
                                     stepsize,
                                     values,
                                     variances,
                                     col,
                                     row,
                                     nSigmaSq)
    if not spatialFit.isValid:
        # throw exception
        pass

    return spatialFit
    

def fitKriging(values, variances, col, row, policy):
    return


def fitPerPixel(differenceImageFootprintInformationList, policy):
    kernelSpatialOrder        = policy.get('kernelSpatialOrder')
    kernelSpatialFunction     = afwMath.PolynomialFunction2D(kernelSpatialOrder)
    backgroundSpatialOrder    = policy.get('backgroundSpatialOrder')
    backgroundSpatialFunction = afwMath.PolynomialFunction2D(backgroundSpatialOrder)
    kernelCols                = policy.get('kernelCols')
    kernelRows                = policy.get('kernelRows')

    # how many good footprints are we dealing with here
    goodDifiList = ipDiffim.getGoodFootprints(differenceImageFootprintInformationList)
    nFootprint   = len(goodDifiList)

    # common to all spatial fits
    footprintExposureCol  = numpy.zeros(nFootprint)
    footprintExposureRow  = numpy.zeros(nFootprint)
    footprintVariances    = numpy.zeros(nFootprint)
    
    for i in range(nFootprint):
        footprintExposureCol[i] = goodDifiList[i].getColcNorm()
        footprintExposureRow[i] = goodDifiList[i].getRowcNorm()
        footprintVariances[i]   = goodDifiList[i].getSingleStats().getResidualVariance()


    # fit the background
    backgroundValues    = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        backgroundValues[i] = goodDifiList[i].getSingleBackground()

    backgroundFunctionFit = fitSpatialFunction(backgroundSpatialFunction,
                                               backgroundValues,
                                               footprintVariances,
                                               footprintExposureCol,
                                               footprintExposureRow,
                                               policy)
    backgroundSpatialFunction.setParameters(backgroundFunctionFit.parameterList)

    # plot background!
    if True:
        print footprintExposureCol
        print footprintExposureRow
        print backgroundValues
        print footprintVariances
        
        ip_diffim_plot.plotBackground(backgroundSpatialFunction, footprintExposureCol, footprintExposureRow)
        
    
    if DOKRIGING:
        backgroundKrigingFit  = fitKriging(backgroundValues,
                                           footprintVariances,
                                           footprintExposureCol,
                                           footprintExposureRow,
                                           policy)


    
    # fit each pixel
    functionFitList = []
    if DOKRIGING:
        krigingFitList  = []
    
    # pre-calculate each kernel's image
    # A VECTOR OF IMAGES IS NOT DEFINED IN IMAGELIB.i
    # for now calculate inside the loop below
    
    for kCol in range(kernelCols):
        for kRow in range(kernelRows):

            #######
            # initialize vectors, one per good kernel
            
            kernelValues = numpy.zeros(nFootprint)
            for i in range(nFootprint):
                singleKernelPtr    = goodDifiList[i].getSingleKernelPtr()
                kernelValues[i]    = singleKernelPtr.computeNewImage(False)[0].getVal(kCol, kRow)
                    
            # initialize vectors, one per good kernel
            #######
            # do the various fitting techniques here

            functionFit = fitSpatialFunction(kernelSpatialFunction,
                                             kernelValues,
                                             footprintVariances,
                                             footprintExposureCol,
                                             footprintExposureRow,
                                             policy)
            kernelSpatialFunction.setParameters(functionFit.parameterList)
            functionFitList.append(kernelSpatialFunction)

            if DOKRIGING:
                krigingFit  = fitKriging(kernelValues,
                                         footprintVariances,
                                         footprintExposureCol,
                                         footprintExposureRow,
                                         policy)
                krigingFitList.append(krigingFit)


            # anything else?
            
            # do the various fitting techniques here
            #######

    # evaluate all the fits at the positions of the objects, create a
    # new kernel, then difference image, the calculate difference
    # image stats
    for i in range(nFootprint):
        kFunctionImage  = afwImage.ImageF(kernelCols, kernelRows)
        functionBackgroundValue = backgroundSpatialFunction(footprintExposureCol[i], footprintExposureRow[i])
        
        if DOKRIGING:
            kKrigingImage   = afwImage.ImageF(kernelCols, kernelRows)
            krigingBackgroundValue  = backgroundKrigingFit.eval(footprintExposureCol[i], footprintExposureRow[i])

        # the *very* slow way to do this is to create a
        # LinearCombinationKernel that is of size kernelCols x
        # kernelRows and have a Function for each one of these.  Lets
        # *not* do that for now...
        for kCol in range(kernelCols):
            for kRow in range(kernelRows):
                functionValue   = functionFitList[i](footprintExposureCol[i], footprintExposureRow[i])
                kFunctionImage.set(kCol, kRow, functionValue)
                
                if DOKRIGING:
                    krigingValue    = krigingList[i].eval(footprintExposureCol[i],  footprintExposureRow[i])
                    kKrigingImage.set(kCol, kRow, krigingValue)
                
    
        functionKernelPtr = afwMath.KernelPtr( afwMath.FixedKernel(kFunctionImage) )
        functionStatistics = goodDifiList[i].computeImageStatistics(functionKernelPtr, functionBackgroundValue)
        
        if DOKRIGING:
            krigingKernelPtr  = afwMath.KernelPtr( afwMath.FixedKernel(kKrigingImage) )
            krigingStatistics  = goodDifiList[i].computeImageStatistics(krigingKernelPtr,  krigingBackgroundValue)





def fitPca(differenceImageFootprintInformationList, policy):
    kernelCols = policy.get('kernelCols')
    kernelRows = policy.get('kernelRows')
    
    # how many good footprints are we dealing with here
    goodDifiList = ipDiffim.getGoodFootprints(differenceImageFootprintInformationList)
    nFootprint   = len(goodDifiList)

    # matrix to invert
    M = numpy.zeros(kernelCols*kernelRows, nFootprint)
    for i in range(nFootprint):
        singleKernelImage  = goodDifiList[i].getSingleKernelPtr().computeNewImage(False)[0]
        singleKernelVector = imageToVector(singleKernelImage)
        M[:,i]             = singleKernelVector
        
    # do the PCA
    meanM = numpy.zeros(kernelCols*kernelRows)
    eVal  = numpy.zeros(nFootprint)
    eVec  = numpy.zeros(kernelCols*kernelRows, nFootprint)
    ipDiffim.computePca(meanM, eVal, eVec, M, True)

    # the values of the PCA-approximated kernels
    approxVectorList = VectorListF(nFootprint)
    for i in range(nFootprint):
        # start with the mean kernel
        approxVectorList[i] = meanM.copy()
        
        # calculate approximate statistics
        approxImage         = vectorToImage(approxImageList[i],
                                            kernelCols, kernelRows)
        approxKernelPtr     = afwMath.KernelPtr( afwMath.FixedKernel(approxImage) )
        approxStatistics    = goodDifiList[i].computeImageStatistics(approxKernelPtr,
                                                                     goodDifiList[i].getSingleBackground())
        

    # the first index runs over the number of eigencoefficients
    # the second index runs over the footprints
    nCoefficients = nFootprint
    eCoefficients = numpy.zeros(nCoefficients, nFootprint)

    # now iterate over all footprints and increment the approximate
    # kernel with eigenKernels
    for i in range(nFootprint):
        singleKernelImage   = goodDifiList[i].getSingleKernelPtr().computeNewImage(False)[0]
        singleKernelVector  = imageToVector(singleKernelImage)

        # subtract off mean for dot products        
        singleKernelVector -= meanM  
            
        # approximate the Kernel basis by basis
        for nCoeff in range(nCoefficients):
            eCoefficient              = numpy.dot(singleKernelVector, eVec[nCoeff,:])
            eCoefficients[nCoeff, i]  = eCoefficient
            
            eContribution             = eCoefficient * eVec[nCoeff,:]
            approxVectorList[i]      += eContribution
            
            # calculate approximate statistics
            approxImage          = vectorToImage(approxImageList[i],
                                                 kernelCols, kernelRows)
            approxKernel         = afwMath.KernelPtr( afwMath.FixedKernel(approxImage) )
            approxStatistics     = goodDifiList[i].computeImageStatistics(approxKernelPtr,
                                                                          goodDifiList[i].getSingleBackground())


    return eVec, eCoefficients



def fitSpatialPca(differenceImageFootprintInformationList, eVec, eCoefficients, policy):
    kernelSpatialOrder     = policy.get('kernelSpatialOrder')
    kernelSpatialFunction  = afwMath.PolynomialFunction2D(kernelSpatialOrder)
    backgroundSpatialOrder = policy.get('backgroundSpatialOrder')
    backgroundFunction     = afwMath.PolynomialFunction2D(backgroundSpatialOrder)
    kernelCols = policy.get('kernelCols')
    kernelRows = policy.get('kernelRows')

    # how many good footprints are we dealing with here
    goodDifiList  = ipDiffim.getGoodFootprints(differenceImageFootprintInformationList)
    nFootprint    = len(goodDifiList)
    nCoefficients = nFootprint

    # common to all spatial fits
    footprintExposureCol  = numpy.zeros(nFootprint)
    footprintExposureRow  = numpy.zeros(nFootprint)
    footprintVariances    = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        footprintExposureCol[i] = goodDifiList[i].getColcNorm()
        footprintExposureRow[i] = goodDifiList[i].getRowcNorm()
        footprintVariances[i]   = goodDifiList[i].getSingleStats().getResidualVariance()

    # we need to fit for the background first
    backgroundValues = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        backgroundValues[i] = goodDifiList[i].getSingleBackground()

    backgroundFunctionFit = fitSpatialFunction(backgroundSpatialFunction,
                                               backgroundValues,
                                               footprintVariances,
                                               footprintExposureCol,
                                               footprintExposureRow)
    backgroundSpatialFunction.setParameters(backgroundFunctionFit.parameterList)
    
    if DOKRIGING:
        backgroundKrigingFit  = fitKriging(backgroundValues,
                                           footprintVariances,
                                           footprintExposureCol,
                                           footprintExposureRow)

    # the values of the spatially-approximated kernels
    approxFunctionVectorList = VectorListF(nFootprint)
    if DOKRIGING:
        approxKrigingVectorList  = VectorListF(nFootprint)
    for i in range(nFootprint):
        # start with the mean kernel
        approxFunctionVectorList[i] = meanM.copy()
        if DOKRIGING:
            approxKrigingVectorList[i]  = meanM.copy()

    
    # lets first fit for spatial variation in each coefficient/kernel
    for nCoeff in range(nCoefficients):
        coefficients = eCoefficients[nCoeff,:]

        coeffFunctionFit = fitSpatialFunction(kernelSpatialFunction,
                                              coefficients,
                                              footprintVariances,
                                              footprintExposureCol,
                                              footprintExposureRow)
        kernelSpatialFunction.setParameters(coeffFunctionFit.parameterList)
        
        if DOKRIGING:
            coeffKrigingFit  = fitKriging(coefficients,
                                          footprintVariances,
                                          footprintExposureCol,
                                          footprintExposureRow)
        
    
        for i in range(nFootprints):
            backgroundFunctionValue      = backgroundSpatialFunction(footprintExposureCol[i], footprintExposureRow[i])
            coeffFunctionValue           = kernelSpatialFunction(footprintExposureCol[i], footprintExposureRow[i])
            approxFunctionVectorList[i] += coeffFunctionValue * eVec[nCoeff,:]
            
            # calculate approximate statistics
            approxFunctionImage          = vectorToImage(approxFunctionVectorList[i],
                                                         kernelCols, kernelRows)
            approxFunctionKernelPtr      = afwMath.KernelPtr( afwMath.FixedKernel(approxFunctionImage) )
            approxFunctionStatistics     = goodDifiList[i].computeImageStatistics(approxFunctionKernelPtr,
                                                                                  backgroundFunctionValue)
            if DOKRIGING:
                backgroundKrigingValue       = backgroundKrigingFit.eval(footprintExposureCol[i], footprintExposureRow[i])
                coeffKrigingValue            = coeffKrigingFit.eval(footprintExposureCol[i], footprintExposureRow[i])
                approxKrigingVectorList[i]  += coeffKrigingValue  * eVec[nCoeff,:]

                # calculate approximate statistics
                approxKrigingImage           = vectorToImage(approxKrigingVectorList[i],
                                                             kernelCols, kernelRows)
                approxKrigingKernelPtr       = afwMath.KernelPtr( afwMath.FixedKernel(approxKrigingImage) )
                approxKrigingStatistics      = goodDifiList[i].computeImageStatistics(approxKrigingKernel,
                                                                                      backgroundKrigingValue)
            

##################
##################
##################

def main():
    defDataDir = eups.productDir('afwdata') or ''
    imageProcDir = eups.productDir('ip_diffim')
    if imageProcDir == None:
        print 'Error: could not set up ip_diffim'
        sys.exit(1)

    defSciencePath = os.path.join(defDataDir, 'CFHT', 'D4', 'cal-53535-i-797722_1')
    defTemplatePath = os.path.join(defDataDir, 'CFHT', 'D4', 'cal-53535-i-797722_1_tmpl')
    defPolicyPath = os.path.join(imageProcDir, 'pipeline', 'ImageSubtractStageDictionary.paf')
    defOutputPath = 'diffImage'
    defVerbosity = 0
    
    usage = """usage: %%prog [options] [scienceImage [templateImage [outputImage]]]]

Notes:
- image arguments are paths to MaskedImage fits files
- image arguments must NOT include the final _img.fits
- the result is science image - template image
- the template image is convolved, the science image is not
- default scienceMaskedImage=%s
- default templateMaskedImage=%s
- default outputImage=%s 
- default --policy=%s
""" % (defSciencePath, defTemplatePath, defOutputPath, defPolicyPath)
    
    parser = optparse.OptionParser(usage)
    parser.add_option('-p', '--policy', default=defPolicyPath, help='policy file')
    parser.add_option('-d', '--debugIO', action='store_true', default=False,
        help='write diagnostic intermediate files')
    parser.add_option('-v', '--verbosity', type=int, default=defVerbosity,
        help='verbosity of diagnostic trace messages; 1 for just warnings, more for more information')
    (options, args) = parser.parse_args()
    
    def getArg(ind, defValue):
        if ind < len(args):
            return args[ind]
        return defValue
    
    sciencePath = getArg(0, defSciencePath)
    templatePath = getArg(1, defTemplatePath)
    outputPath = getArg(2, defOutputPath)
    policyPath = options.policy
    
    print 'Science image: ', sciencePath
    print 'Template image:', templatePath
    print 'Output image:  ', outputPath
    print 'Policy file:   ', policyPath
    
    templateMaskedImage = afwImage.MaskedImageF()
    templateMaskedImage.readFits(templatePath)
    
    scienceMaskedImage  = afwImage.MaskedImageF()
    scienceMaskedImage.readFits(sciencePath)

    policy = Policy.createPolicy(policyPath)
    if options.debugIO:
        policy.set('debugIO', True)

    kernelCols = policy.get('kernelCols')
    kernelRows = policy.get('kernelRows')
    
    if options.verbosity > 0:
        print 'Verbosity =', options.verbosity
        Trace.setVerbosity('lsst.ip.diffim', options.verbosity)

    kernelBasisList = ipDiffim.generateDeltaFunctionKernelSet(kernelCols,
                                                              kernelRows)
    # lets just get a couple for debugging and speed
    policy.set('getCollectionOfFootprintsForPsfMatching.minimumCleanFootprints', 5)
    policy.set('getCollectionOfFootprintsForPsfMatching.footprintDetectionThreshold', 6350.)

    footprintList = ipDiffim.getCollectionOfFootprintsForPsfMatching(templateMaskedImage,
                                                                     scienceMaskedImage,
                                                                     policy)
    kImage = afwImage.ImageD(kernelCols, kernelRows)
    
    convolveDifiPtrList   = ipDiffim.DifiPtrListF()
    deconvolveDifiPtrList = ipDiffim.DifiPtrListF()
    for footprintID, iFootprintPtr in enumerate(footprintList):
        footprintBBox = iFootprintPtr.getBBox()
        fpMin = footprintBBox.min()
        fpMax = footprintBBox.max()
        
        Trace('lsst.ip.diffim.computePsfMatchingKernelForFootprint', 5,
              'Footprint %d = %d,%d -> %d,%d' % (footprintID,
                                                 footprintBBox.min().x(), footprintBBox.min().y(),
                                                 footprintBBox.max().x(), footprintBBox.max().y()))
        templateStampPtr = templateMaskedImage.getSubImage(footprintBBox)
        imageStampPtr = scienceMaskedImage.getSubImage(footprintBBox)

        # Assuming the template is better seeing, find convolution kernel
        convKernelCoeffList, convBackground = ipDiffim.computePsfMatchingKernelForFootprint(
            templateStampPtr.get(),
            imageStampPtr.get(),
            kernelBasisList,
            policy,
        )
        convKernelPtr = afwMath.LinearCombinationKernelPtr(
            afwMath.LinearCombinationKernel(kernelBasisList, convKernelCoeffList))

        # Assuming the template is better seeing, find deconvolution kernel
        deconvKernelCoeffList, deconvBackground = ipDiffim.computePsfMatchingKernelForFootprint(
            imageStampPtr.get(),
            templateStampPtr.get(),
            kernelBasisList,
            policy,
        )
        deconvKernelPtr = afwMath.LinearCombinationKernelPtr(
            afwMath.LinearCombinationKernel(kernelBasisList, deconvKernelCoeffList))

        convDiffIm   = ipDiffim.convolveAndSubtract(templateStampPtr.get(), imageStampPtr.get(), convKernelPtr, convBackground)
        
        convData     = imageToVector(convDiffIm.getImage())
        convVariance = imageToVector(convDiffIm.getVariance())
        convMask     = imageToVector(convDiffIm.getMask())
        convIdx      = numpy.where(convMask == 0)
        convSigma    = convData[convIdx] / numpy.sqrt(convVariance[convIdx])
        
        deconvDiffIm   = ipDiffim.convolveAndSubtract(imageStampPtr.get(), templateStampPtr.get(), deconvKernelPtr, deconvBackground)
        deconvData     = imageToVector(deconvDiffIm.getImage())
        deconvVariance = imageToVector(deconvDiffIm.getVariance())
        deconvMask     = imageToVector(deconvDiffIm.getMask())
        deconvIdx      = numpy.where(deconvMask == 0)
        deconvSigma    = deconvData[deconvIdx] / numpy.sqrt(deconvVariance[deconvIdx])

        convInfo     = (imageToMatrix(templateStampPtr.getImage()),
                        imageToMatrix(imageStampPtr.getImage()),
                        imageToMatrix(convDiffIm.getImage()),
                        imageToMatrix(convKernelPtr.computeNewImage(False)[0]),
                        convSigma)
        deconvInfo   = (imageToMatrix(imageStampPtr.getImage()),
                        imageToMatrix(templateStampPtr.getImage()),
                        imageToMatrix(deconvDiffIm.getImage()),
                        imageToMatrix(deconvKernelPtr.computeNewImage(False)[0]),
                        deconvSigma)
        if False:
            ip_diffim_plot.sigmaHistograms(convInfo, deconvInfo)

        convDifiPtr   = ipDiffim.DifiPtrF(
            ipDiffim.DifferenceImageFootprintInformationF(iFootprintPtr, templateStampPtr, imageStampPtr)
            ) 
        convDifiPtr.setID(footprintID)
        convDifiPtr.setColcNorm( float(fpMin.x() + fpMax.x()) / templateMaskedImage.getCols() - 1.0 )
        convDifiPtr.setRowcNorm( float(fpMin.y() + fpMax.y()) / templateMaskedImage.getRows() - 1.0 ) 
        convDifiPtr.setSingleKernelPtr( convKernelPtr )
        convDifiPtr.setSingleBackground( convBackground )
        convDifiPtr.setSingleStats(
            convDifiPtr.computeImageStatistics(convKernelPtr, convBackground)
            )
        convDifiPtr.setStatus( True )
        convolveDifiPtrList.append(convDifiPtr)

        deconvDifiPtr = ipDiffim.DifiPtrF(
            ipDiffim.DifferenceImageFootprintInformationF(iFootprintPtr, imageStampPtr, templateStampPtr)
            )
        deconvDifiPtr.setID(footprintID)
        deconvDifiPtr.setColcNorm( float(fpMin.x() + fpMax.x()) / templateMaskedImage.getCols() - 1.0 )
        deconvDifiPtr.setRowcNorm( float(fpMin.y() + fpMax.y()) / templateMaskedImage.getRows() - 1.0 ) 
        deconvDifiPtr.setSingleKernelPtr( deconvKernelPtr )
        deconvDifiPtr.setSingleBackground( deconvBackground )
        deconvDifiPtr.setSingleStats(
            deconvDifiPtr.computeImageStatistics(deconvKernelPtr, deconvBackground)
            )
        deconvDifiPtr.setStatus( True )
        deconvolveDifiPtrList.append(deconvDifiPtr)

    rejectKernelOutliers(convolveDifiPtrList, policy)
    rejectKernelOutliers(deconvolveDifiPtrList, policy)

    # now we get to the good stuff!
    convResiduals1   = fitPerPixel(convolveDifiPtrList, policy)
    deconvResiduals1 = fitPerPixel(deconvolveDifiPtrList, policy)

    convEVec, convECoeffs     = fitPca(convolveDifiPtrList, policy)
    deconvEVec, deconvECoeffs = fitPca(deconvolveDifiPtrList, policy)

    convResiduals2   = fitSpatialPca(convolveDifiPtrList, convEVec, convECoeffs, policy)
    deconvResiduals2 = fitSpatialPca(deconvolveDifiPtrList, deconvEVec, deconvECoeffs, policy)

def run():
    Log.getDefaultLog()                 # leaks a DataProperty
    
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), 'Objects leaked:'
        print dafBase.Citizen_census(dafBase.cout, memId0)

if __name__ == '__main__':
    run()
    
    
