import sys, os, optparse
import numpy
import eups

# python
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.daf.base as dafBase
import lsst.ip.diffim as ipDiffim

from lsst.pex.logging import Log
from lsst.pex.policy import Policy

def fitSpatialFunction(spatialFunction, values, variances, col, row, policy):
    nSigmaSq = policy.get('afwMath.minimize.nSigmaSq')
    stepsize = policy.get('afwMath.minimize.stepsize')

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
    return spatialFit
    

def fitKriging(values, variances, col, row, policy):
    import sgems
    
    return

def fitPerPixel(differenceImageFootprintInformationList, policy):
    kernelSpatialOrder     = policy.get('kernelSpatialOrder')
    kernelSpatialFunction  = afwMath.PolynomialFunction2D(kernelSpatialOrder)
    backgroundSpatialOrder = policy.get('backgroundSpatialOrder')
    backgroundFunction     = afwMath.PolynomialFunction2D(backgroundSpatialOrder)
    kernelCols             = policy.get('kernelCols')
    kernelRows             = policy.get('kernelRows')

    # how many good footprints are we dealing with here
    goodDifiList = differenceImageFootprintInformationList.getGoodFootprints()
    nFootprint   = len(goodDifiList)

    # common to all spatial fits
    footprintExposureCol  = numpy.zeros(nFootprint)
    footprintExposureRow  = numpy.zeros(nFootprint)
    footprintVariances    = numpy.zeros(nFootprint)
    
    for i in range(nFootprint):
        footprintExposureCol[i] = goodDifiList[i].getColcNorm()
        footprintExposureRow[i] = goodDifiList[i].getRowcNorm()
        footprintVariances[i]   = goodDifiList[i].getSingleStats().getFootprintResidualVariance()
    
    # fit each pixel
    functionFitList = []
    krigingFitList  = []
    
    for kCol in range(kernelCols):
        for kRow in range(kernelRows):

            #######
            # initialize vectors, one per good kernel
            
            kernelValues = numpy.zeros(nFootprint)
            for i in range(nFootprint):
                singleKernel       = goodDifiList[i].getSingleKernel()
                kernelValues[i]    = singleKernel.getImage().getValue(kCol, kRow)
                    
            # initialize vectors, one per good kernel
            #######
            # do the various fitting techniques here

            functionFit = fitSpatialFunction(kernelSpatialFunction,
                                             kernelValues,
                                             footprintVariances,
                                             footprintExposureCol,
                                             footprintExposureRow)
            
            krigingFit  = fitKriging(kernelValues,
                                     footprintVariances,
                                     footprintExposureCol,
                                     footprintExposureRow)

            functionFitList.append(functionFit)
            krigingFitList.append(krigingFit)

            # anything else?
            
            # do the various fitting techniques here
            #######

    # fit the background
    backgroundValues    = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        backgroundValues[i] = goodDifiList[i].getSingleBackground()

    backgroundFunctionFit = fitSpatialFunction(backgroundSpatialFunction,
                                               backgroundValues,
                                               footprintVariances,
                                               footprintExposureCol,
                                               footprintExposureRow)
    backgroundKrigingFit  = fitKriging(backgroundValues,
                                       footprintVariances,
                                       footprintExposureCol,
                                       footprintExposureRow)


    # evaluate all the fits at the positions of the objects, create a
    # new kernel, then difference image, the calculate difference
    # image stats
    for i in range(nFootprint):
        kFunctionImage  = afwImage.ImageD(kernelCols, kernelRows)
        kKrigingImage   = afwImage.ImageD(kernelCols, kernelRows)

        functionBackgroundValue = backgroundFunctionFit.eval(footprintExposureCol[i], footprintExposureRow[i])
        krigingBackgroundValue  = backgroundKrigingFit.eval(footprintExposureCol[i], footprintExposureRow[i])

        # the *very* slow way to do this is to create a
        # LinearCombinationKernel that is of size kernelCols x
        # kernelRows and have a Function for each one of these.  Lets
        # *not* do that for now...
        for kCol in range(kernelCols):
            for kRow in range(kernelRows):
                functionValue   = functionList[i].eval(footprintExposureCol[i], footprintExposureRow[i])
                krigingValue    = krigingList[i].eval(footprintExposureCol[i],  footprintExposureRow[i])
                
                kFunctionImage.setValue(kCol, kRow, functionValue)
                kKrigingImage.setValue(kCol, kRow, krigingValue)
    
        functionKernelPtr = afwMath.KernelPtr( afwMath.Kernel(kFunctionImage) )
        krigingKernelPtr  = afwMath.KernelPtr( afwMath.Kernel(kKrigingImage) )

        functionStatistics = goodDifiList[i].computeImageStatistics(functionKernelPtr, functionBackgroundValue)
        krigingStatistics  = goodDifiList[i].computeImageStatistics(krigingKernelPtr,  krigingBackgroundValue)




def fitPca(differenceImageFootprintInformationList, policy):
    kernelCols = policy.get('kernelCols')
    kernelRows = policy.get('kernelRows')
    
    # how many good footprints are we dealing with here
    goodDifiList = differenceImageFootprintInformationList.getGoodFootprints()
    nFootprint   = len(goodDifiList)

    # matrix to invert
    M = numpy.zeros(kernelCols*kernelRows, nFootprint)
    for i in range(nFootprint):
        singleKernelImage  = goodDifiList[i].getSingleKernel().getImage()
        singleKernelVector = afwMath.imageToVector(singleKernelImage, kernelCols, kernelRows)
        M[:,i]             = singleKernelVector
        
    # do the PCA
    meanM = numpy.zeros(kernelCols*kernelRows)
    eVal  = numpy.zeros(nFootprint)
    eVec  = numpy.zeros(kernelCols*kernelRows, nFootprint)
    ipDiffim.computePca(meanM, eVal, eVec, M, True)

    # the values of the PCA-approximated kernels
    approxVectorList = VectorListD(nFootprint)
    for i in range(nFootprint):
        # start with the mean kernel
        approxVectorList[i] = meanM.copy()
        
        # calculate approximate statistics
        approxImage         = afwMath.vectorToImage(approxImageList[i],
                                                    kernelCols, kernelRows)
        approxKernelPtr     = afwMath.KernelPtr( afwMath.Kernel(approxImage) )
        approxStatistics    = goodDifiList[i].computeImageStatistics(approxKernelPtr,
                                                                     goodDifiList[i].getSingleBackground())
        

    # the first index runs over the number of eigencoefficients
    # the second index runs over the footprints
    nCoefficients = nFootprint
    eCoefficients = numpy.zeros(nCoefficients, nFootprint)

    # now iterate over all footprints and increment the approximate
    # kernel with eigenKernels
    for i in range(nFootprint):
        singleKernelImage   = goodDifiList[i].getSingleKernel().getImage()
        singleKernelVector  = afwMath.imageToVector(singleKernelImage, kernelCols, kernelRows)

        # subtract off mean for dot products        
        singleKernelVector -= meanM  
            
        # approximate the Kernel basis by basis
        for nCoeff in range(nCoefficients):
            eCoefficient              = numpy.dot(singleKernelVector, eVec[nCoeff,:])
            eCoefficients[nCoeff, i]  = eCoefficient
            
            eContribution             = eCoefficient * eVec[nCoeff,:]
            approxVectorList[i]      += eContribution
            
            # calculate approximate statistics
            approxImage          = afwMath.vectorToImage(approxImageList[i],
                                                         kernelCols, kernelRows)
            approxKernel         = afwMath.KernelPtr( afwMath.Kernel(approxImage) )
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
    goodDifiList  = differenceImageFootprintInformationList.getGoodFootprints()
    nFootprint    = len(goodDifiList)
    nCoefficients = nFootprint

    # common to all spatial fits
    footprintExposureCol  = numpy.zeros(nFootprint)
    footprintExposureRow  = numpy.zeros(nFootprint)
    footprintVariances    = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        footprintExposureCol[i] = goodDifiList[i].getColcNorm()
        footprintExposureRow[i] = goodDifiList[i].getRowcNorm()
        footprintVariances[i]   = goodDifiList[i].getSingleStats().getFootprintResidualVariance()

    # we need to fit for the background first
    backgroundValues = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        backgroundValues[i] = goodDifiList[i].getSingleBackground()

    backgroundFunctionFit = fitSpatialFunction(backgroundSpatialFunction,
                                               backgroundValues,
                                               footprintVariances,
                                               footprintExposureCol,
                                               footprintExposureRow)
    backgroundKrigingFit  = fitKriging(backgroundValues,
                                       footprintVariances,
                                       footprintExposureCol,
                                       footprintExposureRow)

    # the values of the spatially-approximated kernels
    approxFunctionVectorList = VectorListD(nFootprint)
    approxKrigingVectorList  = VectorListD(nFootprint)
    for i in range(nFootprint):
        # start with the mean kernel
        approxFunctionVectorList[i] = meanM.copy()
        approxKrigingVectorList[i]  = meanM.copy()

    
    # lets first fit for spatial variation in each coefficient/kernel
    for nCoeff in range(nCoefficients):
        coefficients = eCoefficients[nCoeff,:]

        coeffFunctionFit = fitSpatialFunction(kernelSpatialFunction,
                                              coefficients,
                                              footprintVariances,
                                              footprintExposureCol,
                                              footprintExposureRow)
        coeffKrigingFit  = fitKriging(coefficients,
                                      footprintVariances,
                                      footprintExposureCol,
                                      footprintExposureRow)
        
    
        for i in range(nFootprints):
            backgroundFunctionValue = backgroundFunctionFit.eval(footprintExposureCol[i], footprintExposureRow[i])
            backgroundKrigingValue  = backgroundKrigingFit.eval(footprintExposureCol[i], footprintExposureRow[i])
            coeffFunctionValue      = coeffFunctionFit.eval(footprintExposureCol[i], footprintExposureRow[i])
            coeffKrigingValue       = coeffKrigingFit.eval(footprintExposureCol[i], footprintExposureRow[i])

            approxFunctionVectorList[i] += coeffFunctionValue * eVec[nCoeff,:]
            approxKrigingVectorList[i]  += coeffKrigingValue  * eVec[nCoeff,:]

            # calculate approximate statistics
            approxFunctionImage          = afwMath.vectorToImage(approxFunctionVectorList[i],
                                                                 kernelCols, kernelRows)
            approxFunctionKernelPtr      = afwMath.KernelPtr( afwMath.Kernel(approxFunctionImage) )
            approxFunctionStatistics     = goodDifiList[i].computeImageStatistics(approxFunctionKernelPtr,
                                                                                  backgroundFunctionValue)
            
            
            approxKrigingImage           = afwMath.vectorToImage(approxKrigingVectorList[i],
                                                                 kernelCols, kernelRows)
            approxKrigingKernelPtr       = afwMath.KernelPtr( afwMath.Kernel(approxKrigingImage) )
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
    parser.add_option('-s', '--switchConvolve', action='store_true', default=False,
        help='switch which image is convolved; still detect on template')
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
    if options.switchConvolve:
        policy.set('switchConvolve', True)

    kernelCols = policy.get('kernelCols')
    kernelRows = policy.get('kernelRows')
    
    if options.verbosity > 0:
        print 'Verbosity =', options.verbosity
        lsst.pex.logging.Trace_setVerbosity('lsst.ip.diffim', options.verbosity)

    print 'foo'
    kernelBasisList = ipDiffim.generateDeltaFunctionKernelSet(kernelCols,
                                                              kernelRows)
    print 'foo'
    footprintList = ipDiffim.getCollectionOfFootprintsForPsfMatching(templateMaskedImage,
                                                                     scienceMaskedImage,
                                                                     policy)
    print 'foo'

    kImage = afwImage.ImageD(kernelCols, kernelRows)
    
    convolveDifiList   = ipDiffim.vectorDiffImContainerD()
    deconvolveDifiList = ipDiffim.vectorDiffImContainerD()
    for footprintID, iFootprintPtr in enumerate(footprintList):
        footprintBBox = iFootprintPtr.getBBox()
        fpMin = footprintBBox.min()
        fpMax = footprintBBox.max()
        
        logging.Trace('lsst.ip.diffim.computePsfMatchingKernelForMaskedImage', 5,
                      'Footprint %d = %d,%d -> %d,%d' % (footprintID,
                                                         footprintBBox.min().x(), footprintBBox.min().y(),
                                                         footprintBBox.max().x(), footprintBBox.max().y()))
        templateStampPtr = templateMaskedImage.getSubImage(footprintBBox)
        imageStampPtr = scienceMaskedImage.getSubImage(footprintBBox)

        # Assuming the template is better seeing, find convolution kernel
        convKernelCoeffList, convBackground = ipDiffim.computePsfMatchingKernelForPostageStamp(
            templateStampPtr.get(),
            imageStampPtr.get(),
            kernelBasisList,
            policy,
        )
        convKernelPtr = afwMath.LinearCombinationKernelPtr(
            afwMath.LinearCombinationKernel(kernelBasisList, convKernelCoeffList))

        # Assuming the template is better seeing, find deconvolution kernel
        deconvKernelCoeffList, deconvBackground = ipDiffim.computePsfMatchingKernelForPostageStamp(
            imageStampPtr.get(),
            templateStampPtr.get(),
            kernelBasisList,
            policy,
        )
        deconvKernelPtr = afwMath.LinearCombinationKernelPtr(
            afwMath.LinearCombinationKernel(kernelBasisList, deconvKernelCoeffList))
        
        convDifi   = ipDiffim.DiffImContainerD()
        convDifi.setId(footprintID)
        convDifi.setConvolveMI(templateStampPtr.get())
        convDifi.setNotConvolveMI(imageStampPtr.get())
        convDifi.setColcNorm( 0.5 * (fpMin.x() + fpMax.x()) ) # or can i use footprint.center or something?
        convDifi.setRowcNorm( 0.5 * (fpMin.y() + fpMax.y()) ) 
        convDifi.setSingleKernel( convKernelPtr )
        convDifi.setSingleBackground( convBackground )
        convolveDifiList.append(convDifi)

        deconvDifi = ipDiffim.DiffImContainerD()
        deconvDifi.setId(footprintID)
        deconvDifi.setConvolveMI(imageStampPtr.get())
        deconvDifi.setNotConvolveMI(templateStampPtr.get())
        deconvDifi.setColcNorm( 0.5 * (fpMin.x() + fpMax.x()) )
        deconvDifi.setRowcNorm( 0.5 * (fpMin.y() + fpMax.y()) ) 
        deconvDifi.setSingleKernel( deconvKernelPtr )
        deconvDifi.setSingleBackground( deconvBackground )
        deconvolveDifiList.append(deconvDifi)

    ipDiffim.rejectKernelOutliers(convolveDifiList, policy)
    ipDiffim.rejectKernelOutliers(deconvolveDifiList, policy)

    # now we get to the good stuff!
    convResiduals1   = fitPerPixel(convolveDifiList, policy)
    deconvResiduals1 = fitPerPixel(deconvolveDifiList, policy)

    convEVec, convECoeffs     = fitPca(convolveDifiList, policy)
    deconvEVec, deconvECoeffs = fitPca(deconvolveDifiList, policy)

    convResiduals2   = fitSpatialPca(convolveDifiList, convEVec, convECoeffs, policy)
    deconvResiduals2 = fitSpatialPca(deconvolveDifiList, deconvEVec, deconvECoeffs, policy)

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
    
    
