import sys, os, optparse
import numpy
import ip_diffim_plot
import eups
import pylab

# python
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.daf.base as dafBase
import lsst.ip.diffim as ipDiffim

from lsst.pex.logging import Log
from lsst.pex.logging import Trace
from lsst.pex.policy import Policy
import lsst.pex.exceptions as Exceptions 

DOKRIGING  = False

# Temporary functions until we formalize this in the build system somewhere

### NOTE THESE ARE ALSO IN lsst.afw.image.testUtils
### Which to use?
def imageToMatrix(im, dtype=float):
    arr = numpy.zeros([im.getCols(), im.getRows()], dtype=dtype)
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            arr[col, row] = im.getVal(col, row)
    return arr

def imageToVector(im, dtype=float):
    arr = numpy.zeros([im.getCols()*im.getRows()], dtype=dtype)
    n   = 0
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            arr[n] = im.getVal(col, row)
            n += 1
    return arr

def matrixToImage(arr):
    im = afwImage.ImageF(arr.shape[0], arr.shape[1])
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            im.set(col, row, arr[col, row])
    return im

def vectorToImage(arr, nCols, nRows):
    assert len(arr) == nCols * nRows
    im = afwImage.ImageF(nCols, nRows)
    n  = 0
    for row in range(im.getRows()):
        for col in range(im.getCols()):
            im.set(col, row, arr[n])
            n += 1
    return im

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
            raise Exceptions.LsstOutOfRange('No good kernels found')

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

def fitSpatialFunction(spatialFunction, values, errors, col, row, policy):
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
                                     errors,
                                     col,
                                     row,
                                     nSigmaSq)
    if not spatialFit.isValid:
        # throw exception
        pass

    return spatialFit
    

def fitKriging(values, errors, col, row, policy):
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
    footprintStd          = numpy.zeros(nFootprint)
    
    for i in range(nFootprint):
        footprintExposureCol[i] = goodDifiList[i].getColcNorm()
        footprintExposureRow[i] = goodDifiList[i].getRowcNorm()
        footprintStd[i]         = goodDifiList[i].getSingleStats().getResidualStd()


    # fit the background
    backgroundValues    = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        backgroundValues[i] = goodDifiList[i].getSingleBackground()

    backgroundFunctionFit = fitSpatialFunction(backgroundSpatialFunction,
                                               backgroundValues,
                                               footprintStd,
                                               footprintExposureCol,
                                               footprintExposureRow,
                                               policy)
    backgroundSpatialFunction.setParameters(backgroundFunctionFit.parameterList)

#    # plot background!
#    if policy.get('debugPlot') == True:
#        #print footprintExposureCol
#        #print footprintExposureRow
#        #print backgroundValues
#        #print footprintStd
#        ip_diffim_plot.plotBackground(backgroundSpatialFunction, footprintExposureCol, footprintExposureRow)
        
    return backgroundSpatialFunction, footprintExposureCol, footprintExposureRow

    if DOKRIGING:
        backgroundKrigingFit  = fitKriging(backgroundValues,
                                           footprintStd,
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
                                             footprintStd,
                                             footprintExposureCol,
                                             footprintExposureRow,
                                             policy)
            kernelSpatialFunction.setParameters(functionFit.parameterList)
            functionFitList.append(kernelSpatialFunction)

            if DOKRIGING:
                krigingFit  = fitKriging(kernelValues,
                                         footprintStd,
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
    M = numpy.zeros((kernelCols*kernelRows, nFootprint))
    for i in range(nFootprint):
        singleKernelImage  = goodDifiList[i].getSingleKernelPtr().computeNewImage(False)[0]
        singleKernelVector = imageToVector(singleKernelImage)
        M[:,i]             = singleKernelVector
        
    # do the PCA
    #meanM = numpy.zeros((kernelCols*kernelRows))
    #eVal  = numpy.zeros((nFootprint))
    #eVec  = numpy.zeros((kernelCols*kernelRows, nFootprint))
    # unfortunately, this depends on vw::math.
    # ipDiffim.computePca(meanM, eVal, eVec, M, True)
    # do it in python for now...

    # Structure of numpy arrays :
    #   numpy.zeros( (2,3) )
    #   array([[ 0.,  0.,  0.],
    #          [ 0.,  0.,  0.]])
    # The first dimension is rows, second is columns

    # We are going to put the features down a given row; the instances
    # one in each column.  We need to subtract off the Mean feature -
    # array.mean(0) returns the mean for each column; array.mean(1)
    # returns the mean for each row.  Therefore the mean Kernel is
    # M.mean(1).
    meanM = M.mean(1)

    # Subtract off the mean Kernel from each column.
    M    -= meanM[:,numpy.newaxis]

    # def svd(a, full_matrices=1, compute_uv=1):
    #
    #    """Singular Value Decomposition.
    #
    #    Factorizes the matrix a into two unitary matrices U and Vh and
    #    an 1d-array s of singular values (real, non-negative) such that
    #    a == U S Vh  if S is an suitably shaped matrix of zeros whose
    #    main diagonal is s.
    #
    #    Parameters
    #    ----------
    #    a : array-like, shape (M, N)
    #        Matrix to decompose
    #    full_matrices : boolean
    #        If true,  U, Vh are shaped  (M,M), (N,N)
    #        If false, the shapes are    (M,K), (K,N) where K = min(M,N)
    #    compute_uv : boolean
    #        Whether to compute U and Vh in addition to s
    #
    #    Returns
    #    -------
    #    U:  array, shape (M,M) or (M,K) depending on full_matrices
    #    s:  array, shape (K,)
    #        The singular values, sorted so that s[i] >= s[i+1]
    #        K = min(M, N)
    #    Vh: array, shape (N,N) or (K,N) depending on full_matrices
    #
    # 

    # A given eigenfeature of M will be down a row of U; the different
    # eigenfeatures are in the columns of U.  I.e the primary
    # eigenfeature is U[:,0].
    #
    # The eigenvalues correspond to s**2 and are already sorted
    U,s,Vh = numpy.linalg.svd( M, full_matrices=0 )
    eVal   = s**2

    # Find the contribution of each eigenKernel to each Kernel.
    # Simple dot product, transpose of M dot the eigenCoeff matrix.
    # The contribution of eigenKernel X to Kernel Y is in eCoeff[Y,X].
    #
    # I.e. M[:,X] = numpy.sum(U * eCoeff[X], 1)
    eCoeff = numpy.dot(M.T, U)
    for i in range(nFootprint):
        residual = numpy.sum(U * eCoeff[i], 1) - M[:,i]
        assert(numpy.sum(residual) < 1e-10)
    
    # Now make PCA-approximations of each Kernel, and watch how the
    # diffim residuals decrease as a function of # of eigenKernels.
    for i in range(nFootprint):
        kModel   = meanM.copy()

        # Calculate statistics using approximated Kernel
        approxImage         = vectorToImage(kModel, kernelCols, kernelRows)
        #approxKernelPtr     = afwMath.KernelPtr( afwMath.FixedKernel(approxImage) )
        #approxStatistics    = goodDifiList[i].computeImageStatistics(approxKernelPtr,
        #                                                             goodDifiList[i].getSingleBackground())

        eContrib = U * eCoeff[i]
        for nKernel in range(eContrib.shape[1]):
            kModel += eContrib[:,nKernel]

            # Calculate statistics ...
            approxImage         = vectorToImage(kModel, kernelCols, kernelRows)
            #approxKernelPtr     = afwMath.KernelPtr( afwMath.FixedKernel(approxImage) )
            #approxStatistics    = goodDifiList[i].computeImageStatistics(approxKernelPtr,
            #                                                             goodDifiList[i].getSingleBackground())
            

    # Turn into Kernels
    #eKernelPtrVector = afwMath.VectorKernel()
    #mKernelImage     = vectorToImage( meanM, kernelCols, kernelRows )
    #mKernelPtr       = afwMath.KernelPtr( afwMath.FixedKernel(mKernelImage) )
    #eKernelPtrVector.append(mKernelPtr)
    #for i in range(U.shape[1]):
    #    eKernelImage = vectorToImage( U[:,i], kernelCols, kernelRows )
    #    eKernelPtr   = afwMath.KernelPtr( afwMath.FixedKernel(eKernelImage) )
    #    #eKernelPtrVector.append(eKernelPtr)

    # return eKernelPtrVector, eVal, eCoeff

    # HACK UNTIL I CAN MAKE FIXEDKERNELS
    info = ( imageToMatrix( vectorToImage( meanM, kernelCols, kernelRows ) ),
             imageToMatrix( vectorToImage( U[:,0], kernelCols, kernelRows ) ),
             imageToMatrix( vectorToImage( U[:,1], kernelCols, kernelRows ) ),
             imageToMatrix( vectorToImage( U[:,2], kernelCols, kernelRows ) ),
             imageToMatrix( vectorToImage( U[:,3], kernelCols, kernelRows ) ),
             imageToMatrix( vectorToImage( U[:,4], kernelCols, kernelRows ) ),
             eVal )

    return info
    



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
    footprintStd          = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        footprintExposureCol[i] = goodDifiList[i].getColcNorm()
        footprintExposureRow[i] = goodDifiList[i].getRowcNorm()
        footprintStd[i]         = goodDifiList[i].getSingleStats().getResidualStd()

    # we need to fit for the background first
    backgroundValues = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        backgroundValues[i] = goodDifiList[i].getSingleBackground()

    backgroundFunctionFit = fitSpatialFunction(backgroundSpatialFunction,
                                               backgroundValues,
                                               footprintStd,
                                               footprintExposureCol,
                                               footprintExposureRow)
    backgroundSpatialFunction.setParameters(backgroundFunctionFit.parameterList)
    
    if DOKRIGING:
        backgroundKrigingFit  = fitKriging(backgroundValues,
                                           footprintStd,
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
                                              footprintStd,
                                              footprintExposureCol,
                                              footprintExposureRow)
        kernelSpatialFunction.setParameters(coeffFunctionFit.parameterList)
        
        if DOKRIGING:
            coeffKrigingFit  = fitKriging(coefficients,
                                          footprintStd,
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
    #policy.set('getCollectionOfFootprintsForPsfMatching.footprintDetectionThreshold', 6350.)
    policy.set('getCollectionOfFootprintsForPsfMatching.footprintDetectionThreshold', 1000.)

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
        
        Trace('lsst.ip.diffim.computePsfMatchingKernelForFootprint', 2,
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
            afwMath.LinearCombinationKernel(kernelBasisList, convKernelCoeffList)
            )
        convDiffIm    = ipDiffim.convolveAndSubtract(templateStampPtr.get(), imageStampPtr.get(), convKernelPtr, convBackground)

        # Assuming the template is better seeing, find deconvolution kernel
        deconvKernelCoeffList, deconvBackground = ipDiffim.computePsfMatchingKernelForFootprint(
            imageStampPtr.get(),
            templateStampPtr.get(),
            kernelBasisList,
            policy,
        )
        deconvKernelPtr = afwMath.LinearCombinationKernelPtr(
            afwMath.LinearCombinationKernel(kernelBasisList, deconvKernelCoeffList)
            )
        deconvDiffIm    = ipDiffim.convolveAndSubtract(imageStampPtr.get(), templateStampPtr.get(), deconvKernelPtr, deconvBackground)

        if policy.get('debugIO') == True:
            imageStampPtr.writeFits('iFoot_%d' % (footprintID))
            templateStampPtr.writeFits('tFoot_%d' % (footprintID))

            convKernelPtr.computeNewImage(False)[0].writeFits('cKernel_%d.fits' % (footprintID))
            deconvKernelPtr.computeNewImage(False)[0].writeFits('dcKernel_%d.fits' % (footprintID))

            convDiffIm.writeFits('cDiff_%d' % (footprintID))
            deconvDiffIm.writeFits('dcDiff_%d' % (footprintID))

            cSigma = imageToMatrix(convDiffIm.getImage()) / numpy.sqrt(imageToMatrix(convDiffIm.getVariance()))
            cSigma = matrixToImage(cSigma)
            cSigma.writeFits('cSig_%d.fits' % (footprintID))

            dcSigma = imageToMatrix(deconvDiffIm.getImage()) / numpy.sqrt(imageToMatrix(deconvDiffIm.getVariance()))
            dcSigma = matrixToImage(dcSigma)
            dcSigma.writeFits('dcSig_%d.fits' % (footprintID))
            
        
        if policy.get('debugPlot') == True:
            convData     = imageToVector(convDiffIm.getImage())
            convVariance = imageToVector(convDiffIm.getVariance())
            convMask     = imageToVector(convDiffIm.getMask())
            convIdx      = numpy.where(convMask == 0)
            convSigma    = convData[convIdx] / numpy.sqrt(convVariance[convIdx])
            Trace('lsst.ip.diffim', 7,
                  'Kernel %d : Python diffim residuals = %.2f +/- %.2f sigma' % (footprintID, convSigma.mean(), convSigma.std()))
           
            deconvData     = imageToVector(deconvDiffIm.getImage())
            deconvVariance = imageToVector(deconvDiffIm.getVariance())
            deconvMask     = imageToVector(deconvDiffIm.getMask())
            deconvIdx      = numpy.where(deconvMask == 0)
            deconvSigma    = deconvData[deconvIdx] / numpy.sqrt(deconvVariance[deconvIdx])
            Trace('lsst.ip.diffim', 7,
                  'Kernel %d : Python diffim residuals = %.2f +/- %.2f sigma' % (footprintID, deconvSigma.mean(), deconvSigma.std()))
            
            convInfo     = (imageToMatrix(templateStampPtr.getImage()),
                            imageToMatrix(imageStampPtr.getImage()),
                            imageToMatrix(convDiffIm.getImage()) / numpy.sqrt(imageToMatrix(convDiffIm.getVariance())),
                            imageToMatrix(convKernelPtr.computeNewImage(False)[0]),
                            convSigma)
            deconvInfo   = (imageToMatrix(imageStampPtr.getImage()),
                            imageToMatrix(templateStampPtr.getImage()),
                            imageToMatrix(deconvDiffIm.getImage()) / numpy.sqrt(imageToMatrix(deconvDiffIm.getVariance())),
                            imageToMatrix(deconvKernelPtr.computeNewImage(False)[0]),
                            deconvSigma)
            ip_diffim_plot.sigmaHistograms(convInfo, deconvInfo, title='Kernel %d' % (footprintID), outfile='Kernel_%d.ps' % (footprintID) )
            

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
        convDifiPtr.setStatus( convDifiPtr.getSingleStats().evaluateQuality(policy) )
        convolveDifiPtrList.append(convDifiPtr)
        if convDifiPtr.getStatus() == True:
            prefix = ''
        else:
            prefix = '#'
        Trace('lsst.ip.diffim', 5,
              '%sKernel %d : Kernel Sum = %.2f, Diffim residuals = %.2f +/- %.2f sigma' % (
            prefix,
            convDifiPtr.getID(),
            convDifiPtr.getSingleKernelPtr().computeNewImage(False)[1],
            convDifiPtr.getSingleStats().getResidualMean(),
            convDifiPtr.getSingleStats().getResidualStd()
            ))

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
        deconvDifiPtr.setStatus( deconvDifiPtr.getSingleStats().evaluateQuality(policy) )
        deconvolveDifiPtrList.append(deconvDifiPtr)
        if deconvDifiPtr.getStatus() == True:
            prefix = ''
        else:
            prefix = '#'
        Trace('lsst.ip.diffim', 5,
              '%sKernel %d : Kernel Sum = %.2f, Diffim residuals = %.2f +/- %.2f sigma' % (
            prefix,
            deconvDifiPtr.getID(),
            deconvDifiPtr.getSingleKernelPtr().computeNewImage(False)[1],
            deconvDifiPtr.getSingleStats().getResidualMean(),
            deconvDifiPtr.getSingleStats().getResidualStd()
            ))

    rejectKernelOutliers(convolveDifiPtrList, policy)
    rejectKernelOutliers(deconvolveDifiPtrList, policy)

    # now we get to the good stuff!
    # HACK UNTIL I CAN MAKE FIXEDKERNELS
    #convEKernel,   convEVals,   convECoeffs   = fitPca(convolveDifiPtrList, policy)
    #deconvEKernel, deconvEVals, deconvECoeffs = fitPca(deconvolveDifiPtrList, policy)

    convInfo   = fitPca(convolveDifiPtrList, policy)
    deconvInfo = fitPca(deconvolveDifiPtrList, policy)

    if policy.get('debugPlot') == True:
        ip_diffim_plot.eigenKernelPlot(convInfo, deconvInfo, outfile='EKernel.ps')

    # foo for plotting
    fooconvResiduals1   = fitPerPixel(convolveDifiPtrList, policy)
    foodeconvResiduals1 = fitPerPixel(deconvolveDifiPtrList, policy)
    if policy.get('debugPlot') == True:
        ip_diffim_plot.plotBackground(fooconvResiduals1, outfile='cBackground.ps')
        ip_diffim_plot.plotBackground(foodeconvResiduals1, outfile='dcBackground.ps')

    #convResiduals1   = fitPerPixel(convolveDifiPtrList, policy)
    #deconvResiduals1 = fitPerPixel(deconvolveDifiPtrList, policy)

    #convResiduals2   = fitSpatialPca(convolveDifiPtrList, convEVec, convECoeffs, policy)
    #deconvResiduals2 = fitSpatialPca(deconvolveDifiPtrList, deconvEVec, deconvECoeffs, policy)

    if policy.get('debugPlot') == True:
        pylab.show()


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
    
    
