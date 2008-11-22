import sys, os, optparse
import numpy
import eups
import pylab

import pdb

# python
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.daf.base as dafBase
import lsst.ip.diffim as ipDiffim

import lsst.ip.diffim.diffimPlot as ipDiffimPlot
import lsst.ip.diffim.diffimTools as ipDiffimTools
import lsst.ip.diffim.diffimDebug as ipDiffimDebug

from lsst.pex.logging import Log
from lsst.pex.logging import Trace
from lsst.pex.policy import Policy
import lsst.pex.exceptions as Exceptions 

DOKRIGING  = False

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
                goodDifiList[i].setStatus(False)
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


def runPca(differenceImageFootprintInformationList, policy):
    kernelCols = policy.get('kernelCols')
    kernelRows = policy.get('kernelRows')
    
    # how many good footprints are we dealing with here
    goodDifiList = ipDiffim.getGoodFootprints(differenceImageFootprintInformationList)
    nFootprint   = len(goodDifiList)

    # matrix to invert
    M = numpy.zeros((kernelCols*kernelRows, nFootprint))
    for i in range(nFootprint):
        singleKernelImage  = goodDifiList[i].getSingleKernelPtr().computeNewImage(False)[0]
        singleKernelVector = ipDiffimTools.imageToVector(singleKernelImage)
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

    # Turn into Kernels
    meanKernelPtr    = ipDiffimTools.vectorToKernelPtr( meanM, kernelCols, kernelRows )
    eKernelPtrVector = afwMath.VectorKernel()
    for i in range(U.shape[1]):
        eKernelPtr   = ipDiffimTools.vectorToKernelPtr( U[:,i], kernelCols, kernelRows )
        eKernelPtrVector.append(eKernelPtr)

    return meanKernelPtr, eKernelPtrVector, eVal, eCoeff


def fitPerPixel(differenceImageFootprintInformationList, label, policy):
    kernelSpatialOrder        = policy.get('kernelSpatialOrder')
    backgroundSpatialOrder    = policy.get('backgroundSpatialOrder')
    kernelCols                = policy.get('kernelCols')
    kernelRows                = policy.get('kernelRows')

    # how many good footprints are we dealing with here
    goodDifiList = ipDiffim.getGoodFootprints(differenceImageFootprintInformationList)
    nFootprint   = len(goodDifiList)

    # common to all spatial fits; position of constraints
    footprintExposureCol = numpy.zeros(nFootprint)
    footprintExposureRow = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        footprintExposureCol[i] = goodDifiList[i].getColcNorm()
        footprintExposureRow[i] = goodDifiList[i].getRowcNorm()

    # fit the background
    backgroundValues = numpy.zeros(nFootprint)
    backgroundErrors = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        backgroundValues[i] = goodDifiList[i].getSingleBackground()
        backgroundErrors[i] = goodDifiList[i].getSingleBackgroundError()

    backgroundSpatialFunction = afwMath.PolynomialFunction2D(backgroundSpatialOrder)
    backgroundFunctionFit = fitSpatialFunction(backgroundSpatialFunction,
                                               backgroundValues,
                                               backgroundErrors,
                                               footprintExposureCol,
                                               footprintExposureRow,
                                               policy)
    Trace('lsst.ip.diffim', 5,
          'Spatial %d background fit parameters : %s' % (label, ' '.join([('%10.3e' % (x)) for x in backgroundFunctionFit.parameterList])))
    backgroundSpatialFunction.setParameters(backgroundFunctionFit.parameterList)
    
    # fit each pixel
    j = 0
    functionFitList = []
    for kCol in range(kernelCols):
        for kRow in range(kernelRows):

            #######
            # initialize vectors, one per good kernel
            
            kernelValues = numpy.zeros(nFootprint)
            kernelErrors = numpy.zeros(nFootprint)
            for i in range(nFootprint):
                singleKernelPtr      = goodDifiList[i].getSingleKernelPtr()
                singleKernelErrorPtr = goodDifiList[i].getSingleKernelErrorPtr()
                kernelValues[i]      = singleKernelPtr.computeNewImage(False)[0].getVal(kCol, kRow)
                kernelErrors[i]      = singleKernelErrorPtr.computeNewImage(False)[0].getVal(kCol, kRow)
                    
            # initialize vectors, one per good kernel
            #######
            # do the various fitting techniques here
            kernelSpatialFunction = afwMath.PolynomialFunction2D(kernelSpatialOrder)
            functionFit = fitSpatialFunction(kernelSpatialFunction,
                                             kernelValues,
                                             kernelErrors,
                                             footprintExposureCol,
                                             footprintExposureRow,
                                             policy)
            Trace('lsst.ip.diffim', 5,
                  'Kernel %s spatial pixel %d fit parameters : %s' % (label, j, ' '.join([('%10.3e' % (x)) for x in functionFit.parameterList])))
            kernelSpatialFunction.setParameters(functionFit.parameterList)
            functionFitList.append(kernelSpatialFunction)

            j += 1
            # do the various fitting techniques here
            #######

    # Evaluate all the fits at the positions of the objects, create a
    # new kernel, then difference image, the calculate difference
    # image stats
    for i in range(nFootprint):
        kFunctionImage          = afwImage.ImageD(kernelCols, kernelRows)
        functionBackgroundValue = backgroundSpatialFunction(footprintExposureCol[i], footprintExposureRow[i])

        j = 0
        for kCol in range(kernelCols):
            for kRow in range(kernelRows):
                functionValue   = functionFitList[j](footprintExposureCol[i], footprintExposureRow[i])
                kFunctionImage.set(kCol, kRow, functionValue)
                j += 1

        if policy.get('debugIO') == True:
            kFunctionImage.writeFits('SKernel_%s%d.fits' % (label, goodDifiList[i].getID()))
    
        functionKernelPtr  = afwMath.KernelPtr( afwMath.FixedKernel(kFunctionImage) )
        spatialDiffIm      = ipDiffim.convolveAndSubtract(goodDifiList[i].getImageToConvolvePtr().get(),
                                                          goodDifiList[i].getImageToNotConvolvePtr().get(),
                                                          functionKernelPtr,
                                                          functionBackgroundValue)
        diffImStatistics   = ipDiffim.DifferenceImageStatisticsF(spatialDiffIm)
                                                     
        if policy.get('debugPlot') == True:
            ipDiffimDebug.plotDiffImQuality1(spatialDiffIm,
                                             functionKernelPtr,
                                             goodDifiList[i].getImageToConvolvePtr().get(),
                                             goodDifiList[i].getImageToNotConvolvePtr().get(),
                                             'Spatial %s kernel %d' % (label, goodDifiList[i].getID()),
                                             'SKernel_%s%d.ps' % (label, goodDifiList[i].getID())
                                             )
    



def fitSpatialPca(differenceImageFootprintInformationList, label, mKernel, eKernels, eCoeffs, policy):
    kernelSpatialOrder     = policy.get('kernelSpatialOrder')
    backgroundSpatialOrder = policy.get('backgroundSpatialOrder')
    kernelCols             = policy.get('kernelCols')
    kernelRows             = policy.get('kernelRows')

    # how many good footprints are we dealing with here
    goodDifiList  = ipDiffim.getGoodFootprints(differenceImageFootprintInformationList)
    nFootprint    = len(goodDifiList)
    nCoefficient  = nFootprint

    # common to all spatial fits
    footprintExposureCol  = numpy.zeros(nFootprint)
    footprintExposureRow  = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        footprintExposureCol[i] = goodDifiList[i].getColcNorm()
        footprintExposureRow[i] = goodDifiList[i].getRowcNorm()

    # fit the background
    backgroundValues = numpy.zeros(nFootprint)
    backgroundErrors = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        backgroundValues[i] = goodDifiList[i].getSingleBackground()
        backgroundErrors[i] = goodDifiList[i].getSingleBackgroundError()

    backgroundSpatialFunction = afwMath.PolynomialFunction2D(backgroundSpatialOrder)
    backgroundFunctionFit = fitSpatialFunction(backgroundSpatialFunction,
                                               backgroundValues,
                                               backgroundErrors,
                                               footprintExposureCol,
                                               footprintExposureRow,
                                               policy)
    Trace('lsst.ip.diffim', 5,
          'Spatial PCA %s background fit parameters : %s' % (label, ' '.join([('%10.3e' % (x)) for x in backgroundFunctionFit.parameterList])))
    backgroundSpatialFunction.setParameters(backgroundFunctionFit.parameterList)

    # fit spatial variation of each eigenkernel
    functionFitList = []
    for i in range(nCoefficient):
        coefficients = eCoeffs[i,:]
        # NOTE TO SELF : WHAT TO DO ABOUT THEIR UNCERTAINTY?
        # THIS IS A HACK TO GET THINGS MOVING
        uncertainties = numpy.sqrt( numpy.abs(coefficients) )
        
        kernelSpatialFunction = afwMath.PolynomialFunction2D(kernelSpatialOrder)
        coeffFunctionFit = fitSpatialFunction(kernelSpatialFunction,
                                              coefficients,
                                              uncertainties,
                                              footprintExposureCol,
                                              footprintExposureRow,
                                              policy)
        Trace('lsst.ip.diffim', 5,
              'Kernel %s spatial PCA %d fit parameters : %s' % (label, i, ' '.join([('%10.3e' % (x)) for x in coeffFunctionFit.parameterList])))
        kernelSpatialFunction.setParameters(coeffFunctionFit.parameterList)
        functionFitList.append(kernelSpatialFunction)
        
    # now see how well this approximates the diffim
    for i in range(nFootprint):
        backgroundValue   = backgroundSpatialFunction(footprintExposureCol[i], footprintExposureRow[i])
        approxKernelImage = mKernel.copy()
        
        for j in range(nCoefficient):
            coeffFunctionValue = kernelSpatialFunction[j](footprintExposureCol[i], footprintExposureRow[i])
            approxKernelImage += coeffFunctionValue * eKernels[j,:]
            approxKernelPtr    = ipDiffimTools.vectorToKernelPtr(approxKernelImage, kernelCols, kernelRows)

            spatialDiffIm      = ipDiffim.convolveAndSubtract(goodDifiList[i].getImageToConvolvePtr().get(),
                                                              goodDifiList[i].getImageToNotConvolvePtr().get(),
                                                              approxKernelPtr,
                                                              backgroundValue)
            diffImStatistics   = ipDiffim.DifferenceImageStatisticsF(spatialDiffIm)
            
            if policy.get('debugPlot') == True:
                ipDiffimDebug.plotDiffImQuality1(spatialDiffIm,
                                                 approxKernelPtr,
                                                 goodDifiList[i].getImageToConvolvePtr().get(),
                                                 goodDifiList[i].getImageToNotConvolvePtr().get(),
                                                 'PCA %s kernel_%d ek_%d' % (label, goodDifiList[i].getID(), j),
                                                 'PKernel_%s%d_%d.ps' % (label, goodDifiList[i].getID(), j)
                                                 )
                
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
    #policy.set('getCollectionOfFootprintsForPsfMatching.footprintDetectionThreshold', 6350.)  # gets 5
    policy.set('getCollectionOfFootprintsForPsfMatching.footprintDetectionThreshold', 5000.)  # gets 12
    #policy.set('getCollectionOfFootprintsForPsfMatching.footprintDetectionThreshold', 1000.) # gets full test suite
    #policy.set('getCollectionOfFootprintsForPsfMatching.footprintDetectionThreshold', 250.) # gets a lot
    
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
        imageStampPtr    = scienceMaskedImage.getSubImage(footprintBBox)

        # initial estimate of the variance per pixel, straight subtraction
        subtractedStamp  = afwImage.MaskedImageF(templateStampPtr.getCols(), templateStampPtr.getRows())
        subtractedStamp += imageStampPtr.get()
        subtractedStamp -= templateStampPtr.get()

        # Assuming the template is better seeing, find convolution kernel
        convKernelVectorPair1 = ipDiffim.computePsfMatchingKernelForFootprint2(
            templateStampPtr.get(),
            imageStampPtr.get(),
            subtractedStamp,
            kernelBasisList,
            policy
        )
        convKernelVector1, convKernelErrorVector1, convBackground1, convBackgroundError1 = ipDiffimTools.vectorPairToVectors(convKernelVectorPair1)
        convKernelPtr1 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kernelBasisList, convKernelVector1)
            )
        convKernelErrorPtr1 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kernelBasisList, convKernelErrorVector1)
            )
        convDiffIm1    = ipDiffim.convolveAndSubtract(templateStampPtr.get(), imageStampPtr.get(), convKernelPtr1, convBackground1)
        #
        #### NOW REDO THIS USING A BETTER ESTIMATE OF THE VARIANCE FROM THE SUBTRACTED IMAGE!
        #
        convKernelVectorPair2 = ipDiffim.computePsfMatchingKernelForFootprint2(
            templateStampPtr.get(),
            imageStampPtr.get(),
            convDiffIm1,
            kernelBasisList,
            policy
        )
        convKernelVector2, convKernelErrorVector2, convBackground2, convBackgroundError2 = ipDiffimTools.vectorPairToVectors(convKernelVectorPair2)
        convKernelPtr2 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kernelBasisList, convKernelVector2)
            )
        convKernelErrorPtr2 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kernelBasisList, convKernelErrorVector2)
            )
        convDiffIm2    = ipDiffim.convolveAndSubtract(templateStampPtr.get(), imageStampPtr.get(), convKernelPtr2, convBackground2)
        # Things tend to converge after an single iteration, so just do this once.

        

        # Assuming the template is better seeing, find deconvolution kernel
        deconvKernelVectorPair1 = ipDiffim.computePsfMatchingKernelForFootprint2(
            imageStampPtr.get(),
            templateStampPtr.get(),
            subtractedStamp,
            kernelBasisList,
            policy
        )
        deconvKernelVector1, deconvKernelErrorVector1, deconvBackground1, deconvBackgroundError1 = ipDiffimTools.vectorPairToVectors(deconvKernelVectorPair1)
        deconvKernelPtr1 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kernelBasisList, deconvKernelVector1)
            )
        deconvKernelErrorPtr1 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kernelBasisList, deconvKernelErrorVector1)
            )
        deconvDiffIm1    = ipDiffim.convolveAndSubtract(imageStampPtr.get(), templateStampPtr.get(), deconvKernelPtr1, deconvBackground1)
        #
        #### NOW REDO THIS USING A BETTER ESTIMATE OF THE VARIANCE FROM THE SUBTRACTED IMAGE!
        #
        deconvKernelVectorPair2 = ipDiffim.computePsfMatchingKernelForFootprint2(
            imageStampPtr.get(),
            templateStampPtr.get(),
            deconvDiffIm1,
            kernelBasisList,
            policy,
        )
        deconvKernelVector2, deconvKernelErrorVector2, deconvBackground2, deconvBackgroundError2 = ipDiffimTools.vectorPairToVectors(deconvKernelVectorPair2)
        deconvKernelPtr2 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kernelBasisList, deconvKernelVector2)
            )
        deconvKernelErrorPtr2 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kernelBasisList, deconvKernelErrorVector2)
            )
        deconvDiffIm2    = ipDiffim.convolveAndSubtract(imageStampPtr.get(), templateStampPtr.get(), deconvKernelPtr2, deconvBackground2)
        # Things tend to converge after an single iteration, so just do this once.


        # Create and fill Difi
        convDifiPtr   = ipDiffim.DifiPtrF(
            ipDiffim.DifferenceImageFootprintInformationF(iFootprintPtr, templateStampPtr, imageStampPtr)
            ) 
        convDifiPtr.setID(footprintID)
        convDifiPtr.setColcNorm( float(fpMin.x() + fpMax.x()) / templateMaskedImage.getCols() - 1.0 )
        convDifiPtr.setRowcNorm( float(fpMin.y() + fpMax.y()) / templateMaskedImage.getRows() - 1.0 ) 
        convDifiPtr.setSingleKernelPtr( convKernelPtr2 )
        convDifiPtr.setSingleKernelErrorPtr( convKernelErrorPtr2 )
        convDifiPtr.setSingleBackground( convBackground2  )
        convDifiPtr.setSingleBackgroundError( convBackgroundError2 )
        convDifiPtr.setSingleStats(
            convDifiPtr.computeImageStatistics(convKernelPtr2, convBackground2)
            )
        convDifiPtr.setStatus( convDifiPtr.getSingleStats().evaluateQuality(policy) )
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
        convolveDifiPtrList.append(convDifiPtr)

        #

        deconvDifiPtr = ipDiffim.DifiPtrF(
            ipDiffim.DifferenceImageFootprintInformationF(iFootprintPtr, imageStampPtr, templateStampPtr)
            )
        deconvDifiPtr.setID(footprintID)
        deconvDifiPtr.setColcNorm( float(fpMin.x() + fpMax.x()) / templateMaskedImage.getCols() - 1.0 )
        deconvDifiPtr.setRowcNorm( float(fpMin.y() + fpMax.y()) / templateMaskedImage.getRows() - 1.0 ) 
        deconvDifiPtr.setSingleKernelPtr( deconvKernelPtr2 )
        deconvDifiPtr.setSingleKernelErrorPtr( deconvKernelErrorPtr2 )
        deconvDifiPtr.setSingleBackground( deconvBackground2 )
        deconvDifiPtr.setSingleBackgroundError( deconvBackgroundError2 )
        deconvDifiPtr.setSingleStats(
            deconvDifiPtr.computeImageStatistics(deconvKernelPtr2, deconvBackground2)
            )
        deconvDifiPtr.setStatus( deconvDifiPtr.getSingleStats().evaluateQuality(policy) )
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
        deconvolveDifiPtrList.append(deconvDifiPtr)

        # Debugging information that is once per kernel
        if policy.get('debugPlot') == True:
            ipDiffimDebug.plotDiffImQuality2(footprintID, 1,
                                             convDiffIm1, convKernelPtr1, templateStampPtr, imageStampPtr, 
                                             deconvDiffIm1, deconvKernelPtr1, imageStampPtr, templateStampPtr)
            ipDiffimDebug.plotDiffImQuality2(footprintID, 2,
                                             convDiffIm2, convKernelPtr2, templateStampPtr, imageStampPtr, 
                                             deconvDiffIm2, deconvKernelPtr2, imageStampPtr, templateStampPtr)

        if policy.get('debugIO') == True:
            ipDiffimDebug.writeDiffImages(footprintID,
                                          templateStampPtr, imageStampPtr,
                                          convDifiPtr, convDiffIm2, convKernelPtr2,
                                          deconvDifiPtr, deconvDiffIm2, deconvKernelPtr2)
            
    rejectKernelOutliers(convolveDifiPtrList, policy)
    rejectKernelOutliers(deconvolveDifiPtrList, policy)

    # now we get to the good stuff!
    convMKernel,   convEKernels,   convEVals,   convECoeffs   = runPca(convolveDifiPtrList, policy)
    deconvMKernel, deconvEKernels, deconvEVals, deconvECoeffs = runPca(deconvolveDifiPtrList, policy)

    Trace('lsst.ip.diffim', 5,
          'EigenValues 1 : %s' % (' '.join([str(x) for x in convEVals])))
    Trace('lsst.ip.diffim', 5,
          'EigenValues 2 : %s' % (' '.join([str(x) for x in deconvEVals])))
    
    if policy.get('debugPlot') == True:
        ipDiffimPlot.eigenKernelPlot( (convMKernel,   convEKernels,   convEVals),
                                      (deconvMKernel, deconvEKernels, deconvEVals),
                                      outfile='EKernel.ps'
                                      )

        goodDifiList = ipDiffim.getGoodFootprints(convolveDifiPtrList)
        nFootprint   = len(goodDifiList)
        for i in range(nFootprint):
            ipDiffimPlot.approxKernelPlot(goodDifiList[i],
                                          (convMKernel, convEKernels),
                                          title='Kernel %d conv' % (goodDifiList[i].getID()),
                                          outfile='PCA_Kc_%d.ps' % (goodDifiList[i].getID())
                                          )

        goodDifiList = ipDiffim.getGoodFootprints(deconvolveDifiPtrList)
        nFootprint   = len(goodDifiList)
        for i in range(nFootprint):
            ipDiffimPlot.approxKernelPlot(goodDifiList[i],
                                          (deconvMKernel, deconvEKernels),
                                          title='Kernel %d dconv' % (goodDifiList[i].getID()),
                                          outfile='PCA_Kdc_%d.ps' % (goodDifiList[i].getID())
                                          )
        

    convResiduals2   = fitSpatialPca(convolveDifiPtrList, 'conv', convMKernel, convEKernels, convECoeffs, policy)
    deconvResiduals2 = fitSpatialPca(deconvolveDifiPtrList, 'deconv', deconvMKernel, deconvEKernels, deconvECoeffs, policy)

    convResiduals1   = fitPerPixel(convolveDifiPtrList, 'conv', policy)
    deconvResiduals1 = fitPerPixel(deconvolveDifiPtrList, 'deconc', policy)


    return
    
    if policy.get('debugPlot') == True:
        ipDiffimPlot.plotBackground(convResiduals1, outfile='cBackground.ps')
        ipDiffimPlot.plotBackground(deconvResiduals1, outfile='dcBackground.ps')


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
    
    
