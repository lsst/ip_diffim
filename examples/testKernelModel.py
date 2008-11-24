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
    kCols = policy.get('kernelCols')
    kRows = policy.get('kernelRows')
    
    # how many good footprints are we dealing with here
    goodDifiList = ipDiffim.getGoodFootprints(differenceImageFootprintInformationList)
    nFootprint   = len(goodDifiList)

    # matrix to invert
    M = numpy.zeros((kCols*kRows, nFootprint))
    for i in range(nFootprint):
        kernelImage  = goodDifiList[i].getSingleKernelPtr().computeNewImage(False)[0]
        M[:,i]       = ipDiffimTools.imageToVector(kernelImage)
        
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
    meanKernelPtr    = ipDiffimTools.vectorToKernelPtr( meanM, kCols, kRows )
    eKernelPtrVector = afwMath.VectorKernel()
    for i in range(U.shape[1]):
        eKernelPtr   = ipDiffimTools.vectorToKernelPtr( U[:,i], kCols, kRows )
        eKernelPtrVector.append(eKernelPtr)

    return meanKernelPtr, eKernelPtrVector, eVal, eCoeff


def fitPerPixel(differenceImageFootprintInformationList, label, policy):
    kSpatialOrder  = policy.get('kernelSpatialOrder')
    bgSpatialOrder = policy.get('backgroundSpatialOrder')
    kCols          = policy.get('kernelCols')
    kRows          = policy.get('kernelRows')

    # how many good footprints are we dealing with here
    goodDifiList = ipDiffim.getGoodFootprints(differenceImageFootprintInformationList)
    nFootprint   = len(goodDifiList)

    # common to all spatial fits; position of constraints
    fpExposureCol = numpy.zeros(nFootprint)
    fpExposureRow = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        fpExposureCol[i] = goodDifiList[i].getColcNorm()
        fpExposureRow[i] = goodDifiList[i].getRowcNorm()

    # fit the background
    bgValues = numpy.zeros(nFootprint)
    bgErrors = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        bgValues[i] = goodDifiList[i].getSingleBackground()
        bgErrors[i] = goodDifiList[i].getSingleBackgroundError()

    bgFunction = afwMath.PolynomialFunction2D(bgSpatialOrder)
    bgFit = fitSpatialFunction(bgFunction,
                               bgValues,
                               bgErrors,
                               fpExposureCol,
                               fpExposureRow,
                               policy)
    Trace('lsst.ip.diffim', 5,
          'Spatial %s background fit parameters : %s' % (label, ' '.join([('%10.3e' % (x)) for x in bgFit.parameterList])))
    bgFunction.setParameters(bgFit.parameterList)
    
    # fit each pixel
    np = 0
    fitList = []
    for kCol in range(kCols):
        for kRow in range(kRows):

            #######
            # initialize vectors, one per good kernel
            
            kValues = numpy.zeros(nFootprint)
            kErrors = numpy.zeros(nFootprint)
            for i in range(nFootprint):
                singleKernelPtr      = goodDifiList[i].getSingleKernelPtr()
                singleKernelErrorPtr = goodDifiList[i].getSingleKernelErrorPtr()
                kValues[i]           = singleKernelPtr.computeNewImage(False)[0].getVal(kCol, kRow)
                kErrors[i]           = singleKernelErrorPtr.computeNewImage(False)[0].getVal(kCol, kRow)
                    
            # initialize vectors, one per good kernel
            #######
            # do the various fitting techniques here
            kFunction = afwMath.PolynomialFunction2D(kSpatialOrder)
            kFit = fitSpatialFunction(kFunction,
                                      kValues,
                                      kErrors,
                                      fpExposureCol,
                                      fpExposureRow,
                                      policy)
            
            Trace('lsst.ip.diffim', 5,
                  'Kernel %s spatial pixel %d fit parameters : %s' % (label, np, ' '.join([('%10.3e' % (x)) for x in kFit.parameterList])))
            kFunction.setParameters(kFit.parameterList)
            fitList.append(kFunction)

            np += 1
            # do the various fitting techniques here
            #######

    # Evaluate all the fits at the positions of the objects, create a
    # new kernel, then difference image, the calculate difference
    # image stats
    for i in range(nFootprint):
        kImage  = afwImage.ImageD(kCols, kRows)
        bgValue = bgFunction(fpExposureCol[i], fpExposureRow[i])

        np = 0
        for kCol in range(kCols):
            for kRow in range(kRows):
                pixelValue = fitList[np](fpExposureCol[i], fpExposureRow[i])
                kImage.set(kCol, kRow, pixelValue)
                np += 1

        if policy.get('debugIO') == True:
            kImage.writeFits('SKernel_%s%d.fits' % (label, goodDifiList[i].getID()))
    
        kernelPtr = afwMath.KernelPtr( afwMath.FixedKernel(kImage) )
        diffIm = ipDiffim.convolveAndSubtract(goodDifiList[i].getImageToConvolvePtr().get(),
                                              goodDifiList[i].getImageToNotConvolvePtr().get(),
                                              kernelPtr, bgValue)
        diffImStats = ipDiffim.DifferenceImageStatisticsF(diffIm)
                                                     
        if policy.get('debugPlot') == True:
            ipDiffimDebug.plotDiffImQuality1(goodDifiList[i],
                                             diffIm,
                                             kernelPtr,
                                             label='Spatial %s kernel %d' % (label, goodDifiList[i].getID()),
                                             outfile='SKernel_%s%d.ps' % (label, goodDifiList[i].getID())
                                             )
    



def fitSpatialPca(differenceImageFootprintInformationList, label, mKernel, eKernels, eCoeffs, policy):
    kSpatialOrder  = policy.get('kernelSpatialOrder')
    bgSpatialOrder = policy.get('backgroundSpatialOrder')
    kCols          = policy.get('kernelCols')
    kRows          = policy.get('kernelRows')

    # how many good footprints are we dealing with here
    goodDifiList  = ipDiffim.getGoodFootprints(differenceImageFootprintInformationList)
    nFootprint    = len(goodDifiList)
    nCoefficient  = nFootprint

    # common to all spatial fits
    fpExposureCol  = numpy.zeros(nFootprint)
    fpExposureRow  = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        fpExposureCol[i] = goodDifiList[i].getColcNorm()
        fpExposureRow[i] = goodDifiList[i].getRowcNorm()

    # fit the background
    bgValues = numpy.zeros(nFootprint)
    bgErrors = numpy.zeros(nFootprint)
    for i in range(nFootprint):
        bgValues[i] = goodDifiList[i].getSingleBackground()
        bgErrors[i] = goodDifiList[i].getSingleBackgroundError()

    bgFunction = afwMath.PolynomialFunction2D(bgSpatialOrder)
    bgFit = fitSpatialFunction(bgFunction,
                               bgValues,
                               bgErrors,
                               fpExposureCol,
                               fpExposureRow,
                               policy)
    Trace('lsst.ip.diffim', 5,
          'Spatial PCA %s background fit parameters : %s' % (label, ' '.join([('%10.3e' % (x)) for x in bgFit.parameterList])))
    bgFunction.setParameters(bgFit.parameterList)

    # fit spatial variation of each eigenkernel
    fitList = []
    for npc in range(nCoefficient):
        coefficients = eCoeffs[npc,:]
        # NOTE TO SELF : WHAT TO DO ABOUT THEIR UNCERTAINTY?
        # THIS IS A HACK TO GET THINGS MOVING
        uncertainties = numpy.sqrt( numpy.abs(coefficients) )
        
        kFunction = afwMath.PolynomialFunction2D(kSpatialOrder)
        kFit = fitSpatialFunction(kFunction,
                                  coefficients,
                                  uncertainties,
                                  fpExposureCol,
                                  fpExposureRow,
                                  policy)
        Trace('lsst.ip.diffim', 5,
              'Kernel %s spatial PCA %d fit parameters : %s' % (label, npc, ' '.join([('%10.3e' % (x)) for x in kFit.parameterList])))
        kFunction.setParameters(kFit.parameterList)
        fitList.append(kFunction)
        
    # now see how well this approximates the diffim
    for i in range(nFootprint):
        bgValue = bgFunction(fpExposureCol[i], fpExposureRow[i])
        kernelImage = mKernel.computeNewImage(False)[0]
        
        for npc in range(nCoefficient):
            coeff        = fitList[npc](fpExposureCol[i], fpExposureRow[i])
            eigenImage   = eKernels[npc].computeNewImage(False)[0]
            eigenImage  *= coeff
            kernelImage += eigenImage
            kernelPtr    = afwMath.KernelPtr( afwMath.FixedKernel(kernelImage) ) 

            diffIm       = ipDiffim.convolveAndSubtract(goodDifiList[i].getImageToConvolvePtr().get(),
                                                        goodDifiList[i].getImageToNotConvolvePtr().get(),
                                                        kernelPtr, bgValue)
            diffImStats  = ipDiffim.DifferenceImageStatisticsF(diffIm)
            
            if policy.get('debugPlot') == True:
                ipDiffimDebug.plotDiffImQuality1(goodDifiList[i],
                                                 diffIm,
                                                 kernelPtr,
                                                 label = 'PCA %s kernel_%d ek_%d' % (label, goodDifiList[i].getID(), npc),
                                                 outfile = 'PKernel_%s%d_%d.ps' % (label, goodDifiList[i].getID(), npc)
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

    kCols = policy.get('kernelCols')
    kRows = policy.get('kernelRows')
    
    if options.verbosity > 0:
        print 'Verbosity =', options.verbosity
        Trace.setVerbosity('lsst.ip.diffim', options.verbosity)

    kBasisList = ipDiffim.generateDeltaFunctionKernelSet(kCols, kRows)
    
    # lets just get a couple for debugging and speed
    policy.set('getCollectionOfFootprintsForPsfMatching.minimumCleanFootprints', 5)
    #policy.set('getCollectionOfFootprintsForPsfMatching.footprintDetectionThreshold', 6350.)  # gets 5
    #policy.set('getCollectionOfFootprintsForPsfMatching.footprintDetectionThreshold', 5000.)  # gets 12
    policy.set('getCollectionOfFootprintsForPsfMatching.footprintDetectionThreshold', 1000.) # gets full test suite
    #policy.set('getCollectionOfFootprintsForPsfMatching.footprintDetectionThreshold', 250.) # gets a lot
    
    fpList = ipDiffim.getCollectionOfFootprintsForPsfMatching(templateMaskedImage,
                                                              scienceMaskedImage,
                                                              policy)
    cDifiPtrList = ipDiffim.DifiPtrListF()
    dDifiPtrList = ipDiffim.DifiPtrListF()
    for fpID, fpPtr in enumerate(fpList):
        fpBBox = fpPtr.getBBox()
        fpMin  = fpBBox.min()
        fpMax  = fpBBox.max()
        
        Trace('lsst.ip.diffim.computePsfMatchingKernelForFootprint', 2,
              'Footprint %d = %d,%d -> %d,%d' % (fpID,
                                                 fpBBox.min().x(), fpBBox.min().y(),
                                                 fpBBox.max().x(), fpBBox.max().y()))
        templatePtr = templateMaskedImage.getSubImage(fpBBox)
        imagePtr    = scienceMaskedImage.getSubImage(fpBBox)

        # initial estimate of the variance per pixel, straight subtraction
        diffIm0  = afwImage.MaskedImageF(templatePtr.getCols(), templatePtr.getRows())
        diffIm0 += imagePtr.get()
        diffIm0 -= templatePtr.get()

        # Assuming the template is better seeing, find convolution kernel
        convVP1 = ipDiffim.computePsfMatchingKernelForFootprint2(
            templatePtr.get(),
            imagePtr.get(),
            diffIm0,
            kBasisList,
            policy
        )
        cKernelVec1, cKernelErrorVec1, cBg1, cBgError1 = ipDiffimTools.vectorPairToVectors(convVP1)
        cKernelPtr1 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kBasisList, cKernelVec1)
            )
        cKernelErrorPtr1 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kBasisList, cKernelErrorVec1)
            )
        cDiffIm1 = ipDiffim.convolveAndSubtract(templatePtr.get(), imagePtr.get(), cKernelPtr1, cBg1)
        #
        #### NOW REDO THIS USING A BETTER ESTIMATE OF THE VARIANCE FROM THE SUBTRACTED IMAGE!
        #
        convVP2 = ipDiffim.computePsfMatchingKernelForFootprint2(
            templatePtr.get(),
            imagePtr.get(),
            cDiffIm1,
            kBasisList,
            policy
        )
        cKernelVec2, cKernelErrorVec2, cBg2, cBgError2 = ipDiffimTools.vectorPairToVectors(convVP2)
        cKernelPtr2 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kBasisList, cKernelVec2)
            )
        cKernelErrorPtr2 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kBasisList, cKernelErrorVec2)
            )
        cDiffIm2 = ipDiffim.convolveAndSubtract(templatePtr.get(), imagePtr.get(), cKernelPtr2, cBg2)
        # Things tend to converge after an single iteration, so just do this once.

        

        # Assuming the template is better seeing, find deconvolution kernel
        deconvVP1 = ipDiffim.computePsfMatchingKernelForFootprint2(
            imagePtr.get(),
            templatePtr.get(),
            diffIm0,
            kBasisList,
            policy
        )
        dKernelVec1, dKernelErrorVec1, dBg1, dBgError1 = ipDiffimTools.vectorPairToVectors(deconvVP1)
        dKernelPtr1 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kBasisList, dKernelVec1)
            )
        dKernelErrorPtr1 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kBasisList, dKernelErrorVec1)
            )
        dDiffIm1 = ipDiffim.convolveAndSubtract(imagePtr.get(), templatePtr.get(), dKernelPtr1, dBg1)
        #
        #### NOW REDO THIS USING A BETTER ESTIMATE OF THE VARIANCE FROM THE SUBTRACTED IMAGE!
        #
        deconvVP2 = ipDiffim.computePsfMatchingKernelForFootprint2(
            imagePtr.get(),
            templatePtr.get(),
            dDiffIm1,
            kBasisList,
            policy,
        )
        dKernelVec2, dKernelErrorVec2, dBg2, dBgError2 = ipDiffimTools.vectorPairToVectors(deconvVP2)
        dKernelPtr2 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kBasisList, dKernelVec2)
            )
        dKernelErrorPtr2 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(kBasisList, dKernelErrorVec2)
            )
        dDiffIm2 = ipDiffim.convolveAndSubtract(imagePtr.get(), templatePtr.get(), dKernelPtr2, dBg2)
        # Things tend to converge after an single iteration, so just do this once.


        # Create and fill Difi
        cDifiPtr   = ipDiffim.DifiPtrF(
            ipDiffim.DifferenceImageFootprintInformationF(fpPtr, templatePtr, imagePtr)
            ) 
        cDifiPtr.setID(fpID)
        cDifiPtr.setColcNorm( float(fpMin.x() + fpMax.x()) / templateMaskedImage.getCols() - 1.0 )
        cDifiPtr.setRowcNorm( float(fpMin.y() + fpMax.y()) / templateMaskedImage.getRows() - 1.0 ) 
        cDifiPtr.setSingleKernelPtr( cKernelPtr2 )
        cDifiPtr.setSingleKernelErrorPtr( cKernelErrorPtr2 )
        cDifiPtr.setSingleBackground( cBg2  )
        cDifiPtr.setSingleBackgroundError( cBgError2 )
        cDifiPtr.setSingleStats(
            cDifiPtr.computeImageStatistics(cKernelPtr2, cBg2)
            )
        cDifiPtr.setStatus( cDifiPtr.getSingleStats().evaluateQuality(policy) )
        if cDifiPtr.getStatus() == True:
            prefix = ''
        else:
            prefix = '#'
        Trace('lsst.ip.diffim', 5,
              '%sKernel %d : Kernel Sum = %.2f, Diffim residuals = %.2f +/- %.2f sigma' % (
            prefix,
            cDifiPtr.getID(),
            cDifiPtr.getSingleKernelPtr().computeNewImage(False)[1],
            cDifiPtr.getSingleStats().getResidualMean(),
            cDifiPtr.getSingleStats().getResidualStd()
            ))
        cDifiPtrList.append(cDifiPtr)

        #

        dDifiPtr = ipDiffim.DifiPtrF(
            ipDiffim.DifferenceImageFootprintInformationF(fpPtr, imagePtr, templatePtr)
            )
        dDifiPtr.setID(fpID)
        dDifiPtr.setColcNorm( float(fpMin.x() + fpMax.x()) / templateMaskedImage.getCols() - 1.0 )
        dDifiPtr.setRowcNorm( float(fpMin.y() + fpMax.y()) / templateMaskedImage.getRows() - 1.0 ) 
        dDifiPtr.setSingleKernelPtr( dKernelPtr2 )
        dDifiPtr.setSingleKernelErrorPtr( dKernelErrorPtr2 )
        dDifiPtr.setSingleBackground( dBg2 )
        dDifiPtr.setSingleBackgroundError( dBgError2 )
        dDifiPtr.setSingleStats(
            dDifiPtr.computeImageStatistics(dKernelPtr2, dBg2)
            )
        dDifiPtr.setStatus( dDifiPtr.getSingleStats().evaluateQuality(policy) )
        if dDifiPtr.getStatus() == True:
            prefix = ''
        else:
            prefix = '#'
        Trace('lsst.ip.diffim', 5,
              '%sKernel %d : Kernel Sum = %.2f, Diffim residuals = %.2f +/- %.2f sigma' % (
            prefix,
            dDifiPtr.getID(),
            dDifiPtr.getSingleKernelPtr().computeNewImage(False)[1],
            dDifiPtr.getSingleStats().getResidualMean(),
            dDifiPtr.getSingleStats().getResidualStd()
            ))
        dDifiPtrList.append(dDifiPtr)

        # Debugging information that is once per kernel
        if policy.get('debugPlot') == True:
            ipDiffimDebug.plotDiffImQuality2(fpID, 1,
                                             cDiffIm1, cKernelPtr1, templatePtr, imagePtr, 
                                             dDiffIm1, dKernelPtr1, imagePtr, templatePtr)
            ipDiffimDebug.plotDiffImQuality2(fpID, 2,
                                             cDiffIm2, cKernelPtr2, templatePtr, imagePtr, 
                                             dDiffIm2, dKernelPtr2, imagePtr, templatePtr)

        if policy.get('debugIO') == True:
            ipDiffimDebug.writeDiffImages(fpID,
                                          templatePtr, imagePtr,
                                          cDifiPtr, cDiffIm2, cKernelPtr2,
                                          dDifiPtr, dDiffIm2, dKernelPtr2)
            
    rejectKernelOutliers(cDifiPtrList, policy)
    rejectKernelOutliers(dDifiPtrList, policy)

    # now we get to the good stuff!
    cMKernel, cEKernels, cEVals, cECoeffs = runPca(cDifiPtrList, policy)
    dMKernel, dEKernels, dEVals, dECoeffs = runPca(dDifiPtrList, policy)

    Trace('lsst.ip.diffim', 5,
          'EigenValues 1 : %s' % (' '.join([str(x) for x in cEVals])))
    Trace('lsst.ip.diffim', 5,
          'EigenValues 2 : %s' % (' '.join([str(x) for x in dEVals])))
    
    if policy.get('debugPlot') == True:
        ipDiffimPlot.eigenKernelPlot( (cMKernel, cEKernels, cEVals),
                                      (dMKernel, dEKernels, dEVals),
                                      outfile='EKernel.ps'
                                      )

        goodDifiList = ipDiffim.getGoodFootprints(cDifiPtrList)
        nFootprint   = len(goodDifiList)
        for i in range(nFootprint):
            ipDiffimPlot.approxKernelPlot(goodDifiList[i],
                                          (cMKernel, cEKernels),
                                          title='Kernel %d conv' % (goodDifiList[i].getID()),
                                          outfile='PCA_Kc_%d.ps' % (goodDifiList[i].getID())
                                          )

        goodDifiList = ipDiffim.getGoodFootprints(dDifiPtrList)
        nFootprint   = len(goodDifiList)
        for i in range(nFootprint):
            ipDiffimPlot.approxKernelPlot(goodDifiList[i],
                                          (dMKernel, dEKernels),
                                          title='Kernel %d dconv' % (goodDifiList[i].getID()),
                                          outfile='PCA_Kdc_%d.ps' % (goodDifiList[i].getID())
                                          )
        

    cResiduals2 = fitSpatialPca(cDifiPtrList, 'conv',   cMKernel, cEKernels, cECoeffs, policy)
    dResiduals2 = fitSpatialPca(dDifiPtrList, 'deconv', dMKernel, dEKernels, dECoeffs, policy)

    cResiduals1 = fitPerPixel(cDifiPtrList, 'conv', policy)
    dResiduals1 = fitPerPixel(dDifiPtrList, 'deconv', policy)


    return
    
    if policy.get('debugPlot') == True:
        ipDiffimPlot.plotBackground(cResiduals1, outfile='cBackground.ps')
        ipDiffimPlot.plotBackground(dResiduals1, outfile='dcBackground.ps')


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
    
    
