import eups
import sys, os, optparse
import numpy

import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.detection as detection

from lsst.pex.logging import Log
from lsst.pex.logging import Trace
from lsst.pex.policy import Policy
import lsst.pex.exceptions as Exceptions

import lsst.ip.diffim.diffimTools as ipDiffimTools

import pdb

################
################
#####  Pixel-based spatial analysis
################
################

def spatialModelByPixel(spatialCells, policy):
    kSpatialOrder  = policy.get('kernelSpatialOrder')
    bgSpatialOrder = policy.get('backgroundSpatialOrder')
    kCols          = policy.get('kernelCols')
    kRows          = policy.get('kernelRows')

    idList   = []
    for scID, scPtr in enumerate(spatialCells):
        # Is the cell usable at all?
        if not scPtr.isUsable():
            continue
        # Is the contained model usable?
        if scPtr.getCurrentModel().getQaStatus():
            idList.append( scID )
            
    nCells = len(idList)

    # common to all spatial fits; position of constraints
    cCol = numpy.zeros(nCells)
    cRow = numpy.zeros(nCells)
    for idx in range(nCells):
        cCol[idx] = spatialCells[ idList[idx] ].getCurrentModel().getColcNorm()
        cRow[idx] = spatialCells[ idList[idx] ].getCurrentModel().getRowcNorm()

    # fit the background
    bgValues = numpy.zeros(nCells)
    bgErrors = numpy.zeros(nCells)
    for idx in range(nCells):
        bgValues[idx] = spatialCells[ idList[idx] ].getCurrentModel().getBg()
        bgErrors[idx] = spatialCells[ idList[idx] ].getCurrentModel().getBgErr()

    bgFunction = afwMath.PolynomialFunction2D(bgSpatialOrder)
    bgFit      = fitFunction(bgFunction, bgValues, bgErrors,
                             cCol, cRow, policy)
    bgFunction.setParameters(bgFit.parameterList)
    
    Trace('lsst.ip.diffim.spatialModelByPixel', 5,
          'Background fit parameters : %s' % (' '.join([('%10.3e' % (x)) for x in bgFit.parameterList])))
    
    # Fit function to each pixel
    #
    # First, make an array of Images so you don't have to do this for each pixel
    # There are no swigged arrays of Images!
    # Lists seem to work
    kImages    = []
    kErrImages = []
    for idx in range(nCells):
        kPtr    = spatialCells[ idList[idx] ].getCurrentModel().getKernelPtr()
        kErrPtr = spatialCells[ idList[idx] ].getCurrentModel().getKernelErrPtr()
        kImages.append( kPtr.computeNewImage(False)[0] )
        kErrImages.append( kErrPtr.computeNewImage(False)[0] )
    
    
    np = 0
    # Is there a vector of Functions?
    pFunctionList = []
    for kRow in range(kRows):
        for kCol in range(kCols):

            #######
            # initialize vectors, one per good kernel
            
            pValues = numpy.zeros(nCells)
            pErrors = numpy.zeros(nCells)
            for idx in range(nCells):
                pValues[idx] = kImages[idx].getVal(kCol, kRow)
                pErrors[idx] = kErrImages[idx].getVal(kCol, kRow)
                    
            # initialize vectors, one per good kernel
            #######
            # do the various fitting techniques here

            pFunction = afwMath.PolynomialFunction2D(kSpatialOrder)
            pFit      = fitFunction(pFunction, pValues, pErrors,
                                    cCol, cRow, policy)
            pFunction.setParameters(pFit.parameterList)
            
            Trace('lsst.ip.diffim.spatialModelByPixel', 5,
                  'Pixel %d fit parameters : %s' % (np, ' '.join([('%10.3e' % (x)) for x in pFit.parameterList])))

            pFunctionList.append(pFunction)
            np += 1
            
            # do the various fitting techniques here
            #######

    return bgFunction, pFunctionList


def evaluateModelByPixel_deprecated(spatialCells, bgFunction, pFunctionList, policy, reject=True):
    # This version makes a kernel image pixel-by-pixel.  Instead make
    # a spatially varying LinearCombinationKernel.
    
    kCols          = policy.get('kernelCols')
    kRows          = policy.get('kernelRows')

    nRejected = 0
    
    idList   = []
    for scID, scPtr in enumerate(spatialCells):
        # Is the cell usable at all?
        if not scPtr.isUsable():
            continue
        # Is the contained model usable?
        if scPtr.getCurrentModel().getQaStatus():
            idList.append( scID )
            
    nCells = len(idList)

    # common to all spatial fits; position of constraints
    cCol = numpy.zeros(nCells)
    cRow = numpy.zeros(nCells)
    for idx in range(nCells):
        cCol[idx] = spatialCells[ idList[idx] ].getCurrentModel().getColcNorm()
        cRow[idx] = spatialCells[ idList[idx] ].getCurrentModel().getRowcNorm()
    
    # Evaluate all the fits at the positions of the objects, create a
    # new kernel, then difference image, the calculate difference
    # image stats
    for idx in range(nCells):
        bgValue = bgFunction(cCol[idx], cRow[idx])

        # Create image representing Kernel, and then a Fixed Kernel from it
        np = 0
        kImage  = afwImage.ImageD(kCols, kRows)
        for kRow in range(kRows):
            for kCol in range(kCols):
                pValue = pFunctionList[np](cCol[idx], cRow[idx])
                kImage.set(kCol, kRow, pValue)
                np += 1
        kPtr = afwMath.KernelPtr( afwMath.FixedKernel(kImage) )

        #if policy.get('debugIO'):
        #    kImage.writeFits('skImage1_%d.fits' % (idx))
        #    spatialCells[ idList[idx] ].getCurrentModel().getKernelPtr().computeNewImage(False)[0].writeFits('kImage_%d.fits' % (idx))
        
        # Create difference image using Kernel model
        diffIm = ipDiffim.convolveAndSubtract(spatialCells[ idList[idx] ].getCurrentModel().getMiToConvolvePtr().get(),
                                              spatialCells[ idList[idx] ].getCurrentModel().getMiToNotConvolvePtr().get(),
                                              kPtr, bgValue)

        # Find quality of difference image
        diffImStats = ipDiffim.DifferenceImageStatisticsF(diffIm)
        if (diffImStats.evaluateQuality(policy) == False):
            if (reject == True):
                # This is only bad in the context of the spatial model
                # May be needed in the future
                #
                # spatialCells[ idList[idx] ].getCurrentModel().setQaStatus(False)                
                spatialCells[ idList[idx] ].increment()
                nRejected += 1

                label = 'Rejected:'
            else:
                label = 'Poor'
        else:
            label = 'OK'

        Trace('lsst.ip.diffim.evaluateModelByPixel', 5,
              '%s Kernel %d : %s Spatial residuals = %.2f +/- %.2f sigma' %
              (spatialCells[ idList[idx] ].getLabel(),
               spatialCells[ idList[idx] ].getCurrentModel().getID(),
               label, diffImStats.getResidualMean(), diffImStats.getResidualStd()))

#        if policy.get('debugPlot') == True:
#            ipDiffimDebug.plotDiffImQuality1(goodDifiList[i],
#                                             diffIm,
#                                             kernelPtr,
#                                             label='Spatial %s kernel %d' % (label, goodDifiList[i].getID()),
#                                             outfile='SKernel_%s%d.ps' % (label, goodDifiList[i].getID())
#                                             )

    return nRejected

def evaluateModelByPixel(spatialCells, bgFunction, sKernel, policy, reject=True):
    kCols          = policy.get('kernelCols')
    kRows          = policy.get('kernelRows')

    nRejected = 0
    
    idList   = []
    for scID, scPtr in enumerate(spatialCells):
        # Is the cell usable at all?
        if not scPtr.isUsable():
            continue
        # Is the contained model usable?
        if scPtr.getCurrentModel().getQaStatus():
            idList.append( scID )
            
    nCells = len(idList)

    # common to all spatial fits; position of constraints
    cCol = numpy.zeros(nCells)
    cRow = numpy.zeros(nCells)
    for idx in range(nCells):
        cCol[idx] = spatialCells[ idList[idx] ].getCurrentModel().getColcNorm()
        cRow[idx] = spatialCells[ idList[idx] ].getCurrentModel().getRowcNorm()
    
    # Evaluate all the fits at the positions of the objects, create a
    # new kernel, then difference image, the calculate difference
    # image stats
    for idx in range(nCells):
        bgValue = bgFunction(cCol[idx], cRow[idx])
        sImage  = sKernel.computeNewImage(False, cCol[idx], cRow[idx])[0]

        #if policy.get('debugIO'):
        #    sImage.writeFits('skImage_%d.fits' % (idx))
            
        kPtr    = afwMath.KernelPtr( afwMath.FixedKernel(sImage) )
        # Create difference image using Kernel model
        diffIm  = ipDiffim.convolveAndSubtract(spatialCells[ idList[idx] ].getCurrentModel().getMiToConvolvePtr().get(),
                                               spatialCells[ idList[idx] ].getCurrentModel().getMiToNotConvolvePtr().get(),
                                               kPtr, bgValue)

        # Find quality of difference image
        diffImStats = ipDiffim.DifferenceImageStatisticsF(diffIm)
        if (diffImStats.evaluateQuality(policy) == False):
            if (reject == True):
                # This is only bad in the context of the spatial model
                # May be needed in the future
                #
                # spatialCells[ idList[idx] ].getCurrentModel().setQaStatus(False)                
                spatialCells[ idList[idx] ].increment()
                nRejected += 1

                label = 'Rejected:'
            else:
                label = 'Poor'
        else:
            label = 'OK'

        Trace('lsst.ip.diffim.evaluateModelByPixel2', 5,
              '%s Kernel %d : %s Spatial residuals = %.2f +/- %.2f sigma' %
              (spatialCells[ idList[idx] ].getLabel(),
               spatialCells[ idList[idx] ].getCurrentModel().getID(),
               label, diffImStats.getResidualMean(), diffImStats.getResidualStd()))

#        if policy.get('debugPlot') == True:
#            ipDiffimDebug.plotDiffImQuality1(goodDifiList[i],
#                                             diffIm,
#                                             kernelPtr,
#                                             label='Spatial %s kernel %d' % (label, goodDifiList[i].getID()),
#                                             outfile='SKernel_%s%d.ps' % (label, goodDifiList[i].getID())
#                                             )

    return nRejected

################
################
#####  Pca-based spatial analysis
################
################


def spatialModelByPca(spatialCells, eCoeffs, neVal, policy):
    kSpatialOrder  = policy.get('kernelSpatialOrder')
    bgSpatialOrder = policy.get('backgroundSpatialOrder')
    kCols          = policy.get('kernelCols')
    kRows          = policy.get('kernelRows')

    idList   = []
    for scID, scPtr in enumerate(spatialCells):
        # Is the cell usable at all?
        if not scPtr.isUsable():
            continue
        # Is the contained model usable?
        if scPtr.getCurrentModel().getQaStatus():
            idList.append( scID )
            
    nCells = len(idList)
    nCoeff = neVal

    # common to all spatial fits; position of constraints
    cCol = numpy.zeros(nCells)
    cRow = numpy.zeros(nCells)
    for idx in range(nCells):
        cCol[idx] = spatialCells[ idList[idx] ].getCurrentModel().getColcNorm()
        cRow[idx] = spatialCells[ idList[idx] ].getCurrentModel().getRowcNorm()

    # fit the background
    bgValues = numpy.zeros(nCells)
    bgErrors = numpy.zeros(nCells)
    for idx in range(nCells):
        bgValues[idx] = spatialCells[ idList[idx] ].getCurrentModel().getBg()
        bgErrors[idx] = spatialCells[ idList[idx] ].getCurrentModel().getBgErr()

    bgFunction = afwMath.PolynomialFunction2D(bgSpatialOrder)
    bgFit      = fitFunction(bgFunction, bgValues, bgErrors,
                             cCol, cRow, policy)
    bgFunction.setParameters(bgFit.parameterList)
    
    Trace('lsst.ip.diffim.spatialModelByPca', 5,
          'Background fit parameters : %s' % (' '.join([('%10.3e' % (x)) for x in bgFit.parameterList])))

    # Fit spatial variation of each eigenkernel
    eFunctionList = []
    for nc in range(nCoeff):
        # The contribution of eigenKernel X to Kernel Y is in eCoeff[Y,X].
        coeffs = eCoeffs[:,nc]

        # The uncertainties on these coefficients are not known
        # This is an approximation 
        errors = numpy.sqrt( numpy.abs(coeffs) )
        
        eFunction = afwMath.PolynomialFunction2D(kSpatialOrder)
        eFit = fitFunction(eFunction,
                           coeffs,
                           errors,
                           cCol,
                           cRow,
                           policy)
        eFunction.setParameters(eFit.parameterList)
        
        Trace('lsst.ip.diffim.spatialModelByPca', 5,
              'eigenKernel %d fit parameters : %s' % (nc, ' '.join([('%10.3e' % (x)) for x in eFit.parameterList])))
        eFunctionList.append(eFunction)

    return bgFunction, eFunctionList


def evaluateModelByPca(spatialCells, bgFunction, eKernel, policy, reject=True):
    nRejected = 0

    idList   = []
    for scID, scPtr in enumerate(spatialCells):
        # Is the cell usable at all?
        if not scPtr.isUsable():
            continue
        # Is the contained model usable?
        if scPtr.getCurrentModel().getQaStatus():
            idList.append( scID )
            
    nCells = len(idList)

    # common to all spatial fits; position of constraints
    cCol = numpy.zeros(nCells)
    cRow = numpy.zeros(nCells)
    for idx in range(nCells):
        cCol[idx] = spatialCells[ idList[idx] ].getCurrentModel().getColcNorm()
        cRow[idx] = spatialCells[ idList[idx] ].getCurrentModel().getRowcNorm()

    # Evaluate all the fits at the positions of the objects, create a
    # new kernel, then difference image, the calculate difference
    # image stats
    for idx in range(nCells):
        bgValue = bgFunction(cCol[idx], cRow[idx])
        sImage  = eKernel.computeNewImage(False, cCol[idx], cRow[idx])[0]
        kPtr    = afwMath.KernelPtr( afwMath.FixedKernel(sImage) )
        diffIm  = ipDiffim.convolveAndSubtract(spatialCells[ idList[idx] ].getCurrentModel().getMiToConvolvePtr().get(),
                                               spatialCells[ idList[idx] ].getCurrentModel().getMiToNotConvolvePtr().get(),
                                               kPtr, bgValue)

        # Find quality of difference image
        diffImStats = ipDiffim.DifferenceImageStatisticsF(diffIm)
        if (diffImStats.evaluateQuality(policy) == False):
            if (reject == True):
                # This is only bad in the context of the spatial model
                # May be needed in the future
                #
                # spatialCells[ idList[idx] ].getCurrentModel().setQaStatus(False)                
                spatialCells[ idList[idx] ].increment()
                nRejected += 1

                label = 'Rejected:'
            else:
                label = 'Poor'
        else:
            label = 'OK'

        Trace('lsst.ip.diffim.evaluateModelByPca', 5,
              '%s Kernel %d : %s Pca residuals = %.2f +/- %.2f sigma' %
              (spatialCells[ idList[idx] ].getLabel(),
               spatialCells[ idList[idx] ].getCurrentModel().getID(),
               label, diffImStats.getResidualMean(), diffImStats.getResidualStd()))
        
    return nRejected

################
################
#####  Pca
################
################

def modelPca(spatialCells, policy):
    kCols = policy.get('kernelCols')
    kRows = policy.get('kernelRows')
    
    idList   = []
    for scID, scPtr in enumerate(spatialCells):
        # Is the cell usable at all?
        if not scPtr.isUsable():
            continue
        # Is the contained model usable?
        if scPtr.getCurrentModel().getQaStatus():
            idList.append( scID )
            
    nCells = len(idList)
    
    # matrix to invert
    M = numpy.zeros((kCols*kRows, nCells))
    for idx in range(nCells):
        # Ideally, we want a "ravel" method on each Model to use both PSF/Kernels here
        kernelImage  = spatialCells[ idList[idx] ].getCurrentModel().getKernelPtr().computeNewImage(False)[0]
        M[:,idx]     = ipDiffimTools.imageToVector(kernelImage)
        
    # do the PCA
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

    Trace('lsst.ip.diffim.modelPca', 5,
          'EigenValues : %s' % (' '.join([('%10.3e' % (x)) for x in eVal])))

    # Find the contribution of each eigenKernel to each Kernel.
    # Simple dot product, transpose of M dot the eigenCoeff matrix.
    # The contribution of eigenKernel X to Kernel Y is in eCoeff[Y,X].
    #
    # I.e. M[:,X] = numpy.sum(U * eCoeff[X], 1)
    eCoeff = numpy.dot(M.T, U)
    for i in range(nCells):
        residual = numpy.sum(U * eCoeff[i], 1) - M[:,i]
        assert(numpy.sum(residual) < 1e-10)

    # Again, I specialize this here for Kernels; need to generalize
    # more if we want to use PSFs
    #
    # Turn into Kernels
    mKernelPtr       = ipDiffimTools.vectorToKernelPtr( meanM, kCols, kRows )
    if policy.get('debugIO'):
        ipDiffimTools.vectorToImage(meanM, kCols, kRows).writeFits('mKernel.fits')

    eKernelPtrVector = afwMath.VectorKernel()
    for i in range(U.shape[1]):
        eKernelPtr   = ipDiffimTools.vectorToKernelPtr( U[:,i], kCols, kRows )
        eKernelPtrVector.append(eKernelPtr)

        if policy.get('debugIO'):
            ipDiffimTools.vectorToImage(U[:,i], kCols, kRows).writeFits('eKernel_%d.fits' % (i))
            

    return mKernelPtr, eKernelPtrVector, eVal, eCoeff

################
################
#####  Pca
################
################

def fitFunction(function, values, errors, cols, rows, policy):
    nSigmaSq = policy.get('nSigmaSq')
    stepsize = policy.get('stepsize')

    # initialize fit parameters
    nPar = function.getNParameters()
    pars = numpy.zeros(nPar)       # start with no spatial variation
    pars[0]  = numpy.mean(values)  # except for the constant term
    stepsize = stepsize * numpy.ones(nPar)
    fit      = afwMath.minimize(function,
                                pars,
                                stepsize,
                                values,
                                errors,
                                cols,
                                rows,
                                nSigmaSq)
    if not fit.isValid:
        # throw exception
        pass

    return fit

################
################
#####  Building the initial set of Kernels
################
################


def rejectKernelOutliers(spatialCells, policy):
    # Compare kernel sums; things that are bad (e.g. variable stars, moving objects, cosmic rays)
    # will have a kernel sum far from the mean; reject these.
    
    maxOutlierIterations = policy.get('maxOutlierIterations')
    maxOutlierSigma      = policy.get('maxOutlierSigma')

    # So are we using pexLog.Trace or pexLog.Log
    Trace('lsst.ip.diffim.rejectKernelOutliers', 4,
          'Rejecting kernels with deviant kSums');
    
    for nIter in xrange(maxOutlierIterations):
        kSumList = []
        idList   = []
        for scID, scPtr in enumerate(spatialCells):
            # Is the cell usable at all?
            if not scPtr.isUsable():
                continue
            # Is the contained model usable?
            if scPtr.getCurrentModel().getQaStatus():
                kSumList.append( scPtr.getCurrentModel().getKernelSum() )
                idList.append( scID )

        nCells = len(kSumList)

        if nCells == 0:
            raise Exceptions.LsstOutOfRange('No good cells found')

        kSumArray = numpy.array(kSumList)
        kSumMean  = kSumArray.mean()
        kSumStd   = kSumArray.std()

        # Reject kernels with aberrent statistics
        nRejected = 0
        for idx in range(nCells):
            if numpy.fabs( (kSumArray[idx]-kSumMean)/kSumStd ) > maxOutlierSigma:
                Trace('lsst.ip.diffim.rejectKernelOutliers', 5,
                      '# %s Kernel %d (kSum=%.3f) REJECTED due to bad kernel sum (mean=%.3f, std=%.3f)' %
                      (spatialCells[ idList[idx] ].getLabel(),
                       spatialCells[ idList[idx] ].getCurrentModel().getID(),
                       kSumArray[idx], kSumMean, kSumStd))
                
                # Move to the next footprint in the cell
                # Here we actually set it as bad, since its original model is aberrant
                spatialCells[ idList[idx] ].getCurrentModel().setQaStatus(False)
                spatialCells[ idList[idx] ].increment()
                nRejected += 1
                
        Trace('lsst.ip.diffim.rejectKernelOutliers', 4,
              'Kernel Sum Iteration %d, rejected %d kernels : Kernel Sum = %0.3f +/- %0.3f' %
              (nIter, nRejected, kSumMean, kSumStd))

        if nRejected == 0:
            break

    if nIter == (maxOutlierIterations-1) and nRejected != 0:
        Trace('lsst.ip.diffim.rejectKernelOutliers', 3,
              'Detection of kernels with deviant kSums reached its limit of %d iterations' %
              (maxOutlierIterations))

    Trace('lsst.ip.diffim.rejectKernelOutliers', 3,
          'Kernel Sum : %0.3f +/- %0.3f, %d kernels' % (kSumMean, kSumStd, nCells))

def segmentMaskedImage(fpList,
                       templateMaskedImage,
                       scienceMaskedImage,
                       kBasisList,
                       policy):

    
    templateMaskedImagePtr = afwImage.MaskedImageFPtr(templateMaskedImage)
    scienceMaskedImagePtr  = afwImage.MaskedImageFPtr(scienceMaskedImage)
    
    nSegmentCol = policy.get('nSegmentCol')
    nSegmentRow = policy.get('nSegmentRow')

    nSegmentColPix = int( templateMaskedImage.getCols() / nSegmentCol )
    nSegmentRowPix = int( templateMaskedImage.getRows() / nSegmentRow )

    spatialCellPtrs = ipDiffim.SpatialModelCellPtrListFK()
    cellCount = 0
    for col in range(nSegmentCol):
        colMin    = max(0, col*nSegmentColPix)
        colMax    = min(templateMaskedImage.getCols(), (col+1)*nSegmentColPix)
        colCenter = int( 0.5 * (colMin + colMax) )
        
        for row in range(nSegmentRow):
            rowMin    = max(0, row*nSegmentRowPix)
            rowMax    = min(templateMaskedImage.getRows(), (row+1)*nSegmentRowPix)
            rowCenter = int( 0.5 * (rowMin + rowMax) )

            fpPtrList         = detection.FootprintContainerT()
            modelPtrList      = ipDiffim.KernelModelQaPtrListF()

            # This is a bit blunt and could be more clever
            # Should never really have a loop within a loop within a loop
            # But we will not have that many Footprints...
            for fpID, fpPtr in enumerate(fpList):
                fpBBox = fpPtr.getBBox()
                fpMin  = fpBBox.min()
                fpMax  = fpBBox.max()
                fpColC = 0.5 * (fpMin.x() + fpMax.x())
                fpRowC = 0.5 * (fpMin.y() + fpMax.y())

                if (fpColC >= colMin) and (fpColC < colMax) and (fpRowC >= rowMin) and (fpRowC < rowMax):
                    fpPtrList.push_back(fpPtr)
                    modelPtrList.push_back( ipDiffim.KernelModelQaPtrF (
                        ipDiffim.KernelModelQaF(fpPtr, templateMaskedImagePtr, scienceMaskedImagePtr, kBasisList, policy, False)
                        ))

            label = 'Cell %d' % cellCount
            spatialCellPtr    = ipDiffim.SpatialModelCellPtrFK(
                ipDiffim.SpatialModelCellFK(label, colCenter, rowCenter, fpPtrList, modelPtrList)
                )
            spatialCellPtrs.push_back(spatialCellPtr)

            # Formatting to the screen 
            Trace('lsst.ip.diffim.segmentMaskedImage', 2, '')

            cellCount += 1

    return spatialCellPtrs

################
################
#####  Main fitting loop
################
################

def spatialTesting(spatialCells, kBasisList, policy, scID):
    kCols = policy.get('kernelCols')
    kRows = policy.get('kernelRows')
    
    try:
        rejectKernelOutliers(spatialCells, policy)
    except:
        Trace('lsst.ip.diffim', 2,
              'LOOP %d FAILED; no good kernels' % (scID))
        return
    
    # LOOP 2 : spatial order
    for order in range(3):
        policy.set('kernelSpatialOrder', order)
        policy.set('backgroundSpatialOrder', order)
        
        kSpatialOrder  = policy.get('kernelSpatialOrder')
        bgSpatialOrder = policy.get('backgroundSpatialOrder')
        
        Trace('lsst.ip.diffim', 1,
              'LOOP %d %d' % (scID, order))
        
        #############
        # Pixel-by-pixel fit
        maxSpatialIterations = policy.get('maxSpatialIterations')
        nRejected = 1
        nIter     = 0
        
        # LOOP 3 : spatial pixel sigma clipping
        while (nRejected != 0) and (nIter < maxSpatialIterations):
        
            # do spatial fit here pixel by pixel
            bgFunction, pFunctionList = spatialModelByPixel(spatialCells, policy)

            # ideally...
            # sKernelPtr = afwMath.LinearCombinationKernel(kBasisList, pFunctionList)
            #
            # instead
            # try to build a spatial kernel with a function per pixel
            sKernelFunc  = afwMath.PolynomialFunction2D(kSpatialOrder)
            kParams      = numpy.zeros( (kCols*kRows, sKernelFunc.getNParameters()) )
            for p in range(kCols*kRows):
                kParams[p] = pFunctionList[p].getParameters()
            # Create spatially varying kernel pointer
            sKernelPtr = afwMath.LinearCombinationKernel(kBasisList, sKernelFunc)
            sKernelPtr.setSpatialParameters(kParams)
            
            nRejected  = evaluateModelByPixel(spatialCells,
                                              bgFunction, sKernelPtr, 
                                              policy, reject=False)
            nIter += 1

        #############
        # PCA fit
        nRejected = 1
        nIter     = 0
        
        # LOOP 4a : PCA sigma clipping
        while (nRejected != 0) and (nIter < maxSpatialIterations):
            
            # Run the PCA
            mKernelPtr, eKernelPtrVector, eVal, eCoeff = modelPca(spatialCells, policy)
            
            # Here we make a decision on how many eigenComponents to use based
            # on eVal, etc
            #
            # While we are testing, check them all
            
            # Find spatial variation of only those components
            # Remove this line after being done testing
            # We fit them all first
            bgFunction, eFunctionList = spatialModelByPca(spatialCells, eCoeff, len(eVal), policy)
    
            # LOOP 4b : Number of principal components
            for neVal in range( len(eVal) ):
            
                Trace('lsst.ip.diffim', 3,
                      'Using varying eigenkernels : N = %d' % (neVal))
                
                # Find spatial variation of only those components
                # Comment this back in when we are not testing
                #
                # bgFunction, eFunctionList = spatialModelByPca(spatialCells, eCoeff, neVal, policy)
            
                # Build LinearCombinationKernel for only neVal components
                # Start with mean Kernel image
                eKernelBases = afwMath.KernelListD()
                eKernelBases.push_back(mKernelPtr)
                # Append eigenKernel images
                for ek in range(neVal):
                    eKernelBases.push_back(eKernelPtrVector[ek])
                    
                # Mean kernel has no spatial variation
                eKernelFunc   = afwMath.PolynomialFunction2D(kSpatialOrder)
                kParams       = numpy.zeros( (neVal+1, eKernelFunc.getNParameters()) )
                kParams[0][0] = 1.0
                # Add already-fit-for spatial variation of eigenKernels
                for ek in range(neVal):
                    kParams[ek+1] = eFunctionList[ek].getParameters()
    
                # Create spatially varying eigenKernel pointer
                eKernelPtr = afwMath.LinearCombinationKernel(eKernelBases, eKernelFunc)
                eKernelPtr.setSpatialParameters(kParams)
        
                # Evaluate quality of spatial fit
                nRejected = evaluateModelByPca(spatialCells, bgFunction, eKernelPtr,
                                               policy, reject=False)
            
            nIter += 1

  
################
################
#####  o  ######
################
################

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
    policy.set('minimumCleanFootprints', 5)
    policy.set('footprintDetectionThreshold', 250.)

    # if you are convolving the template
    #policy.set('iterateKernel', False)
    # if you are convolving the image
    # policy.set('iterateKernel', True)
    
    fpList = ipDiffim.getCollectionOfFootprintsForPsfMatching(templateMaskedImage,
                                                              scienceMaskedImage,
                                                              policy)

    # LOOP 1 : convolution vs deconvolution
    Trace('lsst.ip.diffim', 1, 'SC List 1')
    spatialCellsC = segmentMaskedImage(fpList, templateMaskedImage, scienceMaskedImage, kBasisList, policy)
    Trace('lsst.ip.diffim', 1, 'SC List 2')
    spatialCellsD = segmentMaskedImage(fpList, scienceMaskedImage, templateMaskedImage, kBasisList, policy)

    spatialTesting(spatialCellsC, kBasisList, policy, 0)
    spatialTesting(spatialCellsD, kBasisList, policy, 1)

    return




    scList = (spatialCellsC, spatialCellsD)
    for scID in range(len(scList)):
        spatialCells = scList[scID]
        
        try:
            rejectKernelOutliers(spatialCells, policy)
        except:
            continue

        # LOOP 2 : spatial order
        for order in range(3):
            policy.set('kernelSpatialOrder', order)
            policy.set('backgroundSpatialOrder', order)

            kSpatialOrder  = policy.get('kernelSpatialOrder')
            bgSpatialOrder = policy.get('backgroundSpatialOrder')
            
            Trace('lsst.ip.diffim', 1,
                  'LOOP %d %d' % (scID, order))

            #############
            # Pixel-by-pixel fit
            maxSpatialIterations = policy.get('maxSpatialIterations')
            nRejected = 1
            nIter     = 0

            # LOOP 3 : spatial pixel sigma clipping
            while (nRejected != 0) and (nIter < maxSpatialIterations):
            
                # do spatial fit here pixel by pixel
                bgFunction, pFunctionList = spatialModelByPixel(spatialCells, policy)

                # ideally...
                # sKernelPtr = afwMath.LinearCombinationKernel(kBasisList, pFunctionList)
                #
                # instead
                # try to build a spatial kernel with a function per pixel
                sKernelFunc  = afwMath.PolynomialFunction2D(kSpatialOrder)
                kParams      = numpy.zeros( (kCols*kRows, sKernelFunc.getNParameters()) )
                for p in range(kCols*kRows):
                    kParams[p] = pFunctionList[p].getParameters()
                # Create spatially varying kernel pointer
                sKernelPtr = afwMath.LinearCombinationKernel(kBasisList, sKernelFunc)
                sKernelPtr.setSpatialParameters(kParams)
                
                nRejected  = evaluateModelByPixel(spatialCells,
                                                  bgFunction, sKernelPtr, 
                                                  policy, reject=False)
                nIter += 1
    
            #############
            # PCA fit
            nRejected = 1
            nIter     = 0

            # LOOP 4a : PCA sigma clipping
            while (nRejected != 0) and (nIter < maxSpatialIterations):

                # Run the PCA
                mKernelPtr, eKernelPtrVector, eVal, eCoeff = modelPca(spatialCells, policy)
                
                # Here we make a decision on how many eigenComponents to use based
                # on eVal, etc
                #
                # While we are testing, check them all
                
                # Find spatial variation of only those components
                # Remove this line after being done testing
                # We fit them all first
                bgFunction, eFunctionList = spatialModelByPca(spatialCells, eCoeff, len(eVal), policy)
        
                # LOOP 4b : Number of principal components
                for neVal in range( len(eVal) ):
                
                    Trace('lsst.ip.diffim', 3,
                          'Using varying eigenkernels : N = %d' % (neVal))
                    
                    # Find spatial variation of only those components
                    # Comment this back in when we are not testing
                    #
                    # bgFunction, eFunctionList = spatialModelByPca(spatialCells, eCoeff, neVal, policy)
                
                    # Build LinearCombinationKernel for only neVal components
                    # Start with mean Kernel image
                    eKernelBases = afwMath.KernelListD()
                    eKernelBases.push_back(mKernelPtr)
                    # Append eigenKernel images
                    for ek in range(neVal):
                        eKernelBases.push_back(eKernelPtrVector[ek])

                    # Mean kernel has no spatial variation
                    eKernelFunc   = afwMath.PolynomialFunction2D(kSpatialOrder)
                    kParams       = numpy.zeros( (neVal+1, eKernelFunc.getNParameters()) )
                    kParams[0][0] = 1.0
                    # Add already-fit-for spatial variation of eigenKernels
                    for ek in range(neVal):
                        kParams[ek+1] = eFunctionList[ek].getParameters()
        
                    # Create spatially varying eigenKernel pointer
                    eKernelPtr = afwMath.LinearCombinationKernel(eKernelBases, eKernelFunc)
                    eKernelPtr.setSpatialParameters(kParams)
            
                    # Evaluate quality of spatial fit
                    nRejected = evaluateModelByPca(spatialCells, bgFunction, eKernelPtr,
                                                   policy, reject=False)
                
                nIter += 1
    
def run():
    Log.getDefaultLog()
    memId0 = dafBase.Citizen_getNextMemId()
    main()
    # check for memory leaks
    if dafBase.Citizen_census(0, memId0) != 0:
        print dafBase.Citizen_census(0, memId0), 'Objects leaked:'
        print dafBase.Citizen_census(dafBase.cout, memId0)

if __name__ == '__main__':
    run()
