import numpy
import lsst.afw.detection as detection
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.pex.exceptions as Exceptions
from   lsst.pex.logging import Trace

# for DM <-> numpy conversions
import lsst.afw.image.testUtils as imTestUtils

# relative imports, since these are in __init__.py
import diffimLib 
import diffimDebug

import pdb

#######
# Helper functions used by all ip_diffim packages
#######

def fitFunction(function, values, errors, cols, rows, policy):
    # Fit a function to data
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
        raise Exceptions.LsstRuntime('Spatial fit fails')

    return fit

def rejectKernelSumOutliers(spatialCells, policy):
    # Compare kernel sums; things that are bad (e.g. variable stars,
    # moving objects, cosmic rays) will have a kernel sum far from the
    # mean; reject these.

    maxOutlierIterations = policy.get('maxOutlierIterations')
    maxOutlierSigma      = policy.get('maxOutlierSigma')

    # So are we using pexLog.Trace or pexLog.Log?
    Trace('lsst.ip.diffim.rejectKernelSumOutliers', 4,
          'Rejecting kernels with deviant kSums');
    
    for nIter in xrange(maxOutlierIterations):
        kSumList = []
        idList   = []
        for scID, scPtr in enumerate(spatialCells):
            # Is the cell usable at all?
            if not scPtr.isUsable():
                continue
            # Is the contained model usable?
            if scPtr.getCurrentModel().getSdqaStatus():
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
                Trace('lsst.ip.diffim.rejectKernelSumOutliers', 5,
                      '# %s Kernel %d (kSum=%.3f) REJECTED due to bad kernel sum (mean=%.3f, std=%.3f)' %
                      (spatialCells[ idList[idx] ].getLabel(),
                       spatialCells[ idList[idx] ].getCurrentModel().getID(),
                       kSumArray[idx], kSumMean, kSumStd))
                
                # Set it as fully bad, since its base model is aberrant
                spatialCells[ idList[idx] ].getCurrentModel().setSdqaStatus(False)

                # Move to the next footprint in the cell
                spatialCells[ idList[idx] ].incrementModel()
                nRejected += 1
                
        Trace('lsst.ip.diffim.rejectKernelSumOutliers', 4,
              'Kernel Sum Iteration %d, rejected %d kernels : Kernel Sum = %0.3f +/- %0.3f' %
              (nIter, nRejected, kSumMean, kSumStd))

        if nRejected == 0:
            break

    if nIter == (maxOutlierIterations-1) and nRejected != 0:
        Trace('lsst.ip.diffim.rejectKernelSumOutliers', 3,
              'Detection of kernels with deviant kSums reached its limit of %d iterations' %
              (maxOutlierIterations))

    Trace('lsst.ip.diffim.rejectKernelSumOutliers', 3,
          'Kernel Sum : %0.3f +/- %0.3f, %d kernels' % (kSumMean, kSumStd, nCells))


def createSpatialModelKernelCells(templateMaskedImage,
                                  scienceMaskedImage,
                                  fpInList,
                                  kFunctor,
                                  policy,
                                  cFlag='c'):

    nSegmentCol = policy.get('nSegmentCol')
    nSegmentRow = policy.get('nSegmentRow')

    nSegmentColPix = int( templateMaskedImage.getWidth() / nSegmentCol )
    nSegmentRowPix = int( templateMaskedImage.getHeight() / nSegmentRow )

    spatialCells   = diffimLib.VectorSpatialModelCellF()
    
    cellCount = 0
    for col in range(nSegmentCol):
        colMin    = max(0, col*nSegmentColPix)
        colMax    = min(templateMaskedImage.getWidth(), (col+1)*nSegmentColPix)
        colCenter = int( 0.5 * (colMin + colMax) )
        
        for row in range(nSegmentRow):
            rowMin     = max(0, row*nSegmentRowPix)
            rowMax     = min(templateMaskedImage.getHeight(), (row+1)*nSegmentRowPix)
            rowCenter  = int( 0.5 * (rowMin + rowMax) )
            label      = 'c%d' % cellCount

            modelList = diffimLib.VectorSpatialModelKernelF()

            # This is a bit dumb and could be more clever
            # Should never really have a loop within a loop within a loop
            # But we will not have *that many* Footprints...
            for fpID, fpPtr in enumerate(fpInList):
                
                fpBBox = fpPtr.getBBox()
                fpColC = 0.5 * (fpBBox.getX0() + fpBBox.getX1())
                fpRowC = 0.5 * (fpBBox.getY0() + fpBBox.getY1())

                if (fpColC >= colMin) and (fpColC < colMax) and (fpRowC >= rowMin) and (fpRowC < rowMax):

                    tSubImage = afwImage.MaskedImageF(templateMaskedImage, fpBBox)
                    iSubimage = afwImage.MaskedImageF(scienceMaskedImage, fpBBox)
                    model = diffimLib.SpatialModelKernelF(fpPtr,
                                                          tSubImage,
                                                          iSubimage,
                                                          kFunctor,
                                                          policy,
                                                          False)
                    if policy.get('debugIO'):
                        diffimDebug.writeDiffImages(cFlag, '%s_%d' % (label, fpID), model)
                        
                    modelList.push_back( model )

            spatialCell = diffimLib.SpatialModelCellF(label, colCenter, rowCenter, modelList)
            spatialCells.push_back(spatialCell)

            # Formatting to the screen 
            Trace('lsst.ip.diffim.createSpatialModelKernelCells', 2, '')

            cellCount += 1

    return spatialCells


def runPca(M, policy):
    if type(M) != numpy.ndarray:
        raise Exceptions.LsstInvalidParameter('Invalid matrix type sent to runPca')
    
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

    Trace('lsst.ip.diffim.runPca', 5,
          'EigenValues : %s' % (' '.join([('%10.3e' % (x)) for x in eVal])))
    eSum  = numpy.cumsum(eVal)
    eSum /= eSum[-1]
    Trace('lsst.ip.diffim.runPca', 5,
          'EigenFraction : %s' % (' '.join([('%10.3e' % (x)) for x in eSum])))

    # Find the contribution of each eigenKernel to each Kernel.
    # Simple dot product, transpose of M dot the eigenCoeff matrix.
    # The contribution of eigenKernel X to Kernel Y is in eCoeff[Y,X].
    #
    # I.e. M[:,X] = numpy.sum(U * eCoeff[X], 1)
    eCoeff = numpy.dot(M.T, U)

    # This should probably be moved to a unit test
    for i in range(M.shape[1]):
        residual = numpy.sum(U * eCoeff[i], 1) - M[:,i]
        assert(numpy.sum(residual) < 1e-10)

    return meanM, U, eVal, eCoeff


#######
# Expansions of functionality found in imTestUtils
#######

def vectorFromImage(im, dtype=float):
    vec = numpy.zeros(im.getWidth()*im.getHeight(), dtype=dtype)
    idx = 0
    for row in range(im.getHeight()):
        for col in range(im.getWidth()):
            vec[idx] = im.get(col, row)
            idx     += 1
    return vec

def imageFromVector(vec, width, height, retType=afwImage.ImageF):
    im  = retType(width, height)
    idx = 0
    for row in range(height):
        for col in range(width):
            # need to cast numpy.float64 as float
            im.set(col, row, float(vec[idx]))
            idx     += 1
    return im

