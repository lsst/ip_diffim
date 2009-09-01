import numpy
import lsst.pex.exceptions as pexExceptions

def runPca(M, policy):
    if type(M) != numpy.ndarray:
        raise pexExceptions.LsstInvalidParameter('Invalid matrix type sent to runPca')
    
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
        if (numpy.sum(residual) > 1e-10):
            raise pexExceptions.RangeError("lsst.ip.diffim.runPca PCA invalid")
        assert

    return meanM, U, eVal, eCoeff


