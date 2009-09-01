import numpy

from diffimTools import findIqr

from lsst.pex.logging import Trace


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

            # print nIter, scID, scPtr.getLabel(), scPtr.getCurrentModel().getStatus(), scPtr.getCurrentModel().getKernelSum()
            
            # Is the contained model usable?
            if scPtr.getCurrentModel().getStatus():
                kSumList.append( scPtr.getCurrentModel().getKernelSum() )
                idList.append( scID )

        nCells = len(kSumList)

        if nCells == 0:
            raise RuntimeError('No good cells found')

        # kSumArray and kSumList need to stay synchronized with spatialCells[]
        kSumArray  = numpy.array(kSumList)
        kSumMedian, kSumIqr, kSumSigMean, kSumSigMedian = findIqr(kSumArray)

        # Reject kernels with aberrent statistics
        nRejected = 0
        for idx in range(nCells):
            if numpy.fabs( (kSumArray[idx]-kSumMedian) ) > (kSumSigMedian*maxOutlierSigma):
                Trace('lsst.ip.diffim.rejectKernelSumOutliers', 4,
                      '# %s Kernel %d (kSum=%.3f) REJECTED due to bad kernel sum (median=%.3f, std=%.3f)' %
                      (spatialCells[ idList[idx] ].getLabel(),
                       spatialCells[ idList[idx] ].getCurrentId(),
                       kSumArray[idx], kSumMedian, kSumSigMedian))
                
                # Set it as fully bad, since its base model is aberrant
                spatialCells[ idList[idx] ].getCurrentModel().setStatus(False)

                # Move to the next footprint in the cell
                spatialCells[ idList[idx] ].incrementModel()
                nRejected += 1
                
        Trace('lsst.ip.diffim.rejectKernelSumOutliers', 4,
              'Kernel Sum Iteration %d, rejected %d kernels : Kernel Sum = %0.3f +/- %0.3f' %
              (nIter, nRejected, kSumMedian, kSumSigMedian))

        if nRejected == 0:
            break

    if nIter == (maxOutlierIterations-1) and nRejected != 0:
        Trace('lsst.ip.diffim.rejectKernelSumOutliers', 3,
              'Detection of kernels with deviant kSums reached its limit of %d iterations' %
              (maxOutlierIterations))

    Trace('lsst.ip.diffim.rejectKernelSumOutliers', 3,
          'Found Kernel Sum : %0.3f +/- %0.3f, %d kernels' % (kSumMedian, kSumSigMedian, nCells))


