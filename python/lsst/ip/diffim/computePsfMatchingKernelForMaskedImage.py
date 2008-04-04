#!/usr/bin/env python
import numpy

import lsst.pex.policy
import lsst.pex.logging as logging
import lsst.pex.exceptions as pex_ex
from lsst.pex.logging import Log

import lsst.afw.Core.afwLib as afw
import lsst.detection.detectionLib as detection
import ip_diffimLib as ip_diffim # relative import, since this is in __init__.py

__all__ = ['computePsfMatchingKernelForMaskedImage']

GetMaskPlanesFromPolicy = False

def computePsfMatchingKernelForMaskedImage(kernelFunctionPtr, backgroundFunctionPtr,
    imageToConvolve, imageToNotConvolve, kernelBasisList, footprintList, policy):
    """Compute spatially varying PSF matching kernel for image subtraction.

    Implements the main use case of Image Subtraction.
    Solves for kernel K and background B in the model:
 
        Ic.conv.K = Inc + B
 
    The difference image is (Inc + B - Ic.conv.K)
    
    In/Out:
    - kernelFunctionPtr: shared pointer to function for spatial variation of kernel;
        the parameters are computed by this function
    - backgroundFunctionPtr: shared pointer to function for spatial variation of background;
        the parameters are computed by this function
    
    Inputs:
    - imageToConvolve: the MaskedImage to be convolved (Inc in the equation above)
    - imageToNotConvolve: the MaskedImage to not be convolved (Ic in the equation above)
    - kernelBasisList: list of basis kernels
    - footprintList: list of footprints on which to do difference imaging
    - policy: items used include:
        - computePsfMatchingKernelForMaskedImage
            - maximumFootprintResidualMean
            - maximumFootprintResidualVariance
            - maximumKernelSumNSigma
            - maximumIterationsKernelSum
        - kernelCols
        - kernelRows
        - debugIO (optional; if True various intermediate stage files are written out)
    
    Returns:
    - kernelPtr, a pointer to the kernel K in the equation above
    """
    maximumFootprintResidualMean = policy.get('computePsfMatchingKernelForMaskedImage').get('maximumFootprintResidualMean')
    maximumFootprintResidualVariance = policy.get('computePsfMatchingKernelForMaskedImage').get('maximumFootprintResidualVariance')
    maximumKernelSumNSigma = policy.get('computePsfMatchingKernelForMaskedImage').get('maximumKernelSumNSigma')
    maximumIterationsKernelSum = policy.get('computePsfMatchingKernelForMaskedImage').get('maximumIterationsKernelSum')
    kernelCols = policy.get('kernelCols')
    kernelRows = policy.get('kernelRows')
    debugIO = policy.get('debugIO', False)

    diffImLog = Log(Log.getDefaultLog(), "ip.diffim.computePsfMatchingKernelForMaskedImage")
    
    if GetMaskPlanesFromPolicy:
        badMaskBit = policy.get('badMaskBit')
        edgeMaskBit = policy.get('edgeMaskBit')
    else:
        badMaskBit = imageToConvolve.getMask().getMaskPlane('BAD')
        edgeMaskBit = imageToConvolve.getMask().getMaskPlane('EDGE')
        #assert(badMaskBit == imageToNotConvolve.getMask().getMaskPlane('BAD'))
        #assert(edgeMaskBit == imageToNotConvolve.getMask().getMaskPlane('EDGE'))
    badPixelMask = 0
    if badMaskBit > 0:
        badPixelMask |= 1 << badMaskBit
    if edgeMaskBit > 0:
        badPixelMask |= 1 << edgeMaskBit

    logging.Trace('lsst.ip.diffim.computePsfMatchingKernelForMaskedImage', 3, 
               """Policy:
   * maximumFootprintResidualMean=%r
   * maximumFootprintResidualVariance=%r
   * maximumKernelSumNSigma=%r
   * maximumIterationsKernelSum=%r
   * debugIO=%r
   Mask planes (GetMaskPlanesFromPolicy = %r):
   * badMaskBit=%r
   * edgeMaskBit=%r
   """ % (
        maximumFootprintResidualMean,
        maximumFootprintResidualVariance,
        maximumKernelSumNSigma,
        maximumIterationsKernelSum,
        debugIO,
        GetMaskPlanesFromPolicy,
        badMaskBit,
        edgeMaskBit,
        ))
    
    kImage = afw.ImageD(kernelCols, kernelRows)
    diffImContainerList = ip_diffim.vectorDiffImContainerD()
    for footprintID, iFootprintPtr in enumerate(footprintList):
        footprintBBox = iFootprintPtr.getBBox()
        logging.Trace('lsst.ip.diffim.computePsfMatchingKernelForMaskedImage', 5,
            'Footprint %d = %d,%d -> %d,%d' % (footprintID,
            footprintBBox.min().x(), footprintBBox.min().y(),
            footprintBBox.max().x(), footprintBBox.max().y()))

        imageToConvolveStampPtr    = imageToConvolve.getSubImage(footprintBBox)
        imageToNotConvolveStampPtr = imageToNotConvolve.getSubImage(footprintBBox)

        if debugIO:
            imageToConvolveStampPtr.writeFits('tFits_%d' % (footprintID,))
            imageToNotConvolveStampPtr.writeFits('sFits_%d' % (footprintID,))
        
        # Find best single kernel for this stamp
        kernelCoeffList, background = ip_diffim.computePsfMatchingKernelForPostageStamp(
            imageToConvolveStampPtr.get(),
            imageToNotConvolveStampPtr.get(),
            kernelBasisList,
            policy,
        )
        footprintKernelPtr = afw.LinearCombinationKernelPtr(
            afw.LinearCombinationKernel(kernelBasisList, kernelCoeffList))
        
        # Structure holding information about this footprint and its fit to a kernel
        diffImFootprintContainer                    = ip_diffim.DiffImContainerD()
        diffImFootprintContainer.id                 = footprintID
        
        diffImFootprintContainer.setImageToNotConvolve(imageToNotConvolveStampPtr.get())
        diffImFootprintContainer.setImageToConvolve(imageToConvolveStampPtr.get())
        diffImFootprintContainer.isGood             = True
        diffImFootprintContainer.diffImFootprintPtr = iFootprintPtr

        fpMin = footprintBBox.min()
        fpMax = footprintBBox.max()
        diffImFootprintContainer.colcNorm = (fpMin.x() + fpMax.x()) / 2.0
        diffImFootprintContainer.rowcNorm = (fpMin.y() + fpMax.y()) / 2.0

        # Append to vectors collectively; might want to make this a method on a class to synchronize this
        diffImFootprintContainer.addKernel(footprintKernelPtr,
                                           ip_diffim.MaskedImageDiffimStats(),
                                           background,
                                           footprintKernelPtr.computeImage(kImage, 0.0, 0.0, False))

        logging.Trace('lsst.ip.diffim.computePsfMatchingKernelForMaskedImage', 4,
                   'Developing kernel %d at %d,%d -> %d,%d.  Kernel Sum = %.3f, Background = %.3f' % (
            diffImFootprintContainer.id,
            footprintBBox.min().x(), footprintBBox.min().y(),
            footprintBBox.max().x(), footprintBBox.max().y(),
            diffImFootprintContainer.kernelSums[diffImFootprintContainer.nKernel],
            diffImFootprintContainer.backgrounds[diffImFootprintContainer.nKernel],
            ))

        # Calculate and fill in difference image stats
        ip_diffim.computeDiffImStats(diffImFootprintContainer, diffImFootprintContainer.nKernel, policy)

        diffImContainerList.append(diffImFootprintContainer)
        

    # Compare kernel sums; things that are bad (e.g. variable stars, moving objects, cosmic rays)
    # will have a kernel sum far from the mean; reject these.

    logging.Trace('lsst.ip.diffim.computePsfMatchingKernelForMaskedImage', 3,
               'Rejecting kernels with deviant kSums');

    anyRejected = False
    kernelSumVector = numpy.zeros(len(diffImContainerList))
    for iterNum in xrange(maximumIterationsKernelSum):
        # compute statistics on current list of good kernels
        nGoodKernels = 0
        kMeanMean = 0.0
        kMeanVariance = 0.0

        for containerID in xrange(len(diffImContainerList)):
            # for diffImContainer in diffImContainerList: # THIS DOES NOT WORK; NEITHER DOES ENUMERATE
            
            if not diffImContainerList[containerID].isGood:
                continue
            whichKernel = diffImContainerList[containerID].nKernel
            
            kernelSumVector[nGoodKernels] = diffImContainerList[containerID].kernelSums[whichKernel]
            kMeanMean     += diffImContainerList[containerID].getFootprintResidualMean()
            kMeanVariance += diffImContainerList[containerID].getFootprintResidualVariance()
            nGoodKernels  += 1

        if nGoodKernels == 0:
            raise pex_ex.LsstOutOfRange("No good kernels found")
            
        kSumMean       = kernelSumVector[:nGoodKernels+1].mean()
        kSumStdDev     = kernelSumVector[:nGoodKernels+1].std()
        kMeanMean     /= float(nGoodKernels)
        kMeanVariance /= float(nGoodKernels)
        
        # reject kernels with aberrent statistics
        numRejected = 0
        for containerID in xrange(len(diffImContainerList)):
            if not diffImContainerList[containerID].isGood:
                continue
            whichKernel = diffImContainerList[containerID].nKernel
            
            if numpy.fabs(diffImContainerList[containerID].kernelSums[whichKernel] - kSumMean) > (kSumStdDev * maximumKernelSumNSigma):
                diffImContainerList[containerID].isGood = False
                numRejected += 1
                logging.Trace('lsst.ip.diffim.computePsfMatchingKernelForMaskedImage', 5, 
                           '# Kernel %d (kSum=%.3f) REJECTED due to bad kernel sum (mean=%.3f, std=%.3f)' %
                           (diffImContainerList[containerID].id,
                            diffImContainerList[containerID].kernelSums[whichKernel],
                            kSumMean, kSumStdDev))
                diffImLog.log(Log.INFO,
                              '#Kernel %d (kSum=%.3f) REJECTED due to bad kernel sum (mean=%.3f, std=%.3f)' %
                              (diffImContainerList[containerID].id,
                               diffImContainerList[containerID].kernelSums[whichKernel],
                               kSumMean, kSumStdDev))

        logging.Trace('lsst.ip.diffim.computePsfMatchingKernelForMaskedImage', 4, 
                   'Iteration %d, rejected %d kernels : Kernel Sum = %0.3f (%0.3f); Average Mean Residual = %0.3f; Average Residual Variance = %0.3f' %
                   (iterNum, numRejected, kSumMean, kSumStdDev, kMeanMean, kMeanVariance))
        diffImLog.log(Log.INFO,
                      'Iteration %d, rejected %d kernels : Kernel Sum = %0.3f (%0.3f); Average Mean Residual = %0.3f; Average Residual Variance = %0.3f' %
                      (iterNum, numRejected, kSumMean, kSumStdDev, kMeanMean, kMeanVariance))

        if numRejected == 0:
            # rejected all bad kernels; time to quit!
            break
    else:
        logging.Trace('lsst.ip.diffim.computePsfMatchingKernelForMaskedImage', 3, 
                   'Detection of kernels with deviant kSums reached its limit of %d iterations' % (
            maximumIterationsKernelSum,
            ))
      
    # In the end we want to test if the kernelInBasisList is Delta Functions; if so, do PCA
    # For DC2 we just do it
    kernelOutBasisList = ip_diffim.computePcaKernelBasis(diffImContainerList, policy)
    
    kernelPtr = ip_diffim.computeSpatiallyVaryingPsfMatchingKernel(
        kernelFunctionPtr,
        backgroundFunctionPtr,
        diffImContainerList,
        kernelOutBasisList,
        policy,
    )

    # Lets do some debugging here of the kernel sum; first from the good Footprints
    kernelSums = []
    for containerID in xrange(len(diffImContainerList)):
        if not diffImContainerList[containerID].isGood:
            continue
        whichKernel = diffImContainerList[containerID].nKernel
        kernelSums.append( diffImContainerList[containerID].kernelSums[whichKernel] )
    kernelSumArray = numpy.array(kernelSums)
    diffImLog.log(Log.INFO,
                  'Final Kernel Sum from %d Footprints : %0.3f (%0.3f)' % 
                  (len(kernelSumArray), kernelSumArray.mean(), kernelSumArray.std()))
        
        
    # Next from the corners
    kernelSums = []
    for nCol in [0, imageToConvolve.getCols()]:
        for nRow in [0, imageToConvolve.getRows()]:
            kernelSums.append( kernelPtr.computeImage(kImage, nCol, nRow, False) )
    kernelSumArray = numpy.array(kernelSums)
    diffImLog.log(Log.INFO,
                  'Final Kernel Sum from Image Corners : %0.3f (%0.3f)' % 
                  (kernelSumArray.mean(), kernelSumArray.std()))
    
    return kernelPtr
