#!/usr/bin/env python
import os
import pdb
import sys
import numpy
import math
import unittest

import lsst.mwi.data as mwiData
import lsst.mwi.policy
import lsst.fw.Core.fwLib as fw
import lsst.detection.detectionLib as detection
import imageprocLib as imageproc
import lsst.mwi.utils as mwiu

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

    mwiu.Trace('lsst.imageproc.computePsfMatchingKernelForMaskedImage', 3, 
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
    
    kImage = fw.ImageD(kernelCols, kernelRows)
    diffImContainerList = imageproc.vectorDiffImContainer_D()
    for footprintID, iFootprintPtr in enumerate(footprintList):
        footprintBBox = iFootprintPtr.getBBox()
        mwiu.Trace('lsst.imageproc.computePsfMatchingKernelForMaskedImage', 5,
            'Footprint %d = %d,%d -> %d,%d' % (footprintID,
            footprintBBox.min().x(), footprintBBox.min().y(),
            footprintBBox.max().x(), footprintBBox.max().y()))

        imageToConvolveStampPtr = imageToConvolve.getSubImage(footprintBBox)
        imageToNotConvolveStampPtr  = imageToNotConvolve.getSubImage(footprintBBox)

        if debugIO:
            imageToConvolveStampPtr.writeFits('tFits_%d' % (footprintID,))
            imageToNotConvolveStampPtr.writeFits('sFits_%d' % (footprintID,))
        
        # Find best single kernel for this stamp
        kernelCoeffList, background = imageproc.computePsfMatchingKernelForPostageStamp(
            imageToConvolveStampPtr.get(),
            imageToNotConvolveStampPtr.get(),
            kernelBasisList,
            policy,
        )
        footprintKernelPtr = fw.LinearCombinationKernelDPtr(
            fw.LinearCombinationKernelD(kernelBasisList, kernelCoeffList))
        
        # Structure holding information about this footprint and its fit to a kernel
        diffImFootprintContainer = imageproc.DiffImContainer_D()
        diffImFootprintContainer.id = footprintID
        diffImFootprintContainer.isGood = True
        diffImFootprintContainer.diffImFootprintPtr = iFootprintPtr
        diffImFootprintContainer.diffImKernelPtr = footprintKernelPtr
        diffImFootprintContainer.background = background
        diffImFootprintContainer.kSum = footprintKernelPtr.computeImage(kImage, 0.0, 0.0, False)

        fpMin = footprintBBox.min()
        fpMax = footprintBBox.max()
        diffImFootprintContainer.colcNorm = (fpMin.x() + fpMax.x()) / 2.0
        diffImFootprintContainer.rowcNorm = (fpMin.y() + fpMax.y()) / 2.0

        mwiu.Trace('lsst.imageproc.computePsfMatchingKernelForMaskedImage', 4,
                   'Developing kernel %d at %d,%d -> %d,%d.  Background = %.3f' % (
            diffImFootprintContainer.id,
            footprintBBox.min().x(), footprintBBox.min().y(),
            footprintBBox.max().x(), footprintBBox.max().y(),
            diffImFootprintContainer.background,
            ))

        # Calculate the residual of the subtracted image
        convolvedImageStamp = fw.convolve(
            imageToConvolveStampPtr.get(),
            footprintKernelPtr.get(),
            edgeMaskBit,
            False,
        )

        # Add in background
        convolvedImageStamp += background
        if debugIO:
            kImage.writeFits('kFits_%d.fits' % (footprintID,))
            convolvedImageStamp.writeFits('cFits_%d' %  (footprintID,))
        
        
        # Do actual subtraction
        convolvedImageStamp -= imageToNotConvolveStampPtr.get()
        convolvedImageStamp *= -1.0

        if debugIO:
            convolvedImageStamp.writeFits('dFits_%d' % (footprintID,))

        # Calculate stats of the difference image
        nGoodPixels, meanOfResiduals, varianceOfResiduals = imageproc.calculateMaskedImageResiduals(
            convolvedImageStamp, badPixelMask)

        diffImFootprintContainer.footprintResidualMean = meanOfResiduals
        diffImFootprintContainer.footprintResidualVariance = varianceOfResiduals/nGoodPixels # per pixel
        
        if not numpy.isfinite(diffImFootprintContainer.footprintResidualMean) \
            or not numpy.isfinite(diffImFootprintContainer.footprintResidualVariance):
            diffImFootprintContainer.isGood = False
            mwiu.Trace('lsst.imageproc.computePsfMatchingKernelForMaskedImage', 4, 
                       '# Kernel %d (kSum=%.3f), analysis of footprint failed : %.3f %.3f (%d pixels)' % (
                diffImFootprintContainer.id,
                diffImFootprintContainer.kSum,
                diffImFootprintContainer.footprintResidualMean,
                diffImFootprintContainer.footprintResidualVariance,
                nGoodPixels,
                ))
            
        elif math.fabs(diffImFootprintContainer.footprintResidualMean) > maximumFootprintResidualMean:
            diffImFootprintContainer.isGood = False
            mwiu.Trace('lsst.imageproc.computePsfMatchingKernelForMaskedImage', 5, 
                       '# Kernel %d (kSum=%.3f), bad mean residual of footprint : %.3f (%d pixels)' % (
                diffImFootprintContainer.id,
                diffImFootprintContainer.kSum,
                diffImFootprintContainer.footprintResidualMean,
                nGoodPixels,
                ))
            
        elif diffImFootprintContainer.footprintResidualVariance > maximumFootprintResidualVariance:
            diffImFootprintContainer.isGood = False
            mwiu.Trace('lsst.imageproc.computePsfMatchingKernelForMaskedImage', 4, 
                       '# Kernel %d (kSum=%.3f), bad residual variance of footprint : %.3f (%d pixels)' % (
                diffImFootprintContainer.id,
                diffImFootprintContainer.kSum,
                diffImFootprintContainer.footprintResidualVariance,
                nGoodPixels,
                ))

        else:
            mwiu.Trace('lsst.imageproc.computePsfMatchingKernelForMaskedImage', 5, 
                       'Kernel %d (kSum=%.3f), mean and variance of residuals in footprint : %.3f %.3f (%d pixels)' % (
                diffImFootprintContainer.id,
                diffImFootprintContainer.kSum,
                diffImFootprintContainer.footprintResidualMean,
                diffImFootprintContainer.footprintResidualVariance,
                nGoodPixels,
                ))
            
        diffImContainerList.append(diffImFootprintContainer)
        

    # Compare kernel sums; things that are bad (e.g. variable stars, moving objects, cosmic rays)
    # will have a kernel sum far from the mean; reject these.

    mwiu.Trace('lsst.imageproc.computePsfMatchingKernelForMaskedImage', 3,
               'Rejecting kernels with deviant kSums');

    anyRejected = False
    kernelSumVector = numpy.zeros(len(diffImContainerList))
    for iterNum in xrange(maximumIterationsKernelSum):
        # compute statistics on current list of good kernels
        nGoodKernels = 0
        kMeanMean = 0.0
        kMeanVariance = 0.0
        for ind, diffImContainer in enumerate(diffImContainerList):
            if not diffImContainer.isGood:
                continue
            kernelSumVector[nGoodKernels] = diffImContainer.kSum
            kMeanMean += diffImContainer.footprintResidualMean
            kMeanVariance += diffImContainer.footprintResidualVariance
            nGoodKernels += 1
            
        kSumMean = kernelSumVector[0:nGoodKernels].mean()
        kSumStdDev = kernelSumVector[0:nGoodKernels].std()
        kMeanMean /= float(nGoodKernels)
        kMeanVariance /= float(nGoodKernels)
        
        # reject kernels with aberrent statistics
        numRejected = 0
        for diffImContainer in diffImContainerList:
            if not diffImContainer.isGood:
                continue
            if math.fabs(diffImContainer.kSum - kSumMean) > (kSumStdDev * maximumKernelSumNSigma):
                diffImContainer.isGood = False
                numRejected += 1
                mwiu.Trace('lsst.imageproc.computePsfMatchingKernelForMaskedImage', 5, 
                           '# Kernel %d (kSum=%.3f) REJECTED due to bad kernel sum (mean=%.3f, std=%.3f)' %
                           (diffImContainer.id, diffImContainer.kSum, kSumMean, kSumStdDev))

        mwiu.Trace('lsst.imageproc.computePsfMatchingKernelForMaskedImage', 4, 
                   'Iteration %d, rejected %d kernels : Kernel Sum = %0.3f (%0.3f); Average Mean Residual = %0.3f; Average Residual Variance = %0.3f' %
                   (iterNum, numRejected, kSumMean, kSumStdDev, kMeanMean, kMeanVariance))

        if numRejected == 0:
            # rejected all bad kernels; time to quit!
            break;
    else:
        mwiu.Trace('lsst.imageproc.computePsfMatchingKernelForMaskedImage', 3, 
                   'Detection of kernels with deviant kSums reached its limit of %d iterations' % (
            maximumIterationsKernelSum,
            ))
      
    # In the end we want to test if the kernelInBasisList is Delta Functions; if so, do PCA
    # For DC2 we just do it
    kernelOutBasisList = imageproc.computePcaKernelBasis(diffImContainerList, policy)
    
    kernelPtr = imageproc.computeSpatiallyVaryingPsfMatchingKernel(
        kernelFunctionPtr,
        backgroundFunctionPtr,
        diffImContainerList,
        kernelOutBasisList,
        policy,
    )
    
    return kernelPtr
