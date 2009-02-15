import eups
import sys, os, optparse
import numpy

import lsst.afw.image as afwImage
import lsst.daf.base  as dafBase
import lsst.afw.math  as afwMath

from   lsst.pex.logging import Log
from   lsst.pex.logging import Trace
from   lsst.pex.policy  import Policy
import lsst.pex.exceptions as Exceptions

# relative imports, since these are in __init__.py
import diffimLib 
import diffimTools

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
        if scPtr.getCurrentModel().getSdqaStatus():
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
    bgFit      = diffimTools.fitFunction(bgFunction, bgValues, bgErrors,
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
            pFit      = diffimTools.fitFunction(pFunction, pValues, pErrors,
                                                cCol, cRow, policy)
            pFunction.setParameters(pFit.parameterList)
            
            Trace('lsst.ip.diffim.spatialModelByPixel', 5,
                  'Pixel %d fit parameters : %s' % (np, ' '.join([('%10.3e' % (x)) for x in pFit.parameterList])))

            pFunctionList.append(pFunction)
            np += 1
            
            # do the various fitting techniques here
            #######

    return bgFunction, pFunctionList


def evaluateModelByPixel(spatialCells, bgFunction, sKernel, policy, reject=True):
    kCols          = policy.get('kernelCols')
    kRows          = policy.get('kernelRows')

    nRejected = 0
    nGood     = 0
    
    idList   = []
    for scID, scPtr in enumerate(spatialCells):
        # Is the cell usable at all?
        if not scPtr.isUsable():
            continue
        # Is the contained model usable?
        if scPtr.getCurrentModel().getSdqaStatus():
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
            
        kernel  = afwMath.FixedKernel(sImage)
        # Create difference image using Kernel model
        diffIm  = diffim.convolveAndSubtract(spatialCells[ idList[idx] ].getCurrentModel().getMiToConvolvePtr().get(),
                                             spatialCells[ idList[idx] ].getCurrentModel().getMiToNotConvolvePtr().get(),
                                             kernel, bgValue)

        # Find quality of difference image
        diffImStats = diffimLib.DifferenceImageStatisticsF(diffIm)
        if (diffImStats.evaluateQuality(policy) == False):
            if (reject == True):
                # This is only bad in the context of the spatial model
                # May be needed in the future
                #
                # spatialCells[ idList[idx] ].getCurrentModel().setSdqaStatus(False)
                spatialCells[ idList[idx] ].incrementModel()
                nRejected += 1

                label = 'Rejected:'
            else:
                label = 'Poor'
        else:
            nGood += 1
            label  = 'OK'

        Trace('lsst.ip.diffim.evaluateModelByPixel', 5,
              '%s Kernel %d : %s Spatial residuals = %.2f +/- %.2f sigma' %
              (spatialCells[ idList[idx] ].getLabel(),
               spatialCells[ idList[idx] ].getCurrentModel().getID(),
               label, diffImStats.getResidualMean(), diffImStats.getResidualStd()))

    Trace('lsst.ip.diffim.evaluateModelByPixel', 3,
          'Spatial model by pixel : %d / %d Kernels acceptable' % (nGood, nCells))

#        if policy.get('debugPlot') == True:
#            ipDiffimDebug.plotDiffImQuality1(goodDifiList[i],
#                                             diffIm,
#                                             kernelPtr,
#                                             label='Spatial %s kernel %d' % (label, goodDifiList[i].getID()),
#                                             outfile='SKernel_%s%d.ps' % (label, goodDifiList[i].getID())
#                                             )

    return nRejected

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
        if scPtr.getCurrentModel().getSdqaStatus():
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
        kernel = afwMath.FixedKernel(kImage)

        #if policy.get('debugIO'):
        #    kImage.writeFits('skImage1_%d.fits' % (idx))
        #    spatialCells[ idList[idx] ].getCurrentModel().getKernelPtr().computeNewImage(False)[0].writeFits('kImage_%d.fits' % (idx))
        
        # Create difference image using Kernel model
        diffIm = diffimLib.convolveAndSubtract(spatialCells[ idList[idx] ].getCurrentModel().getMiToConvolvePtr().get(),
                                               spatialCells[ idList[idx] ].getCurrentModel().getMiToNotConvolvePtr().get(),
                                               kernel, bgValue)

        # Find quality of difference image
        diffImStats = diffimLib.DifferenceImageStatisticsF(diffIm)
        if (diffImStats.evaluateQuality(policy) == False):
            if (reject == True):
                # This is only bad in the context of the spatial model
                # May be needed in the future
                #
                # spatialCells[ idList[idx] ].getCurrentModel().setSdqaStatus(False)                
                spatialCells[ idList[idx] ].incrementModel()
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


################
################
#####  Pca-based spatial analysis
################
################

def spatialModelKernelPca(spatialCells, policy, id):
    from lsst.ip.diffim.runPca import runPca

    kCols = policy.get('kernelCols')
    kRows = policy.get('kernelRows')
    
    idList   = []
    for scID, scPtr in enumerate(spatialCells):
        # Is the cell usable at all?
        if not scPtr.isUsable():
            continue
        # Is the contained model usable?
        if scPtr.getCurrentModel().getSdqaStatus():
            idList.append( scID )
            
    nCells = len(idList)
    
    # matrix to invert
    M = numpy.zeros((kCols*kRows, nCells))
    for idx in range(nCells):
        kernelImage  = spatialCells[ idList[idx] ].getCurrentModel().getKernelPtr().computeNewImage(False)[0]
        M[:,idx]     = diffimTools.imageToVector(kernelImage)

    # Call numpy Pca
    meanM, U, eVal, eCoeff = runPca(M, policy)

    # Turn principal components into Kernels
    mKernel       = diffimTools.vectorToKernelPtr( meanM, kCols, kRows )
    if policy.get('debugIO'):
        diffimTools.vectorToImage(meanM, kCols, kRows).writeFits('mKernel%s.fits' % (id))

    eKernelVector = afwMath.VectorKernel()
    for i in range(U.shape[1]):
        eKernel   = diffimTools.vectorToKernelPtr( U[:,i], kCols, kRows )
        eKernelVector.append(eKernel)

        if policy.get('debugIO'):
            diffimTools.vectorToImage(U[:,i], kCols, kRows).writeFits('eKernel%s_%d.fits' % (id, i))
            
    return mKernel, eKernelVector, eVal, eCoeff
    

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
        if scPtr.getCurrentModel().getSdqaStatus():
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
    bgFit      = diffimTools.fitFunction(bgFunction, bgValues, bgErrors,
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
        eFit = diffimTools.fitFunction(eFunction,
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
    nGood     = 0
    
    idList   = []
    for scID, scPtr in enumerate(spatialCells):
        # Is the cell usable at all?
        if not scPtr.isUsable():
            continue
        # Is the contained model usable?
        if scPtr.getCurrentModel().getSdqaStatus():
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
        kernel  = afwMath.FixedKernel(sImage)
        diffIm  = diffimLib.convolveAndSubtract(spatialCells[ idList[idx] ].getCurrentModel().getMiToConvolvePtr().get(),
                                                spatialCells[ idList[idx] ].getCurrentModel().getMiToNotConvolvePtr().get(),
                                                kernel, bgValue)

        # Find quality of difference image
        diffImStats = diffimLib.DifferenceImageStatisticsF(diffIm)
        if (diffImStats.evaluateQuality(policy) == False):
            if (reject == True):
                # This is only bad in the context of the spatial model
                # May be needed in the future
                #
                # spatialCells[ idList[idx] ].getCurrentModel().setSdqaStatus(False)                
                spatialCells[ idList[idx] ].incrementModel()
                nRejected += 1

                label = 'Rejected:'
            else:
                label = 'Poor'
        else:
            nGood += 1
            label = 'OK'

        Trace('lsst.ip.diffim.evaluateModelByPca', 5,
              '%s Kernel %d : %s Pca residuals = %.2f +/- %.2f sigma' %
              (spatialCells[ idList[idx] ].getLabel(),
               spatialCells[ idList[idx] ].getCurrentModel().getID(),
               label, diffImStats.getResidualMean(), diffImStats.getResidualStd()))
        
    Trace('lsst.ip.diffim.evaluateModelByPca', 3,
          'Spatial model by PCA : %d / %d Kernels acceptable' % (nGood, nCells))

    return nRejected
