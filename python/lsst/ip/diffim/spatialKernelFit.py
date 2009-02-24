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
        cCol[idx] = spatialCells[ idList[idx] ].getCurrentModel().getColc()
        cRow[idx] = spatialCells[ idList[idx] ].getCurrentModel().getRowc()

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

        # pretty sure kernels are double
        kImage    = afwImage.ImageD(kPtr.getDimensions())
        kErrImage = afwImage.ImageD(kErrPtr.getDimensions())

        kPtr.computeImage(kImage, False)
        kErrPtr.computeImage(kErrImage, False)
        
        kImages.append( kImage )
        kErrImages.append( kErrImage )
    
    
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
            
    # Turn list of pixel functions into spatially varying kernel
    sKernelFunc  = afwMath.PolynomialFunction2D(kSpatialOrder)
    kParams      = numpy.zeros( (kCols*kRows, sKernelFunc.getNParameters()) )
    for p in range(kCols*kRows):
        kParams[p] = pFunctionList[p].getParameters()
    # Create spatially varying kernel pointer
    sKernel = afwMath.LinearCombinationKernel(kBasisList, sKernelFunc)
    sKernel.setSpatialParameters(kParams)
    
    return sKernel, bgFunction, 


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
        cCol[idx] = spatialCells[ idList[idx] ].getCurrentModel().getColc()
        cRow[idx] = spatialCells[ idList[idx] ].getCurrentModel().getRowc()
    
    # Evaluate all the fits at the positions of the objects, create a
    # new kernel, then difference image, the calculate difference
    # image stats
    imStats = diffimLib.ImageStatisticsF()

    for idx in range(nCells):
        bgValue = bgFunction(cCol[idx], cRow[idx])
        sImage  = afwImage.ImageD(sKernel.getDimensions())
        sKernel.computeImage(sImage, False, cCol[idx], cRow[idx])
        kernel  = afwMath.FixedKernel(sImage)
        diffIm  = diffim.convolveAndSubtract(spatialCells[ idList[idx] ].getCurrentModel().getMiToConvolvePtr(),
                                             spatialCells[ idList[idx] ].getCurrentModel().getMiToNotConvolvePtr(),
                                             kernel, bgValue)

        # Find quality of difference image
        imStats.apply(diffIm)
        if (imStats.evaluateQuality(policy) == False):
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
               label, imStats.getResidualMean(), imStats.getResidualStd()))

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
        cCol[idx] = spatialCells[ idList[idx] ].getCurrentModel().getColc()
        cRow[idx] = spatialCells[ idList[idx] ].getCurrentModel().getRowc()
    
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
        diffIm = diffimLib.convolveAndSubtract(spatialCells[ idList[idx] ].getCurrentModel().getMiToConvolvePtr(),
                                               spatialCells[ idList[idx] ].getCurrentModel().getMiToNotConvolvePtr(),
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

def spatialModelKernelPca(spatialCells, policy, id=''):

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
        kPtr         = spatialCells[ idList[idx] ].getCurrentModel().getKernelPtr()
        kImage       = afwImage.ImageD(kPtr.getDimensions())
        kPtr.computeImage(kImage, False)
        M[:,idx]     = diffimTools.vectorFromImage(kImage)

    # Call numpy Pca
    meanM, U, eVal, eCoeff = diffimTools.runPca(M, policy)

    # Turn principal components into Kernels (which are double)
    mImage  = diffimTools.imageFromVector(meanM, kCols, kRows, retType=afwImage.ImageD)
    mKernel = afwMath.FixedKernel(mImage)
    if policy.get('debugIO'):
        mImage.writeFits('mKernel%s.fits' % (id))

    eKernelVector = afwMath.VectorKernel()
    for i in range(U.shape[1]):
        eImage  = diffimTools.imageFromVector(U[:,i], kCols, kRows, retType=afwImage.ImageD)
        eKernel = afwMath.FixedKernel(eImage)
        eKernelVector.append(eKernel)

        if policy.get('debugIO'):
            eImage.writeFits('eKernel%s_%d.fits' % (id, i))
            
    return mKernel, eKernelVector, eVal, eCoeff
    

def spatialModelByPca(spatialCells, mKernel, eKernelVector, eCoeffs, nEVal, policy, makeKernel=True):
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
        cCol[idx] = spatialCells[ idList[idx] ].getCurrentModel().getColc()
        cRow[idx] = spatialCells[ idList[idx] ].getCurrentModel().getRowc()

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
    for nc in range(nEVal):
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

    if makeKernel:
        # Build LinearCombinationKernel
        eKernelBases = afwMath.KernelListD()

        # Start with mean Kernel image
        eKernelBases.push_back(mKernel)
        # And append eigenKernel images
        for ek in range(nEVal):
            eKernelBases.push_back(eKernelVector[ek])
                    
        # Mean kernel has no spatial variation
        eKernelFunc   = afwMath.PolynomialFunction2D(kSpatialOrder)
        kParams       = numpy.zeros( (nEVal+1, eKernelFunc.getNParameters()) )
        kParams[0][0] = 1.0
        # Add already-fit-for spatial variation of eigenKernels
        for ek in range(nEVal):
            kParams[ek+1] = eFunctionList[ek].getParameters()
    
        # Create spatially varying eigenKernel pointer
        eKernel = afwMath.LinearCombinationKernel(eKernelBases, eKernelFunc)
        eKernel.setSpatialParameters(kParams)
        
        return eKernel, bgFunction
    else:
        # for spatial testing
        return eFunctionList, bgFunction


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
        cCol[idx] = spatialCells[ idList[idx] ].getCurrentModel().getColc()
        cRow[idx] = spatialCells[ idList[idx] ].getCurrentModel().getRowc()

    # Evaluate all the fits at the positions of the objects, create a
    # new kernel, then difference image, the calculate difference
    # image stats
    imStats = diffimLib.ImageStatisticsF()

    for idx in range(nCells):
        bgValue = bgFunction(cCol[idx], cRow[idx])
        eImage  = afwImage.ImageD(eKernel.getDimensions())
        eKernel.computeImage(eImage, False, cCol[idx], cRow[idx])
        kernel  = afwMath.FixedKernel(eImage)
        diffIm  = diffimLib.convolveAndSubtract(spatialCells[ idList[idx] ].getCurrentModel().getMiToConvolvePtr(),
                                                spatialCells[ idList[idx] ].getCurrentModel().getMiToNotConvolvePtr(),
                                                kernel, bgValue)

        # Find quality of difference image
        imStats.apply(diffIm)
        if (imStats.evaluateQuality(policy) == False):
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
               label, imStats.getMean(), imStats.getRms()))
        
    Trace('lsst.ip.diffim.evaluateModelByPca', 3,
          'Spatial model by PCA : %d / %d Kernels acceptable' % (nGood, nCells))

    return nRejected

################
################
#####  Testing loop
################
################

def spatialKernelTesting(spatialCells, kBasisList, policy, scID):
    kCols = policy.get('kernelCols')
    kRows = policy.get('kernelRows')
    rejectKernels = policy.get('spatialKernelRejection')
    
    try:
        ipDiffim.rejectKernelSumOutliers(spatialCells, policy)
    except:
        Trace('lsst.ip.diffim', 2,
              'LOOP %d FAILED; no good kernels' % (scID))
        return

    bgList            = {}
    bgList['spatial'] = []
    bgList['pca']     = []

    kList             = {}
    kList['spatial']  = []
    kList['pca']      = []
    
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
            sKernel, bgFunction = spatialModelByPixel(spatialCells, policy)

            # and check quality
            nRejected  = evaluateModelByPixel(spatialCells,
                                              bgFunction, sKernel, 
                                              policy, reject=rejectKernels)
            nIter += 1

        # In future versions of the code, there will be no need to make pointers like this.
        kList['spatial'].append( sKernel )
        bgList['spatial'].append( bgFunction )
        
        #############
        # PCA fit
        nRejected = 1
        nIter     = 0
        
        # LOOP 4a : PCA sigma clipping
        while (nRejected != 0) and (nIter < maxSpatialIterations):
            
            # Run the PCA
            mKernel, eKernelVector, eVal, eCoeff = spatialModelKernelPca(spatialCells, policy, scID)
            
            # Here we make a decision on how many eigenComponents to use based
            # on eVal, etc
            #
            # While we are testing, check them all
            
            # Find spatial variation of only those components
            # Remove this line after being done testing
            # We fit them all first
            bgFunction, eFunctionList = spatialModelByPca(spatialCells, eCoeff, len(eVal), policy, False)
    
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
                eKernelBases.push_back(mKernel)
                # Append eigenKernel images
                for ek in range(neVal):
                    eKernelBases.push_back(eKernelVector[ek])
                    
                # Mean kernel has no spatial variation
                eKernelFunc   = afwMath.PolynomialFunction2D(kSpatialOrder)
                kParams       = numpy.zeros( (neVal+1, eKernelFunc.getNParameters()) )
                kParams[0][0] = 1.0
                # Add already-fit-for spatial variation of eigenKernels
                for ek in range(neVal):
                    kParams[ek+1] = eFunctionList[ek].getParameters()
    
                # Create spatially varying eigenKernel pointer
                eKernel = afwMath.LinearCombinationKernel(eKernelBases, eKernelFunc)
                eKernel.setSpatialParameters(kParams)
        
                # Evaluate quality of spatial fit
                nRejected = evaluateModelByPca(spatialCells, bgFunction, eKernel,
                                               policy, reject=rejectKernels)
                
            nIter += 1

        # NOTE to self : you get the "final" PCA kernel here with all elements
        # In future versions of the code, there will be no need to make pointers like this.
        kList['pca'].append( eKernel )
        bgList['pca'].append( bgFunction )
                
    return bgList, kList
