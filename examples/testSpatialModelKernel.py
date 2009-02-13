import eups
import sys, os, optparse
import numpy

# SWIGged interfaces
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim

from lsst.pex.logging import Log
from lsst.pex.logging import Trace
from lsst.pex.policy import Policy
import lsst.pex.exceptions as Exceptions

# For degugging needs
import pdb

################
################
#####  Main fitting loop
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
            sKernel = afwMath.LinearCombinationKernel(kBasisList, sKernelFunc)
            sKernel.setSpatialParameters(kParams)
            
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
    
    templateMaskedImage = afwImage.MaskedImageF(templatePath)
    scienceMaskedImage  = afwImage.MaskedImageF(sciencePath)

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
    #policy.set('minimumCleanFootprints', 5)
    #policy.set('footprintDetectionThreshold', 5.)

    # if you are convolving the template
    # policy.set('iterateKernel', False)
    # if you are convolving the image
    # policy.set('iterateKernel', True)

    # NOTE : if you get a runtime error like
    # Wrong number of arguments for overloaded function 'getCollectionOfFootprintsForPsfMatching'.
    # this is because Swig gets confused over the type of mask
    # In particular this line in diffimLib.i causes the problems :
    # typedef unsigned short boost::uint16_t;
    fpList = ipDiffim.getCollectionOfFootprintsForPsfMatching(templateMaskedImage,
                                                              scienceMaskedImage,
                                                              policy)

    # LOOP 1 : convolution vs deconvolution
    Trace('lsst.ip.diffim', 1, 'SC List 1')
    spatialCellsC = ipDiffim.createSpatialModelKernelCells(fpList,
                                                           templateMaskedImage,
                                                           scienceMaskedImage,
                                                           kBasisList,
                                                           policy,
                                                           cFlag='c')
    if policy.get('spatialKernelTesting') == True:
        resultsC = spatialKernelTesting(spatialCellsC, kBasisList, policy, 0)
        if resultsC == None:
            Trace('lsst.ip.diffim', 3, 'Spatial testing failed')
        else:
            bgListC, kListC = resultsC
            for key1 in bgListC.keys():
                for key2 in range(len(bgListC[key1])):
                    background = bgListC[key1][key2]
                    kernel     = kListC[key1][key2]
    
                    diffIm = ipDiffim.convolveAndSubtract(
                        templateMaskedImage,
                        scienceMaskedImage,
                        kernel,
                        background)
                    diffIm.writeFits('diffImC_%s_%d' % (key1, key2))

    #
    ###
    #

    Trace('lsst.ip.diffim', 1, 'SC List 2')
    spatialCellsD = ipDiffim.createSpatialModelKernelCells(fpList,
                                                           scienceMaskedImage,
                                                           templateMaskedImage,
                                                           kBasisList,
                                                           policy,
                                                           cFlag='d')
    if policy.get('spatialKernelTesting') == True:
        resultsD = spatialKernelTesting(spatialCellsD, kBasisList, policy, 1)
        if resultsD == None:
            Trace('lsst.ip.diffim', 3, 'Spatial testing failed')
        else:
            bgListD, kListD = resultsD
            for key1 in fitListD.keys():
                for key2 in range(len(bgListD[key1])):
                    background = bgListD[key1][key2]
                    kernel     = kListD[key1][key2]
    
                    diffIm = ipDiffim.convolveAndSubtract(
                        scienceMaskedImage,
                        templateMaskedImage,
                        kernel,
                        background)
                    diffIm.writeFits('diffImD_%s_%d' % (key1, key2))
    
    return


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
