import eups
import sys, os, optparse
import numpy

import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim

from lsst.pex.logging import Trace
from lsst.pex.logging import Log
from lsst.pex.policy import Policy

# For degugging needs
import pdb

def subtractMaskedImage(templateMaskedImage, scienceMaskedImage, policy, fpList=None):
    # Make sure they are the EXACT same dimensions in pixels
    # This is non-negotiable
    assert (templateMaskedImage.getDimensions() == scienceMaskedImage.getDimensions())
    
    kCols = policy.get('kernelCols')
    kRows = policy.get('kernelRows')

    kBasisList = ipDiffim.generateDeltaFunctionKernelSet(kCols, kRows)
    kFunctor   = ipDiffim.PsfMatchingFunctorF(kBasisList)

    if fpList == None:
        # Need to find own footprints
        fpList = ipDiffim.getCollectionOfFootprintsForPsfMatching(templateMaskedImage,
                                                                  scienceMaskedImage,
                                                                  policy)

    # Set up grid for spatial model
    spatialCells = ipDiffim.createSpatialModelKernelCells(templateMaskedImage,
                                                          scienceMaskedImage,
                                                          fpList,
                                                          kFunctor,
                                                          policy)

    # Set up fitting loop 
    maxSpatialIterations = policy.getInt('maxSpatialIterations')
    rejectKernels        = policy.getBool('spatialKernelRejection')
    nRejected = -1
    nIter     =  0
    
    # And fit spatial kernel model
    if policy.get('spatialKernelModel') == 'pca':
        # Fit spatial variation of principal components

        minPrincipalComponents = policy.getInt('minPrincipalComponents')
        maxPrincipalComponents = policy.getInt('maxPrincipalComponents')
        fracEigenVal           = policy.getFloat('fracEigenVal')
        
        while (nRejected != 0) and (nIter < maxSpatialIterations):
            # Run the PCA
            mKernel, eKernelVector, eVal, eCoeff = spatialModelKernelPca(spatialCells, policy)

            # Make the decision on how many components to use
            eFrac  = num.cumsum(eVal)
            eFrac /= eFrac[-1]
            nEval  = len(num.where(eFrac > fracEigenVal))
            nEval  = min(minPrincipalComponents, len(eKernelVector))
            nEval  = min(nEval, maxPrincipalComponents)
            nEval  = max(nEval, minPrincipalComponents)

            # do spatial fit here by Principal Component
            sKernel, bgFunction = spatialModelByPca(spatialCells, eCoeff, nEval, policy)

            # Evaluate quality of spatial fit
            nRejected = evaluateModelByPca(spatialCells, bgFunction, sKernel,
                                           policy, reject=rejectKernels)
                
            nIter += 1

    elif policy.get('spatialKernelModel') == 'pixel':
        # Fit function to each pixel

        while (nRejected != 0) and (nIter < maxSpatialIterations):
            # do spatial fit here pixel by pixel
            sKernel, bgFunction = spatialModelByPixel(spatialCells, policy)
            # and check quality
            nRejected  = evaluateModelByPixel(spatialCells,
                                              bgFunction, sKernel, 
                                              policy, reject=rejectKernels)
            nIter += 1
        
    else:
        # All that is supported
        # Throw exception!
        pass

    differenceMaskedImage = ipDiffim.convolveAndSubtract(templateMaskedimage,
                                                         scienceMaskedImage,
                                                         sKernel,
                                                         bgFunction)
    # Do something about this eventually
    Sdqa = None

    # What kind of metadata do we add here to MaskedImage?
    
    return differenceMaskedImage, sKernel, bgFunction, SdqaList

def testingLoop(templateMaskedImage, scienceMaskedImage, policy, fpList=None):

    kCols = policy.get('kernelCols')
    kRows = policy.get('kernelRows')

    kBasisList = ipDiffim.generateDeltaFunctionKernelSet(kCols, kRows)
    kFunctor   = ipDiffim.PsfMatchingFunctorF(kBasisList)

    if fpList == None:
        # Need to find own footprints
        fpList = ipDiffim.getCollectionOfFootprintsForPsfMatching(templateMaskedImage,
                                                                  scienceMaskedImage,
                                                                  policy)
    # LOOP 1 : convolution vs deconvolution
    Trace('lsst.ip.diffim', 1, 'SC List 1')
    spatialCellsC = ipDiffim.createSpatialModelKernelCells(templateMaskedImage,
                                                           scienceMaskedImage,
                                                           fpList,
                                                           kFunctor,
                                                           policy,
                                                           cFlag='c')
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
    spatialCellsD = ipDiffim.createSpatialModelKernelCells(scienceMaskedImage,
                                                           templateMaskedImage,
                                                           fpList,
                                                           kFunctor,
                                                           policy,
                                                           cFlag='d')
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

    
def main():
    defDataDir = eups.productDir('afwdata') or ''
    imageProcDir = eups.productDir('ip_diffim')
    if imageProcDir == None:
        print 'Error: could not set up ip_diffim'
        sys.exit(1)

    defSciencePath  = os.path.join(defDataDir, 'CFHT', 'D4', 'cal-53535-i-797722_1')
    defTemplatePath = os.path.join(defDataDir, 'CFHT', 'D4', 'cal-53535-i-797722_1_tmpl')
    defPolicyPath   = os.path.join(imageProcDir, 'pipeline', 'ImageSubtractStageDictionary.paf')
    defOutputPath   = 'diffImage'
    defVerbosity    = 0
    
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
    
    sciencePath  = getArg(0, defSciencePath)
    templatePath = getArg(1, defTemplatePath)
    outputPath   = getArg(2, defOutputPath)
    policyPath   = options.policy
    
    print 'Science image: ', sciencePath
    print 'Template image:', templatePath
    print 'Output image:  ', outputPath
    print 'Policy file:   ', policyPath
    
    templateMaskedImage = afwImage.MaskedImageF(templatePath)
    scienceMaskedImage  = afwImage.MaskedImageF(sciencePath)
    policy              = Policy.createPolicy(policyPath)
    
    if options.debugIO:
        policy.set('debugIO', True)

    if options.verbosity > 0:
        print 'Verbosity =', options.verbosity
        Trace.setVerbosity('lsst.ip.diffim', options.verbosity)

    if policy.get('spatialKernelTesting') == True:
        testingLoop(templateMaskedImage, scienceMaskedImage, policy)
    else:
        subtractMaskedImage(templateMaskedImage, scienceMaskedImage, policy)
    


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
