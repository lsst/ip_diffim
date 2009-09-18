import eups
import sys, os, optparse
import numpy

import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.image.testUtils as imTestUtils
import lsst.ip.diffim as ipDiffim
import lsst.sdqa as sdqa

from lsst.pex.logging import Trace
from lsst.pex.logging import Log
from lsst.pex.policy import Policy

from lsst.ip.diffim import subtractMaskedImage

# For degugging needs
import pdb

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
    defDataDir = eups.productDir('afwdata') 
    if defDataDir == None:
        print 'Error: afwdata not set up'
        sys.exit(1)
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
    parser.add_option('-a', '--alard', action='store_true', default=False, help='use alard basis')
    parser.add_option('-p', '--policy', default=defPolicyPath, help='policy file')
    parser.add_option('-d', '--debugIO', action='store_true', default=False,
        help='write diagnostic intermediate files')
    parser.add_option('-s', '--display', action='store_true', default=False,
        help='display diagnostic images')
    parser.add_option('-v', '--verbosity', type=int, default=defVerbosity,
        help='verbosity of diagnostic trace messages; 1 for just warnings, more for more information')
    parser.add_option('-i', '--invert', action='store_true', default=False,
        help='invert the image to convolve')
                      
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

    # Need to conform mask planes since these data are relatively old
    templateMaskedImage = afwImage.MaskedImageF(templatePath)
    templateMaskedImage.getMask().conformMaskPlanes(afwImage.MaskU(0,0).getMaskPlaneDict())
    scienceMaskedImage  = afwImage.MaskedImageF(sciencePath)
    scienceMaskedImage.getMask().conformMaskPlanes(afwImage.MaskU(0,0).getMaskPlaneDict())
    policy              = Policy.createPolicy(policyPath)
    
    if options.debugIO:
        print 'DebugIO =', options.debugIO
        policy.set('debugIO', True)

    invert = False
    if options.invert:
        print 'Invert =', options.invert
        invert = True

    alard = False
    if options.alard:
        print 'Using Alard basis'
        alard = True

    display = False
    if options.display:
        print 'Displaying Images'
        display = True
        
    if options.verbosity > 0:
        print 'Verbosity =', options.verbosity
        Trace.setVerbosity('lsst.ip.diffim', options.verbosity)
        
    log = Log(Log.getDefaultLog(),
              "ip.diffim.subtractMaskedImage")

    if policy.get('spatialKernelTesting') == True:
        testingLoop(templateMaskedImage, scienceMaskedImage, policy)
    else:
        differenceMaskedImage, sKernel, bgFunction, sdqaList =  subtractMaskedImage(templateMaskedImage,
                                                                                    scienceMaskedImage,
                                                                                    policy, log,
                                                                                    useAlard=alard,
                                                                                    invert=invert,
                                                                                    display=display)
    differenceMaskedImage.writeFits(outputPath)
    


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
