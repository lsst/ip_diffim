import eups
import sys, os, optparse
import numpy

import lsst.daf.base as dafBase
import lsst.afw.image as afwImage

from lsst.pex.logging import Trace
from lsst.pex.logging import Log
from lsst.pex.policy import Policy

from lsst.ip.diffim import subtractMaskedImage

# For degugging needs
import pdb

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
    parser.add_option('-p', '--policy', default=defPolicyPath, help='policy file')
    parser.add_option('-v', '--verbosity', type=int, default=defVerbosity,
                      help='verbosity of Trace messages')
    parser.add_option('-i', '--invert', action='store_true', default=False,
                      help='invert the image to convolve')
    parser.add_option('-d', '--display', action='store_true', default=False,
                      help='display the images')
                      
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
    
    invert = False
    if options.invert:
        print 'Invert =', options.invert
        invert = True

    display = False
    if options.display:
        print 'Display =', options.display
        display = True

    if options.verbosity > 0:
        print 'Verbosity =', options.verbosity
        Trace.setVerbosity('lsst.ip.diffim', options.verbosity)
        
    log = Log(Log.getDefaultLog(),
              "ip.diffim.subtractMaskedImage")

    differenceMaskedImage, sKernel, bgFunction, kernelCellSet =  subtractMaskedImage(templateMaskedImage,
                                                                                     scienceMaskedImage,
                                                                                     policy)
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
