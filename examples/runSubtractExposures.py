import eups
import sys
import os
import optparse
import re

import lsst.daf.base as dafBase
import lsst.afw.image as afwImage

from lsst.pex.logging import Trace
from lsst.pex.logging import Log

from lsst.ip.diffim import psfMatch, makeDefaultPolicy, modifyForImagePsfMatch, makeSdqaRatingVector
import lsst.ip.diffim.diffimTools as diffimTools

def main():
    defDataDir   = eups.productDir('afwdata')
    if defDataDir == None:
        print 'Error: afwdata not set up'
        sys.exit(1)
    imageProcDir = eups.productDir('ip_diffim')
    if imageProcDir == None:
        print 'Error: could not set up ip_diffim'
        sys.exit(1)

    defSciencePath  = os.path.join(defDataDir, "DC3a-Sim", "sci", "v26-e0", "v26-e0-c011-a10.sci")
    defTemplatePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v5-e0", "v5-e0-c011-a10.sci")

    mergePolicyPath = None
    defOutputPath   = 'diffExposure.fits'
    defVerbosity    = 5
    defFwhm         = 3.5
    
    usage = """usage: %%prog [options] [scienceExposure [templateExposure [outputExposure]]]]

Notes:
- image arguments are paths to Expoure fits files
- image arguments must NOT include the final _img.fits
- the result is science image - template image
- the template exposure is convolved, the science exposure is not
- default scienceExposure=%s
- default templateExposure=%s
- default outputExposure=%s 
- default --policy=%s
""" % (defSciencePath, defTemplatePath, defOutputPath, mergePolicyPath)
    
    parser = optparse.OptionParser(usage)
    parser.add_option('-p', '--policy', default=mergePolicyPath, help='policy file to merge with defaults')
    parser.add_option('-v', '--verbosity', type=int, default=defVerbosity,
                      help='verbosity of Trace messages')
    parser.add_option('-d', '--display', action='store_true', default=False,
                      help='display the images')
    parser.add_option('-b', '--bg', action='store_true', default=False,
                      help='subtract backgrounds using afw')

    parser.add_option('-t', '--fwhmTemplate', type=float,
                      help='Template Psf Fwhm (pixel)')
    parser.add_option('-i', '--fwhmImage', type=float,
                      help='Image Psf Fwhm (pixel)')

    (options, args) = parser.parse_args()
    
    def getArg(ind, defValue):
        if ind < len(args):
            return args[ind]
        return defValue
    
    sciencePath     = getArg(0, defSciencePath)
    templatePath    = getArg(1, defTemplatePath)
    outputPath      = getArg(2, defOutputPath)
    mergePolicyPath = options.policy
    
    print 'Science exposure: ', sciencePath
    print 'Template exposure:', templatePath
    print 'Output exposure:  ', outputPath
    print 'Policy file:      ', mergePolicyPath

    if options.fwhmTemplate:
        fwhmTemplate = options.fwhmTemplate
        print 'Fwhm Template =', fwhmTemplate
    else:
        fwhmTemplate = defFwhm

    if options.fwhmImage:
        fwhmImage = options.fwhmImage
        print 'Fwhm Image =', fwhmImage
    else:
        fwhmImage = defFwhm

    display = False
    if options.display:
        print 'Display =', options.display
        display = True

    bgSub = False
    if options.bg:
        print 'Background subtract =', options.bg
        bgSub = True

    if options.verbosity > 0:
        print 'Verbosity =', options.verbosity
        Trace.setVerbosity('lsst.ip.diffim', options.verbosity)

    ####
        
    templateExposure = afwImage.ExposureF(templatePath)
    scienceExposure  = afwImage.ExposureF(sciencePath)
    policy           = makeDefaultPolicy(mergePolicy = mergePolicyPath)
    policy           = modifyForImagePsfMatch(policy, fwhmTemplate, fwhmImage)
    print policy
    
    if bgSub:
        diffimTools.backgroundSubtract(policy.getPolicy("afwBackgroundPolicy"),
                                       [templateExposure.getMaskedImage(),
                                        scienceExposure.getMaskedImage()])
    else:
        policy.set('fitForBackground', True)
        if policy.get('fitForBackground') == False:
            print 'NOTE: no background subtraction at all is requested'

        
    psfmatch = psfMatch.ImagePsfMatch(policy)
    results  = psfmatch.subtractExposures(templateExposure, scienceExposure)

    differenceExposure = results[0]
    differenceExposure.writeFits(outputPath)

    psfMatchingKernel = results[1]
    backgroundModel   = results[2]
    kernelCellSet     = results[3]
    makeSdqaRatingVector(kernelCellSet, psfMatchingKernel, backgroundModel)

    diffimTools.writeKernelCellSet(kernelCellSet, psfMatchingKernel, backgroundModel,
                                   re.sub('.fits', '', outputPath))

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
