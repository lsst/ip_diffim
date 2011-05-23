import eups
import sys
import os
import optparse

import lsst.daf.base as dafBase
import lsst.afw.image as afwImage

from lsst.pex.logging import Trace
from lsst.pex.logging import Log

from lsst.ip.diffim import subtractExposures, makeDefaultPolicy
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
    defOutputPath   = 'diffExposure'
    defVerbosity    = 0
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
    parser.add_option('-b', '--bg', action='store_true', default=True,
                      help='subtract backgrounds')
    parser.add_option('-f', '--fwhm', type=float,
                      help='Psf Fwhm (pixel)')

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

    fwhm = defFwhm
    if options.fwhm:
        print 'Fwhm =', options.fwhm
        fwhm = options.fwhm

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
    policy           = makeDefaultPolicy(mergePolicyPath = mergePolicyPath, fwhm=fwhm)

    if bgSub:
        diffimTools.backgroundSubtract(policy.getPolicy("afwBackgroundPolicy"),
                                       [templateExposure.getMaskedImage(),
                                        scienceExposure.getMaskedImage()])

    results = subtractExposures(templateExposure,
                                scienceExposure,
                                policy,
                                display = display)
    differenceExposure = results[0]
    differenceExposure.writeFits(outputPath)

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
