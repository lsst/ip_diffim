import eups
import sys
import os
import optparse

import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9

from lsst.pex.logging import Trace
from lsst.pex.logging import Log

from lsst.ip.diffim import subtractMaskedImages, makeDefaultPolicy
import lsst.ip.diffim.diffimTools as diffimTools

def main():
    defDataDir = eups.productDir('afwdata') 
    if defDataDir == None:
        print 'Error: afwdata not set up'
        sys.exit(1)
    imageProcDir = eups.productDir('ip_diffim')
    if imageProcDir == None:
        print 'Error: could not set up ip_diffim'
        sys.exit(1)

    defSciencePath  = None
    defTemplatePath = None
    mergePolicyPath = None
    defOutputPath   = 'diffImage.fits'
    defVerbosity    = 0
    defFwhm         = 3.5
    
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
""" % (defSciencePath, defTemplatePath, defOutputPath, mergePolicyPath)
    
    parser = optparse.OptionParser(usage)
    parser.add_option('-p', '--policy', default=mergePolicyPath, help='policy file')
    parser.add_option('-v', '--verbosity', type=int, default=defVerbosity,
                      help='verbosity of Trace messages')
    parser.add_option('-d', '--display', action='store_true', default=False,
                      help='display the images')
    parser.add_option('-b', '--bg', action='store_true', default=False,
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

    if sciencePath == None or templatePath == None:
        parser.print_help()
        sys.exit(1)
    
    print 'Science image: ', sciencePath
    print 'Template image:', templatePath
    print 'Output image:  ', outputPath
    print 'Policy file:   ', mergePolicyPath

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
        
    templateMaskedImage = afwImage.MaskedImageF(templatePath)
    scienceMaskedImage  = afwImage.MaskedImageF(sciencePath)
    policy              = makeDefaultPolicy(mergePolicyPath = mergePolicyPath, fwhm=fwhm)
    
    if bgSub:
        diffimTools.backgroundSubtract(policy.getPolicy("afwBackgroundPolicy"),
                                       [templateMaskedImage, scienceMaskedImage])
    else:
        if policy.get('fitForBackground') == False:
            print 'NOTE: no background subtraction at all is requested'

    psfmatch = PsfMatch.ImagePsfMatch(policy)
    results  = psfmatch.subtractMaskedImages(templateMaskedImage, scienceMaskedImage)

    differenceMaskedImage = results[0]
    differenceMaskedImage.writeFits(outputPath)

    spatialKernel = results[1]
    print spatialKernel.getSpatialParameters()
    
    if display:
        ds9.mtv(differenceMaskedImage)

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
