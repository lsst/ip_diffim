import eups
import sys, os, optparse
import numpy

from subtractMaskedImage import subtractMaskedImage

def subtractExposure(templateExposure, scienceExposure, policy):
    # Make sure they end up the same dimensions on the sky
    templateWcs      = templateExposure.getWcs() 
    scienceWcs       = scienceExposure.getWcs()

    templateOrign    = templateWcs.xyToRaDec(0,0)
    scienceOrign     = scienceWcs.xyToRaDec(0,0)
    assert(templateOrigin[0] == scienceOrigin[0])
    assert(templateOrigin[1] == scienceOrigin[1])

    templateLimit    = templateWcs.xyToRaDec(templateExposure.getHeight()
                                             templateExposure.getWidth())
    scienceLimit     = scienceWcs.xyToRaDec(scienceExposure.getHeight()
                                            scienceExposure.getWidth())
    # Within some tolerance; do we have sky distance methods?
    assert(templateLimit[0]  == scienceLimit[0])
    assert(templateLimit[1]  == scienceLimit[1])

    # Make sure they end up the EXACT same dimensions in pixels
    # This is non-negotiable
    assert (templateExposure.getDimensions() == scienceExposure.getDimensions())

    # Subtract their MaskedImages
    differenceMaskedImage, spatialKernel, backgroundModel, SdqaList = subtractMaskedImage(templateExposure.getMaskedImage(),
                                                                                          scienceExposure.getMaskedImage(),
                                                                                          policy)
    # Note : we assume that the Template is warped to the science image's WCS
    #      : meaning that the scienceWcs is the correct one to store in the diffim
    differenceExposure = afwImage.ExposureF(differenceMaskedImage, scienceWcs)

    # What kind of metadata do we add here to Exposure?
    
    return differenceExposure, spatialKernel, backgroundModel, SdqaList


def main():
    defDataDir   = eups.productDir('afwdata') or ''
    imageProcDir = eups.productDir('ip_diffim')
    if imageProcDir == None:
        print 'Error: could not set up ip_diffim'
        sys.exit(1)

    defSciencePath  = os.path.join(defDataDir, 'CFHT', 'D4', 'cal-53535-i-797722_1')
    defTemplatePath = os.path.join(defDataDir, 'CFHT', 'D4', 'cal-53535-i-797722_1_tmpl')
    defPolicyPath   = os.path.join(imageProcDir, 'pipeline', 'ImageSubtractStageDictionary.paf')
    defOutputPath   = 'diffExposure'
    defVerbosity    = 0
    
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
    
    print 'Science exposure: ', sciencePath
    print 'Template exposure:', templatePath
    print 'Output exposure:  ', outputPath
    print 'Policy file:      ', policyPath
    
    templateExposure = afwImage.ExposureF(templatePath)
    scienceExposure  = afwImage.ExposureF(sciencePath)
    policy           = Policy.createPolicy(policyPath)

    if options.debugIO:
        policy.set('debugIO', True)

    if options.verbosity > 0:
        print 'Verbosity =', options.verbosity
        Trace.setVerbosity('lsst.ip.diffim', options.verbosity)

    differenceExposure, spatialKernel, backgroundModel, SdqaList = subtractExpoure(templateExposure,
                                                                                   scienceExposure,
                                                                                   policy)
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
