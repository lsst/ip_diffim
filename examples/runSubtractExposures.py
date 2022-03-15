import lsst.utils
import sys
import os
import optparse
import re
import numpy as np

import lsst.afw.image as afwImage
import lsst.utils.logging as logUtils
from lsst.meas.algorithms import SingleGaussianPsf

import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools


def main():
    defDataDir = lsst.utils.getPackageDir('afwdata')

    defSciencePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v26-e0", "v26-e0-c011-a10.sci.fits")
    defTemplatePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v5-e0", "v5-e0-c011-a10.sci.fits")

    defOutputPath = 'diffExposure.fits'
    defVerbosity = 0
    defFwhm = 3.5
    sigma2fwhm = 2. * np.sqrt(2. * np.log(2.))

    usage = """usage: %%prog [options] [scienceExposure [templateExposure [outputExposure]]]]

Notes:
- image arguments are paths to Expoure fits files
- image arguments must NOT include the final _img.fits
- the result is science image - template image
- the template exposure is convolved, the science exposure is not
- default scienceExposure=%s
- default templateExposure=%s
- default outputExposure=%s
""" % (defSciencePath, defTemplatePath, defOutputPath)

    parser = optparse.OptionParser(usage)
    parser.add_option('-v', '--verbosity', type=int, default=defVerbosity,
                      help='verbosity of Trace messages')
    parser.add_option('-d', '--display', action='store_true', default=False,
                      help='display the images')
    parser.add_option('-b', '--bg', action='store_true', default=False,
                      help='subtract backgrounds using afw')
    parser.add_option('--fwhmS', type=float,
                      help='Science Image Psf Fwhm (pixel)')
    parser.add_option('--fwhmT', type=float,
                      help='Template Image Psf Fwhm (pixel)')

    (options, args) = parser.parse_args()

    def getArg(ind, defValue):
        if ind < len(args):
            return args[ind]
        return defValue

    sciencePath = getArg(0, defSciencePath)
    templatePath = getArg(1, defTemplatePath)
    outputPath = getArg(2, defOutputPath)

    if sciencePath is None or templatePath is None:
        parser.print_help()
        sys.exit(1)

    print('Science exposure: ', sciencePath)
    print('Template exposure:', templatePath)
    print('Output exposure:  ', outputPath)

    templateExposure = afwImage.ExposureF(templatePath)
    scienceExposure = afwImage.ExposureF(sciencePath)

    config = ipDiffim.ImagePsfMatchTask.ConfigClass()
    config.kernel.name = "AL"
    subconfig = config.kernel.active

    fwhmS = defFwhm
    if options.fwhmS:
        if scienceExposure.hasPsf():
            sciPsf = scienceExposure.getPsf()
            fwhm = sciPsf.computeShape(sciPsf.getAveragePosition()).getDeterminantRadius() * sigma2fwhm
            print('NOTE: Embedded Psf has FwhmS =', fwhm)
        print('USING: FwhmS =', options.fwhmS)
        fwhmS = options.fwhmS

    fwhmT = defFwhm
    if options.fwhmT:
        if templateExposure.hasPsf():
            templPsf = templateExposure.getPsf()
            fwhm = templPsf.computeShape(templPsf.getAveragePosition()).getDeterminantRadius() * sigma2fwhm
            print('NOTE: Embedded Psf has FwhmT =', fwhm)
        print('USING: FwhmT =', options.fwhmT)
        fwhmT = options.fwhmT

    bgSub = False
    if options.bg:
        print('Background subtract =', options.bg)
        bgSub = True

    if options.verbosity > 0:
        print('Verbosity =', options.verbosity)
        logUtils.trace_set_at("lsst.ip.diffim", options.verbosity)

    ####

    if bgSub:
        diffimTools.backgroundSubtract(subconfig.afwBackgroundConfig,
                                       [templateExposure.getMaskedImage(),
                                        scienceExposure.getMaskedImage()])
    else:
        if not subconfig.fitForBackground:
            print('NOTE: no background subtraction at all is requested')

    # ImagePsfMatchTask requires a candidateList or that the exposoure have a PSF for detection
    if not scienceExposure.hasPsf():
        sigmaS = fwhmS/sigma2fwhm
        kSize = int(6 * sigmaS)  # minimum kernel size for FWHM
        oddKSize = kSize + 1 if kSize % 2 == 0 else kSize
        scienceExposure.setPsf(SingleGaussianPsf(oddKSize, oddKSize, sigmaS))

    psfmatch = ipDiffim.ImagePsfMatchTask(config)
    results = psfmatch.subtractExposures(templateExposure, scienceExposure,
                                         templateFwhmPix=fwhmT, scienceFwhmPix=fwhmS)

    differenceExposure = results.subtractedExposure
    differenceExposure.writeFits(outputPath)

    if False:
        psfMatchingKernel = results.psfMatchingKernel
        backgroundModel = results.backgroundModel
        kernelCellSet = results.kernelCellSet

        diffimTools.writeKernelCellSet(kernelCellSet, psfMatchingKernel, backgroundModel,
                                       re.sub('.fits', '', outputPath))


def run():
    main()


if __name__ == '__main__':
    run()
