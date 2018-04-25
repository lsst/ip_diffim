import sys
import optparse

import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9
import lsst.log.utils as logUtils

import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools


def main():
    defSciencePath = None
    defTemplatePath = None
    defOutputPath = 'matchedImage.fits'
    defVerbosity = 5
    defFwhm = 3.5

    usage = """usage: %%prog [options] [scienceImage [templateImage [outputImage]]]]

Notes:
- image arguments are paths to MaskedImage fits files
- image arguments must NOT include the final _img.fits
- the result is science image matched template image
- the template image is convolved, the science image is not
- default scienceMaskedImage=%s
- default templateMaskedImage=%s
- default outputImage=%s
""" % (defSciencePath, defTemplatePath, defOutputPath)

    parser = optparse.OptionParser(usage)
    parser.add_option('-v', '--verbosity', type=int, default=defVerbosity,
                      help='verbosity of Trace messages')
    parser.add_option('-d', '--display', action='store_true', default=False,
                      help='display the images')
    parser.add_option('-b', '--bg', action='store_true', default=False,
                      help='subtract backgrounds')
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

    print('Science image: ', sciencePath)
    print('Template image:', templatePath)
    print('Output image:  ', outputPath)

    fwhmS = defFwhm
    if options.fwhmS:
        print('FwhmS =', options.fwhmS)
        fwhmS = options.fwhmS

    fwhmT = defFwhm
    if options.fwhmT:
        print('Fwhmt =', options.fwhmT)
        fwhmT = options.fwhmT

    display = False
    if options.display:
        print('Display =', options.display)
        display = True

    bgSub = False
    if options.bg:
        print('Background subtract =', options.bg)
        bgSub = True

    if options.verbosity > 0:
        print('Verbosity =', options.verbosity)
        logUtils.traceSetAt("ip.diffim", options.verbosity)

    ####

    templateMaskedImage = afwImage.MaskedImageF(templatePath)
    scienceMaskedImage = afwImage.MaskedImageF(sciencePath)

    config = ipDiffim.ImagePsfMatch.ConfigClass()
    config.kernel.name = "AL"
    subconfig = config.kernel.active

    if bgSub:
        diffimTools.backgroundSubtract(subconfig.afwBackgroundConfig,
                                       [templateMaskedImage, scienceMaskedImage])
    else:
        if not subconfig.fitForBackground:
            print('NOTE: no background subtraction at all is requested')

    psfmatch = ipDiffim.ImagePsfMatch(subconfig)
    results = psfmatch.matchMaskedImages(templateMaskedImage, scienceMaskedImage,
                                         templateFwhmPix=fwhmT, scienceFwhmPix=fwhmS)

    matchMaskedImage = results[0]
    matchMaskedImage.writeFits(outputPath)

    if False:
        spatialKernel = results[1]
        print(spatialKernel.getSpatialParameters())

    if display:
        ds9.mtv(differenceMaskedImage)


def run():
    main()


if __name__ == '__main__':
    run()
