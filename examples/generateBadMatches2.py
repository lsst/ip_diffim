from builtins import range
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import numpy as num
import lsst.log.utils as logUtils
import lsst.pex.config as pexConfig
import lsst.afw.display.ds9 as ds9
verbosity = 4
logUtils.traceSetAt("ip.diffim", verbosity)

imSize = 2**7
kSize = 2**5 + 1
rdm = afwMath.Random(afwMath.Random.MT19937, 10101)
scaling = 10000
doNorm = True
# gScale         = 1.0 # FFT works!
gScale = 5.0  # FFT fails!?

doAddNoise = False
writeFits = False


def makeTest1(doAddNoise):
    gaussian1 = afwMath.GaussianFunction2D(1. * gScale, 1. * gScale, 0.)
    kernel1 = afwMath.AnalyticKernel(imSize, imSize, gaussian1)
    image1 = afwImage.ImageD(kernel1.getDimensions())
    kernel1.computeImage(image1, doNorm)
    image1 *= scaling  # total counts = scaling
    image1 = image1.convertF()
    mask1 = afwImage.Mask(kernel1.getDimensions())
    var1 = afwImage.ImageF(image1, True)
    mi1 = afwImage.MaskedImageF(image1, mask1, var1)
    if doAddNoise:
        addNoise(mi1)

    gaussian2 = afwMath.GaussianFunction2D(2. * gScale, 1.5 * gScale, 0.5 * num.pi)
    kernel2 = afwMath.AnalyticKernel(imSize, imSize, gaussian2)
    image2 = afwImage.ImageD(kernel2.getDimensions())
    kernel2.computeImage(image2, doNorm)
    image2 *= scaling  # total counts = scaling
    image2 = image2.convertF()
    mask2 = afwImage.Mask(kernel2.getDimensions())
    var2 = afwImage.ImageF(image2, True)
    mi2 = afwImage.MaskedImageF(image2, mask2, var2)
    if doAddNoise:
        addNoise(mi2)
    return mi1, mi2


def makeTest2(doAddNoise, shiftX=int(2.0 * gScale), shiftY=int(1.0 * gScale)):
    gaussian1 = afwMath.GaussianFunction2D(1. * gScale, 1. * gScale, 0.)
    kernel1 = afwMath.AnalyticKernel(imSize, imSize, gaussian1)
    image1 = afwImage.ImageD(kernel1.getDimensions())
    kernel1.computeImage(image1, doNorm)
    image1 = image1.convertF()
    ####
    boxA = afwGeom.Box2I(afwGeom.PointI(0, 0),
                         afwGeom.ExtentI(imSize-shiftX, imSize - shiftY))
    boxB = afwGeom.Box2I(afwGeom.PointI(shiftX, shiftY),
                         afwGeom.ExtentI(imSize-shiftX, imSize - shiftY))
    subregA = afwImage.ImageF(image1, boxA, afwImage.PARENT)
    subregB = afwImage.ImageF(image1, boxB, afwImage.PARENT, True)
    subregA += subregB
    # this messes up the total counts so rescale
    counts = afwMath.makeStatistics(image1, afwMath.SUM).getValue()
    image1 /= counts
    image1 *= scaling
    ###
    mask1 = afwImage.Mask(kernel1.getDimensions())
    var1 = afwImage.ImageF(image1, True)
    mi1 = afwImage.MaskedImageF(image1, mask1, var1)
    if doAddNoise:
        addNoise(mi1)

    gaussian2 = afwMath.GaussianFunction2D(2. * gScale, 1.5 * gScale, 0.5 * num.pi)
    kernel2 = afwMath.AnalyticKernel(imSize, imSize, gaussian2)
    image2 = afwImage.ImageD(kernel2.getDimensions())
    kernel2.computeImage(image2, doNorm)
    image2 *= scaling  # total counts = scaling
    image2 = image2.convertF()
    mask2 = afwImage.Mask(kernel2.getDimensions())
    var2 = afwImage.ImageF(image2, True)
    mi2 = afwImage.MaskedImageF(image2, mask2, var2)
    if doAddNoise:
        addNoise(mi2)

    return mi1, mi2


def makeTest3(doAddNoise):
    gaussian1 = afwMath.GaussianFunction2D(1. * gScale, 1. * gScale, 0.)
    kernel1 = afwMath.AnalyticKernel(imSize, imSize, gaussian1)
    image1 = afwImage.ImageD(kernel1.getDimensions())
    kernel1.computeImage(image1, doNorm)
    image1 *= scaling  # total counts = scaling
    image1 = image1.convertF()
    mask1 = afwImage.Mask(kernel1.getDimensions())
    var1 = afwImage.ImageF(image1, True)
    mi1 = afwImage.MaskedImageF(image1, mask1, var1)

    gaussian2 = afwMath.GaussianFunction2D(2. * gScale, 1.5 * gScale, 0.5 * num.pi)
    kernel2 = afwMath.AnalyticKernel(imSize, imSize, gaussian2)
    image2 = afwImage.ImageD(kernel2.getDimensions())
    kernel2.computeImage(image2, doNorm)
    image2 *= scaling  # total counts = scaling
    image2 = image2.convertF()
    mask2 = afwImage.Mask(kernel2.getDimensions())
    var2 = afwImage.ImageF(image2, True)
    mi2 = afwImage.MaskedImageF(image2, mask2, var2)

    image3 = afwImage.ImageF(image1, True)
    for y in range(imSize):
        for x in range(imSize//2):
            image3.set(x, y, image2.get(x, y))
    counts = afwMath.makeStatistics(image3, afwMath.SUM).getValue()
    image3 /= counts
    image3 *= scaling

    mask3 = afwImage.Mask(image3.getDimensions())
    var3 = afwImage.ImageF(image3, True)
    mi3 = afwImage.MaskedImageF(image3, mask3, var3)

    if doAddNoise:
        addNoise(mi1)
        addNoise(mi2)
        addNoise(mi3)

    return mi1, mi2, mi3


def fft(im1, im2, fftSize):
    arr1 = im1.getArray()
    arr2 = im2.getArray()

    fft1 = num.fft.rfft2(arr1)
    fft2 = num.fft.rfft2(arr2)
    rat = fft2 / fft1

    kfft = num.fft.irfft2(rat, s=fftSize)
    kfft = num.fft.fftshift(kfft)
    kim = afwImage.ImageF(fftSize)
    kim.getArray()[:] = kfft

    ds9.mtv(kim, frame=5)


# If we don't add noise, the edges of the Gaussian images go to zero,
# and that boundary causes artificial artefacts in the kernels
def addNoise(mi):
    sfac = 1.0
    img = mi.getImage()
    rdmImage = img.Factory(img.getDimensions())
    afwMath.randomGaussianImage(rdmImage, rdm)
    rdmImage *= sfac
    img += rdmImage

    # and don't forget to add to the variance
    var = mi.getVariance()
    var += sfac


if __name__ == '__main__':

    configAL = ipDiffim.ImagePsfMatchTask.ConfigClass()
    configAL.kernel.name = "AL"
    subconfigAL = configAL.kernel.active

    configDF = ipDiffim.ImagePsfMatchTask.ConfigClass()
    configDF.kernel.name = "DF"
    subconfigDF = configDF.kernel.active

    subconfigAL.fitForBackground = False
    subconfigDF.fitForBackground = False

    # Super-important for these faked-up kernels...
    subconfigAL.constantVarianceWeighting = True
    subconfigDF.constantVarianceWeighting = True

    subconfigAL.kernelSize = kSize
    subconfigDF.kernelSize = kSize

    alardSigGauss = subconfigAL.alardSigGauss
    subconfigAL.alardSigGauss = [x * gScale for x in alardSigGauss]

    fnum = 1

    for switch in ['A', 'B', 'C']:
        if switch == 'A':
            # Default Alard Lupton
            config = subconfigAL
        elif switch == 'B':
            # Add more AL bases (typically 4 3 2)
            config = subconfigAL
            config.alardDegGauss = (8, 6, 4)
        elif switch == 'C':
            # Delta function
            config = subconfigDF
            config.useRegularization = False

        kList = ipDiffim.makeKernelBasisList(config)

        policy = pexConfig.makePolicy(config)
        bskv = ipDiffim.BuildSingleKernelVisitorF(kList, policy)

        # TEST 1
        tmi, smi = makeTest1(doAddNoise)
        kc = ipDiffim.makeKernelCandidate(0.0, 0.0, tmi, smi, policy)
        bskv.processCandidate(kc)

        kernel = kc.getKernel(ipDiffim.KernelCandidateF.ORIG)
        kimage = afwImage.ImageD(kernel.getDimensions())
        kernel.computeImage(kimage, False)
        diffim = kc.getDifferenceImage(ipDiffim.KernelCandidateF.ORIG)

        ds9.mtv(tmi, frame=fnum)
        fnum += 1
        ds9.mtv(smi, frame=fnum)
        fnum += 1
        ds9.mtv(kimage, frame=fnum)
        fnum += 1
        ds9.mtv(diffim, frame=fnum)
        fnum += 1

        if writeFits:
            tmi.writeFits("template1.fits")
            smi.writeFits("science1.fits")
            kimage.writeFits("kernel1.fits")
            diffim.writeFits("diffim1.fits")

        # TEST 2
        tmi, smi = makeTest2(doAddNoise)
        kc = ipDiffim.makeKernelCandidate(0.0, 0.0, tmi, smi, policy)
        bskv.processCandidate(kc)

        kernel = kc.getKernel(ipDiffim.KernelCandidateF.ORIG)
        kimage = afwImage.ImageD(kernel.getDimensions())
        kernel.computeImage(kimage, False)
        diffim = kc.getDifferenceImage(ipDiffim.KernelCandidateF.ORIG)

        ds9.mtv(tmi, frame=fnum)
        fnum += 1
        ds9.mtv(smi, frame=fnum)
        fnum += 1
        ds9.mtv(kimage, frame=fnum)
        fnum += 1
        ds9.mtv(diffim, frame=fnum)
        fnum += 1

        if writeFits:
            tmi.writeFits("template2.fits")
            smi.writeFits("science2.fits")
            kimage.writeFits("kernel2.fits")
            diffim.writeFits("diffim2.fits")

        # TEST 3
        smi1, smi2, tmi = makeTest3(doAddNoise)
        kc1 = ipDiffim.makeKernelCandidate(0.0, 0.0, tmi, smi1, policy)
        kc2 = ipDiffim.makeKernelCandidate(0.0, 0.0, tmi, smi2, policy)
        bskv.processCandidate(kc1)
        bskv.processCandidate(kc2)

        kernel1 = kc1.getKernel(ipDiffim.KernelCandidateF.ORIG)
        kimage1 = afwImage.ImageD(kernel1.getDimensions())
        kernel1.computeImage(kimage1, False)
        diffim1 = kc1.getDifferenceImage(ipDiffim.KernelCandidateF.ORIG)

        kernel2 = kc2.getKernel(ipDiffim.KernelCandidateF.ORIG)
        kimage2 = afwImage.ImageD(kernel2.getDimensions())
        kernel2.computeImage(kimage2, False)
        diffim2 = kc2.getDifferenceImage(ipDiffim.KernelCandidateF.ORIG)

        ds9.mtv(tmi, frame=fnum)
        fnum += 1
        ds9.mtv(smi1, frame=fnum)
        fnum += 1
        ds9.mtv(kimage1, frame=fnum)
        fnum += 1
        ds9.mtv(diffim1, frame=fnum)
        fnum += 1

        ds9.mtv(tmi, frame=fnum)
        fnum += 1
        ds9.mtv(smi2, frame=fnum)
        fnum += 1
        ds9.mtv(kimage2, frame=fnum)
        fnum += 1
        ds9.mtv(diffim2, frame=fnum)
        fnum += 1

        if writeFits:
            tmi.writeFits("template3a.fits")
            smi1.writeFits("science3a.fits")
            kimage1.writeFits("kernel3a.fits")
            diffim1.writeFits("diffim3a.fits")

            tmi.writeFits("template3b.fits")
            smi2.writeFits("science3b.fits")
            kimage2.writeFits("kernel3b.fits")
            diffim2.writeFits("diffim3b.fits")
