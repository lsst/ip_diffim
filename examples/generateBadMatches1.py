# This file is part of ip_diffim.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import numpy as num
import lsst.log.utils as logUtils
import lsst.pex.config as pexConfig
import lsst.afw.display as afwDisplay
verbosity = 4
logUtils.traceSetAt("ip.diffim", verbosity)

imSize = 75
rdm = afwMath.Random(afwMath.Random.MT19937, 10101)
writeFits = False


def makeTest1(doAddNoise):
    gaussian1 = afwMath.GaussianFunction2D(1., 1., 0.)
    kernel1 = afwMath.AnalyticKernel(imSize, imSize, gaussian1)
    image1 = afwImage.ImageD(kernel1.getDimensions())
    kernel1.computeImage(image1, False)
    image1 *= 10000
    image1 = image1.convertF()
    mask1 = afwImage.Mask(kernel1.getDimensions())
    var1 = afwImage.ImageF(image1, True)
    mi1 = afwImage.MaskedImageF(image1, mask1, var1)
    if doAddNoise:
        addNoise(mi1)

    gaussian2 = afwMath.GaussianFunction2D(2., 1.5, 0.5 * num.pi)
    kernel2 = afwMath.AnalyticKernel(imSize, imSize, gaussian2)
    image2 = afwImage.ImageD(kernel2.getDimensions())
    kernel2.computeImage(image2, False)
    image2 *= 10000
    image2 = image2.convertF()
    mask2 = afwImage.Mask(kernel2.getDimensions())
    var2 = afwImage.ImageF(image2, True)
    mi2 = afwImage.MaskedImageF(image2, mask2, var2)
    if doAddNoise:
        addNoise(mi2)
    return mi1, mi2


def makeTest2(doAddNoise, shiftX=5, shiftY=3):
    gaussian1 = afwMath.GaussianFunction2D(1., 1., 0.)
    kernel1 = afwMath.AnalyticKernel(imSize, imSize, gaussian1)
    image1 = afwImage.ImageD(kernel1.getDimensions())
    kernel1.computeImage(image1, False)
    image1 *= 10000
    image1 = image1.convertF()
    ####
    boxA = afwGeom.Box2I(afwGeom.PointI(imSize//2, imSize//2),
                         afwGeom.ExtentI(imSize//2, imSize//2))
    boxB = afwGeom.Box2I(afwGeom.PointI(imSize//2 - shiftX, imSize//2 - shiftY),
                         afwGeom.ExtentI(imSize//2, imSize//2))
    subregA = afwImage.ImageF(image1, boxA, afwImage.PARENT)
    subregB = afwImage.ImageF(image1, boxB, afwImage.PARENT, True)
    subregA += subregB
    ###
    mask1 = afwImage.Mask(kernel1.getDimensions())
    var1 = afwImage.ImageF(image1, True)
    mi1 = afwImage.MaskedImageF(image1, mask1, var1)
    if doAddNoise:
        addNoise(mi1)

    gaussian2 = afwMath.GaussianFunction2D(2., 1.5, 0.5 * num.pi)
    kernel2 = afwMath.AnalyticKernel(imSize, imSize, gaussian2)
    image2 = afwImage.ImageD(kernel2.getDimensions())
    kernel2.computeImage(image2, False)
    image2 *= 10000
    image2 = image2.convertF()
    mask2 = afwImage.Mask(kernel2.getDimensions())
    var2 = afwImage.ImageF(image2, True)
    mi2 = afwImage.MaskedImageF(image2, mask2, var2)
    if doAddNoise:
        addNoise(mi2)

    return mi1, mi2


def fft(im1, im2, fftSize):
    arr1 = im1.getArray()
    arr2 = im2.getArray()

    fft1 = num.fft.rfft2(arr1)
    fft2 = num.fft.rfft2(arr2)
    rat = fft2 / fft1

    kfft = num.fft.irfft2(rat, s=fftSize)
    kim = afwImage.ImageF(fftSize)
    kim.getArray()[:] = kfft

    afwDisplay.Display(frame=5).mtv(kim)


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
    doAddNoise = True

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

    fnum = 1

    for switch in ['A', 'B', 'C']:
        if switch == 'A':
            # AL
            config = subconfigAL
        elif switch == 'B':
            # AL with ~320 bases
            config = subconfigAL
            config.alardDegGauss = (15, 10, 5)
        elif switch == 'C':
            config = subconfigDF
            config.useRegularization = False

        kList = ipDiffim.makeKernelBasisList(config)

        ps = pexConfig.makePropertySet(config)
        bskv = ipDiffim.BuildSingleKernelVisitorF(kList, ps)

        # TEST 1
        tmi, smi = makeTest1(doAddNoise)
        kc = ipDiffim.makeKernelCandidate(0.0, 0.0, tmi, smi, ps)
        bskv.processCandidate(kc)

        kernel = kc.getKernel(ipDiffim.KernelCandidateF.ORIG)
        kimage = afwImage.ImageD(kernel.getDimensions())
        kernel.computeImage(kimage, False)
        diffim = kc.getDifferenceImage(ipDiffim.KernelCandidateF.ORIG)

        afwDisplay.Display(frame=fnum).mtv(tmi)
        fnum += 1
        afwDisplay.Display(frame=fnum).mtv(smi)
        fnum += 1
        afwDisplay.Display(frame=fnum).mtv(kimage)
        fnum += 1
        afwDisplay.Display(frame=fnum).mtv(diffim)
        fnum += 1

        if writeFits:
            tmi.writeFits("template1.fits")
            smi.writeFits("science1.fits")
            kimage.writeFits("kernel1.fits")
            diffim.writeFits("diffim1.fits")

        # TEST 2
        tmi, smi = makeTest2(doAddNoise, shiftX=2, shiftY=2)
        kc = ipDiffim.makeKernelCandidate(0.0, 0.0, tmi, smi, ps)
        bskv.processCandidate(kc)

        kernel = kc.getKernel(ipDiffim.KernelCandidateF.ORIG)
        kimage = afwImage.ImageD(kernel.getDimensions())
        kernel.computeImage(kimage, False)
        diffim = kc.getDifferenceImage(ipDiffim.KernelCandidateF.ORIG)

        afwDisplay.Display(frame=fnum).mtv(tmi)
        fnum += 1
        afwDisplay.Display(frame=fnum).mtv(smi)
        fnum += 1
        afwDisplay.Display(frame=fnum).mtv(kimage)
        fnum += 1
        afwDisplay.Display(frame=fnum).mtv(diffim)
        fnum += 1

        if writeFits:
            tmi.writeFits("template2.fits")
            smi.writeFits("science2.fits")
            kimage.writeFits("kernel2.fits")
            diffim.writeFits("diffim2.fits")
