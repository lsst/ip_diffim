# This file is part of ip_diffim.
#
# LSST Data Management System
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
# See COPYRIGHT file at the top of the source tree.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#

import numpy as np
import unittest

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.daf.base as dafBase
import lsst.geom as geom
from lsst.meas.algorithms.testUtils import plantSources
import lsst.utils.tests

from lsst.ip.diffim.imageDecorrelation import DecorrelateALKernelTask, DecorrelateALKernelConfig
from lsst.ip.diffim.imagePsfMatch import ImagePsfMatchTask, ImagePsfMatchConfig
from lsst.ip.diffim.zogy import ZogyTask, ZogyConfig


class ImageDifferenceTestCase(lsst.utils.tests.TestCase):
    """A test case for comparing image differencing algorithms.

    Attributes
    ----------
    bbox : `lsst.afw.geom.Box2I`
        Bounding box of the test model.
    bufferSize : `int`
        Distance from the inner edge of the bounding box
        to avoid placing test sources in the model images.
    nRandIter : int
        Description
    rng : TYPE
        Description
    statsCtrl : TYPE
        Description
    """

    def setUp(self):
        """Define the filter, DCR parameters, and the bounding box for the tests.
        """
        self.nRandIter = 10  # Number of iterations to repeat each test with random numbers.
        self.bufferSize = 5
        xSize = 250
        ySize = 260
        x0 = 12345
        y0 = 67890
        self.bbox = geom.Box2I(geom.Point2I(x0, y0), geom.Extent2I(xSize, ySize))
        config = ImagePsfMatchConfig()
        self.matcher = ImagePsfMatchTask(config=config)

        self.statsCtrl = afwMath.StatisticsControl()
        self.statsCtrl.setNumSigmaClip(3.)
        self.statsCtrl.setNumIter(3)
        self.statsCtrl.setAndMask(afwImage.Mask
                                  .getPlaneBitMask(["INTRP", "EDGE", "SAT", "CR",
                                                    "BAD", "NO_DATA"]))

    def makeTestImages(self, seed=5, nSrc=5, psfSize=2., noiseLevel=5.,
                       fluxLevel=500., fluxRange=2.):
        """Make reproduceable PSF-convolved masked images for testing.

        Parameters
        ----------
        seed : `int`, optional
            Seed value to initialize the random number generator.
        nSrc : `int`, optional
            Number of sources to simulate.
        psfSize : `float`, optional
            Width of the PSF of the simulated sources, in pixels.
        noiseLevel : `float`, optional
            Standard deviation of the noise to add to each pixel.
        fluxLevel : float, optional
            Description
        fluxRange : `float`, optional
            Range in flux amplitude of the simulated sources.

        Returns
        -------
        modelImages : `list` of `lsst.afw.image.Image`
            A list of images, each containing the model for one subfilter

        Deleted Parameters
        ------------------
        sourceSigma : `float`, optional
            Average amplitude of the simulated sources,
            relative to ``noiseLevel``
        detectionSigma : `float`, optional
            Threshold amplitude of the image to set the "DETECTED" mask.
        """
        rng = np.random.RandomState(seed)
        x0, y0 = self.bbox.getBegin()
        xSize, ySize = self.bbox.getDimensions()
        xLoc = rng.rand(nSrc)*(xSize - 2*self.bufferSize) + self.bufferSize + x0
        yLoc = rng.rand(nSrc)*(ySize - 2*self.bufferSize) + self.bufferSize + y0

        flux = (rng.rand(nSrc)*(fluxRange - 1.) + 1.)*fluxLevel
        sigmas = [psfSize for src in range(nSrc)]
        coordList = list(zip(xLoc, yLoc, flux, sigmas))
        kernelSize = int(xSize/2)  # Need a careful explanation of this kernel size choice
        skyLevel = 0
        # Don't use the built in poisson noise: it modifies the global state of numpy random
        model = plantSources(self.bbox, kernelSize, skyLevel, coordList, addPoissonNoise=False)
        # model.variance.array += noiseLevel
        noise = rng.rand(ySize, xSize)*noiseLevel
        model.image.array += noise
        model.variance.array = (np.sqrt(np.abs(model.image.array)) + noiseLevel
                                - np.mean(np.sqrt(np.abs(noise))))

        # Run source detection to set up the mask plane
        sourceCat = self.matcher.getSelectSources(model)

        model.setWcs(self._makeWcs())
        return model, sourceCat

    @staticmethod
    def _makeWcs(offset=0):
        """ Make a fake Wcs

        Parameters
        ----------
        offset : float
          offset the Wcs by this many pixels.
        """
        # taken from $AFW_DIR/tests/testMakeWcs.py
        metadata = dafBase.PropertySet()
        metadata.set("SIMPLE", "T")
        metadata.set("BITPIX", -32)
        metadata.set("NAXIS", 2)
        metadata.set("NAXIS1", 1024)
        metadata.set("NAXIS2", 1153)
        metadata.set("RADESYS", 'FK5')
        metadata.set("EQUINOX", 2000.)
        metadata.setDouble("CRVAL1", 215.604025685476)
        metadata.setDouble("CRVAL2", 53.1595451514076)
        metadata.setDouble("CRPIX1", 1109.99981456774 + offset)
        metadata.setDouble("CRPIX2", 560.018167811613 + offset)
        metadata.set("CTYPE1", 'RA---SIN')
        metadata.set("CTYPE2", 'DEC--SIN')
        metadata.setDouble("CD1_1", 5.10808596133527E-05)
        metadata.setDouble("CD1_2", 1.85579539217196E-07)
        metadata.setDouble("CD2_2", -5.10281493481982E-05)
        metadata.setDouble("CD2_1", -8.27440751733828E-07)
        return afwGeom.makeSkyWcs(metadata)

    def _retrieveMaskPixels(self, exposure, maskPlane):
        notBadPix = (exposure.mask.array & self.statsCtrl.getAndMask()) == 0
        finitePix = np.isfinite(exposure.image.array)
        goodPix = (exposure.mask.array & afwImage.Mask.getPlaneBitMask(maskPlane)) > 0
        return goodPix & notBadPix & finitePix

    def diffimMetric0(self, residual):
        goodPix = self._retrieveMaskPixels(residual, 'DETECTED')
        nNeg = np.sum(residual.image.array[goodPix] < 0)
        nPos = np.sum(residual.image.array[goodPix] > 0)
        metric = (nPos - nNeg)/(nPos + nNeg)
        return metric

    def diffimMetric1(self, residual, sourceCat, radius=2):
        nNeg = 0
        nPos = 0
        # print(f"Length of sourceCat: {len(sourceCat)}")
        for src in sourceCat:
            srcX = int(src.getX()) - residual.getBBox().getBeginX()
            srcY = int(src.getY()) - residual.getBBox().getBeginY()
            srcRes = residual.image.array[srcX - radius: srcX + radius + 1, srcY - radius: srcY + radius + 1]
            # print(f"Shape of srcRes: {srcRes.shape}")
            nNeg += np.sum(srcRes > 0)
            nPos += np.sum(srcRes < 0)
        # goodPix = self._retrieveMaskPixels(residual, 'DETECTED')
        # nNeg = np.sum(residual.image.array[goodPix] < 0)
        # nPos = np.sum(residual.image.array[goodPix] > 0)
        print(f"Number of finite pixels: {np.sum(np.isfinite(residual.image.array))}")
        print(f"Pixels used in the metric: {nPos + nNeg}")

        if (nPos + nNeg) == 0:
            metric = 0.
        else:
            metric = (nPos - nNeg)/(nPos + nNeg)
        return metric

    @staticmethod
    def replaceNans(exposure, value=0.):
        imNans = np.isnan(exposure.image.array)
        varNans = np.isnan(exposure.variance.array)
        exposure.image.array[imNans] = value
        exposure.variance.array[varNans] = value

    def computeVarianceMean(self, exposure):
        statObj = afwMath.makeStatistics(exposure.getMaskedImage().getVariance(),
                                         exposure.getMaskedImage().getMask(),
                                         afwMath.MEANCLIP, self.statsCtrl)
        var = statObj.getValue(afwMath.MEANCLIP)
        return var

    @staticmethod
    def wrapZogyDiffim(config, refExposure, sciExposure):
        config.scaleByCalibration = False
        zogyTask = ZogyTask(config=config)

        result = zogyTask.run(sciExposure, refExposure)
        return result.diffExp

    @staticmethod
    def wrapAlDiffim(config, templateExposure, scienceExposure, convolveTemplate=True, returnKernel=False):
        alTask = ImagePsfMatchTask(config=config)
        templateFwhmPix = templateExposure.getPsf().getSigma()
        scienceFwhmPix = scienceExposure.getPsf().getSigma()
        result = alTask.subtractExposures(templateExposure, scienceExposure,
                                          templateFwhmPix=templateFwhmPix,
                                          scienceFwhmPix=scienceFwhmPix,
                                          doWarping=False,
                                          convolveTemplate=convolveTemplate,
                                          )
        if returnKernel:
            return result.psfMatchingKernel
        else:
            return result.subtractedExposure

    def testModelImages(self):
        """Check that the simulated images are useable.
        """
        sciPsf = 2.4
        refPsf = 2.
        sciNoise = 5.
        refNoise = 1.5
        fluxRatio = refPsf**2/sciPsf**2
        sciIm, src = self.makeTestImages(psfSize=sciPsf, noiseLevel=sciNoise)
        sciIm2, _ = self.makeTestImages(psfSize=sciPsf, noiseLevel=sciNoise)
        refIm, _ = self.makeTestImages(psfSize=refPsf, noiseLevel=refNoise)

        # Making the test images should be repeatable
        self.assertFloatsAlmostEqual(sciIm.image.array, sciIm2.image.array)

        diffIm = sciIm.clone()
        diffIm.image.array -= refIm.image.array

        # The "reference" image has a smaller PSF but the same source fluxes, so the peak should be greater.
        self.assertGreater(np.max(refIm.image.array), np.max(sciIm.image.array))
        # The difference image won't be zero since the two images have different PSFs,
        #  but the peak should be much lower.
        sciPeak = np.max(sciIm.image.array)
        residualPeak = np.sqrt(1 - fluxRatio)*sciPeak
        self.assertGreater(residualPeak, np.max(abs(diffIm.image.array)))

        # It should be possible to compute the diffim metric from the science and reference images
        refMetric = self.diffimMetric1(refIm, src)
        sciMetric = self.diffimMetric1(sciIm, src)
        self.assertGreaterEqual(refMetric, -1)
        self.assertGreaterEqual(sciMetric, -1)

    def testSimAlSciNotModified(self):
        "Running AL and convolving the template should not change the science image."
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 37
        rng = np.random.RandomState(seed)
        alConfig = ImagePsfMatchConfig()

        sciPsf = sciPsfBase + rng.random()*2.
        ref, _ = self.makeTestImages(seed=seed, nSrc=20, psfSize=refPsf,
                                     noiseLevel=refNoise, fluxLevel=500)
        sci, src = self.makeTestImages(seed=seed, nSrc=20, psfSize=sciPsf,
                                       noiseLevel=sciNoise, fluxLevel=500)
        # Make a deep copy of the images first
        sci2 = sci.clone()

        # Basic AL, but we don't care about the result.
        self.wrapAlDiffim(alConfig, ref, sci, convolveTemplate=True)
        self.assertMaskedImagesEqual(sci.maskedImage, sci2.maskedImage)

    @unittest.expectedFailure
    def testSimAlSciModified2(self):
        "Running AL and convolving science image should change it."
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 37
        rng = np.random.RandomState(seed)
        alConfig = ImagePsfMatchConfig()

        sciPsf = sciPsfBase + rng.random()*2.
        ref, _ = self.makeTestImages(seed=seed, nSrc=20, psfSize=refPsf,
                                     noiseLevel=refNoise, fluxLevel=500)
        sci, src = self.makeTestImages(seed=seed, nSrc=20, psfSize=sciPsf,
                                       noiseLevel=sciNoise, fluxLevel=500)
        # Make a deep copy of the images first
        sci2 = sci.clone()

        # Basic AL, but we don't care about the result.
        self.wrapAlDiffim(alConfig, ref, sci, convolveTemplate=False)
        self.assertMaskedImagesEqual(sci.maskedImage, sci2.maskedImage)

    @unittest.expectedFailure
    def testSimAlSciModified3(self):
        "The science image should be different if the template is convolved vs the science image."
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 37
        rng = np.random.RandomState(seed)
        alConfig = ImagePsfMatchConfig()

        sciPsf = sciPsfBase + rng.random()*2.
        ref1, _ = self.makeTestImages(seed=seed, nSrc=20, psfSize=refPsf,
                                      noiseLevel=refNoise, fluxLevel=500)
        sci, src = self.makeTestImages(seed=seed, nSrc=20, psfSize=sciPsf,
                                       noiseLevel=sciNoise, fluxLevel=500)
        # Make a deep copy of the images first
        sci2 = sci.clone()
        sci1B = sci.clone()
        sci2B = sci2.clone()
        ref1B = ref1.clone()

        # Basic AL, but we don't care about the result.
        self.wrapAlDiffim(alConfig, ref1B, sci1B, convolveTemplate=True)
        self.wrapAlDiffim(alConfig, ref1, sci, convolveTemplate=False)
        self.assertMaskedImagesEqual(sci2.maskedImage, sci2B.maskedImage)

    @unittest.expectedFailure
    def testSimAlRefModified(self):
        "If we convolve the template, it should be changed."
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 37
        rng = np.random.RandomState(seed)
        alConfig = ImagePsfMatchConfig()

        sciPsf = sciPsfBase + rng.random()*2.
        ref, _ = self.makeTestImages(seed=seed, nSrc=20, psfSize=refPsf,
                                     noiseLevel=refNoise, fluxLevel=500)
        sci, src = self.makeTestImages(seed=seed, nSrc=20, psfSize=sciPsf,
                                       noiseLevel=sciNoise, fluxLevel=500)
        # Make a deep copy of the images first
        ref2 = ref.clone()

        # Basic AL, but we don't care about the result.
        self.wrapAlDiffim(alConfig, ref, sci, convolveTemplate=True)
        self.assertMaskedImagesEqual(ref.maskedImage, ref2.maskedImage)

    def testSimAlRefNotModified(self):
        "If we don't convolve the template, it should not be changed."
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 37
        rng = np.random.RandomState(seed)
        alConfig = ImagePsfMatchConfig()

        sciPsf = sciPsfBase + rng.random()*2.
        ref, _ = self.makeTestImages(seed=seed, nSrc=20, psfSize=refPsf,
                                     noiseLevel=refNoise, fluxLevel=500)
        sci, src = self.makeTestImages(seed=seed, nSrc=20, psfSize=sciPsf,
                                       noiseLevel=sciNoise, fluxLevel=500)
        # Make a deep copy of the images first
        ref2 = ref.clone()

        # Basic AL, but we don't care about the result.
        self.wrapAlDiffim(alConfig, ref, sci, convolveTemplate=False)
        self.assertMaskedImagesEqual(ref.maskedImage, ref2.maskedImage)

    def testSimDiffim(self):
        nIter = 5
        refPsf = 2.4
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        metrics = {"AL": [], "ZOGY": [], "Decorr": []}
        seed = 8
        decorrelate = DecorrelateALKernelTask()
        zogyConfig = ZogyConfig()
        alConfig = ImagePsfMatchConfig()

        for s in range(nIter):
            sciPsf = sciPsfBase + s*0.2
            print(f"Iteration {s} with PSF {sciPsf}")
            ref, _ = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=refPsf,
                                         noiseLevel=refNoise, fluxLevel=500)
            sci, src = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=sciPsf,
                                           noiseLevel=sciNoise, fluxLevel=500)
            # The diffim tasks can modify the images, so make a deep copy to make sure they are independent
            sci2 = sci.clone()
            ref2 = ref.clone()

            resAl = self.wrapAlDiffim(alConfig, ref, sci)
            resZogy = self.wrapZogyDiffim(zogyConfig, ref2, sci2)
            metricZogy = self.diffimMetric1(resZogy, src)
            metrics["ZOGY"].append(metricZogy)
            metricAl = self.diffimMetric1(resAl, src)
            metrics["AL"].append(metricAl)

            mKernel = self.wrapAlDiffim(alConfig, ref, sci, returnKernel=True)
            resDecorr = decorrelate.run(sci, ref, resAl, mKernel).correctedExposure
            metricDecorr = self.diffimMetric1(resDecorr, src)
            metrics["Decorr"].append(metricDecorr)
            self.assertGreaterEqual(metricZogy, -1)
            self.assertGreaterEqual(metricAl, -1)
            self.assertGreaterEqual(metricDecorr, -1)
            print(f"Metrics: {metricAl} (AL), {metricDecorr} (AL-D), {metricZogy} (ZOGY)\n")

    def testSimReverseZogy(self):
        nIter = 5
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 18
        rng = np.random.RandomState(seed)
        zogyConfig = ZogyConfig()

        for s in range(nIter):
            sciPsf = sciPsfBase + rng.random()*2.
            ref, _ = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=refPsf,
                                         noiseLevel=refNoise, fluxLevel=500)
            sci, src = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=sciPsf,
                                           noiseLevel=sciNoise, fluxLevel=500)
            # The diffim tasks can modify the images, so make a deep copy to make sure they are independent
            sci2 = sci.clone()
            ref2 = ref.clone()

            res = self.wrapZogyDiffim(zogyConfig, ref, sci)
            resR = self.wrapZogyDiffim(zogyConfig, sci2, ref2)
            metric = self.diffimMetric1(res, src)
            metricR = self.diffimMetric1(resR, src)
            self.assertFloatsAlmostEqual(metric, -metricR)

    def testSimReverseAlNoDecorrEqualNoise(self):
        nIter = 5
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 5
        seed = 37
        rng = np.random.RandomState(seed)
        alConfig = ImagePsfMatchConfig()

        for s in range(nIter):
            sciPsf = sciPsfBase + rng.random()*2.
            ref, _ = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=refPsf,
                                         noiseLevel=refNoise, fluxLevel=500)
            sci, src = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=sciPsf,
                                           noiseLevel=sciNoise, fluxLevel=500)
            # The diffim tasks can modify the images, so make a deep copy to make sure they are independent
            sci2 = sci.clone()
            ref2 = ref.clone()

            res = self.wrapAlDiffim(alConfig, ref, sci, convolveTemplate=True)
            resR = self.wrapAlDiffim(alConfig, sci2, ref2, convolveTemplate=False)
            metric = self.diffimMetric1(res, src)
            metricR = self.diffimMetric1(resR, src)
            # Alard&Lupton is not fully reversable, but the answers should be close.
            # Partly this needs the decorrelation afterburner
            # It might also be a difference in background subtraction
            self.assertFloatsAlmostEqual(metric, -metricR, atol=.05, rtol=0.1)

    def testSimReverseAlNoDecorrUnequalNoise(self):
        nIter = 5
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 37
        rng = np.random.RandomState(seed)
        alConfig = ImagePsfMatchConfig()

        for s in range(nIter):
            sciPsf = sciPsfBase + rng.random()*2.
            ref, _ = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=refPsf,
                                         noiseLevel=refNoise, fluxLevel=500)
            sci, src = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=sciPsf,
                                           noiseLevel=sciNoise, fluxLevel=500)
            # The diffim tasks can modify the images, so make a deep copy to make sure they are independent
            sci2 = sci.clone()
            ref2 = ref.clone()

            res = self.wrapAlDiffim(alConfig, ref, sci, convolveTemplate=True)
            resR = self.wrapAlDiffim(alConfig, sci2, ref2, convolveTemplate=False)
            metric = self.diffimMetric1(res, src)
            metricR = self.diffimMetric1(resR, src)
            # Alard&Lupton is not fully reversable, but the answers should be close.
            # Partly this needs the decorrelation afterburner
            # It might also be a difference in background subtraction
            self.assertFloatsAlmostEqual(metric, -metricR, atol=.05, rtol=0.1)

    def testSimAlDecorr(self):
        nIter = 1
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 37
        rng = np.random.RandomState(seed)
        decorrelateConfig = DecorrelateALKernelConfig()
        decorrelate = DecorrelateALKernelTask(config=decorrelateConfig)
        alConfig = ImagePsfMatchConfig()

        for s in range(nIter):
            sciPsf = sciPsfBase + rng.random()*2.
            ref, _ = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=refPsf,
                                         noiseLevel=refNoise, fluxLevel=500)
            sci, src = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=sciPsf,
                                           noiseLevel=sciNoise, fluxLevel=500)
            # The diffim tasks can modify the images, so make a deep copy to make sure they are independent
            sci2 = sci.clone()
            ref2 = ref.clone()

            # Basic AL
            res = self.wrapAlDiffim(alConfig, ref, sci, convolveTemplate=True)

            # Decorrelated AL
            mKernel = self.wrapAlDiffim(alConfig, ref, sci, convolveTemplate=True, returnKernel=True)
            resD = decorrelate.run(sci, ref, res, mKernel).correctedExposure
            metricD = self.diffimMetric1(resD, src)

            # Swap the "science" and "reference" images, and alse swap which image is convolved.
            # The result is that the same image should be convolved as above
            resR = self.wrapAlDiffim(alConfig, sci2, ref2, convolveTemplate=False)

            # Swap the images as above, and also decorrelate.
            mKernelR = self.wrapAlDiffim(alConfig, sci2, ref2, convolveTemplate=False, returnKernel=True)
            resDR = decorrelate.run(ref2, sci2, resR, mKernelR).correctedExposure
            metricDR = self.diffimMetric1(resDR, src)

            self.assertFloatsAlmostEqual(metricD, -metricDR, atol=.01, rtol=0.05)
