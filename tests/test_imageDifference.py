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

import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.daf.base as dafBase
import lsst.geom as geom
from lsst.meas.algorithms.testUtils import plantSources
import lsst.utils.tests

from lsst.ip.diffim.imageDecorrelation import DecorrelateALKernelTask, DecorrelateALKernelConfig
from lsst.ip.diffim.imagePsfMatch import ImagePsfMatchTask, ImagePsfMatchConfig
from lsst.ip.diffim.zogy import ZogyTask, ZogyConfig


class ImageDifferenceTestBase(lsst.utils.tests.TestCase):
    """A test case for comparing image differencing algorithms.

    Attributes
    ----------
    bbox : `lsst.afw.geom.Box2I`
        Bounding box of the test model.
    bufferSize : `int`
        Distance from the inner edge of the bounding box
        to avoid placing test sources in the model images.
    nRandIter : `int`
        Number of iterations to repeat each test with random numbers.
    statsCtrl : `lsst.afw.math.StatisticsControl`
        Statistics control object.
    """

    def setUp(self):
        """Define the filter, DCR parameters, and the bounding box for the tests.
        """
        self.nRandIter = 5  # Number of iterations to repeat each test with random numbers.
        self.bufferSize = 5
        xSize = 250
        ySize = 260
        x0 = 12345
        y0 = 67890
        self.bbox = geom.Box2I(geom.Point2I(x0, y0), geom.Extent2I(xSize, ySize))
        self.statsCtrl = afwMath.StatisticsControl()
        self.statsCtrl.setNumSigmaClip(3.)
        self.statsCtrl.setNumIter(3)

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
        fluxLevel : `float`, optional
            Reference flux of the simulated sources.
        fluxRange : `float`, optional
            Range in flux amplitude of the simulated sources.

        Returns
        -------
        modelImages : `lsst.afw.image.ExposureF`
            The model image, with the mask and variance planes.
        sourceCat : `lsst.afw.table.SourceCatalog`
            Catalog of sources detected on the model image.
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
        noise = rng.rand(ySize, xSize)*noiseLevel
        model.image.array += noise
        model.variance.array = (np.sqrt(np.abs(model.image.array)) + noiseLevel
                                - np.mean(np.sqrt(np.abs(noise))))

        # Run source detection to set up the mask plane
        psfMatchTask = ImagePsfMatchTask(config=ImagePsfMatchConfig())
        sourceCat = psfMatchTask.getSelectSources(model)

        model.setWcs(self._makeWcs())
        return model, sourceCat

    @staticmethod
    def _makeWcs(offset=0):
        """Make a fake Wcs.

        Parameters
        ----------
        offset : `float`
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

    def diffimMetricBasic(self, residual, sourceCat, radius=2, sigma=0.):
        """Compute a basic metric based on the total number of positive and
        negative pixels in a residual image.

        Parameters
        ----------
        residual : `lsst.afw.image.ExposureF`
            A residual image resulting from image differencing.
        sourceCat : `lsst.afw.table.SourceCatalog`
            Source catalog containing the locations to calculate the metric.
        radius : `int`, optional
            Radius in pixels to use around each source location for the metric.
        sigma : `float`, optional
            Threshold to include pixel values in the metric.

        Returns
        -------
        `float`
            Metric assessing the image differencing residual.
        """
        nNeg = 0
        nPos = 0
        threshold = sigma*self.computeExposureStddev(residual)
        for src in sourceCat:
            srcX = int(src.getX()) - residual.getBBox().getBeginX()
            srcY = int(src.getY()) - residual.getBBox().getBeginY()
            srcRes = residual.image.array[srcY - radius: srcY + radius + 1, srcX - radius: srcX + radius + 1]
            nPos += np.sum(srcRes > threshold)
            nNeg += np.sum(srcRes < -threshold)

        if (nPos + nNeg) == 0:
            metric = 0.
        else:
            metric = (nPos - nNeg)/(nPos + nNeg)
        return metric

    def computeExposureStddev(self, exposure):
        """Compute the standard deviation of an exposure, using the mask plane.

        Parameters
        ----------
        exposure : `lsst.afw.image.ExposureF`
            The input exposure.

        Returns
        -------
        `float`
            The standard deviation of the unmasked pixels of the input image.
        """
        statObj = afwMath.makeStatistics(exposure.maskedImage.image,
                                         exposure.maskedImage.mask,
                                         afwMath.STDEVCLIP, self.statsCtrl)
        var = statObj.getValue(afwMath.STDEVCLIP)
        return var

    @staticmethod
    def wrapZogyDiffim(config, templateExposure, scienceExposure):
        """Prepare and run ZOGY-style image differencing.

        Parameters
        ----------
        config : `lsst.pex.config.Config`
            The image differencing Task configuration settings.
        templateExposure : `lsst.afw.image.ExposureF`
            The reference image to subtract from the science image.
        scienceExposure : `lsst.afw.image.ExposureF`
            The science image.

        Returns
        -------
        `lsst.afw.image.ExposureF`
            The image difference.
        """
        config.scaleByCalibration = False
        zogyTask = ZogyTask(config=config)

        result = zogyTask.run(scienceExposure, templateExposure)
        return result.diffExp

    @staticmethod
    def wrapAlDiffim(config, templateExposure, scienceExposure, convolveTemplate=True, returnKernel=False,
                     precomputeKernelCandidates=False):
        """Prepare and run Alard&Lupton-style image differencing.

        Parameters
        ----------
        config : `lsst.pex.config.Config`
            The image differencing Task configuration settings.
        templateExposure : `lsst.afw.image.ExposureF`
            The reference image to subtract from the science image.
        scienceExposure : `lsst.afw.image.ExposureF`
            The science image.
        convolveTemplate : `bool`, optional
            Option to convolve the template or the science image.
        returnKernel : `bool`, optional
            Option to return the residual image or the matching kernel.

        Returns
        -------
        `lsst.afw.image.ExposureF` or `lsst.afw.math.LinearCombinationKernel`
            The image difference, or the PSF matching kernel.
        """
        alTask = ImagePsfMatchTask(config=config)
        candidateList = None
        if precomputeKernelCandidates:
            if convolveTemplate:
                candidateList = alTask.getSelectSources(scienceExposure.clone())
            else:
                candidateList = alTask.getSelectSources(templateExposure.clone())
        templateFwhmPix = templateExposure.getPsf().getSigma()
        scienceFwhmPix = scienceExposure.getPsf().getSigma()
        result = alTask.subtractExposures(templateExposure, scienceExposure,
                                          templateFwhmPix=templateFwhmPix,
                                          scienceFwhmPix=scienceFwhmPix,
                                          doWarping=False,
                                          convolveTemplate=convolveTemplate,
                                          candidateList=candidateList,
                                          )
        if returnKernel:
            return result.psfMatchingKernel
        else:
            return result.subtractedExposure


class ImageDifferenceTestVerification(ImageDifferenceTestBase):

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
        refMetric = self.diffimMetricBasic(refIm, src, sigma=3)
        sciMetric = self.diffimMetricBasic(sciIm, src, sigma=3)
        self.assertGreaterEqual(refMetric, -1)
        self.assertGreaterEqual(sciMetric, -1)

    def testSimDiffim(self):
        "Basic smoke test to verify that the test code itself can run."
        refPsf = 2.4
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 8
        fluxLevel = 500
        decorrelate = DecorrelateALKernelTask()
        zogyConfig = ZogyConfig()
        alConfig = ImagePsfMatchConfig()

        for s in range(self.nRandIter):
            sciPsf = sciPsfBase + s*0.2
            ref, _ = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=refPsf,
                                         noiseLevel=refNoise, fluxLevel=fluxLevel)
            sci, src = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=sciPsf,
                                           noiseLevel=sciNoise, fluxLevel=fluxLevel)
            # The diffim tasks might modify the images,
            # so make a deep copy to make sure they are independent
            sci2 = sci.clone()
            ref2 = ref.clone()

            resAl = self.wrapAlDiffim(alConfig, ref, sci)
            resZogy = self.wrapZogyDiffim(zogyConfig, ref2, sci2)
            metricZogy = self.diffimMetricBasic(resZogy, src, sigma=3)
            metricAl = self.diffimMetricBasic(resAl, src, sigma=3)
            mKernel = self.wrapAlDiffim(alConfig, ref, sci, returnKernel=True)
            resDecorr = decorrelate.run(sci, ref, resAl, mKernel).correctedExposure
            metricDecorr = self.diffimMetricBasic(resDecorr, src, sigma=3)
            self.assertGreaterEqual(metricZogy, -1)
            self.assertGreaterEqual(metricAl, -1)
            self.assertGreaterEqual(metricDecorr, -1)


class ImageDifferenceTestAlardLupton(ImageDifferenceTestBase):

    def testSimAlRefNotModified(self):
        "Image differencing should not modify the original template image."
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 37
        fluxLevel = 500
        rng = np.random.RandomState(seed)
        alConfig = ImagePsfMatchConfig()

        sciPsf = sciPsfBase + rng.random()*2.
        refOriginal, _ = self.makeTestImages(seed=seed, nSrc=20, psfSize=refPsf,
                                             noiseLevel=refNoise, fluxLevel=fluxLevel)
        sciOriginal, src = self.makeTestImages(seed=seed, nSrc=20, psfSize=sciPsf,
                                               noiseLevel=sciNoise, fluxLevel=fluxLevel)
        # Make a deep copy of the images first
        sciTest1 = sciOriginal.clone()
        refTest1 = refOriginal.clone()

        # Basic AL, but we don't care about the result.
        self.wrapAlDiffim(alConfig, refTest1, sciTest1, convolveTemplate=False)
        self.assertMaskedImagesEqual(refOriginal.maskedImage, refTest1.maskedImage)

        # Basic AL, but we don't care about the result.
        self.wrapAlDiffim(alConfig, refTest1, sciTest1, convolveTemplate=True)
        self.assertMaskedImagesEqual(refOriginal.maskedImage, refTest1.maskedImage)

    def testSimAlSciNotModified(self):
        "Image differencing should not modify the original science image."
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 37
        fluxLevel = 500
        rng = np.random.RandomState(seed)
        alConfig = ImagePsfMatchConfig()

        sciPsf = sciPsfBase + rng.random()*2.
        refOriginal, _ = self.makeTestImages(seed=seed, nSrc=20, psfSize=refPsf,
                                             noiseLevel=refNoise, fluxLevel=fluxLevel)
        sciOriginal, src = self.makeTestImages(seed=seed, nSrc=20, psfSize=sciPsf,
                                               noiseLevel=sciNoise, fluxLevel=fluxLevel)
        # Make a deep copy of the images first
        sciTest1 = sciOriginal.clone()
        refTest1 = refOriginal.clone()

        # Basic AL, but we don't care about the result.
        # Note that selecting KernelCandidates *does* change the science image slightly
        # because a background is subtracted before detection, then added back in.
        # For this test, we separate out that known modification by precomputing the
        # kernel candidates in wrapAlDiffim and using a deep copy of the science image.
        # This test is therefore checking that there are no other, unknown, modifications
        # of the science image.
        self.wrapAlDiffim(alConfig, refTest1, sciTest1, convolveTemplate=True,
                          precomputeKernelCandidates=True)

        self.assertMaskedImagesEqual(sciOriginal.maskedImage, sciTest1.maskedImage)

        # Basic AL, but we don't care about the result.
        self.wrapAlDiffim(alConfig, refTest1, sciTest1, convolveTemplate=False,
                          precomputeKernelCandidates=True)

        self.assertMaskedImagesEqual(sciOriginal.maskedImage, sciTest1.maskedImage)

    def testSimReverseAlNoDecorrEqualNoise(self):
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 5
        seed = 37
        metricSigma = 0
        fluxLevel = 500
        rng = np.random.RandomState(seed)
        alConfig = ImagePsfMatchConfig()

        for s in range(self.nRandIter):
            sciPsf = sciPsfBase + rng.random()*2.
            ref, _ = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=refPsf,
                                         noiseLevel=refNoise, fluxLevel=fluxLevel)
            sci, src = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=sciPsf,
                                           noiseLevel=sciNoise, fluxLevel=fluxLevel)

            res = self.wrapAlDiffim(alConfig, ref, sci, convolveTemplate=True)
            resR = self.wrapAlDiffim(alConfig, sci, ref, convolveTemplate=False)

            metric = self.diffimMetricBasic(res, src, sigma=metricSigma)
            metricR = self.diffimMetricBasic(resR, src, sigma=metricSigma)
            # Alard&Lupton is not fully reversable, but the answers should be close.
            # Partly this needs the decorrelation afterburner
            # It might also be a difference in background subtraction
            self.assertFloatsAlmostEqual(metric, -metricR, atol=.1, rtol=.1)

    def testSimReverseAlNoDecorrUnequalNoise(self):
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 37
        metricSigma = 0
        fluxLevel = 500
        rng = np.random.RandomState(seed)
        alConfig = ImagePsfMatchConfig()

        for s in range(self.nRandIter):
            sciPsf = sciPsfBase + rng.random()*2.
            ref, _ = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=refPsf,
                                         noiseLevel=refNoise, fluxLevel=fluxLevel)
            sci, src = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=sciPsf,
                                           noiseLevel=sciNoise, fluxLevel=fluxLevel)

            res = self.wrapAlDiffim(alConfig, ref, sci, convolveTemplate=True)
            resR = self.wrapAlDiffim(alConfig, sci, ref, convolveTemplate=False)

            metric = self.diffimMetricBasic(res, src, sigma=metricSigma)
            metricR = self.diffimMetricBasic(resR, src, sigma=metricSigma)
            # Alard&Lupton is not fully reversable, but the answers should be close.
            # Partly this needs the decorrelation afterburner
            # It might also be a difference in background subtraction
            self.assertFloatsAlmostEqual(metric, -metricR, atol=.1, rtol=.1)


class ImageDifferenceTestZogy(ImageDifferenceTestBase):

    def testSimZogySciRefNotModified(self):
        "Image differencing should not modify the original images."
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 37
        fluxLevel = 500
        rng = np.random.RandomState(seed)
        zogyConfig = ZogyConfig()

        sciPsf = sciPsfBase + rng.random()*2.
        refOriginal, _ = self.makeTestImages(seed=seed, nSrc=20, psfSize=refPsf,
                                             noiseLevel=refNoise, fluxLevel=fluxLevel)
        sciOriginal, src = self.makeTestImages(seed=seed, nSrc=20, psfSize=sciPsf,
                                               noiseLevel=sciNoise, fluxLevel=fluxLevel)
        # Make a deep copy of the images first
        sciTest1 = sciOriginal.clone()
        refTest1 = refOriginal.clone()

        # Basic ZOGY, but we don't care about the result.
        self.wrapZogyDiffim(zogyConfig, refTest1, sciTest1)
        self.assertMaskedImagesEqual(refOriginal.maskedImage, refTest1.maskedImage)
        self.assertMaskedImagesEqual(sciOriginal.maskedImage, sciTest1.maskedImage)

    def testSimReverseZogy(self):
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 18
        fluxLevel = 500
        rng = np.random.RandomState(seed)
        zogyConfig = ZogyConfig()

        for s in range(self.nRandIter):
            sciPsf = sciPsfBase + rng.random()*2.
            ref, _ = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=refPsf,
                                         noiseLevel=refNoise, fluxLevel=fluxLevel)
            sci, src = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=sciPsf,
                                           noiseLevel=sciNoise, fluxLevel=fluxLevel)

            res = self.wrapZogyDiffim(zogyConfig, ref, sci)
            resR = self.wrapZogyDiffim(zogyConfig, sci, ref)
            metric = self.diffimMetricBasic(res, src, sigma=3)
            metricR = self.diffimMetricBasic(resR, src, sigma=3)
            self.assertFloatsAlmostEqual(metric, -metricR)


class ImageDifferenceTestDecorrelation(ImageDifferenceTestBase):

    def testSimAlDecorr(self):
        refPsf = 2.
        sciPsfBase = 2.
        sciNoise = 5.
        refNoise = 1.5
        seed = 37
        metricSigma = 0
        fluxLevel = 500
        rng = np.random.RandomState(seed)
        decorrelateConfig = DecorrelateALKernelConfig()
        decorrelate = DecorrelateALKernelTask(config=decorrelateConfig)
        alConfig = ImagePsfMatchConfig()

        for s in range(self.nRandIter):
            sciPsf = sciPsfBase + rng.random()*2.
            ref, _ = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=refPsf,
                                         noiseLevel=refNoise, fluxLevel=fluxLevel)
            sci, src = self.makeTestImages(seed=seed + s, nSrc=20, psfSize=sciPsf,
                                           noiseLevel=sciNoise, fluxLevel=fluxLevel)
            # The diffim tasks can modify the images, so make a deep copy to make sure they are independent
            sci2 = sci.clone()
            ref2 = ref.clone()

            # Basic AL
            res = self.wrapAlDiffim(alConfig, ref, sci, convolveTemplate=True)

            # Decorrelated AL
            mKernel = self.wrapAlDiffim(alConfig, ref, sci, convolveTemplate=True, returnKernel=True)
            resD = decorrelate.run(sci, ref, res, mKernel).correctedExposure
            metricD = self.diffimMetricBasic(resD, src, sigma=metricSigma)

            # Swap the "science" and "reference" images, and alse swap which image is convolved.
            # The result is that the same image should be convolved as above
            resR = self.wrapAlDiffim(alConfig, sci2, ref2, convolveTemplate=False)

            # Swap the images as above, and also decorrelate.
            mKernelR = self.wrapAlDiffim(alConfig, sci2, ref2, convolveTemplate=False, returnKernel=True)
            resDR = decorrelate.run(ref2, sci2, resR, mKernelR).correctedExposure
            metricDR = self.diffimMetricBasic(resDR, src, sigma=metricSigma)

            self.assertFloatsAlmostEqual(metricD, -metricDR, atol=.1, rtol=0.1)
