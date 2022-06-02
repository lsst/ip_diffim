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

import unittest
import numpy as np

import lsst.afw.geom as afwGeom
from lsst.afw.image import PhotoCalib
from lsst.afw.math import Chebyshev1Function2D
import lsst.geom
from lsst.meas.algorithms.testUtils import plantSources
import lsst.ip.diffim.imagePsfMatch
from lsst.ip.diffim import subtractImages
from lsst.ip.diffim.utils import getPsfFwhm
from lsst.pex.config import FieldValidationError
import lsst.utils.tests


def makeFakeWcs():
    """Make a fake, affine Wcs.
    """
    crpix = lsst.geom.Point2D(123.45, 678.9)
    crval = lsst.geom.SpherePoint(0.1, 0.1, lsst.geom.degrees)
    cdMatrix = np.array([[5.19513851e-05, -2.81124812e-07],
                        [-3.25186974e-07, -5.19112119e-05]])
    return afwGeom.makeSkyWcs(crpix, crval, cdMatrix)


def makeTestImage(seed=5, nSrc=20, psfSize=2., noiseLevel=5.,
                  noiseSeed=6, fluxLevel=500., fluxRange=2.,
                  kernelSize=32, templateBorderSize=0,
                  background=None,
                  xSize=256,
                  ySize=256,
                  x0=12345,
                  y0=67890,
                  calibration=1.,
                  doApplyCalibration=False,
                  ):
    """Make a reproduceable PSF-convolved exposure for testing.

    Parameters
    ----------
    seed : `int`, optional
        Seed value to initialize the random number generator for sources.
    nSrc : `int`, optional
        Number of sources to simulate.
    psfSize : `float`, optional
        Width of the PSF of the simulated sources, in pixels.
    noiseLevel : `float`, optional
        Standard deviation of the noise to add to each pixel.
    noiseSeed : `int`, optional
        Seed value to initialize the random number generator for noise.
    fluxLevel : `float`, optional
        Reference flux of the simulated sources.
    fluxRange : `float`, optional
        Range in flux amplitude of the simulated sources.
    kernelSize : `int`, optional
        Size in pixels of the kernel for simulating sources.
    templateBorderSize : `int`, optional
        Size in pixels of the image border used to pad the image.
    background : `lsst.afw.math.Chebyshev1Function2D`, optional
        Optional background to add to the output image.
    xSize, ySize : `int`, optional
        Size in pixels of the simulated image.
    x0, y0 : `int`, optional
        Origin of the image.
    calibration : `float`, optional
        Conversion factor between instFlux and nJy.
    doApplyCalibration : `bool`, optional
        Apply the photometric calibration and return the image in nJy?

    Returns
    -------
    modelExposure : `lsst.afw.image.Exposure`
        The model image, with the mask and variance planes.
    sourceCat : `lsst.afw.table.SourceCatalog`
        Catalog of sources detected on the model image.
    """
    # Distance from the inner edge of the bounding box to avoid placing test
    # sources in the model images.
    bufferSize = kernelSize/2 + templateBorderSize + 1

    bbox = lsst.geom.Box2I(lsst.geom.Point2I(x0, y0), lsst.geom.Extent2I(xSize, ySize))
    if templateBorderSize > 0:
        bbox.grow(templateBorderSize)

    rng = np.random.RandomState(seed)
    rngNoise = np.random.RandomState(noiseSeed)
    x0, y0 = bbox.getBegin()
    xSize, ySize = bbox.getDimensions()
    xLoc = rng.rand(nSrc)*(xSize - 2*bufferSize) + bufferSize + x0
    yLoc = rng.rand(nSrc)*(ySize - 2*bufferSize) + bufferSize + y0

    flux = (rng.rand(nSrc)*(fluxRange - 1.) + 1.)*fluxLevel
    sigmas = [psfSize for src in range(nSrc)]
    coordList = list(zip(xLoc, yLoc, flux, sigmas))
    skyLevel = 0
    # Don't use the built in poisson noise: it modifies the global state of numpy random
    modelExposure = plantSources(bbox, kernelSize, skyLevel, coordList, addPoissonNoise=False)
    modelExposure.setWcs(makeFakeWcs())
    noise = rngNoise.randn(ySize, xSize)*noiseLevel
    modelExposure.image.array += noise
    modelExposure.variance.array = (np.sqrt(np.abs(modelExposure.image.array)) + noiseLevel
                                    - np.mean(np.sqrt(np.abs(noise))))

    # Run source detection to set up the mask plane
    psfMatchTask = lsst.ip.diffim.imagePsfMatch.ImagePsfMatchTask()
    sourceCat = psfMatchTask.getSelectSources(modelExposure)
    modelExposure.setPhotoCalib(PhotoCalib(calibration, 0., bbox))
    if background is not None:
        modelExposure.image += background
    modelExposure.maskedImage /= calibration
    if doApplyCalibration:
        modelExposure.maskedImage = modelExposure.photoCalib.calibrateImage(modelExposure.maskedImage)

    return modelExposure, sourceCat


class AlardLuptonSubtractTest(lsst.utils.tests.TestCase):

    def test_allowed_config_modes(self):
        """Verify the allowable modes for convolution.
        """
        config = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        config.mode = 'auto'
        config.mode = 'convolveScience'
        config.mode = 'convolveTemplate'

        with self.assertRaises(FieldValidationError):
            config.mode = 'aotu'

    def test_mismatched_template(self):
        """Test that an error is raised if the template
        does not fully contain the science image.
        """
        xSize = 200
        ySize = 200
        science, sources = makeTestImage(psfSize=2.4, xSize=xSize + 20, ySize=ySize + 20)
        template, _ = makeTestImage(psfSize=2.4, xSize=xSize, ySize=ySize, doApplyCalibration=True)
        config = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        task = subtractImages.AlardLuptonSubtractTask(config=config)
        with self.assertRaises(AssertionError):
            task.run(template, science, sources)

    def test_equal_images(self):
        """Test that running with enough sources produces reasonable output,
        with the same size psf in the template and science.
        """
        noiseLevel = 1.
        science, sources = makeTestImage(psfSize=2.4, noiseLevel=noiseLevel, noiseSeed=6)
        template, _ = makeTestImage(psfSize=2.4, noiseLevel=noiseLevel, noiseSeed=7,
                                    templateBorderSize=20, doApplyCalibration=True)
        config = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        config.doSubtractBackground = False
        task = subtractImages.AlardLuptonSubtractTask(config=config)
        output = task.run(template, science, sources)
        # There shoud be no NaN values in the difference image
        self.assertTrue(np.all(np.isfinite(output.difference.image.array)))
        # Mean of difference image should be close to zero.
        mean_error = noiseLevel/np.sqrt(output.difference.image.array.size)
        self.assertFloatsAlmostEqual(np.mean(output.difference.image.array), 0, atol=5*mean_error)
        # stddev of difference image should be close to expected value.
        self.assertFloatsAlmostEqual(np.std(output.difference.image.array), np.sqrt(2)*noiseLevel, rtol=0.1)

    def test_auto_convolveTemplate(self):
        """Test that auto mode gives the same result as convolveTemplate when
        the template psf is the smaller.
        """
        noiseLevel = 1.
        science, sources = makeTestImage(psfSize=3.0, noiseLevel=noiseLevel, noiseSeed=6)
        scienceOrig = science.clone()
        template, _ = makeTestImage(psfSize=2.0, noiseLevel=noiseLevel, noiseSeed=7,
                                    templateBorderSize=20, doApplyCalibration=True)
        config = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        config.doSubtractBackground = False
        config.mode = "convolveTemplate"

        task = subtractImages.AlardLuptonSubtractTask(config=config)
        output = task.run(template, science, sources)

        config.mode = "auto"
        task = subtractImages.AlardLuptonSubtractTask(config=config)
        outputAuto = task.run(template, scienceOrig, sources)
        self.assertMaskedImagesEqual(output.difference.maskedImage, outputAuto.difference.maskedImage)

    def test_auto_convolveScience(self):
        """Test that auto mode gives the same result as convolveScience when
        the science psf is the smaller.
        """
        noiseLevel = 1.
        science, sources = makeTestImage(psfSize=2.0, noiseLevel=noiseLevel, noiseSeed=6)
        scienceOrig = science.clone()
        template, _ = makeTestImage(psfSize=3.0, noiseLevel=noiseLevel, noiseSeed=7,
                                    templateBorderSize=20, doApplyCalibration=True)
        config = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        config.doSubtractBackground = False
        config.mode = "convolveScience"

        task = subtractImages.AlardLuptonSubtractTask(config=config)
        output = task.run(template, science, sources)

        config.mode = "auto"
        task = subtractImages.AlardLuptonSubtractTask(config=config)
        outputAuto = task.run(template, scienceOrig, sources)
        self.assertMaskedImagesEqual(output.difference.maskedImage, outputAuto.difference.maskedImage)

    def test_science_better(self):
        """Test that running with enough sources produces reasonable output,
        with the science psf being smaller than the template.
        """
        noiseLevel = 1.
        science, sources = makeTestImage(psfSize=2.0, noiseLevel=noiseLevel, noiseSeed=6)
        template, _ = makeTestImage(psfSize=3.0, noiseLevel=noiseLevel, noiseSeed=7,
                                    templateBorderSize=20, doApplyCalibration=True)
        config = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        config.doSubtractBackground = False
        task = subtractImages.AlardLuptonSubtractTask(config=config)
        output = task.run(template, science, sources)
        # Mean of difference image should be close to zero.
        nGoodPix = np.sum(np.isfinite(output.difference.image.array))
        mean_error = noiseLevel/np.sqrt(nGoodPix)
        self.assertFloatsAlmostEqual(np.nanmean(output.difference.image.array), 0, atol=5*mean_error)
        # stddev of difference image should be close to expected value.
        self.assertFloatsAlmostEqual(np.nanstd(output.difference.image.array),
                                     np.sqrt(2)*noiseLevel, rtol=0.1)

    def test_template_better(self):
        """Test that running with enough sources produces reasonable output,
        with the template psf being smaller than the science.
        """
        noiseLevel = 1.
        science, sources = makeTestImage(psfSize=3.0, noiseLevel=noiseLevel, noiseSeed=6)
        template, _ = makeTestImage(psfSize=2.0, noiseLevel=noiseLevel, noiseSeed=7,
                                    templateBorderSize=20, doApplyCalibration=True)
        config = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        config.doSubtractBackground = False
        task = subtractImages.AlardLuptonSubtractTask(config=config)
        output = task.run(template, science, sources)
        # There should be no NaNs in the image if we convolve the template with a buffer
        self.assertTrue(np.all(np.isfinite(output.difference.image.array)))
        # Mean of difference image should be close to zero.
        mean_error = noiseLevel/np.sqrt(output.difference.image.array.size)
        self.assertFloatsAlmostEqual(np.mean(output.difference.image.array), 0, atol=5*mean_error)
        # stddev of difference image should be close to expected value.
        self.assertFloatsAlmostEqual(np.std(output.difference.image.array), np.sqrt(2)*noiseLevel, rtol=0.1)

    def test_symmetry(self):
        """Test that convolving the science and convolving the template are
        symmetric: if the psfs are switched between them, the difference image
        should be nearly the same.
        """
        noiseLevel = 1.
        # Don't include a border for the template, in order to make the results
        #  comparable when we swap which image is treated as the "science" image.
        science, sources = makeTestImage(psfSize=2.0, noiseLevel=noiseLevel,
                                         noiseSeed=6, templateBorderSize=0)
        template, _ = makeTestImage(psfSize=3.0, noiseLevel=noiseLevel,
                                    noiseSeed=7, templateBorderSize=0, doApplyCalibration=True)
        scienceOrig = science.clone()
        config = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        config.mode = 'auto'
        config.doSubtractBackground = False
        task = subtractImages.AlardLuptonSubtractTask(config=config)

        # The science image will be modified in place, so use a copy for the second run.
        science_better = task.run(template, science, sources)
        template_better = task.run(scienceOrig, template, sources)

        delta = template_better.difference.clone()
        delta.image -= science_better.difference.image
        delta.variance -= science_better.difference.variance
        delta.mask.array -= science_better.difference.mask.array

        # Mean of delta should be very close to zero.
        nGoodPix = np.sum(np.isfinite(delta.image.array))
        mean_error = 2*noiseLevel/np.sqrt(nGoodPix)
        self.assertFloatsAlmostEqual(np.nanmean(delta.image.array), 0, atol=5*mean_error)
        # stddev of difference image should be close to expected value
        self.assertFloatsAlmostEqual(np.nanstd(delta.image.array), 2*np.sqrt(2)*noiseLevel, rtol=.1)

    def test_few_sources(self):
        """Test with only 1 source, to check that we get a useful error.
        """
        xSize = 256
        ySize = 256
        science, sources = makeTestImage(psfSize=2.4, nSrc=1, xSize=xSize, ySize=ySize)
        template, _ = makeTestImage(psfSize=2.0, nSrc=1, xSize=xSize, ySize=ySize, doApplyCalibration=True)
        config = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        task = subtractImages.AlardLuptonSubtractTask(config=config)
        with self.assertRaisesRegex(lsst.pex.exceptions.Exception,
                                    'Unable to determine kernel sum; 0 candidates'):
            task.run(template, science, sources)

    def test_images_unmodified(self):
        """Verify that image subtraction does not change the input images.
        """
        noiseLevel = 1.
        science, sources = makeTestImage(psfSize=2.0, noiseLevel=noiseLevel, noiseSeed=6)
        template, _ = makeTestImage(psfSize=2.0, noiseLevel=noiseLevel, noiseSeed=7,
                                    templateBorderSize=20, doApplyCalibration=True)
        scienceOrig = science.clone()
        templateOrig = template.clone()
        config = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        config.mode = "convolveTemplate"
        task = subtractImages.AlardLuptonSubtractTask(config=config)
        task.run(template, science, sources)
        self.assertMaskedImagesEqual(template.maskedImage, templateOrig.maskedImage)
        # The science image will have its photometric calibration applied.
        self.assertMaskedImagesEqual(science.maskedImage, scienceOrig.maskedImage)

        # Reset the science image for the second run
        science = scienceOrig.clone()
        config.mode = "convolveScience"
        task = subtractImages.AlardLuptonSubtractTask(config=config)
        task.run(template, science, sources)
        self.assertMaskedImagesEqual(template.maskedImage, templateOrig.maskedImage)

        self.assertMaskedImagesEqual(science.maskedImage,
                                     scienceOrig.maskedImage)

    def test_background_subtraction(self):
        """Check that we can recover the background,
        and that it is subtracted correctly in the difference image.
        """
        noiseLevel = 1.
        xSize = 512
        ySize = 512
        x0 = 123
        y0 = 456
        template, _ = makeTestImage(psfSize=2.0, noiseLevel=noiseLevel, noiseSeed=7,
                                    templateBorderSize=20,
                                    xSize=xSize, ySize=ySize, x0=x0, y0=y0,
                                    doApplyCalibration=True)
        params = [2.2, 2.1, 2.0, 1.2, 1.1, 1.0]

        bbox2D = lsst.geom.Box2D(lsst.geom.Point2D(x0, y0), lsst.geom.Extent2D(xSize, ySize))
        background_model = Chebyshev1Function2D(params, bbox2D)
        science, sources = makeTestImage(psfSize=2.0, noiseLevel=noiseLevel, noiseSeed=6,
                                         background=background_model,
                                         xSize=xSize, ySize=ySize, x0=x0, y0=y0)
        config = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        config.doSubtractBackground = True

        config.makeKernel.kernel.name = "AL"
        config.makeKernel.kernel.active.fitForBackground = True
        config.makeKernel.kernel.active.spatialKernelOrder = 1
        config.makeKernel.kernel.active.spatialBgOrder = 2

        def _run_and_check_images(config, mode):
            """Check that the fit background matches the input model.
            """
            config.mode = mode
            task = subtractImages.AlardLuptonSubtractTask(config=config)
            output = task.run(template, science, sources)

            # We should be fitting the same number of parameters as were in the input model
            self.assertEqual(output.backgroundModel.getNParameters(), background_model.getNParameters())

            # The parameters of the background fit should be close to the input model
            self.assertFloatsAlmostEqual(np.array(output.backgroundModel.getParameters()),
                                         np.array(params), rtol=0.3)

            # stddev of difference image should be close to expected value.
            # This will fail if we have mis-subtracted the background.
            self.assertFloatsAlmostEqual(np.nanstd(output.difference.image.array),
                                         np.sqrt(2)*noiseLevel, rtol=0.1)

        _run_and_check_images(config, "convolveTemplate")
        _run_and_check_images(config, "convolveScience")

    def test_scale_variance(self):
        """Check variance scaling of the image difference.
        """
        noiseLevel = 1.
        scaleFactor = 2.345
        science, sources = makeTestImage(psfSize=3.0, noiseLevel=noiseLevel, noiseSeed=6)
        template, _ = makeTestImage(psfSize=2.0, noiseLevel=noiseLevel, noiseSeed=7,
                                    templateBorderSize=20, doApplyCalibration=True)
        config = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        config.doSubtractBackground = False
        scienceVarianceOrig = science.variance.array[...]
        templateVarianceOrig = template.variance.array[...]

        def _run_and_check_images(doDecorrelation, doScaleVariance, scaleFactor=1.):
            """Check that the variance plane matches the expected value for
            different configurations of ``doDecorrelation`` and ``doScaleVariance``.
            """
            scienceVarMean = np.mean(science.variance.array)
            templateVarMean = np.mean(template.variance.array)
            config.doDecorrelation = doDecorrelation
            config.doScaleVariance = doScaleVariance
            task = subtractImages.AlardLuptonSubtractTask(config=config)
            output = task.run(template, science, sources)
            if doDecorrelation:
                if doScaleVariance:
                    templateNoise = templateVarMean*scaleFactor
                    scienceNoise = scienceVarMean*scaleFactor
                else:
                    templateNoise = templateVarMean
                    scienceNoise = scienceVarMean
            else:
                if doScaleVariance:
                    templateNoise = np.mean(output.matchedTemplate.variance.array)*scaleFactor
                    scienceNoise = scienceVarMean*scaleFactor
                else:
                    templateNoise = np.mean(output.matchedTemplate.variance.array)*scaleFactor
                    scienceNoise = scienceVarMean
            self.assertFloatsAlmostEqual(np.mean(output.difference.variance.array),
                                         scienceNoise + templateNoise, rtol=0.1)

        # Verify that the variance plane of the difference image is correct
        #  when the template and science variance planes are correct
        _run_and_check_images(doDecorrelation=True, doScaleVariance=True)
        _run_and_check_images(doDecorrelation=True, doScaleVariance=False)
        _run_and_check_images(doDecorrelation=False, doScaleVariance=True)
        _run_and_check_images(doDecorrelation=False, doScaleVariance=False)

        # Verify that the variance plane of the difference image is correct
        #  when the template and science variance planes are incorrect
        science.variance.array[...] = scienceVarianceOrig/scaleFactor
        template.variance.array[...] = templateVarianceOrig/scaleFactor
        _run_and_check_images(doDecorrelation=True, doScaleVariance=True, scaleFactor=scaleFactor)
        _run_and_check_images(doDecorrelation=True, doScaleVariance=False, scaleFactor=scaleFactor)
        _run_and_check_images(doDecorrelation=False, doScaleVariance=True, scaleFactor=scaleFactor)
        _run_and_check_images(doDecorrelation=False, doScaleVariance=False, scaleFactor=scaleFactor)

    def test_exposure_properties_convolve_template(self):
        """Check that all necessary exposure metadata is included
        when the template is convolved.
        """
        noiseLevel = 1.
        seed = 37
        rng = np.random.RandomState(seed)
        science, sources = makeTestImage(psfSize=3.0, noiseLevel=noiseLevel, noiseSeed=6)
        psf = science.psf
        psfAvgPos = psf.getAveragePosition()
        psfSize = getPsfFwhm(science.psf)
        psfImg = psf.computeKernelImage(psfAvgPos)
        template, _ = makeTestImage(psfSize=2.0, noiseLevel=noiseLevel, noiseSeed=7,
                                    templateBorderSize=20, doApplyCalibration=True)

        # Generate a random aperture correction map
        apCorrMap = lsst.afw.image.ApCorrMap()
        for name in ("a", "b", "c"):
            apCorrMap.set(name, lsst.afw.math.ChebyshevBoundedField(science.getBBox(), rng.randn(3, 3)))
        science.info.setApCorrMap(apCorrMap)

        config = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        config.mode = "convolveTemplate"

        def _run_and_check_images(doDecorrelation):
            """Check that the metadata is correct with or without decorrelation.
            """
            config.doDecorrelation = doDecorrelation
            task = subtractImages.AlardLuptonSubtractTask(config=config)
            output = task.run(template, science, sources)
            psfOut = output.difference.psf
            psfAvgPos = psfOut.getAveragePosition()
            if doDecorrelation:
                # Decorrelation requires recalculating the PSF,
                #  so it will not be the same as the input
                psfOutSize = getPsfFwhm(science.psf)
                self.assertFloatsAlmostEqual(psfSize, psfOutSize)
            else:
                psfOutImg = psfOut.computeKernelImage(psfAvgPos)
                self.assertImagesAlmostEqual(psfImg, psfOutImg)

            # check PSF, WCS, bbox, filterLabel, photoCalib, aperture correction
            self._compare_apCorrMaps(apCorrMap, output.difference.info.getApCorrMap())
            self.assertWcsAlmostEqualOverBBox(science.getWcs(), output.difference.getWcs(), science.getBBox())
            self.assertEqual(science.getFilterLabel(), output.difference.getFilterLabel())
            self.assertEqual(science.getPhotoCalib(), output.difference.getPhotoCalib())
        _run_and_check_images(doDecorrelation=True)
        _run_and_check_images(doDecorrelation=False)

    def test_exposure_properties_convolve_science(self):
        """Check that all necessary exposure metadata is included
        when the science image is convolved.
        """
        noiseLevel = 1.
        seed = 37
        rng = np.random.RandomState(seed)
        science, sources = makeTestImage(psfSize=2.0, noiseLevel=noiseLevel, noiseSeed=6)
        template, _ = makeTestImage(psfSize=3.0, noiseLevel=noiseLevel, noiseSeed=7,
                                    templateBorderSize=20, doApplyCalibration=True)
        psf = template.psf
        psfAvgPos = psf.getAveragePosition()
        psfSize = getPsfFwhm(template.psf)
        psfImg = psf.computeKernelImage(psfAvgPos)

        # Generate a random aperture correction map
        apCorrMap = lsst.afw.image.ApCorrMap()
        for name in ("a", "b", "c"):
            apCorrMap.set(name, lsst.afw.math.ChebyshevBoundedField(science.getBBox(), rng.randn(3, 3)))
        science.info.setApCorrMap(apCorrMap)

        config = subtractImages.AlardLuptonSubtractTask.ConfigClass()
        config.mode = "convolveScience"

        def _run_and_check_images(doDecorrelation):
            """Check that the metadata is correct with or without decorrelation.
            """
            config.doDecorrelation = doDecorrelation
            task = subtractImages.AlardLuptonSubtractTask(config=config)
            output = task.run(template, science, sources)
            if doDecorrelation:
                # Decorrelation requires recalculating the PSF,
                #  so it will not be the same as the input
                psfOutSize = getPsfFwhm(template.psf)
                self.assertFloatsAlmostEqual(psfSize, psfOutSize)
            else:
                psfOut = output.difference.psf
                psfAvgPos = psfOut.getAveragePosition()
                psfOutImg = psfOut.computeKernelImage(psfAvgPos)
                self.assertImagesAlmostEqual(psfImg, psfOutImg)

            # check PSF, WCS, bbox, filterLabel, photoCalib, aperture correction
            self._compare_apCorrMaps(apCorrMap, output.difference.info.getApCorrMap())
            self.assertWcsAlmostEqualOverBBox(science.getWcs(), output.difference.getWcs(), science.getBBox())
            self.assertEqual(science.getFilterLabel(), output.difference.getFilterLabel())
            self.assertEqual(science.getPhotoCalib(), output.difference.getPhotoCalib())

        _run_and_check_images(doDecorrelation=True)
        _run_and_check_images(doDecorrelation=False)

    def _compare_apCorrMaps(self, a, b):
        """Compare two ApCorrMaps for equality, without assuming that their BoundedFields have the
        same addresses (i.e. so we can compare after serialization).

        This function is taken from ``ApCorrMapTestCase`` in afw/tests/.

        Parameters
        ----------
        a, b : `lsst.afw.image.ApCorrMap`
            The two aperture correction maps to compare.
        """
        self.assertEqual(len(a), len(b))
        for name, value in list(a.items()):
            value2 = b.get(name)
            self.assertIsNotNone(value2)
            self.assertEqual(value.getBBox(), value2.getBBox())
            self.assertFloatsAlmostEqual(
                value.getCoefficients(), value2.getCoefficients(), rtol=0.0)


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
