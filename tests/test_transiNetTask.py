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

import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.geom
import lsst.ip.diffim.imagePsfMatch
import lsst.utils.tests
import numpy as np
from lsst.ip.diffim import subtractImagesTransiNet as subtractImages
from lsst.ip.diffim.utils import (computeRobustStatistics,
                                   getPsfFwhm, makeStats, makeTestImage)
from lsst.pex.exceptions import InvalidParameterError


class TransiNetSubtractTest(lsst.utils.tests.TestCase):
    def test_mismatched_template(self):
        """Test that an error is raised if the template
        does not fully contain the science image.
        """
        xSize = 200
        ySize = 200
        science, sources = makeTestImage(psfSize=2.4, xSize=xSize + 20, ySize=ySize + 20)
        template, _ = makeTestImage(psfSize=2.4, xSize=xSize, ySize=ySize, doApplyCalibration=True)
        config = subtractImages.TransiNetSubtractTask.ConfigClass()
        task = subtractImages.TransiNetSubtractTask(config=config)
        with self.assertRaises(AssertionError):
            task.run(template, science)

    def test_equal_images(self):
        """Test that running with enough sources produces reasonable output,
        with the same size psf in the template and science.
        """
        noiseLevel = 1.
        science, sources = makeTestImage(psfSize=2.4, noiseLevel=noiseLevel, noiseSeed=6)
        template, _ = makeTestImage(psfSize=2.4, noiseLevel=noiseLevel, noiseSeed=7,
                                    templateBorderSize=20, doApplyCalibration=True)
        config = subtractImages.TransiNetSubtractTask.ConfigClass()
        task = subtractImages.TransiNetSubtractTask(config=config)
        output = task.run(template, science)
        # There shoud be no NaN values in the difference image
        self.assertTrue(np.all(np.isfinite(output.difference.image.array)))
        # Mean of difference image should be close to zero.
        meanError = noiseLevel/np.sqrt(output.difference.image.array.size)
        # Make sure to include pixels with the DETECTED mask bit set.
        statsCtrl = makeStats(badMaskPlanes=("EDGE", "BAD", "NO_DATA"))
        differenceMean = computeRobustStatistics(output.difference.image, output.difference.mask, statsCtrl)
        self.assertFloatsAlmostEqual(differenceMean, 0, atol=5*meanError)
        # stddev of difference image should be close to expected value.
        differenceStd = computeRobustStatistics(output.difference.image, output.difference.mask,
                                                makeStats(), statistic=afwMath.STDEV)
        self.assertFloatsAlmostEqual(differenceStd, np.sqrt(2)*noiseLevel, rtol=0.1)

    def test_science_better(self):
        """Test that running with enough sources produces reasonable output,
        with the science psf being smaller than the template.
        """
        statsCtrl = makeStats()
        statsCtrlDetect = makeStats(badMaskPlanes=("EDGE", "BAD", "NO_DATA"))

        def _run_and_check_images(statsCtrl, statsCtrlDetect, scienceNoiseLevel, templateNoiseLevel):
            science, sources = makeTestImage(psfSize=2.0, noiseLevel=scienceNoiseLevel, noiseSeed=6)
            template, _ = makeTestImage(psfSize=3.0, noiseLevel=templateNoiseLevel, noiseSeed=7,
                                        templateBorderSize=20, doApplyCalibration=True)
            config = subtractImages.TransiNetSubtractTask.ConfigClass()
            config.doSubtractBackground = False
            config.mode = "convolveScience"
            task = subtractImages.TransiNetSubtractTask(config=config)
            output = task.run(template, science, sources)
            self.assertFloatsAlmostEqual(task.metadata["scaleTemplateVarianceFactor"], 1., atol=.05)
            self.assertFloatsAlmostEqual(task.metadata["scaleScienceVarianceFactor"], 1., atol=.05)
            # Mean of difference image should be close to zero.
            nGoodPix = np.sum(np.isfinite(output.difference.image.array))
            meanError = (scienceNoiseLevel + templateNoiseLevel)/np.sqrt(nGoodPix)
            diffimMean = computeRobustStatistics(output.difference.image, output.difference.mask,
                                                 statsCtrlDetect)

            self.assertFloatsAlmostEqual(diffimMean, 0, atol=5*meanError)
            # stddev of difference image should be close to expected value.
            noiseLevel = np.sqrt(scienceNoiseLevel**2 + templateNoiseLevel**2)
            varianceMean = computeRobustStatistics(output.difference.variance, output.difference.mask,
                                                   statsCtrl)
            diffimStd = computeRobustStatistics(output.difference.image, output.difference.mask,
                                                statsCtrl, statistic=afwMath.STDEV)
            self.assertFloatsAlmostEqual(varianceMean, noiseLevel**2, rtol=0.1)
            self.assertFloatsAlmostEqual(diffimStd, noiseLevel, rtol=0.1)

        _run_and_check_images(statsCtrl, statsCtrlDetect, scienceNoiseLevel=1., templateNoiseLevel=1.)
        _run_and_check_images(statsCtrl, statsCtrlDetect, scienceNoiseLevel=1., templateNoiseLevel=.1)
        _run_and_check_images(statsCtrl, statsCtrlDetect, scienceNoiseLevel=.1, templateNoiseLevel=.1)

    def test_template_better(self):
        """Test that running with enough sources produces reasonable output,
        with the template psf being smaller than the science.
        """
        statsCtrl = makeStats()
        statsCtrlDetect = makeStats(badMaskPlanes=("EDGE", "BAD", "NO_DATA"))

        def _run_and_check_images(statsCtrl, statsCtrlDetect, scienceNoiseLevel, templateNoiseLevel):
            science, sources = makeTestImage(psfSize=3.0, noiseLevel=scienceNoiseLevel, noiseSeed=6)
            template, _ = makeTestImage(psfSize=2.0, noiseLevel=templateNoiseLevel, noiseSeed=7,
                                        templateBorderSize=20, doApplyCalibration=True)
            config = subtractImages.TransiNetSubtractTask.ConfigClass()
            config.doSubtractBackground = False
            task = subtractImages.TransiNetSubtractTask(config=config)
            output = task.run(template, science, sources)
            self.assertFloatsAlmostEqual(task.metadata["scaleTemplateVarianceFactor"], 1., atol=.05)
            self.assertFloatsAlmostEqual(task.metadata["scaleScienceVarianceFactor"], 1., atol=.05)
            # There should be no NaNs in the image if we convolve the template with a buffer
            self.assertTrue(np.all(np.isfinite(output.difference.image.array)))
            # Mean of difference image should be close to zero.
            meanError = (scienceNoiseLevel + templateNoiseLevel)/np.sqrt(output.difference.image.array.size)

            diffimMean = computeRobustStatistics(output.difference.image, output.difference.mask,
                                                 statsCtrlDetect)
            self.assertFloatsAlmostEqual(diffimMean, 0, atol=5*meanError)
            # stddev of difference image should be close to expected value.
            noiseLevel = np.sqrt(scienceNoiseLevel**2 + templateNoiseLevel**2)
            varianceMean = computeRobustStatistics(output.difference.variance, output.difference.mask,
                                                   statsCtrl)
            diffimStd = computeRobustStatistics(output.difference.image, output.difference.mask,
                                                statsCtrl, statistic=afwMath.STDEV)
            self.assertFloatsAlmostEqual(varianceMean, noiseLevel**2, rtol=0.1)
            self.assertFloatsAlmostEqual(diffimStd, noiseLevel, rtol=0.1)

        _run_and_check_images(statsCtrl, statsCtrlDetect, scienceNoiseLevel=1., templateNoiseLevel=1.)
        _run_and_check_images(statsCtrl, statsCtrlDetect, scienceNoiseLevel=1., templateNoiseLevel=.1)
        _run_and_check_images(statsCtrl, statsCtrlDetect, scienceNoiseLevel=.1, templateNoiseLevel=.1)

    def test_few_sources(self):
        """Test with only 1 source, to check that we get a useful error.
        """
        xSize = 256
        ySize = 256
        science, sources = makeTestImage(psfSize=2.4, nSrc=10, xSize=xSize, ySize=ySize)
        template, _ = makeTestImage(psfSize=2.0, nSrc=10, xSize=xSize, ySize=ySize, doApplyCalibration=True)
        config = subtractImages.TransiNetSubtractTask.ConfigClass()
        task = subtractImages.TransiNetSubtractTask(config=config)
        sources = sources[0:1]
        with self.assertRaisesRegex(RuntimeError,
                                    "Cannot compute PSF matching kernel: too few sources selected."):
            task.run(template, science, sources)

    def test_scale_variance_convolve_science(self):
        """Check variance scaling of the image difference.
        """
        scienceNoiseLevel = 4.
        templateNoiseLevel = 2.
        scaleFactor = 1.345
        # Make sure to include pixels with the DETECTED mask bit set.
        statsCtrl = makeStats(badMaskPlanes=("EDGE", "BAD", "NO_DATA"))

        def _run_and_check_images(science, template, sources, statsCtrl,
                                  doDecorrelation, doScaleVariance, scaleFactor=1.):
            """Check that the variance plane matches the expected value for
            different configurations of ``doDecorrelation`` and ``doScaleVariance``.
            """

            config = subtractImages.TransiNetSubtractTask.ConfigClass()
            config.mode = "convolveScience"
            config.doSubtractBackground = False
            config.doDecorrelation = doDecorrelation
            config.doScaleVariance = doScaleVariance
            task = subtractImages.TransiNetSubtractTask(config=config)
            output = task.run(template.clone(), science.clone(), sources)
            if doScaleVariance:
                self.assertFloatsAlmostEqual(task.metadata["scaleTemplateVarianceFactor"],
                                             scaleFactor, atol=0.05)
                self.assertFloatsAlmostEqual(task.metadata["scaleScienceVarianceFactor"],
                                             scaleFactor, atol=0.05)

            templateNoise = computeRobustStatistics(template.variance, template.mask, statsCtrl)
            if doDecorrelation:
                scienceNoise = computeRobustStatistics(science.variance, science.mask, statsCtrl)
            else:
                scienceNoise = computeRobustStatistics(output.matchedScience.variance,
                                                       output.matchedScience.mask,
                                                       statsCtrl)

            if doScaleVariance:
                templateNoise *= scaleFactor
                scienceNoise *= scaleFactor

            varMean = computeRobustStatistics(output.difference.variance, output.difference.mask, statsCtrl)
            self.assertFloatsAlmostEqual(varMean, scienceNoise + templateNoise, rtol=0.1)

        science, sources = makeTestImage(psfSize=2.0, noiseLevel=scienceNoiseLevel, noiseSeed=6)
        template, _ = makeTestImage(psfSize=3.0, noiseLevel=templateNoiseLevel, noiseSeed=7,
                                    templateBorderSize=20, doApplyCalibration=True)
        # Verify that the variance plane of the difference image is correct
        #  when the template and science variance planes are correct
        _run_and_check_images(science, template, sources, statsCtrl,
                              doDecorrelation=True, doScaleVariance=True)
        _run_and_check_images(science, template, sources, statsCtrl,
                              doDecorrelation=True, doScaleVariance=False)
        _run_and_check_images(science, template, sources, statsCtrl,
                              doDecorrelation=False, doScaleVariance=True)
        _run_and_check_images(science, template, sources, statsCtrl,
                              doDecorrelation=False, doScaleVariance=False)

        # Verify that the variance plane of the difference image is correct
        #  when the template and science variance planes are incorrect
        science.variance.array /= scaleFactor
        template.variance.array /= scaleFactor
        _run_and_check_images(science, template, sources, statsCtrl,
                              doDecorrelation=True, doScaleVariance=True, scaleFactor=scaleFactor)
        _run_and_check_images(science, template, sources, statsCtrl,
                              doDecorrelation=True, doScaleVariance=False, scaleFactor=scaleFactor)
        _run_and_check_images(science, template, sources, statsCtrl,
                              doDecorrelation=False, doScaleVariance=True, scaleFactor=scaleFactor)
        _run_and_check_images(science, template, sources, statsCtrl,
                              doDecorrelation=False, doScaleVariance=False, scaleFactor=scaleFactor)

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

        config = subtractImages.TransiNetSubtractTask.ConfigClass()
        config.mode = "convolveTemplate"

        def _run_and_check_images(doDecorrelation):
            """Check that the metadata is correct with or without decorrelation.
            """
            config.doDecorrelation = doDecorrelation
            task = subtractImages.TransiNetSubtractTask(config=config)
            output = task.run(template.clone(), science.clone(), sources)
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
            self.assertWcsAlmostEqualOverBBox(science.wcs, output.difference.wcs, science.getBBox())
            self.assertEqual(science.filter, output.difference.filter)
            self.assertEqual(science.photoCalib, output.difference.photoCalib)
        _run_and_check_images(doDecorrelation=True)
        _run_and_check_images(doDecorrelation=False)


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
