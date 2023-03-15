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
import lsst.geom
import lsst.ip.diffim.imagePsfMatch
import lsst.utils.tests
import numpy as np
from lsst.ip.diffim import subtractImagesTransiNet as subtractImages
from lsst.ip.diffim.utils import (computeRobustStatistics,
                                  makeStats, makeTestImage)


class TransiNetSubtractTest(lsst.utils.tests.TestCase):

    def setUp(self):
        """Set up the test.
        """
        self.config = subtractImages.TransiNetSubtractTask.ConfigClass()
        self.config.modelPackageName = "TN_39b"

    def test_correct_interface_init(self):
        """Test that task is initialized with the correct interface.
        """
        task = subtractImages.TransiNetSubtractTask(config=self.config)
        self.assertEqual(task.transiNetInterface.model_package_name,
                         self.config.modelPackageName)

    def test_mismatched_template(self):
        """Test that an error is raised if the template
        does not fully contain the science image.
        """
        xSize = 200
        ySize = 200
        science, sources = makeTestImage(psfSize=2.4, xSize=xSize + 20, ySize=ySize + 20)
        template, _ = makeTestImage(psfSize=2.4, xSize=xSize, ySize=ySize, doApplyCalibration=True)
        task = subtractImages.TransiNetSubtractTask(config=self.config)
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
        task = subtractImages.TransiNetSubtractTask(config=self.config)
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

        # TransiNet normally suppresses all noise. So just for the sake of compatibility with
        # tests of other algorithms, we set test whether it is lower than the expected value
        # or not.
        self.assertLess(differenceStd, np.sqrt(2)*noiseLevel)

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
            task = subtractImages.TransiNetSubtractTask(config=self.config)
            output = task.run(template, science)

            # Mean of difference image should be close to zero.
            nGoodPix = np.sum(np.isfinite(output.difference.image.array))
            meanError = (scienceNoiseLevel + templateNoiseLevel)/np.sqrt(nGoodPix)
            diffimMean = computeRobustStatistics(output.difference.image, output.difference.mask,
                                                 statsCtrlDetect)

            self.assertFloatsAlmostEqual(diffimMean, 0, atol=5*meanError)

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
            task = subtractImages.TransiNetSubtractTask(config=self.config)
            output = task.run(template, science)
            # There should be no NaNs in the image if we convolve the template with a buffer
            self.assertTrue(np.all(np.isfinite(output.difference.image.array)))
            # Mean of difference image should be close to zero.
            meanError = (scienceNoiseLevel + templateNoiseLevel)/np.sqrt(output.difference.image.array.size)

            diffimMean = computeRobustStatistics(output.difference.image, output.difference.mask,
                                                 statsCtrlDetect)
            self.assertFloatsAlmostEqual(diffimMean, 0, atol=5*meanError)
            # stddev of difference image should be close to expected value.
            diffimStd = computeRobustStatistics(output.difference.image, output.difference.mask,
                                                statsCtrl, statistic=afwMath.STDEV)

            self.assertFloatsAlmostEqual(diffimStd, 0, atol=0.003)

        _run_and_check_images(statsCtrl, statsCtrlDetect, scienceNoiseLevel=1., templateNoiseLevel=1.)
        _run_and_check_images(statsCtrl, statsCtrlDetect, scienceNoiseLevel=1., templateNoiseLevel=.1)
        _run_and_check_images(statsCtrl, statsCtrlDetect, scienceNoiseLevel=.1, templateNoiseLevel=.1)


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
