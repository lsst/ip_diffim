from __future__ import absolute_import, division, print_function
from builtins import range
from past.builtins import basestring
#
# LSST Data Management System
# Copyright 2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.

import unittest
import numpy as np

import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.pex.config as pexConfig
import lsst.meas.algorithms as measAlg
from lsst.ip.diffim.imageGridder import (ImageGridderTask, ImageGridderConfig,
                                         ImageGridSubtask, ImageGridSubtaskConfig)


def setup_module(module):
    lsst.utils.tests.init()


class TestImageGridSubtaskConfig(ImageGridSubtaskConfig):
    """!
    \anchor TestImageGridSubtaskConfig_

    \brief Configuration parameters for the TestImageGridSubtask
    """
    addAmount = pexConfig.Field(
        dtype=float,
        doc="""Amount to add to image""",
        default=10.
    )

class TestImageGridSubtask(ImageGridSubtask):
    ConfigClass = TestImageGridSubtaskConfig
    _DefaultName = "ip_diffim_TestImageGridSubtask"

    def __init__(self, *args, **kwargs):
        """! Create the image gridding subTask
        @param *args arguments to be passed to lsst.pipe.base.task.Task.__init__
        @param **kwargs keyword arguments to be passed to lsst.pipe.base.task.Task.__init__
        """
        ImageGridSubtask.__init__(self, *args, **kwargs)

    def run(self, subExp, expandedSubExp, fullBBox, **kwargs):
        """! Add `addAmount` to given subExposure.

        @param[in] subExp the sub-exposure of `exposure` upon which to operate
        @param[in] expandedSubExp the expanded sub-exposure of `exposure` upon which to operate
        @param[in] fullBBox the bounding box of the original exposure
        @return a `afw.Exp` or `afw.Exposure`
        """
        subExp = subExp.clone()
        img = subExp.getMaskedImage()
        img += self.config.addAmount
        return subExp

class TestImageGridderConfig(ImageGridderConfig):
    """!
    \anchor TestImageGridderConfig_

    \brief Configuration parameters for the TestImageGridderTask
    """
    gridSubtask = pexConfig.ConfigurableField(
        doc="Subtask to run on each subimage",
        target=TestImageGridSubtask,
    )

class TestImageGridSubtask2(ImageGridSubtask):
    ConfigClass = TestImageGridSubtaskConfig
    _DefaultName = "ip_diffim_TestImageGridSubtask2"

    def __init__(self, *args, **kwargs):
        """! Create the image gridding subTask
        @param *args arguments to be passed to lsst.pipe.base.task.Task.__init__
        @param **kwargs keyword arguments to be passed to lsst.pipe.base.task.Task.__init__
        """
        ImageGridSubtask.__init__(self, *args, **kwargs)

    def run(self, subExp, expandedSubExp, fullBBox, **kwargs):
        """!Return the mean of the given subExposure.

        @param[in] subExp the sub-exposure of `exposure` upon which to operate
        @param[in] expandedSubExp the expanded sub-exposure of `exposure` upon which to operate
        @param[in] fullBBox the bounding box of the original exposure
        @return a float, the mean of `subExp.getMaskedImage().getArray()`.
        """
        return subExp.getMaskedImage().getImage().getArray().mean()

class TestImageGridderConfig2(ImageGridderConfig):
    """!
    \anchor TestImageGridderConfig_

    \brief Configuration parameters for the TestImageGridderTask2
    """
    gridSubtask = pexConfig.ConfigurableField(
        doc="Subtask to run on each subimage",
        target=TestImageGridSubtask2,
    )


class ImageGridderTest(lsst.utils.tests.TestCase):
    """!A test case for the image gridded processing task
    """

    def setUp(self):
        self.exposure = afwImage.ExposureF(128, 128)
        self.exposure.setPsf(measAlg.DoubleGaussianPsf(11, 11, 2.0, 3.7))
        mi = self.exposure.getMaskedImage()
        mi[:, :] = 0.

    def tearDown(self):
        del self.exposure

    def testExampleTaskNoOverlaps(self):
        """Test sample grid task that adds 5.0 to input image and uses
        default 'copy' `reduceOperation`.
        """
        config = TestImageGridderConfig()
        task = ImageGridderTask(config)
        config.gridSubtask.addAmount = 5.
        newExp = task.run(self.exposure)
        newMI = newExp.getMaskedImage()
        newArr = newMI.getImage().getArray()
        mi = self.exposure.getMaskedImage()

        self.assertClose(mi.getImage().getArray().mean(), newArr.mean() - 5.)
        self.assertClose(newArr.mean(), 5.)
        self.assertClose(newArr.min(), 5.)
        self.assertClose(newArr.max(), 5.)

    def testExampleTaskWithOverlaps(self):
        """Test sample grid task that adds 5.0 to input image and uses
        'average' `reduceOperation`.
        """
        config = TestImageGridderConfig()
        config.gridStepX = config.gridStepY = 8
        config.reduceOperation = 'average'
        task = ImageGridderTask(config)
        config.gridSubtask.addAmount = 5.
        newExp = task.run(self.exposure)
        newMI = newExp.getMaskedImage()
        newArr = newMI.getImage().getArray()
        mi = self.exposure.getMaskedImage()

        self.assertClose(mi.getImage().getArray().mean(), newArr.mean() - 5.)
        self.assertClose(newArr.mean(), 5.)
        self.assertClose(newArr.min(), 5.)
        self.assertClose(newArr.max(), 5.)

    def testExampleTaskMean(self):
        """Test sample grid task that returns the mean of the subimages and uses
        'none' `reduceOperation`.
        """
        config = TestImageGridderConfig2()
        config.gridStepX = config.gridStepY = 8
        config.reduceOperation = 'none'
        task = ImageGridderTask(config)
        subMeans = task.run(self.exposure)

        self.assertEqual(len(subMeans), len(task.boxes0))
        self.assertClose(subMeans[0], 0.)

    def testGridCentroids(self):
        """Test sample grid task which is provided a set of `gridCentroids` and
        returns the mean of the subimages surrounding those centroids using 'none'
        for `reduceOperation`.
        """
        config = TestImageGridderConfig2()
        config.gridStepX = config.gridStepY = 8
        config.reduceOperation = 'none'
        centroidsX = [i for i in np.linspace(0, 128, 50)]
        centroidsY = centroidsX
        config.gridCentroidsX = centroidsX
        config.gridCentroidsY = centroidsY
        task = ImageGridderTask(config)
        subMeans = task.run(self.exposure)

        self.assertEqual(len(subMeans), len(centroidsX))
        self.assertClose(subMeans[0], 0.)

    ## TBD: make task where 'gridCentroids' is not None; ensure validity and
    ##  if 'reduceOperation' = 'none' that the length of the resulting list
    ##  is the same as the input centroids

    def testGridValidity(self):
        """!Test sample grids with various spacings and sizes and other options.

        Try to see if we can break it.
        """
        config = TestImageGridderConfig()
        config.reduceOperation = 'copy'
        expectedVal = 1.
        config.gridSubtask.addAmount = expectedVal

        for scaleByFwhm in (False, True):
            config.scaleByFwhm = scaleByFwhm
            for gsx in range(51, 1, -10):
                config.gridStepX = config.gridSizeX = gsx
                for gsy in range(51, 1, -10):
                    config.gridStepY = config.gridSizeY = gsy

                    task = ImageGridderTask(config)
                    boxes = task._generateGrid(self.exposure)
                    if len(boxes[0]) > 1000:  # bypass to prevent slow testing
                        continue
                    newExp = task.run(self.exposure)
                    newMI = newExp.getMaskedImage()
                    newArr = newMI.getImage().getArray()
                    mi = self.exposure.getMaskedImage()

                    mn = newArr.mean()
                    self.assertClose(mi.getImage().getArray().mean(), mn - expectedVal)
                    self.assertClose(mn, expectedVal)

class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
