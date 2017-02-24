from __future__ import absolute_import, division, print_function
from builtins import range
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
import lsst.pex.config as pexConfig
import lsst.meas.algorithms as measAlg
from lsst.ip.diffim.imageMapReduce import (ImageMapReduceTask, ImageMapReduceConfig,
                                           ImageMapperSubtask, ImageMapperSubtaskConfig,
                                           ImageReducerSubtask, ImageReducerSubtaskConfig)


def setup_module(module):
    lsst.utils.tests.init()


class AddAmountImageMapperSubtaskConfig(ImageMapperSubtaskConfig):
    """Configuration parameters for the AddAmountImageMapperSubtask
    """
    addAmount = pexConfig.Field(
        dtype=float,
        doc="Amount to add to image",
        default=10.
    )


class AddAmountImageMapperSubtask(ImageMapperSubtask):
    """Image mapper subTask that adds a constant value to the input subexposure
    """
    ConfigClass = AddAmountImageMapperSubtaskConfig
    _DefaultName = "ip_diffim_AddAmountImageMapperSubtask"

    def run(self, subExposure, expandedSubExp, fullBBox, **kwargs):
        """Add `addAmount` to given `subExposure`.

        Parameters
        ----------
        subExposure : `afwImage.ExposureF`
            Input `subExposure` upon which to operate
        expandedSubExp : `afwImage.ExposureF`
            Input expanded subExposure (not used here)
        fullBBox : `afwGeom.BoundingBox`
            Bounding box of original exposure (not used here)
        kwargs
            Arbitrary keyword arguments (ignored)
        -------
        Returns
        -------
        subExp : `afwImage.ExposureF`
           Copy of `subExposure` to which `addAmount` has been added
        """
        subExp = subExposure.clone()
        img = subExp.getMaskedImage()
        img += self.config.addAmount
        return subExp


class AddAmountImageMapReduceConfig(ImageMapReduceConfig):
    """Configuration parameters for the AddAmountImageMapReduceTask
    """
    mapperSubtask = pexConfig.ConfigurableField(
        doc="Mapper subtask to run on each subimage",
        target=AddAmountImageMapperSubtask,
    )


class GetMeanImageMapperSubtask(ImageMapperSubtask):
    """ImageMapper subtask that computes and returns the mean value of the
    input sub-exposure
    """
    ConfigClass = AddAmountImageMapperSubtaskConfig  # Doesn't need its own config
    _DefaultName = "ip_diffim_GetMeanImageMapperSubtask"

    def run(self, subExposure, expandedSubExp, fullBBox, **kwargs):
        """Compute the mean of the given `subExposure`

        Parameters
        ----------
        subExposure : `afwImage.ExposureF`
            Input `subExposure` upon which to operate
        expandedSubExp : `afwImage.ExposureF`
            Input expanded subExposure (not used here)
        fullBBox : `afwGeom.BoundingBox`
            Bounding box of original exposure (not used here)
        kwargs
            Arbitrary keyword arguments (ignored)
        -------
        Returns
        -------
        mean : `float`
           Mean value of `subExposure` image plane
        """
        subArr = subExposure.getMaskedImage().getImage().getArray()
        return subArr.mean()


class GetMeanImageMapReduceConfig(ImageMapReduceConfig):
    """Configuration parameters for the GetMeanImageMapReduceTask
    """
    mapperSubtask = pexConfig.ConfigurableField(
        doc="Mapper subtask to run on each subimage",
        target=GetMeanImageMapperSubtask,
    )


class ImageMapReduceTest(lsst.utils.tests.TestCase):
    """A test case for the image gridded processing task
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
        config = AddAmountImageMapReduceConfig()
        task = ImageMapReduceTask(config)
        config.mapperSubtask.addAmount = 5.
        newExp = task.run(self.exposure)
        newMI = newExp.getMaskedImage()
        newArr = newMI.getImage().getArray()
        mi = self.exposure.getMaskedImage()

        self.assertFloatsAlmostEqual(mi.getImage().getArray(), newArr - 5.)

    def testExampleTaskWithOverlaps(self):
        """Test sample grid task that adds 5.0 to input image and uses
        'average' `reduceOperation`.
        """
        config = AddAmountImageMapReduceConfig()
        config.gridStepX = config.gridStepY = 8
        config.reducerSubtask.reduceOperation = 'average'
        task = ImageMapReduceTask(config)
        config.mapperSubtask.addAmount = 5.
        newExp = task.run(self.exposure)
        newMI = newExp.getMaskedImage()
        newArr = newMI.getImage().getArray()
        mi = self.exposure.getMaskedImage()

        self.assertFloatsAlmostEqual(mi.getImage().getArray(), newArr - 5.)

    def testExampleTaskMean(self):
        """Test sample grid task that returns the mean of the subimages and uses
        'none' `reduceOperation`.
        """
        config = GetMeanImageMapReduceConfig()
        config.gridStepX = config.gridStepY = 8
        config.reducerSubtask.reduceOperation = 'none'
        task = ImageMapReduceTask(config)
        testExposure = self.exposure.clone()
        testExposure.getMaskedImage()[:, :] = 1.234
        subMeans = task.run(testExposure)

        self.assertEqual(len(subMeans), len(task.boxes0))
        self.assertFloatsAlmostEqual(np.array(subMeans), 1.234, rtol=1e-6)

    def testGridCentroids(self):
        """Test sample grid task which is provided a set of `gridCentroids` and
        returns the mean of the subimages surrounding those centroids using 'none'
        for `reduceOperation`.
        """
        config = GetMeanImageMapReduceConfig()
        config.gridStepX = config.gridStepY = 8
        config.reducerSubtask.reduceOperation = 'none'
        centroidsX = [i for i in np.linspace(0, 128, 50)]
        centroidsY = centroidsX
        config.gridCentroidsX = centroidsX
        config.gridCentroidsY = centroidsY
        task = ImageMapReduceTask(config)
        testExposure = self.exposure.clone()
        testExposure.getMaskedImage()[:, :] = 1.234
        subMeans = task.run(testExposure)

        self.assertEqual(len(subMeans), len(centroidsX))
        self.assertFloatsAlmostEqual(np.array(subMeans), 1.234, rtol=1e-6)

    def testGridValidity(self):
        """Test sample grids with various spacings and sizes and other options.

        Try to see if we can break it.
        """
        expectedVal = 1.
        n_tests = 0

        for reduceOp in ('copy', 'average'):
            for gstepx in range(11, 2, -3):
                for gsizex in gstepx + np.array([0, 1, 2]):
                    for gstepy in range(11, 2, -3):
                        for gsizey in gstepy + np.array([0, 1, 2]):
                            config = AddAmountImageMapReduceConfig()
                            config.reducerSubtask.reduceOperation = reduceOp
                            n_tests += 1
                            self._testGridValidity(config, gstepx, gsizex, gstepy, gsizey, expectedVal)
        print("Ran a total of %d grid validity tests." % n_tests)

    def _testGridValidity(self, config, gstepx, gsizex, gstepy, gsizey, expectedVal=1.):
        """Method to test the grid validity given an input config.

        Here we also iterate over scaleByFwhm in (True, False) and
        ensure that we get more `boxes` when `scaleByFwhm=False` than
        vice versa.

        Parameters
        ----------
        config : `ipDiffim.AddAmountImageMapReduceConfig`
            input AddAmountImageMapReduceConfig
        gstepx : `float`
            grid x-direction step size
        gsizex : `float`
            grid x-direction box size
        gstepy : `float`
            grid y-direction step size
        gsizey : `float`
            grid y-direction box size
        expectedVal : `float`
            float to add to exposure (to compare for testing)
        """
        config.mapperSubtask.addAmount = expectedVal
        lenBoxes = [0, 0]
        for scaleByFwhm in (True, False):
            config.scaleByFwhm = scaleByFwhm
            if scaleByFwhm:
                config.gridStepX = gstepx
                config.gridSizeX = gsizex
                config.gridStepY = gstepy
                config.gridSizeY = gsizey
            else:  # otherwise the grid is too fine and elements too small.
                config.gridStepX = gstepx * 3.
                config.gridSizeX = gsizex * 3.
                config.gridStepY = gstepy * 3.
                config.gridSizeY = gsizey * 3.
            task = ImageMapReduceTask(config)
            boxes = task._generateGrid(self.exposure)
            ind = 0 if scaleByFwhm else 1
            lenBoxes[ind] = len(boxes[0])
            if len(boxes[0]) > 1000:  # bypass to prevent slow testing
                continue
            newExp = task.run(self.exposure)
            newMI = newExp.getMaskedImage()
            newArr = newMI.getImage().getArray()
            mi = self.exposure.getMaskedImage()

            self.assertFloatsAlmostEqual(mi.getImage().getArray(), newArr - expectedVal,
                                         msg='Failed on config: %s' % str(config))

        self.assertTrue(lenBoxes[0] < lenBoxes[1], msg='Failed on config: %s' % str(config))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
