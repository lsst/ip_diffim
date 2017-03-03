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
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.meas.algorithms as measAlg
import lsst.pipe.base as pipeBase

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

    def run(self, subExposure, expandedSubExp, fullBBox, addNans=False, **kwargs):
        """Add `addAmount` to given `subExposure`.

        Optionally add NaNs to check the NaN-safe 'copy' operation.

        Parameters
        ----------
        subExposure : `afwImage.ExposureF`
            Input `subExposure` upon which to operate
        expandedSubExp : `afwImage.ExposureF`
            Input expanded subExposure (not used here)
        fullBBox : `afwGeom.BoundingBox`
            Bounding box of original exposure (not used here)
        addNaNs : boolean
            Set a single pixel of `subExposure` to `np.nan`
        kwargs
            Arbitrary keyword arguments (ignored)

        Returns
        -------
        `pipeBase.Struct` containing (with name 'subExposure') the
        copy of `subExposure` to which `addAmount` has been added
        """
        subExp = subExposure.clone()
        img = subExp.getMaskedImage()
        img += self.config.addAmount
        if addNans:
            img.getImage().getArray()[0, 0] = np.nan
        return pipeBase.Struct(subExposure=subExp)


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

        Returns
        -------
        `pipeBase.Struct` containing the mean value of `subExposure`
        image plane. We name it 'subExposure' to enable the correct
        test in `testNotNoneReduceWithNonExposureMapper`. In real
        operations, use something like 'mean' for the name.
        """
        subMI = subExposure.getMaskedImage()
        statObj = afwMath.makeStatistics(subMI, afwMath.MEAN)
        return pipeBase.Struct(subExposure=statObj.getValue())


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
        self.longMessage = True
        self._makeImage()

    def tearDown(self):
        del self.exposure

    def _makeImage(self):
        self.exposure = afwImage.ExposureF(128, 128)
        self.exposure.setPsf(measAlg.DoubleGaussianPsf(11, 11, 2.0, 3.7))
        mi = self.exposure.getMaskedImage()
        mi.set(0.)

    def testCopySumNoOverlaps(self):
        self._testCopySumNoOverlaps(reduceOp='copy', withNaNs=False)
        self._testCopySumNoOverlaps(reduceOp='copy', withNaNs=True)
        self._testCopySumNoOverlaps(reduceOp='sum', withNaNs=False)
        self._testCopySumNoOverlaps(reduceOp='sum', withNaNs=True)

    def _testCopySumNoOverlaps(self, reduceOp='copy', withNaNs=False):
        """Test sample grid task that adds 5.0 to input image and uses
        default 'copy' `reduceOperation`. Optionally add NaNs to subimages.
        """
        config = AddAmountImageMapReduceConfig()
        task = ImageMapReduceTask(config)
        config.mapperSubtask.addAmount = 5.
        config.reducerSubtask.reduceOperation = reduceOp
        newExp = task.run(self.exposure, addNans=withNaNs).exposure
        newMI = newExp.getMaskedImage()
        newArr = newMI.getImage().getArray()
        isnan = np.isnan(newArr)
        if not withNaNs:
            self.assertEqual(np.sum(isnan), 0,
                             msg='Failed on withNaNs: %s' % str(withNaNs))

        mi = self.exposure.getMaskedImage().getImage().getArray()
        if reduceOp != 'sum':
            self.assertFloatsAlmostEqual(mi[~isnan], newArr[~isnan] - 5.,
                                         msg='Failed on withNaNs: %s' % str(withNaNs))

    def testAverageWithOverlaps(self):
        self._testAverageWithOverlaps(withNaNs=False)
        self._testAverageWithOverlaps(withNaNs=True)

    def _testAverageWithOverlaps(self, withNaNs=False):
        """Test sample grid task that adds 5.0 to input image and uses
        'average' `reduceOperation`. Optionally add NaNs to subimages.
        """
        config = AddAmountImageMapReduceConfig()
        config.gridStepX = config.gridStepY = 8
        config.reducerSubtask.reduceOperation = 'average'
        task = ImageMapReduceTask(config)
        config.mapperSubtask.addAmount = 5.
        newExp = task.run(self.exposure, addNans=withNaNs).exposure
        newMI = newExp.getMaskedImage()
        newArr = newMI.getImage().getArray()
        mi = self.exposure.getMaskedImage()
        isnan = np.isnan(newArr)
        if not withNaNs:
            self.assertEqual(np.sum(isnan), 0,
                             msg='Failed on withNaNs: %s' % str(withNaNs))

        mi = self.exposure.getMaskedImage().getImage().getArray()
        self.assertFloatsAlmostEqual(mi[~isnan], newArr[~isnan] - 5.,
                                     msg='Failed on withNaNs: %s' % str(withNaNs))

    def testAverageVersusCopy(self):
        self.testAverageVersusCopy(withNaNs=False)
        self.testAverageVersusCopy(withNaNs=True)

    def _testAverageVersusCopy(self, withNaNs=False):
        """Re-run `testExampleTaskNoOverlaps` and `testExampleTaskWithOverlaps`
        on a more complex image (with random noise). Ensure that the results are
        identical.
        """
        exposure1 = self.exposure.clone()
        img = exposure1.getMaskedImage().getImage()
        afwMath.randomGaussianImage(img, afwMath.Random())
        config = AddAmountImageMapReduceConfig()
        task = ImageMapReduceTask(config)
        config.mapperSubtask.addAmount = 5.
        newExp = task.run(exposure1, addNans=withNaNs).exposure
        newMI1 = newExp.getMaskedImage()

        exposure2 = self.exposure.clone()
        img = exposure2.getMaskedImage().getImage()
        afwMath.randomGaussianImage(img, afwMath.Random())
        config = AddAmountImageMapReduceConfig()
        config.gridStepX = config.gridStepY = 8
        config.reducerSubtask.reduceOperation = 'average'
        task = ImageMapReduceTask(config)
        config.mapperSubtask.addAmount = 5.
        newExp = task.run(exposure2, addNans=withNaNs).exposure
        newMI2 = newExp.getMaskedImage()

        newMA1 = newMI1.getImage().getArray()
        isnan = np.isnan(newMA1)
        if not withNaNs:
            self.assertEqual(np.sum(isnan), 0)
        newMA2 = newMI2.getImage().getArray()

        self.assertFloatsAlmostEqual(newMA1[~isnan], newMA2[~isnan])

    def testMean(self):
        """Test sample grid task that returns the mean of the subimages and uses
        'none' `reduceOperation`.
        """
        config = GetMeanImageMapReduceConfig()
        config.reducerSubtask.reduceOperation = 'none'
        task = ImageMapReduceTask(config)
        testExposure = self.exposure.clone()
        testExposure.getMaskedImage().set(1.234)
        subMeans = task.run(testExposure).result
        subMeans = [x.subExposure for x in subMeans]

        self.assertEqual(len(subMeans), len(task.boxes0))
        firstPixel = testExposure.getMaskedImage().getImage().getArray()[0, 0]
        self.assertFloatsAlmostEqual(np.array(subMeans), firstPixel)

    def testGridCentroids(self):
        """Test sample grid task which is provided a set of `gridCentroids` and
        returns the mean of the subimages surrounding those centroids using 'none'
        for `reduceOperation`.
        """
        config = GetMeanImageMapReduceConfig()
        config.gridStepX = config.gridStepY = 8
        config.reducerSubtask.reduceOperation = 'none'
        config.gridCentroidsX = [i for i in np.linspace(0, 128, 50)]
        config.gridCentroidsY = config.gridCentroidsX
        task = ImageMapReduceTask(config)
        testExposure = self.exposure.clone()
        testExposure.getMaskedImage().set(1.234)
        subMeans = task.run(testExposure).result
        subMeans = [x.subExposure for x in subMeans]

        self.assertEqual(len(subMeans), len(config.gridCentroidsX))
        firstPixel = testExposure.getMaskedImage().getImage().getArray()[0, 0]
        self.assertFloatsAlmostEqual(np.array(subMeans), firstPixel)

    def testGridCentroidsWrongLength(self):
        """Test sample grid task which is provided a set of `gridCentroids` and
        returns the mean of the subimages surrounding those centroids using 'none'
        for `reduceOperation`. In this case, we ensure that len(task.boxes0) !=
        len(task.boxes1) and check for ValueError.
        """
        config = GetMeanImageMapReduceConfig()
        config.reducerSubtask.reduceOperation = 'none'
        config.gridCentroidsX = [i for i in np.linspace(0, 128, 50)]
        config.gridCentroidsY = [i for i in np.linspace(0, 128, 50)]
        task = ImageMapReduceTask(config)
        task._generateGrid(self.exposure)
        del task.boxes0[-1]  # remove the last box
        with self.assertRaises(ValueError):
            task.run(self.exposure)

    def testNotNoneReduceWithNonExposureMapper(self):
        """Test that a combination of a mapperSubtask that returns a non-exposure
        cannot work correctly with a reducerSubtask with reduceOperation='none'.
        Should raise a TypeError.
        """
        config = GetMeanImageMapReduceConfig()  # mapper returns a float (mean)
        config.gridStepX = config.gridStepY = 8
        config.reducerSubtask.reduceOperation = 'average'  # not 'none'!
        task = ImageMapReduceTask(config)
        with self.assertRaises(TypeError):
            task.run(self.exposure)

    def testGridValidity(self):
        """Test sample grids with various spacings and sizes and other options.
        """
        expectedVal = 1.
        n_tests = 0

        for reduceOp in ('copy', 'average'):
            for adjustGridOption in ('spacing', 'size', 'none'):
                for gstepx in range(11, 3, -4):
                    for gsizex in gstepx + np.array([0, 1, 2]):
                        for gstepy in range(11, 3, -4):
                            for gsizey in gstepy + np.array([0, 1, 2]):
                                config = AddAmountImageMapReduceConfig()
                                config.reducerSubtask.reduceOperation = reduceOp
                                n_tests += 1
                                self._runGridValidity(config, gstepx, gsizex,
                                                      gstepy, gsizey, adjustGridOption,
                                                      expectedVal)
        print("Ran a total of %d grid validity tests." % n_tests)

    def _runGridValidity(self, config, gstepx, gsizex, gstepy, gsizey,
                         adjustGridOption, expectedVal=1.):
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
            config.adjustGridOption = adjustGridOption
            task = ImageMapReduceTask(config)
            task._generateGrid(self.exposure)
            ind = 0 if scaleByFwhm else 1
            lenBoxes[ind] = len(task.boxes0)
            newExp = task.run(self.exposure).exposure
            newMI = newExp.getMaskedImage()
            newArr = newMI.getImage().getArray()
            isnan = np.isnan(newArr)
            self.assertEqual(np.sum(isnan), 0, msg='Failed NaN (%d), on config: %s' %
                             (np.sum(isnan), str(config)))

            mi = self.exposure.getMaskedImage().getImage().getArray()
            self.assertFloatsAlmostEqual(mi[~isnan], newArr[~isnan] - expectedVal,
                                         msg='Failed on config: %s' % str(config))

        self.assertLess(lenBoxes[0], lenBoxes[1], msg='Failed lengths on config: %s' %
                        str(config))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
