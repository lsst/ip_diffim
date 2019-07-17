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
import lsst.geom as geom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.daf.base as dafBase
import lsst.meas.algorithms as measAlg
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.ip.diffim.imageMapReduce import (ImageMapReduceTask, ImageMapReduceConfig,
                                           ImageMapper, ImageMapperConfig)


def setup_module(module):
    lsst.utils.tests.init()


def makeWcs(offset=0):
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


def getPsfMoments(psfArray):
    # Borrowed and modified from meas_algorithms/testCoaddPsf
    sumx2 = sumy2 = sumy = sumx = sumf = 0.0
    for x in range(psfArray.shape[0]):
        for y in range(psfArray.shape[1]):
            f = psfArray[x, y]
            sumx2 += x*x*f
            sumy2 += y*y*f
            sumx += x*f
            sumy += y*f
            sumf += f
    xbar = sumx/sumf
    ybar = sumy/sumf
    mxx = sumx2 - 2*xbar*sumx + xbar*xbar*sumf
    myy = sumy2 - 2*ybar*sumy + ybar*ybar*sumf
    return sumf, xbar, ybar, mxx, myy


def getPsfSecondMoments(psfArray):
    sum, xbar, ybar, mxx, myy = getPsfMoments(psfArray)
    return mxx, myy


class AddAmountImageMapperConfig(ImageMapperConfig):
    """Configuration parameters for the AddAmountImageMapper
    """
    addAmount = pexConfig.Field(
        dtype=float,
        doc="Amount to add to image",
        default=10.
    )


class AddAmountImageMapper(ImageMapper):
    """Image mapper subTask that adds a constant value to the input subexposure
    """
    ConfigClass = AddAmountImageMapperConfig
    _DefaultName = "ip_diffim_AddAmountImageMapper"

    def run(self, subExposure, expandedSubExp, fullBBox, addNans=False, **kwargs):
        """Add `addAmount` to given `subExposure`.

        Optionally add NaNs to check the NaN-safe 'copy' operation.

        Parameters
        ----------
        subExposure : `afwImage.Exposure`
            Input `subExposure` upon which to operate
        expandedSubExp : `afwImage.Exposure`
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
    mapper = pexConfig.ConfigurableField(
        doc="Mapper subtask to run on each subimage",
        target=AddAmountImageMapper,
    )


class GetMeanImageMapper(ImageMapper):
    """ImageMapper subtask that computes and returns the mean value of the
    input sub-exposure
    """
    ConfigClass = AddAmountImageMapperConfig  # Doesn't need its own config
    _DefaultName = "ip_diffim_GetMeanImageMapper"

    def run(self, subExposure, expandedSubExp, fullBBox, **kwargs):
        """Compute the mean of the given `subExposure`

        Parameters
        ----------
        subExposure : `afwImage.Exposure`
            Input `subExposure` upon which to operate
        expandedSubExp : `afwImage.Exposure`
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
    mapper = pexConfig.ConfigurableField(
        doc="Mapper subtask to run on each subimage",
        target=GetMeanImageMapper,
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
        self.exposure.setWcs(makeWcs())  # required for PSF construction via CoaddPsf

    def testCopySumNoOverlaps(self):
        self._testCopySumNoOverlaps(reduceOp='copy', withNaNs=False)
        self._testCopySumNoOverlaps(reduceOp='copy', withNaNs=True)
        self._testCopySumNoOverlaps(reduceOp='sum', withNaNs=False)
        self._testCopySumNoOverlaps(reduceOp='sum', withNaNs=True)

    def _testCopySumNoOverlaps(self, reduceOp='copy', withNaNs=False):
        """Test sample grid task that adds 5.0 to input image and uses
        `reduceOperation = 'copy'`. Optionally add NaNs to subimages.
        """
        config = AddAmountImageMapReduceConfig()
        task = ImageMapReduceTask(config)
        config.mapper.addAmount = 5.
        config.reducer.reduceOperation = reduceOp
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
        else:  # We don't construct a new PSF if reduceOperation == 'copy'.
            self._testCoaddPsf(newExp)

    def testAverageWithOverlaps(self):
        self._testAverageWithOverlaps(withNaNs=False)
        self._testAverageWithOverlaps(withNaNs=True)

    def _testAverageWithOverlaps(self, withNaNs=False):
        """Test sample grid task that adds 5.0 to input image and uses
        'average' `reduceOperation`. Optionally add NaNs to subimages.
        """
        config = AddAmountImageMapReduceConfig()
        config.gridStepX = config.gridStepY = 8.
        config.reducer.reduceOperation = 'average'
        task = ImageMapReduceTask(config)
        config.mapper.addAmount = 5.
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
        self._testCoaddPsf(newExp)

    def _testCoaddPsf(self, newExposure):
        """Test that the new CoaddPsf of the `newExposure` returns PSF images
        ~identical to the input PSF of `self.exposure` across a grid
        covering the entire exposure bounding box.
        """
        origPsf = self.exposure.getPsf()
        newPsf = newExposure.getPsf()
        self.assertTrue(isinstance(newPsf, measAlg.CoaddPsf))
        extentX = int(self.exposure.getWidth()*0.05)
        extentY = int(self.exposure.getHeight()*0.05)
        for x in np.linspace(extentX, self.exposure.getWidth()-extentX, 10):
            for y in np.linspace(extentY, self.exposure.getHeight()-extentY, 10):
                point = geom.Point2D(np.rint(x), np.rint(y))
                oPsf = origPsf.computeImage(point).getArray()
                nPsf = newPsf.computeImage(point).getArray()
                if oPsf.shape[0] < nPsf.shape[0]:  # sometimes CoaddPsf does this.
                    oPsf = np.pad(oPsf, ((1, 1), (0, 0)), mode='constant')
                elif oPsf.shape[0] > nPsf.shape[0]:
                    nPsf = np.pad(nPsf, ((1, 1), (0, 0)), mode='constant')
                if oPsf.shape[1] < nPsf.shape[1]:  # sometimes CoaddPsf does this.
                    oPsf = np.pad(oPsf, ((0, 0), (1, 1)), mode='constant')
                elif oPsf.shape[1] > nPsf.shape[1]:
                    nPsf = np.pad(nPsf, ((0, 0), (1, 1)), mode='constant')
                # pixel-wise comparison -- pretty stringent
                self.assertFloatsAlmostEqual(oPsf, nPsf, atol=1e-4, msg='Failed on Psf')

                origMmts = np.array(getPsfSecondMoments(oPsf))
                newMmts = np.array(getPsfSecondMoments(nPsf))
                self.assertFloatsAlmostEqual(origMmts, newMmts, atol=1e-4, msg='Failed on Psf')

    def testAverageVersusCopy(self):
        self._testAverageVersusCopy(withNaNs=False)
        self._testAverageVersusCopy(withNaNs=True)

    def _testAverageVersusCopy(self, withNaNs=False):
        """Re-run `testExampleTaskNoOverlaps` and `testExampleTaskWithOverlaps`
        on a more complex image (with random noise). Ensure that the results are
        identical (within between 'copy' and 'average' reduceOperation.
        """
        exposure1 = self.exposure.clone()
        img = exposure1.getMaskedImage().getImage()
        afwMath.randomGaussianImage(img, afwMath.Random())
        exposure2 = exposure1.clone()

        config = AddAmountImageMapReduceConfig()
        task = ImageMapReduceTask(config)
        config.mapper.addAmount = 5.
        newExp = task.run(exposure1, addNans=withNaNs).exposure
        newMI1 = newExp.getMaskedImage()

        config.gridStepX = config.gridStepY = 8.
        config.reducer.reduceOperation = 'average'
        task = ImageMapReduceTask(config)
        newExp = task.run(exposure2, addNans=withNaNs).exposure
        newMI2 = newExp.getMaskedImage()

        newMA1 = newMI1.getImage().getArray()
        isnan = np.isnan(newMA1)
        if not withNaNs:
            self.assertEqual(np.sum(isnan), 0)
        newMA2 = newMI2.getImage().getArray()

        # Because the average uses a float accumulator, we can have differences, set a tolerance.
        # Turns out (in practice for this test), only 7 pixels seem to have a small difference.
        self.assertFloatsAlmostEqual(newMA1[~isnan], newMA2[~isnan], rtol=1e-7)

    def testMean(self):
        """Test sample grid task that returns the mean of the subimages and uses
        'none' `reduceOperation`.
        """
        config = GetMeanImageMapReduceConfig()
        config.reducer.reduceOperation = 'none'
        task = ImageMapReduceTask(config)
        testExposure = self.exposure.clone()
        testExposure.getMaskedImage().set(1.234)
        subMeans = task.run(testExposure).result
        subMeans = [x.subExposure for x in subMeans]

        self.assertEqual(len(subMeans), len(task.boxes0))
        firstPixel = testExposure.getMaskedImage().getImage().getArray()[0, 0]
        self.assertFloatsAlmostEqual(np.array(subMeans), firstPixel)

    def testCellCentroids(self):
        """Test sample grid task which is provided a set of `cellCentroids` and
        returns the mean of the subimages surrounding those centroids using 'none'
        for `reduceOperation`.
        """
        config = GetMeanImageMapReduceConfig()
        config.gridStepX = config.gridStepY = 8.
        config.reducer.reduceOperation = 'none'
        config.cellCentroidsX = [i for i in np.linspace(0, 128, 50)]
        config.cellCentroidsY = config.cellCentroidsX
        task = ImageMapReduceTask(config)
        testExposure = self.exposure.clone()
        testExposure.getMaskedImage().set(1.234)
        subMeans = task.run(testExposure).result
        subMeans = [x.subExposure for x in subMeans]

        self.assertEqual(len(subMeans), len(config.cellCentroidsX))
        firstPixel = testExposure.getMaskedImage().getImage().getArray()[0, 0]
        self.assertFloatsAlmostEqual(np.array(subMeans), firstPixel)

    def testCellCentroidsWrongLength(self):
        """Test sample grid task which is provided a set of `cellCentroids` and
        returns the mean of the subimages surrounding those centroids using 'none'
        for `reduceOperation`. In this case, we ensure that len(task.boxes0) !=
        len(task.boxes1) and check for ValueError.
        """
        config = GetMeanImageMapReduceConfig()
        config.reducer.reduceOperation = 'none'
        config.cellCentroidsX = [i for i in np.linspace(0, 128, 50)]
        config.cellCentroidsY = [i for i in np.linspace(0, 128, 50)]
        task = ImageMapReduceTask(config)
        task._generateGrid(self.exposure)
        del task.boxes0[-1]  # remove the last box
        with self.assertRaises(ValueError):
            task.run(self.exposure)

    def testMasks(self):
        """Test the mask for an exposure produced by a sample grid task
        where we provide a set of `cellCentroids` and thus should have
        many invalid pixels.
        """
        config = AddAmountImageMapReduceConfig()
        config.gridStepX = config.gridStepY = 8.
        config.cellCentroidsX = [i for i in np.linspace(0, 128, 50)]
        config.cellCentroidsY = config.cellCentroidsX
        config.reducer.reduceOperation = 'average'
        task = ImageMapReduceTask(config)
        config.mapper.addAmount = 5.
        newExp = task.run(self.exposure).exposure
        newMI = newExp.getMaskedImage()
        newArr = newMI.getImage().getArray()
        mi = self.exposure.getMaskedImage()
        isnan = np.isnan(newArr)
        self.assertGreater(np.sum(isnan), 1000)

        mi = self.exposure.getMaskedImage().getImage().getArray()
        self.assertFloatsAlmostEqual(mi[~isnan], newArr[~isnan] - 5.)

        mask = newMI.getMask()  # Now check the mask
        self.assertGreater(mask.getMaskPlane('INVALID_MAPREDUCE'), 0)
        maskBit = mask.getPlaneBitMask('INVALID_MAPREDUCE')
        nMasked = np.sum(np.bitwise_and(mask.getArray(), maskBit) != 0)
        self.assertGreater(nMasked, 1000)
        self.assertEqual(np.sum(np.isnan(newArr)), nMasked)

    def testNotNoneReduceWithNonExposureMapper(self):
        """Test that a combination of a mapper that returns a non-exposure
        cannot work correctly with a reducer with reduceOperation='none'.
        Should raise a TypeError.
        """
        config = GetMeanImageMapReduceConfig()  # mapper returns a float (mean)
        config.gridStepX = config.gridStepY = 8.
        config.reducer.reduceOperation = 'average'  # not 'none'!
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
                                config.reducer.reduceOperation = reduceOp
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
        config.mapper.addAmount = expectedVal
        lenBoxes = [0, 0]
        for scaleByFwhm in (True, False):
            config.scaleByFwhm = scaleByFwhm
            if scaleByFwhm:
                config.gridStepX = float(gstepx)
                config.cellSizeX = float(gsizex)
                config.gridStepY = float(gstepy)
                config.cellSizeY = float(gsizey)
            else:  # otherwise the grid is too fine and elements too small.
                config.gridStepX = gstepx * 3.
                config.cellSizeX = gsizex * 3.
                config.gridStepY = gstepy * 3.
                config.cellSizeY = gsizey * 3.
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
