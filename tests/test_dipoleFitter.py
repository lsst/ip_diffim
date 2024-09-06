#
# LSST Data Management System
# Copyright 2008-2017 AURA/LSST.
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

"""Tests of the DipoleFitAlgorithm and its related tasks and plugins.

Each test generates a fake image with two synthetic dipoles as input data.
"""
import unittest

import numpy as np

import lsst.utils.tests
import lsst.afw.table as afwTable
from lsst.ip.diffim.dipoleFitTask import (DipoleFitAlgorithm, DipoleFitTask)
import utils


class DipoleTestImage:
    """Create a test dipole image and store the parameters used to make it,
    for comparison with the fitted results.

    Parameters
    ----------
    xc, yc : `list` [`float`]
        x, y coordinate (pixels) of center(s) of input dipole(s).
    flux: `list` [`float`]
        Flux(es) of input dipole(s).
    gradientParams : `tuple`
        Tuple with three parameters for linear background gradient.
    offsets : `list` [`float`]
        Pixel coordinates between lobes of dipoles.
    """

    def __init__(self, xc=None, yc=None, flux=None, offsets=None, gradientParams=None, edgeWidth=None):
        self.xc = xc if xc is not None else [65.3, 24.2]
        self.yc = yc if yc is not None else [38.6, 78.5]
        self.offsets = offsets if offsets is not None else np.array([-2., 2.])
        self.flux = flux if flux is not None else [2500., 2345.]
        self.gradientParams = gradientParams if gradientParams is not None else [10., 3., 5.]
        self.edgeWidth = edgeWidth if edgeWidth is not None else 8

        # The default tolerance for comparisons of fitted parameters with input values.
        # Given the noise in the input images (default noise value of 2.), this is a
        # useful test of algorithm robustness, and will guard against future regressions.
        self.rtol = 0.01

        self.generateTestImage()

    def generateTestImage(self):
        self.testImage = utils.DipoleTestImage(
            w=100, h=100,
            xcenPos=self.xc + self.offsets,
            ycenPos=self.yc + self.offsets,
            xcenNeg=self.xc - self.offsets,
            ycenNeg=self.yc - self.offsets,
            flux=self.flux, fluxNeg=self.flux,
            noise=2.,  # Note the input noise - this affects the relative tolerances used.
            gradientParams=self.gradientParams,
            edgeWidth=self.edgeWidth)


class DipoleFitTest(lsst.utils.tests.TestCase):
    """A test case for separately testing the dipole fit algorithm
    directly, and the single frame measurement.

    In each test, create a simulated diffim with two dipoles, noise,
    and a linear background gradient in the pre-sub images then
    compare the input fluxes/centroids with the fitted results.
    """

    def testDipoleAlgorithm(self):
        """Test the dipole fitting algorithm directly (fitDipole()).

        Test that the resulting fluxes/centroids are very close to the
        input values for both dipoles in the image.
        """
        # Display (plot) the output dipole thumbnails with matplotlib.
        display = False
        # Be verbose during fitting, including the lmfit internal details.
        verbose = False

        dipoleTestImage = DipoleTestImage()
        catalog = dipoleTestImage.testImage.detectDipoleSources(minBinSize=32)

        for s in catalog:
            fp = s.getFootprint()
            self.assertTrue(len(fp.getPeaks()) == 2)

        rtol = dipoleTestImage.rtol
        offsets = dipoleTestImage.offsets
        testImage = dipoleTestImage.testImage
        for i, s in enumerate(catalog):
            alg = DipoleFitAlgorithm(testImage.diffim, testImage.posImage, testImage.negImage)
            result, _ = alg.fitDipole(
                s, rel_weight=0.5, separateNegParams=False,
                verbose=verbose, display=display)

            self.assertFloatsAlmostEqual((result.posFlux + abs(result.negFlux))/2.,
                                         dipoleTestImage.flux[i], rtol=rtol)
            self.assertFloatsAlmostEqual(result.posCentroidX, dipoleTestImage.xc[i] + offsets[i], rtol=rtol)
            self.assertFloatsAlmostEqual(result.posCentroidY, dipoleTestImage.yc[i] + offsets[i], rtol=rtol)
            self.assertFloatsAlmostEqual(result.negCentroidX, dipoleTestImage.xc[i] - offsets[i], rtol=rtol)
            self.assertFloatsAlmostEqual(result.negCentroidY, dipoleTestImage.yc[i] - offsets[i], rtol=rtol)

    def _runDetection(self, dipoleTestImage, maxFootprintArea=None):
        """Run 'diaSource' detection on the diffim, including merging of
        positive and negative sources.

        Then run DipoleFitTask on the image and return the resulting catalog.
        """
        # Create the various tasks and schema -- avoid code reuse.
        testImage = dipoleTestImage.testImage
        detectTask, schema = testImage.detectDipoleSources(doMerge=False, minBinSize=32)

        config = DipoleFitTask.ConfigClass()
        # Also run the older C++ DipoleFlux algorithm for comparison purposes.
        config.plugins.names |= ["ip_diffim_PsfDipoleFlux"]
        if maxFootprintArea:
            config.plugins["ip_diffim_DipoleFit"].maxFootprintArea = maxFootprintArea
        measureTask = DipoleFitTask(schema=schema, config=config)

        table = afwTable.SourceTable.make(schema)
        detectResult = detectTask.run(table, testImage.diffim)
        fpSet = detectResult.positive
        fpSet.merge(detectResult.negative, 2, 2, False)
        sources = afwTable.SourceCatalog(table)
        fpSet.makeSources(sources)

        measureTask.run(sources, testImage.diffim, testImage.posImage, testImage.negImage)
        return sources

    def _checkTaskOutput(self, dipoleTestImage, sources, rtol=None):
        """Compare the fluxes/centroids in `sources` are entered
        into the correct slots of the catalog, and have values that
        are very close to the input values for both dipoles in the
        image.

        Also test that the resulting fluxes are close to those
        generated by the existing ip_diffim_DipoleMeasurement task
        (PsfDipoleFit).
        """

        if rtol is None:
            rtol = dipoleTestImage.rtol
        offsets = dipoleTestImage.offsets
        for i, r1 in enumerate(sources):
            result = r1.extract("ip_diffim_DipoleFit*")
            self.assertFloatsAlmostEqual((result['ip_diffim_DipoleFit_pos_instFlux']
                                          + abs(result['ip_diffim_DipoleFit_neg_instFlux']))/2.,
                                         dipoleTestImage.flux[i], rtol=rtol)
            self.assertFloatsAlmostEqual(result['ip_diffim_DipoleFit_pos_x'],
                                         dipoleTestImage.xc[i] + offsets[i], rtol=rtol)
            self.assertFloatsAlmostEqual(result['ip_diffim_DipoleFit_pos_y'],
                                         dipoleTestImage.yc[i] + offsets[i], rtol=rtol)
            self.assertFloatsAlmostEqual(result['ip_diffim_DipoleFit_neg_x'],
                                         dipoleTestImage.xc[i] - offsets[i], rtol=rtol)
            self.assertFloatsAlmostEqual(result['ip_diffim_DipoleFit_neg_y'],
                                         dipoleTestImage.yc[i] - offsets[i], rtol=rtol)
            # Note this is dependent on the noise (variance) being realistic in the image.
            # otherwise it throws off the chi2 estimate, which is used for classification:
            self.assertTrue(result['ip_diffim_DipoleFit_flag_classification'])

            # compare to the original ip_diffim_PsfDipoleFlux measurements
            result2 = r1.extract("ip_diffim_PsfDipoleFlux*")
            self.assertFloatsAlmostEqual((result['ip_diffim_DipoleFit_pos_instFlux']
                                          + abs(result['ip_diffim_DipoleFit_neg_instFlux']))/2.,
                                         (result2['ip_diffim_PsfDipoleFlux_pos_instFlux']
                                          + abs(result2['ip_diffim_PsfDipoleFlux_neg_instFlux']))/2.,
                                         rtol=rtol)
            self.assertFloatsAlmostEqual(result['ip_diffim_DipoleFit_pos_x'],
                                         result2['ip_diffim_PsfDipoleFlux_pos_centroid_x'],
                                         rtol=rtol)
            self.assertFloatsAlmostEqual(result['ip_diffim_DipoleFit_pos_y'],
                                         result2['ip_diffim_PsfDipoleFlux_pos_centroid_y'],
                                         rtol=rtol)
            self.assertFloatsAlmostEqual(result['ip_diffim_DipoleFit_neg_x'],
                                         result2['ip_diffim_PsfDipoleFlux_neg_centroid_x'],
                                         rtol=rtol)
            self.assertFloatsAlmostEqual(result['ip_diffim_DipoleFit_neg_y'],
                                         result2['ip_diffim_PsfDipoleFlux_neg_centroid_y'],
                                         rtol=rtol)

        return result

    def testDipoleTask(self):
        """Test the dipole fitting singleFramePlugin.

        Test that the resulting fluxes/centroids are entered into the
        correct slots of the catalog, and have values that are very
        close to the input values for both dipoles in the image.

        Also test that the resulting fluxes are close to those
        generated by the existing ip_diffim_DipoleMeasurement task
        (PsfDipoleFit).
        """
        dipoleTestImage = DipoleTestImage()
        sources = self._runDetection(dipoleTestImage)
        self._checkTaskOutput(dipoleTestImage, sources)

    def testDipoleTaskNoPosImage(self):
        """Test the dipole fitting singleFramePlugin in the case where no
        `posImage` is provided. It should be the same as above because
        `posImage` can be constructed from `diffim+negImage`.

        Test that the resulting fluxes/centroids are entered into the
        correct slots of the catalog, and have values that are very
        close to the input values for both dipoles in the image.

        Also test that the resulting fluxes are close to those
        generated by the existing ip_diffim_DipoleMeasurement task
        (PsfDipoleFit).
        """
        dipoleTestImage = DipoleTestImage()
        dipoleTestImage.testImage.posImage = None
        sources = self._runDetection(dipoleTestImage)
        self._checkTaskOutput(dipoleTestImage, sources)

    def testDipoleTaskNoNegImage(self):
        """Test the dipole fitting singleFramePlugin in the case where no
        `negImage` is provided. It should be the same as above because
        `negImage` can be constructed from `posImage-diffim`.

        Test that the resulting fluxes/centroids are entered into the
        correct slots of the catalog, and have values that are very
        close to the input values for both dipoles in the image.

        Also test that the resulting fluxes are close to those
        generated by the existing ip_diffim_DipoleMeasurement task
        (PsfDipoleFit).
        """
        dipoleTestImage = DipoleTestImage()
        dipoleTestImage.testImage.negImage = None
        sources = self._runDetection(dipoleTestImage)
        self._checkTaskOutput(dipoleTestImage, sources)

    def testDipoleTaskNoPreSubImages(self):
        """Test the dipole fitting singleFramePlugin in the case where no
        pre-subtraction data (`posImage` or `negImage`) are provided.
        In this case it just fits a dipole model to the diffim
        (dipole) image alone. Note that this test will only pass for
        widely-separated dipoles.

        Test that the resulting fluxes/centroids are entered into the
        correct slots of the catalog, and have values that are very
        close to the input values for both dipoles in the image.

        Also test that the resulting fluxes are close to those
        generated by the existing ip_diffim_DipoleMeasurement task
        (PsfDipoleFit).
        """
        dipoleTestImage = DipoleTestImage()
        dipoleTestImage.testImage.posImage = dipoleTestImage.testImage.negImage = None
        sources = self._runDetection(dipoleTestImage)
        self._checkTaskOutput(dipoleTestImage, sources)

    def testDipoleEdge(self):
        """Test the too-close-to-image-edge scenario for dipole fitting
        singleFramePlugin.

        Test that the dipoles which are too close to the edge are
        not detected.
        """

        # with no edge we should detect both dipole sources
        dipoleTestImage = DipoleTestImage(xc=[5.3, 4.8], yc=[4.6, 86.5], flux=[200, 210], edgeWidth=0)
        sources = self._runDetection(dipoleTestImage)
        self.assertEqual(len(sources), 2)

        # with a wide edge we should not detect any sources
        dipoleTestImage = DipoleTestImage(xc=[5.3, 4.8], yc=[4.6, 86.5], flux=[200, 210], edgeWidth=20)
        sources = self._runDetection(dipoleTestImage)
        self.assertEqual(len(sources), 0)

    def testDipoleFootprintTooLarge(self):
        """Test that the footprint area cut flags sources."""

        dipoleTestImage = DipoleTestImage()
        # This area is smaller than the area of the test sources (~750).
        sources = self._runDetection(dipoleTestImage, maxFootprintArea=500)

        self.assertTrue(np.all(sources["ip_diffim_DipoleFit_flag"]))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
