from __future__ import absolute_import, division, print_function
#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
import lsst.afw.table as afwTable
import lsst.meas.base as measBase
from lsst.ip.diffim.dipoleFitTask import (DipoleFitAlgorithm, DipoleFitTask)
import lsst.ip.diffim.utils as ipUtils

"""This file contains three unit tests of the DipoelFitAlgorithm and its related tasks and plugins.
Each test generates a fake image with two synthetic dipoles as input data.

DipoleFitTest.testDipoleAlgorithm - tests the DipoleFitAlgorithm directly and ensures that the recovered
   fluxes and centroids are close to the input values.
DipoleFitTest.testDipoleTask - tests the DipoleFitTask on identical synthetic data, after running
   dipole detection on the image.
DipoleFitTest.testDipoleEdge - ensures correct handling of dipoles too close to the edge of the image.
"""


class DipoleFitTestGlobalParams(object):
    """Class to initialize and store global parameters used by all tests below.

    Attributes:
    @var display: Display (plot) the output dipole thumbnails (matplotlib)
    @var verbose: be verbose during fitting
    @var w: width (pixels) of test generated exposures
    @var h: height (pixels) of test generated exposures
    @var xc: x coordinate (pixels) of center(s) of input dipole(s)
    @var yc: y coordinate (pixels) of center(s) of input dipole(s)
    @var flux: flux(es) of input dipole(s)
    @var gradientParams: tuple with three parameters for linear background gradient
    @var offsets: pixel coordinates between lobes of dipoles
    """

    def __init__(self, xc=None, yc=None):
        """Initialize the parameters.

        @param xc iterable x coordinate (pixels) of center(s) of input dipole(s)
        @param yc iterable y coordinate (pixels) of center(s) of input dipole(s)
        """
        np.random.seed(666)
        self.display = False
        self.verbose = False
        self.w, self.h = 100, 100  # size of image

        self.xc = xc  # xcenters of two dipoles in image
        if xc is None:
            xc = self.xc = [65.3, 24.2]
        self.yc = yc  # ycenters of two dipoles
        if yc is None:
            yc = self.yc = [38.6, 78.5]
        self.flux = [2500., 2345.]  # fluxes of pos/neg lobes
        self.gradientParams = [10., 3., 5.]

        self.offsets = np.array([-2., 2.])  # pixel coord offsets between lobes of dipoles
        self.testImage = ipUtils.DipoleTestImage(
            xcenPos=self.xc + self.offsets,
            ycenPos=self.yc + self.offsets,
            xcenNeg=self.xc - self.offsets,
            ycenNeg=self.yc - self.offsets,
            flux=self.flux, fluxNeg=self.flux,
            noise=2.,  # Note the input noise - this affects the relative tolerances used.
            gradientParams=self.gradientParams)

        self.catalog = self.testImage.detectDipoleSources()
        self.rtol = 0.01  # This is okay given the default noise of 2. set above.


# First, test the algorithm itself (fitDipole()):
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
        params = DipoleFitTestGlobalParams()
        rtol = params.rtol

        for s in params.catalog:
            fp = s.getFootprint()
            self.assertTrue(len(fp.getPeaks()) == 2)

        offsets = params.offsets
        testImage = params.testImage
        for i, s in enumerate(params.catalog):
            alg = DipoleFitAlgorithm(testImage.diffim, testImage.posImage, testImage.negImage)
            result = alg.fitDipole(
                s, rel_weight=1., separateNegParams=False,
                verbose=params.verbose, display=params.display)

            self.assertClose((result.psfFitPosFlux + abs(result.psfFitNegFlux))/2.,
                             params.flux[i], rtol=rtol)
            self.assertClose(result.psfFitPosCentroidX, params.xc[i] + offsets[i], rtol=rtol)
            self.assertClose(result.psfFitPosCentroidY, params.yc[i] + offsets[i], rtol=rtol)
            self.assertClose(result.psfFitNegCentroidX, params.xc[i] - offsets[i], rtol=rtol)
            self.assertClose(result.psfFitNegCentroidY, params.yc[i] - offsets[i], rtol=rtol)

    def _runDetection(self, testImage):
        """Run 'diaSource' detection on the diffim, including merging of
        positive and negative sources.
        """

        # Create the various tasks and schema -- avoid code reuse.
        detectTask, schema = testImage.detectDipoleSources(doMerge=False)

        measureConfig = measBase.SingleFrameMeasurementConfig()

        measureConfig.slots.calibFlux = None
        measureConfig.slots.modelFlux = None
        measureConfig.slots.instFlux = None
        measureConfig.slots.shape = None
        measureConfig.slots.centroid = "ip_diffim_NaiveDipoleCentroid"
        measureConfig.doReplaceWithNoise = False

        measureConfig.plugins.names = ["base_CircularApertureFlux",
                                       "base_PixelFlags",
                                       "base_SkyCoord",
                                       "base_PsfFlux",
                                       "ip_diffim_NaiveDipoleCentroid",
                                       "ip_diffim_NaiveDipoleFlux",
                                       "ip_diffim_PsfDipoleFlux"]

        # Disable aperture correction, which requires having an ApCorrMap attached to
        # the Exposure (it'll warn if it's not present and we don't explicitly disable it).
        measureConfig.doApplyApCorr = "no"

        # Here is where we make the dipole fitting task. It can run the other measurements as well.
        # This is an example of how to pass it a custom config.
        measureTask = DipoleFitTask(config=measureConfig, schema=schema)

        table = afwTable.SourceTable.make(schema)
        detectResult = detectTask.run(table, testImage.diffim)
        # catalog = detectResult.sources
        # deblendTask.run(self.dipole, catalog, psf=self.dipole.getPsf())

        fpSet = detectResult.fpSets.positive
        fpSet.merge(detectResult.fpSets.negative, 2, 2, False)
        sources = afwTable.SourceCatalog(table)
        fpSet.makeSources(sources)

        measureTask.run(sources, testImage.diffim, testImage.posImage, testImage.negImage)
        return sources

    def testDipoleTask(self):
        """Test the dipole fitting singleFramePlugin.

        Test that the resulting fluxes/centroids are entered into the
        correct slots of the catalog, and have values that are very
        close to the input values for both dipoles in the image.

        Also test that the resulting fluxes are close to those
        generated by the existing ip_diffim_DipoleMeasurement task
        (PsfDipoleFit).
        """
        params = DipoleFitTestGlobalParams()
        rtol = params.rtol
        sources = self._runDetection(params.testImage)

        offsets = params.offsets
        for i, r1 in enumerate(sources):
            result = r1.extract("ip_diffim_DipoleFit*")
            self.assertClose((result['ip_diffim_DipoleFit_pos_flux'] +
                              abs(result['ip_diffim_DipoleFit_neg_flux']))/2.,
                             params.flux[i], rtol=rtol)
            self.assertClose(result['ip_diffim_DipoleFit_pos_centroid_x'],
                             params.xc[i] + offsets[i], rtol=rtol)
            self.assertClose(result['ip_diffim_DipoleFit_pos_centroid_y'],
                             params.yc[i] + offsets[i], rtol=rtol)
            self.assertClose(result['ip_diffim_DipoleFit_neg_centroid_x'],
                             params.xc[i] - offsets[i], rtol=rtol)
            self.assertClose(result['ip_diffim_DipoleFit_neg_centroid_y'],
                             params.yc[i] - offsets[i], rtol=rtol)
            # Note this is dependent on the noise (variance) being realistic in the image.
            # otherwise it throws off the chi2 estimate, which is used for classification:
            self.assertTrue(result['ip_diffim_DipoleFit_flag_classification'])

            # compare to the original ip_diffim_PsfDipoleFlux measurements
            result2 = r1.extract("ip_diffim_PsfDipoleFlux*")
            self.assertClose((result['ip_diffim_DipoleFit_pos_flux'] +
                              abs(result['ip_diffim_DipoleFit_neg_flux']))/2.,
                             (result2['ip_diffim_PsfDipoleFlux_pos_flux'] +
                              abs(result2['ip_diffim_PsfDipoleFlux_neg_flux']))/2.,
                             rtol=rtol)
            self.assertClose(result['ip_diffim_DipoleFit_pos_centroid_x'],
                             result2['ip_diffim_PsfDipoleFlux_pos_centroid_x'],
                             rtol=rtol)
            self.assertClose(result['ip_diffim_DipoleFit_pos_centroid_y'],
                             result2['ip_diffim_PsfDipoleFlux_pos_centroid_y'],
                             rtol=rtol)
            self.assertClose(result['ip_diffim_DipoleFit_neg_centroid_x'],
                             result2['ip_diffim_PsfDipoleFlux_neg_centroid_x'],
                             rtol=rtol)
            self.assertClose(result['ip_diffim_DipoleFit_neg_centroid_y'],
                             result2['ip_diffim_PsfDipoleFlux_neg_centroid_y'],
                             rtol=rtol)

            if params.display:
                params.testImage.displayCutouts(r1)
        if params.display:
            plt = ipUtils.importMatplotlib()
            if not plt:
                return result
            plt.show()

        return result

    def testDipoleEdge(self):
        """Test the edge case for dipole fitting singleFramePlugin.

        Test that the dipoles which are too close to the edge are
        indicated as so in the catalog.
        """

        params = DipoleFitTestGlobalParams(xc=[5.3, 2.2], yc=[2.6, 98.5])
        sources = self._runDetection(params.testImage)

        for i, r1 in enumerate(sources):
            result = r1.extract("ip_diffim_DipoleFit*")
            self.assertTrue(result.get("ip_diffim_DipoleFit_flag"))


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(DipoleFitTest)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)

