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

"""!Tests of the DipoleFitAlgorithm and its related tasks and plugins.
Each test generates a fake image with two synthetic dipoles as input data.
"""

import unittest

import numpy as np

import lsst.utils.tests
import lsst.afw.table as afwTable
import lsst.meas.base as measBase
from lsst.ip.diffim.dipoleFitTask import (DipoleFitAlgorithm)
import lsst.ip.diffim.utils as ipUtils


class DipoleTestImage(object):
    """!Class to initialize test dipole image used by all tests below.

    @var display: Display (plot) the output dipole thumbnails (matplotlib)
    @var verbose: be verbose during fitting
    @var xc: x coordinate (pixels) of center(s) of input dipole(s)
    @var yc: y coordinate (pixels) of center(s) of input dipole(s)
    @var flux: flux(es) of input dipole(s)
    @var gradientParams: tuple with three parameters for linear background gradient
    @var offsets: pixel coordinates between lobes of dipoles

    Also stores all parameters used to generate the test image (to compare to fitting results).
    """

    def __init__(self, xc=None, yc=None, flux=None, offsets=None, gradientParams=None):
        """!Store the parameters, create the test image and run detection on it.

        @param xc iterable x coordinate (pixels) of center(s) of input dipole(s)
        @param yc iterable y coordinate (pixels) of center(s) of input dipole(s)
        @param offsets iterable pixel coord offsets between lobes of dipole(s)
        @param flux iterable fluxes of pos/neg lobes of dipole(s)
        @param gradientParams iterable three parameters for linear background gradient
        """
        np.random.seed(666)
        self.display = False  # Display (plot) the output dipole thumbnails (matplotlib)
        self.verbose = False  # be verbose during fitting

        self.xc = xc if xc is not None else [65.3, 24.2]
        self.yc = yc if yc is not None else [38.6, 78.5]
        self.offsets = offsets if offsets is not None else np.array([-2., 2.])
        self.flux = flux if flux is not None else [2500., 2345.]
        self.gradientParams = gradientParams if gradientParams is not None else [10., 3., 5.]

        # The default tolerance for comparisons of fitted parameters with input values.
        # Given the noise in the input images (default noise value of 2.), this is a
        # useful test of algorithm robustness, and will guard against future regressions.
        self.rtol = 0.01

        self.generateTestImage()

    def generateTestImage(self):
        self.testImage = ipUtils.DipoleTestImage(
            w=100, h=100,
            xcenPos=self.xc + self.offsets,
            ycenPos=self.yc + self.offsets,
            xcenNeg=self.xc - self.offsets,
            ycenNeg=self.yc - self.offsets,
            flux=self.flux, fluxNeg=self.flux,
            noise=2.,  # Note the input noise - this affects the relative tolerances used.
            gradientParams=self.gradientParams)


class DipoleFitTest(lsst.utils.tests.TestCase):
    """!A test case for separately testing the dipole fit algorithm
    directly, and the single frame measurement.

    In each test, create a simulated diffim with two dipoles, noise,
    and a linear background gradient in the pre-sub images then
    compare the input fluxes/centroids with the fitted results.
    """

    def testDipoleAlgorithm(self):
        """!Test the dipole fitting algorithm directly (fitDipole()).

        Test that the resulting fluxes/centroids are very close to the
        input values for both dipoles in the image.
        """
        params = DipoleTestImage()
        catalog = params.testImage.detectDipoleSources()

        for s in catalog:
            fp = s.getFootprint()
            self.assertTrue(len(fp.getPeaks()) == 2)

        rtol = params.rtol
        offsets = params.offsets
        testImage = params.testImage
        for i, s in enumerate(catalog):
            alg = DipoleFitAlgorithm(testImage.diffim, testImage.posImage, testImage.negImage)
            result, _ = alg.fitDipole(
                s, rel_weight=0.5, separateNegParams=False,
                verbose=params.verbose, display=params.display)

            self.assertClose((result.posFlux + abs(result.negFlux))/2.,
                             params.flux[i], rtol=rtol)
            self.assertClose(result.posCentroidX, params.xc[i] + offsets[i], rtol=rtol)
            self.assertClose(result.posCentroidY, params.yc[i] + offsets[i], rtol=rtol)
            self.assertClose(result.negCentroidX, params.xc[i] - offsets[i], rtol=rtol)
            self.assertClose(result.negCentroidY, params.yc[i] - offsets[i], rtol=rtol)


def suite():
    """!Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(DipoleFitTest)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """!Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)

