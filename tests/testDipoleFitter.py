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
"""
Tests for the DipoleFitAlgorithm and its related tasks and plugins.
Each test generates a fake image with two synthetic dipoles as input data.
"""

import unittest

import numpy as np

import lsst.utils.tests
import lsst.afw.table as afwTable
import lsst.meas.base as measBase
from lsst.ip.diffim.dipoleFitTask import (DipoleFitAlgorithm, DipoleFitTask)
import lsst.ip.diffim.utils as ipUtils


class DipoleFitBase(object):
    """
    Derive from this, TestCase.
    Set self.xc, self.yc in setUp before calling super(), to define the dipole locations.
    """
    def setUp(self):
        # Display (plot) the output dipole thumbnails (matplotlib)
        self.display = False
        # be verbose during fitting
        self.verbose = False
        np.random.seed(666)

        # x and y centers of two dipoles the image
        self.flux = [2500., 2345.]  # fluxes of pos/neg lobes
        self.gradientParams = [10., 3., 5.]

        self.offsets = np.array([-2., 2.])  # pixel coord offsets between lobes of dipoles
        self.rtol = 0.01  # This is okay given the default noise of 2. set above.

        self.testImage = ipUtils.DipoleTestImage(
            xcenPos=self.xc + self.offsets,
            ycenPos=self.yc + self.offsets,
            xcenNeg=self.xc - self.offsets,
            ycenNeg=self.yc - self.offsets,
            flux=self.flux, fluxNeg=self.flux,
            noise=2.,  # Note the input noise - this affects the relative tolerances used.
            gradientParams=self.gradientParams)

        self.catalog = self.testImage.detectDipoleSources()

    def _runDetection(self):
        """Run 'diaSource' detection on the diffim, including merging of
        positive and negative sources.
        """

        # Create the various tasks and schema -- avoid code reuse.
        detectTask, schema = self.testImage.detectDipoleSources(doMerge=False)

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
        detectResult = detectTask.run(table, self.testImage.diffim)
        # catalog = detectResult.sources
        # deblendTask.run(self.dipole, catalog, psf=self.dipole.getPsf())

        fpSet = detectResult.fpSets.positive
        fpSet.merge(detectResult.fpSets.negative, 2, 2, False)
        sources = afwTable.SourceCatalog(table)
        fpSet.makeSources(sources)

        measureTask.run(sources, self.testImage.diffim, self.testImage.posImage, self.testImage.negImage)
        return sources


class DipoleFitTest(DipoleFitBase, lsst.utils.tests.TestCase):
    """A test case for separately testing the dipole fit algorithm
    directly, and the single frame measurement.

    In each test, create a simulated diffim with two dipoles, noise,
    and a linear background gradient in the pre-sub images then
    compare the input fluxes/centroids with the fitted results.
    """
    def setUp(self):
        self.xc = [65.3, 24.2]
        self.yc = [38.6, 78.5]
        super(DipoleFitTest, self).setUp()

    def testDipoleAlgorithm(self):
        """Test the dipole fitting algorithm directly (fitDipole()).

        Test that the resulting fluxes/centroids are very close to the
        input values for both dipoles in the image.
        """
        rtol = self.rtol

        for s in self.catalog:
            fp = s.getFootprint()
            self.assertTrue(len(fp.getPeaks()) == 2)

        offsets = self.offsets
        testImage = self.testImage
        for i, s in enumerate(self.catalog):
            alg = DipoleFitAlgorithm(testImage.diffim, testImage.posImage, testImage.negImage)
            result = alg.fitDipole(
                s, rel_weight=1., separateNegParams=False,
                verbose=self.verbose, display=self.display)

            self.assertClose((result.psfFitPosFlux + abs(result.psfFitNegFlux))/2.,
                             self.flux[i], rtol=rtol)
            self.assertClose(result.psfFitPosCentroidX, self.xc[i] + offsets[i], rtol=rtol)
            self.assertClose(result.psfFitPosCentroidY, self.yc[i] + offsets[i], rtol=rtol)
            self.assertClose(result.psfFitNegCentroidX, self.xc[i] - offsets[i], rtol=rtol)
            self.assertClose(result.psfFitNegCentroidY, self.yc[i] - offsets[i], rtol=rtol)

    def testDipoleTask(self):
        """Test the dipole fitting singleFramePlugin.

        Test that the resulting fluxes/centroids are entered into the
        correct slots of the catalog, and have values that are very
        close to the input values for both dipoles in the image.

        Also test that the resulting fluxes are close to those
        generated by the existing ip_diffim_DipoleMeasurement task
        (PsfDipoleFit).
        """
        rtol = self.rtol
        sources = self._runDetection()

        offsets = self.offsets
        for i, r1 in enumerate(sources):
            result = r1.extract("ip_diffim_DipoleFit*")
            self.assertClose((result['ip_diffim_DipoleFit_pos_flux'] +
                              abs(result['ip_diffim_DipoleFit_neg_flux']))/2.,
                             self.flux[i], rtol=rtol)
            self.assertClose(result['ip_diffim_DipoleFit_pos_centroid_x'],
                             self.xc[i] + offsets[i], rtol=rtol)
            self.assertClose(result['ip_diffim_DipoleFit_pos_centroid_y'],
                             self.yc[i] + offsets[i], rtol=rtol)
            self.assertClose(result['ip_diffim_DipoleFit_neg_centroid_x'],
                             self.xc[i] - offsets[i], rtol=rtol)
            self.assertClose(result['ip_diffim_DipoleFit_neg_centroid_y'],
                             self.yc[i] - offsets[i], rtol=rtol)
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

            if self.display:
                self.testImage.displayCutouts(r1)
        if self.display:
            plt = ipUtils.importMatplotlib()
            if not plt:
                return result
            plt.show()

        return result


class DipoleEdgeTest(DipoleFitBase, lsst.utils.tests.TestCase):
    """Test for dipoles on the edge of images."""

    def setUp(self):
        self.xc = [5.3, 2.2]
        self.yc = [2.6, 98.5]
        super(DipoleEdgeTest, self).setUp()

    def testDipoleEdge(self):
        """
        Test that the dipoles which are too close to the edge are
        indicated as so in the catalog.
        """

        sources = self._runDetection()
        for i, r1 in enumerate(sources):
            result = r1.extract("ip_diffim_DipoleFit*")
            self.assertTrue(result.get("ip_diffim_DipoleFit_flag"))


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(DipoleFitTest)
    suites += unittest.makeSuite(DipoleEdgeTest)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)

