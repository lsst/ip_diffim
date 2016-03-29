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

from numpy import (random as np_random,
                   array as np_array,
                   mgrid as np_mgrid)

import lsst.utils.tests as lsst_tests
from lsst.afw.geom import (Box2I, Point2I, Point2D)
from lsst.meas.algorithms import (SourceDetectionConfig, SourceDetectionTask)
from lsst.afw.table import (SourceTable, SourceCatalog)
from lsst.ip.diffim.dipoleFitTask import (DipoleFitAlgorithm, DipoleFitConfig, DipoleFitTask, DipolePlotUtils)
from lsst.meas.base import SingleFrameMeasurementConfig

# Export DipoleTestUtils to expose fake image generating funcs
__all__ = ("DipoleTestUtils")


class DipoleFitTestGlobalParams(object):
    """
    Class to initialize and store global parameters used by all tests below.

    Attributes
    ----------
    display : boolean
        Display (plot) the output dipole thumbnails (matplotlib)
    verbose : boolean
        be verbose during fitting
    w : integer
        width (pixels) of test generated exposures
    h : integer
        height (pixels) of test generated exposures
    xc : array of floats
        x coordinate (pixels) of center(s) of input dipole(s)
    yc : array of floats
        y coordinate (pixels) of center(s) of input dipole(s)
    flux : array of floats
        flux(es) of input dipole(s)
    gradientParams : array of floats
        three parameters for linear background gradient
    offsets : array of floats
        pixel coordinates between lobes of dipoles
    """

    def __init__(self):
        """
        Initialize the parameters.
        """
        np_random.seed(666)
        self.display = False
        self.verbose = False
        self.w, self.h = 100, 100  # size of image

        self.xc = np_array([65.3, 24.2])  # xcenters of two dipoles in image
        self.yc = np_array([38.6, 78.5])  # ycenters of two dipoles
        self.flux = np_array([2500., 2345.])  # fluxes of pos/neg lobes
        self.gradientParams = np_array([10., 3., 5.])

        self.offsets = np_array([-2., 2.])  # pixel coord offsets between lobes of dipoles


# First, test the algorithm itself (fitDipole()):
# Create a simulated diffim (with dipoles) and a linear background gradient in the pre-sub images
#   then compare the input fluxes/centroids with the fitted results.
class DipoleFitAlgorithmTest(lsst_tests.TestCase):
    """
    A test case for dipole fit algorithm.

    Test the dipole fitting algorithm on two dipoles within an image with
    simulated noise.
    """
    def setUp(self):
        self.params = DipoleFitTestGlobalParams()

        offsets = self.params.offsets
        self.dipole, (self.posImage, self.posCatalog), (self.negImage, self.negCatalog) = \
            DipoleTestUtils.makeDipoleImage(
                xcenPos=self.params.xc + offsets,
                ycenPos=self.params.yc + offsets,
                xcenNeg=self.params.xc - offsets,
                ycenNeg=self.params.yc - offsets,
                flux=self.params.flux, fluxNeg=self.params.flux,
                gradientParams=self.params.gradientParams)

        self.catalog = DipoleTestUtils.detectDipoleSources(self.dipole)

    def tearDown(self):
        del self.dipole, self.posImage, self.negImage
        del self.catalog, self.posCatalog, self.negCatalog
        del self.params

    def testDipoleFitter(self):
        """
        Test the dipole fitting algorithm. Test that the resulting fluxes/centroids
        are very close to the input values for both dipoles in the image.
        """
        if self.params.verbose:
            for s in self.catalog:
                fp = s.getFootprint()
                print(fp.getBBox(), fp.getNpix())
                for pk in fp.getPeaks():
                    print('FOOTPRINT CENTER:', pk.getIy(), pk.getIx(), pk.getPeakValue())

        offsets = self.params.offsets
        for i, s in enumerate(self.catalog):
            alg = DipoleFitAlgorithm(self.dipole, self.posImage, self.negImage)
            result = alg.fitDipole(
                s, rel_weight=1., separateNegParams=False,
                verbose=self.params.verbose, display=self.params.display)

            self.assertClose((result.psfFitPosFlux + abs(result.psfFitNegFlux))/2.,
                             self.params.flux[i], rtol=0.02)
            self.assertClose(result.psfFitPosCentroidX, self.params.xc[i] + offsets[i], rtol=0.01)
            self.assertClose(result.psfFitPosCentroidY, self.params.yc[i] + offsets[i], rtol=0.01)
            self.assertClose(result.psfFitNegCentroidX, self.params.xc[i] - offsets[i], rtol=0.01)
            self.assertClose(result.psfFitNegCentroidY, self.params.yc[i] - offsets[i], rtol=0.01)


# Test the task in the same way as the algorithm:
# Also test that it correctly classifies the dipoles.
class DipoleFitTaskTest(DipoleFitAlgorithmTest):
    """ A test case for dipole fit task"""
    def setUp(self):
        DipoleFitAlgorithmTest.setUp(self)

    def tearDown(self):
        DipoleFitAlgorithmTest.tearDown(self)

    def runDetection(self):

        # Create the various tasks and schema -- avoid code reuse.
        detectTask, schema = DipoleTestUtils.detectDipoleSources(self.dipole, doMerge=False)

        measureConfig = SingleFrameMeasurementConfig()

        measureConfig.slots.modelFlux = "ip_diffim_DipoleFit"

        measureConfig.plugins.names |= ["base_CircularApertureFlux",
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
        dpFitConfig = DipoleFitConfig()
        measureTask = DipoleFitTask(config=measureConfig, schema=schema, dpFitConfig=dpFitConfig)

        table = SourceTable.make(schema)
        detectResult = detectTask.run(table, self.dipole)
        # catalog = detectResult.sources
        # deblendTask.run(self.dipole, catalog, psf=self.dipole.getPsf())

        fpSet = detectResult.fpSets.positive
        fpSet.merge(detectResult.fpSets.negative, 2, 2, False)
        sources = SourceCatalog(table)
        fpSet.makeSources(sources)

        measureTask.run(sources, self.dipole, self.posImage, self.negImage)
        return sources

    def testDipoleFitter(self):
        pass

    def testDipoleTask(self):
        """
        Test the dipole fitting singleFramePlugin. Test that the resulting fluxes/centroids
        are entered into the correct slots of the catalog, and have values that are
        very close to the input values for both dipoles in the image.
        """
        sources = self.runDetection()

        offsets = self.params.offsets
        for i, r1 in enumerate(sources):
            result = r1.extract("ip_diffim_DipoleFit*")
            self.assertClose((result['ip_diffim_DipoleFit_pos_flux'] +
                              abs(result['ip_diffim_DipoleFit_neg_flux']))/2.,
                             self.params.flux[i], rtol=0.02)
            self.assertClose(result['ip_diffim_DipoleFit_pos_centroid_x'],
                             self.params.xc[i] + offsets[i], rtol=0.01)
            self.assertClose(result['ip_diffim_DipoleFit_pos_centroid_y'],
                             self.params.yc[i] + offsets[i], rtol=0.01)
            self.assertClose(result['ip_diffim_DipoleFit_neg_centroid_x'],
                             self.params.xc[i] - offsets[i], rtol=0.01)
            self.assertClose(result['ip_diffim_DipoleFit_neg_centroid_y'],
                             self.params.yc[i] - offsets[i], rtol=0.01)
            # Note this is dependent on the noise (variance) being realistic in the image.
            # otherwise it throws off the chi2 estimate, which is used for classification:
            self.assertTrue(result['ip_diffim_DipoleFit_flag_classification'])

            # compare to the original ip_diffim_PsfDipoleFlux measurements
            result2 = r1.extract("ip_diffim_PsfDipoleFlux*")
            self.assertClose((result['ip_diffim_DipoleFit_pos_flux'] +
                              abs(result['ip_diffim_DipoleFit_neg_flux']))/2.,
                             (result2['ip_diffim_PsfDipoleFlux_pos_flux'] +
                              abs(result2['ip_diffim_PsfDipoleFlux_neg_flux']))/2.,
                             rtol=0.02)
            self.assertClose(result['ip_diffim_DipoleFit_pos_centroid_x'],
                             result2['ip_diffim_PsfDipoleFlux_pos_centroid_x'],
                             rtol=0.01)
            self.assertClose(result['ip_diffim_DipoleFit_pos_centroid_y'],
                             result2['ip_diffim_PsfDipoleFlux_pos_centroid_y'],
                             rtol=0.01)
            self.assertClose(result['ip_diffim_DipoleFit_neg_centroid_x'],
                             result2['ip_diffim_PsfDipoleFlux_neg_centroid_x'],
                             rtol=0.01)
            self.assertClose(result['ip_diffim_DipoleFit_neg_centroid_y'],
                             result2['ip_diffim_PsfDipoleFlux_neg_centroid_y'],
                             rtol=0.01)

            if self.params.display:
                DipolePlotUtils.displayCutouts(r1, self.dipole, self.posImage, self.negImage)
        if self.params.display:
            DipolePlotUtils.plt.show()

        return result


# Test the task in the same way as the algorithm:
# Here test that dipoles too close to the edge are raised correctly
class DipoleFitTaskEdgeTest(DipoleFitTaskTest):
    """ A test case for dipole fit task"""
    def setUp(self):
        # Ensure that both dipoles will fail (too close to edge)
        self.params = DipoleFitTestGlobalParams()
        self.params.xc = np_array([5.3, 2.2])  # xcenters of two dipoles in image
        self.params.yc = np_array([2.6, 98.5])  # ycenters of two dipoles

        offsets = self.params.offsets
        self.dipole, (self.posImage, self.posCatalog), (self.negImage, self.negCatalog) = \
            DipoleTestUtils.makeDipoleImage(
                xcenPos=self.params.xc + offsets,
                ycenPos=self.params.yc + offsets,
                xcenNeg=self.params.xc - offsets,
                ycenNeg=self.params.yc - offsets,
                flux=self.params.flux, fluxNeg=self.params.flux,
                gradientParams=self.params.gradientParams)

        self.catalog = DipoleTestUtils.detectDipoleSources(self.dipole)

    def tearDown(self):
        DipoleFitTaskTest.tearDown(self)

    def testDipoleTask(self):
        pass

    def testDipoleEdge(self):
        """
        Test the dipole fitting singleFramePlugin. Test that the dipoles which are too
        close to the edge raise the correct exception.
        """

        sources = DipoleFitTaskTest.runDetection(self)

        for i, r1 in enumerate(sources):
            result = r1.extract("ip_diffim_DipoleFit*")
            self.assertTrue(result.get("ip_diffim_DipoleFit_flag"))


# UTILITY CLASS WITH STATIC METHODS FOR DIPOLE TESTING ###
class DipoleTestUtils(object):

    @staticmethod
    def makeStarImage(w=101, h=101, xc=[15.3], yc=[18.6], flux=[2500], psfSigma=2., noise=10.0,
                      gradientParams=None, schema=None):

        from lsst.meas.base.tests import TestDataset
        bbox = Box2I(Point2I(0, 0), Point2I(w-1, h-1))
        dataset = TestDataset(bbox, psfSigma=psfSigma, threshold=1.)

        for i in xrange(len(xc)):
            dataset.addSource(flux=flux[i], centroid=Point2D(xc[i], yc[i]))

        if schema is None:
            schema = TestDataset.makeMinimalSchema()
        exposure, catalog = dataset.realize(noise=noise, schema=schema)

        if gradientParams is not None:
            y, x = np_mgrid[:w, :h]
            gp = gradientParams
            gradient = gp[0] + gp[1] * x + gp[2] * y
            if len(gradientParams) > 3:  # it includes a set of 2nd-order polynomial params
                gradient += gp[3] * x*y + gp[4] * x*x + gp[5] * y*y
            imgArr = exposure.getMaskedImage().getArrays()[0]
            imgArr += gradient

        return exposure, catalog

    @staticmethod
    def makeDipoleImage(w=101, h=101, xcenPos=[27.], ycenPos=[25.], xcenNeg=[23.], ycenNeg=[25.],
                        psfSigma=2., flux=[30000.], fluxNeg=None, noise=10., gradientParams=None):

        posImage, posCatalog = DipoleTestUtils.makeStarImage(
            w, h, xcenPos, ycenPos, flux=flux, psfSigma=psfSigma,
            gradientParams=gradientParams, noise=noise)

        if fluxNeg is None:
            fluxNeg = flux
        negImage, negCatalog = DipoleTestUtils.makeStarImage(
            w, h, xcenNeg, ycenNeg, flux=fluxNeg, psfSigma=psfSigma,
            gradientParams=gradientParams, noise=noise)

        dipole = posImage.clone()
        di = dipole.getMaskedImage()
        di -= negImage.getMaskedImage()

        # Carry through pos/neg detection masks to new planes in diffim image
        dm = di.getMask()
        posDetectedBits = posImage.getMaskedImage().getMask().getArray() == dm.getPlaneBitMask("DETECTED")
        negDetectedBits = negImage.getMaskedImage().getMask().getArray() == dm.getPlaneBitMask("DETECTED")
        pos_det = dm.addMaskPlane("DETECTED_POS")  # new mask plane -- different from "DETECTED"
        neg_det = dm.addMaskPlane("DETECTED_NEG")  # new mask plane -- different from "DETECTED_NEGATIVE"
        dma = dm.getArray()
        # set the two custom mask planes to these new masks
        dma[:, :] = posDetectedBits * pos_det + negDetectedBits * neg_det
        return dipole, (posImage, posCatalog), (negImage, negCatalog)

    @staticmethod
    def detectDipoleSources(diffim, doMerge=True, detectSigma=5.5, grow=3):
        """
        Utility function for detecting dipoles. Detects pos/neg sources in the diffim,
        then merges them. A bigger "grow" parameter leads to a larger footprint which
        helps with dipole measurement for faint dipoles.
        """

        # Start with a minimal schema - only the fields all SourceCatalogs need
        schema = SourceTable.makeMinimalSchema()

        # Customize the detection task a bit (optional)
        detectConfig = SourceDetectionConfig()
        detectConfig.returnOriginalFootprints = False  # should be the default

        psfSigma = diffim.getPsf().computeShape().getDeterminantRadius()

        # code from imageDifference.py:
        detectConfig.thresholdPolarity = "both"
        detectConfig.thresholdValue = detectSigma
        # detectConfig.nSigmaToGrow = psfSigma
        detectConfig.reEstimateBackground = True  # if False, will fail often for faint sources on gradients?
        detectConfig.thresholdType = "pixel_stdev"

        # Create the detection task. We pass the schema so the task can declare a few flag fields
        detectTask = SourceDetectionTask(schema, config=detectConfig)

        table = SourceTable.make(schema)
        catalog = detectTask.makeSourceCatalog(table, diffim, sigma=psfSigma)

        # Now do the merge.
        if doMerge:
            fpSet = catalog.fpSets.positive
            fpSet.merge(catalog.fpSets.negative, grow, grow, False)
            sources = SourceCatalog(table)
            fpSet.makeSources(sources)

            return sources

        else:
            return detectTask, schema


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst_tests.init()

    suites = []
    suites += unittest.makeSuite(DipoleFitAlgorithmTest)
    suites += unittest.makeSuite(DipoleFitTaskTest)
    suites += unittest.makeSuite(DipoleFitTaskEdgeTest)
    suites += unittest.makeSuite(lsst_tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(shouldExit = False):
    """Run the tests"""
    lsst_tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)

