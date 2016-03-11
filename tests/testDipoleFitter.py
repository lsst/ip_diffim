#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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
                   array as np_array)

## LSST imports:
import lsst.utils.tests as lsst_tests
from lsst.afw.table import (SourceTable, SourceCatalog)
from lsst.ip.diffim import dipoleFitTask as dft
from lsst.meas.base import SingleFrameMeasurementConfig

class DipoleFitTestGlobalParams():
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
        self.display = True
        self.verbose = False
        self.w, self.h = 100, 100 # size of image

        self.xc = np_array([65.3, 24.2]) ## xcenters of two dipoles in image
        self.yc = np_array([38.6, 78.5]) ## ycenters of two dipoles
        self.flux = np_array([2500., 2345.])  ## fluxes of pos/neg lobes
        self.gradientParams = np_array([10., 3., 5.])

        self.offsets = np_array([-2., 2.]) ## pixel coord offsets between lobes of dipoles

## First, test the algorithm itself (fitDipole_new()):
## Create a simulated diffim (with dipoles) and a linear background gradient in the pre-sub images
##   then compare the input fluxes/centroids with the fitted results.
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
            dft.DipoleUtils.makeDipoleImage(
                xcenPos=self.params.xc + offsets,
                ycenPos=self.params.yc + offsets,
                xcenNeg=self.params.xc - offsets,
                ycenNeg=self.params.yc - offsets,
                flux=self.params.flux, fluxNeg=self.params.flux,
                gradientParams=self.params.gradientParams)

        self.catalog = dft.DipoleUtils.detectDipoleSources(
            self.dipole, self.posImage, self.posCatalog,
            self.negImage, self.negCatalog)

    def tearDown(self):
        del self.dipole, self.posImage, self.negImage
        del self.catalog, self.posCatalog, self.negCatalog
        del self.params

    def testDipoleFitter(self):
        """
        Test the dipole fitting algorithm. Test that the resulting fluxes/centroids
        are very close to the input values for both dipoles in the image.
        """
        for s in self.catalog:
            fp = s.getFootprint()
            print fp.getBBox(), fp.getNpix()
            for pk in fp.getPeaks():
                print 'FOOTPRINT CENTER:', pk.getIy(), pk.getIx(), pk.getPeakValue()

        offsets = self.params.offsets
        for i,s in enumerate(self.catalog):
            result = dft.DipoleFitAlgorithm.fitDipole_new(
                self.dipole, s, self.posImage, self.negImage,
                rel_weight=1., separateNegParams=False,
                verbose=self.params.verbose, display=self.params.display)
            print result
            self.assertClose((result.psfFitPosFlux + abs(result.psfFitNegFlux))/2.,
                             self.params.flux[i], rtol=0.02)
            self.assertClose(result.psfFitPosCentroidX, self.params.xc[i] + offsets[i], rtol=0.01)
            self.assertClose(result.psfFitPosCentroidY, self.params.yc[i] + offsets[i], rtol=0.01)
            self.assertClose(result.psfFitNegCentroidX, self.params.xc[i] - offsets[i], rtol=0.01)
            self.assertClose(result.psfFitNegCentroidY, self.params.yc[i] - offsets[i], rtol=0.01)

## Second, test the task in the same way as the algorithm:
## Also test that it correctly classifies the dipoles.
class DipoleFitTaskTest(DipoleFitAlgorithmTest):
    """ A test case for dipole fit task"""
    def setUp(self):
        DipoleFitAlgorithmTest.setUp(self)

    def tearDown(self):
        DipoleFitAlgorithmTest.tearDown(self)

    def testDipoleTask(self):
        """
        Test the dipole fitting singleFramePlugin. Test that the resulting fluxes/centroids
        are entered into the correct slots of the catalog, and have values that are
        very close to the input values for both dipoles in the image.
        """

        ## Create the various tasks and schema -- avoid code reuse.
        detectTask, deblendTask, schema = dft.DipoleUtils.detectDipoleSources(
            self.dipole, self.posImage,
            self.posCatalog, self.negImage,
            self.negCatalog, doMerge=False)

        #measureConfig = dft.DipoleFitConfig()
        measureConfig = SingleFrameMeasurementConfig()

        # Modify the set of active plugins ('.names' behaves like a Python set)
        measureConfig.plugins.names.remove("base_PsfFlux")
        measureConfig.plugins.names.remove("base_GaussianCentroid")
        # If I remove any of the following plugins, I get an 'base_GaussianFlux_flux' not found error
        #measureConfig.plugins.names.remove("base_SdssCentroid")
        #measureConfig.plugins.names.remove("base_GaussianFlux")
        #measureConfig.plugins.names.remove("base_SdssShape")

        measureConfig.slots.modelFlux = "ip_diffim_DipoleFit"

        measureConfig.plugins.names |= ["base_CircularApertureFlux",
                                        "base_PixelFlags",
                                        "base_SkyCoord",
                                        "base_PsfFlux",
                                        "ip_diffim_NaiveDipoleCentroid",
                                        "ip_diffim_NaiveDipoleFlux",
                                        "ip_diffim_PsfDipoleFlux"]
                                        #"ip_diffim_DipoleFit"]

        # Modify the internal configuration of one of the plugins
        #measureConfig.plugins["base_ClassificationExtendedness"].fluxRatio = 0.985

        # Disable aperture correction, which requires having an ApCorrMap attached to
        # the Exposure (it'll warn if it's not present and we don't explicitly disable it).
        measureConfig.doApplyApCorr = "no"

        measureTask = dft.DipoleFitTask(config=measureConfig, schema=schema)

        table = SourceTable.make(schema)
        detectResult = detectTask.run(table, self.dipole)
        catalog = detectResult.sources
        deblendTask.run(self.dipole, catalog, psf=self.dipole.getPsf())

        fpSet = detectResult.fpSets.positive
        fpSet.merge(detectResult.fpSets.negative, 2, 2, False)
        sources = SourceCatalog(table)
        fpSet.makeSources(sources)

        measureTask.run(sources, self.dipole, self.posImage, self.negImage)

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
            self.assertTrue(result['ip_diffim_DipoleFit_flag_classification'])

            ## compare to the original ip_diffim_PsfDipoleFlux measurements
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


        #for r1 in sources:
        #    print r1.extract("ip_diffim_PsfDipoleFlux*")
        #    print r1.extract("base_CircularApertureFlux_25_*")

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst_tests.init()

    suites = []
    suites += unittest.makeSuite(DipoleFitAlgorithmTest)
    suites += unittest.makeSuite(DipoleFitTaskTest)
    suites += unittest.makeSuite(lsst_tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    """Run the tests"""
    lsst_tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)

