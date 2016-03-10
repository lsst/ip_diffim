import unittest

import numpy as np

## LSST imports:
import lsst.utils.tests as tests
#import lsst.afw.image as afwImage
#import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
#import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
#import lsst.ip.diffim as ipDiffim
from lsst.ip.diffim import dipoleFitTask as dft
## need to separately set up meas_modelfit:
import lsst.meas.modelfit as modelfit
import lsst.meas.base as measBase

# import lsst.meas.base
# import lsst.pex.config
#import lsst.meas.deblender
from lsst.meas.deblender import SourceDeblendTask
import lsst.utils.tests

## First, test the algorithm itself (fitDipole_new()):
## Create a simulated diffim (with dipoles) and a linear background gradient in the pre-sub images
##   then compare the input fluxes/centroids with the fitted results.
class DipoleFitAlgorithmTest(lsst.utils.tests.TestCase):
    """ A test case for dipole fit algorithm"""
    def setUp(self):
        np.random.seed(666)
        self.w, self.h = 100, 100 # size of image
        self.xc, self.yc = 50, 50 # location of center of dipole

        self.xc = np.array([65.3, 24.2]) ## xcenters of two dipoles in image
        self.yc = np.array([38.6, 78.5]) ## ycenters of two dipoles
        self.flux = np.array([2500., 2345.])  ## fluxes of pos/neg lobes
        self.gradientParams = (10., 3., 5.)

        offsets = np.array([-2., 2.])
        self.dipole, (self.posImage, self.posCatalog), (self.negImage, self.negCatalog) = \
            dft.makeDipoleImage_lsst(xcenPos=self.xc + offsets,
                                     ycenPos=self.yc + offsets,
                                     xcenNeg=self.xc - offsets,
                                     ycenNeg=self.yc - offsets,
                                     flux=self.flux, fluxNeg=self.flux,
                                     gradientParams=self.gradientParams)

        self.catalog = dft.detectDipoleSources(self.dipole, self.posImage, self.posCatalog,
                                               self.negImage, self.negCatalog)

    def tearDown(self):
        del self.dipole, self.posImage, self.negImage, self.catalog, self.posCatalog, self.negCatalog

    def testDipoleFitter(self):
        for s in self.catalog:
            fp = s.getFootprint()
            print fp.getBBox(), fp.getNpix()
            for pk in fp.getPeaks():
                print 'FOOTPRINT CENTER:', pk.getIy(), pk.getIx(), pk.getPeakValue()

        offsets = np.array([-2., 2.])
        for i,s in enumerate(self.catalog):
            result = dft.fitDipole_new(self.dipole, s, self.posImage, self.negImage,
                                   rel_weight=1., separateNegParams=False,
                                   verbose=False, display=False)
            self.assertClose((result.psfFitPosFlux + abs(result.psfFitNegFlux))/2., self.flux[i], rtol=0.02)
            self.assertClose(result.psfFitPosCentroidX, self.xc[i] + offsets[i], rtol=0.01)
            self.assertClose(result.psfFitPosCentroidY, self.yc[i] + offsets[i], rtol=0.01)
            self.assertClose(result.psfFitNegCentroidX, self.xc[i] - offsets[i], rtol=0.01)
            self.assertClose(result.psfFitNegCentroidY, self.yc[i] - offsets[i], rtol=0.01)

## Second, test the task in the same way as the algorithm:
## Also test that it correctly classifies the dipoles.
class DipoleFitTaskTest(lsst.utils.tests.TestCase):
    """ A test case for dipole fit task"""
    def setUp(self):
        np.random.seed(666)
        self.w, self.h = 100, 100 # size of image
        self.xc, self.yc = 50, 50 # location of center of dipole

        self.xc = np.array([65.3, 24.2]) ## xcenters of two dipoles in image
        self.yc = np.array([38.6, 78.5]) ## ycenters of two dipoles
        self.flux = np.array([2500., 2345.])  ## fluxes of pos/neg lobes
        self.gradientParams = (10., 3., 5.)

        self.dipole, (self.posImage, self.posCatalog), (self.negImage, self.negCatalog) = \
            dft.makeDipoleImage_lsst(xcenPos=self.xc + np.array([-2., +2.]),
                                     ycenPos=self.yc + np.array([-2., +2.]),
                                     xcenNeg=self.xc - np.array([-2., +2.]),
                                     ycenNeg=self.yc - np.array([-2., +2.]),
                                     flux=self.flux, fluxNeg=self.flux,
                                     gradientParams=self.gradientParams)

        self.catalog = dft.detectDipoleSources(self.dipole, self.posImage, self.posCatalog,
                                               self.negImage, self.negCatalog)

    def tearDown(self):
        del self.dipole, self.posImage, self.negImage, self.catalog, self.posCatalog, self.negCatalog

    def testDipoleTask(self):
        ## Create the various tasks and schema -- avoid code reuse.
        detectTask, deblendTask, schema = dft.detectDipoleSources(self.dipole, self.posImage,
                                                                   self.posCatalog, self.negImage,
                                                                   self.negCatalog, doMerge=False)

        #measureConfig = dft.DipoleFitConfig()
        measureConfig = measBase.SingleFrameMeasurementConfig()

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

        table = afwTable.SourceTable.make(schema)
        detectResult = detectTask.run(table, self.dipole)
        catalog = detectResult.sources
        deblendTask.run(self.dipole, catalog, psf=self.dipole.getPsf())

        fpSet = detectResult.fpSets.positive
        fpSet.merge(detectResult.fpSets.negative, 2, 2, False)
        sources = afwTable.SourceCatalog(table)
        fpSet.makeSources(sources)

        measureTask.run(sources, self.dipole, self.posImage, self.negImage)

        offsets = np.array([-2., 2.])
        for i, r1 in enumerate(sources):
            result = r1.extract("ip_diffim_DipoleFit*")
            self.assertClose((result['ip_diffim_DipoleFit_pos_flux'] +
                              abs(result['ip_diffim_DipoleFit_neg_flux']))/2., self.flux[i], rtol=0.02)
            self.assertClose(result['ip_diffim_DipoleFit_pos_centroid_x'], self.xc[i] + offsets[i], rtol=0.01)
            self.assertClose(result['ip_diffim_DipoleFit_pos_centroid_y'], self.yc[i] + offsets[i], rtol=0.01)
            self.assertClose(result['ip_diffim_DipoleFit_neg_centroid_x'], self.xc[i] - offsets[i], rtol=0.01)
            self.assertClose(result['ip_diffim_DipoleFit_neg_centroid_y'], self.yc[i] - offsets[i], rtol=0.01)
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

    tests.init()

    suites = []
    suites += unittest.makeSuite(DipoleFitAlgorithmTest)
    suites += unittest.makeSuite(DipoleFitTaskTest)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    """Run the tests"""
    tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)

