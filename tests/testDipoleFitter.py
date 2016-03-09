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

## First, test the algorithm itself (fitDipole_new()):
class DipoleFitAlgorithmTest(unittest.TestCase):
    """ A test case for dipole fit algorithm"""
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
        pass

    def testDipoleFitter(self):
        for s in self.catalog:
            fp = s.getFootprint()
            print fp.getBBox(), fp.getNpix()
            for pk in fp.getPeaks():
                print 'FOOTPRINT CENTER:', pk.getIy(), pk.getIx(), pk.getPeakValue()

        for s in self.catalog:
            result = dft.fitDipole_new(self.dipole, s, self.posImage, self.negImage,
                                   rel_weight=1., separateNegParams=False,
                                   verbose=False, display=False)
            print result

## Second, test the task:
class DipoleFitTaskTest(unittest.TestCase):
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
        pass

    def testDipoleTask(self):
        # Start with a minimal schema - only the fields all SourceCatalogs need
        schema = afwTable.SourceTable.makeMinimalSchema()

        # Create and run a task for deblending (optional, but almost always a good idea). 
        # Again, the task defines a few flag fields it will later fill. 
        deblendTask = SourceDeblendTask(schema=schema)

        # Customize the detection task a bit (optional)
        detectConfig = measAlg.SourceDetectionConfig() 
        detectConfig.returnOriginalFootprints = False # should be the default 
        detectConfig.thresholdValue = 10 # only 10-sigma detections

        ## code from imageDifference.py:
        detectConfig.thresholdPolarity = "both"
        detectConfig.thresholdValue = 5.5
        detectConfig.reEstimateBackground = False
        detectConfig.thresholdType = "pixel_stdev"

        # Create the detection task. We pass the schema so the task can declare a few flag fields
        detectTask = measAlg.SourceDetectionTask(config=detectConfig, schema=schema)

        table = afwTable.SourceTable.make(schema)
        detectResult = detectTask.run(table, self.dipole)
        catalog = detectResult.sources
        #results = detectTask.makeSourceCatalog(table, exposure, sigma=psfSigma)

        deblendTask.run(self.dipole, catalog, psf=self.dipole.getPsf())

        measureConfig = measBase.SingleFrameMeasurementConfig()

        # Modify the set of active plugins ('.names' behaves like a Python set)
        #measureConfig.plugins.names.remove("base_GaussianCentroid")

        # Enable some plugins - import the Python module first to make them available
        #measureConfig.plugins.names |= ["modelfit_ShapeletPsfApprox", "modelfit_CModel"]

        # Change which plugin's output we "bless" as the "Model Flux"
        #measureConfig.slots.modelFlux = "modelfit_CModel"

        measureConfig.slots.modelFlux = "ip_diffim_DipoleFit"

        #["base_CircularApertureFlux",
    #                     "base_PixelFlags",
    #                     "base_SkyCoord",
    #                     "base_PsfFlux",
    #                     #"ip_diffim_NaiveDipoleCentroid",
    #                     #"ip_diffim_NaiveDipoleFlux",
    #                     #"ip_diffim_PsfDipoleFlux"]
    #                     DipoleFitTask._DefaultName]        

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

        for r1 in sources: #, r2 in zip(catalog, transformedCatalog):
            print r1.extract("ip_diffim_DipoleFit*")

        # Compare to input:
        # ```
        # # xc = npa([65.3, 24.2])
        # # yc = npa([38.6, 78.5])
        # # flux = npa([2500., 2345.])
        # ```
        # and output directly from `fitDipole_new()`:
        # ```
        # resultsOutput(psfFitPosCentroidX=26.132875038347073, psfFitPosCentroidY=80.489634132198631, psfFitNegCentroidX=22.080022018131551, psfFitNegCentroidY=76.444291337973098, psfFitPosFlux=2366.4792132019111, psfFitNegFlux=-2434.333049782319, psfFitCentroidX=24.106448528239312, psfFitCentroidY=78.466962735085872, psfFitOrientation=44.946864198037069, psfFitSignaltoNoise=57.14172173097851)
        # ```

def suite():
    """Returns a suite containing all the test cases in this module."""

    tests.init()

    suites = []
    suites += unittest.makeSuite(DipoleFitAlgorithmTest)
    suites += unittest.makeSuite(DipoleFitTaskTest)
    #suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    """Run the tests"""
    tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)

