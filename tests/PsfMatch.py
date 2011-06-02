#!/usr/bin/env python
import unittest
import lsst.utils.tests as tests
import lsst.afw.image as afwImage
import lsst.afw.image.utils as imageUtils
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.daf.base as dafBase
import lsst.pex.policy as pexPolicy

import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)

class DiffimTestCases(unittest.TestCase):

    def setUp(self):
        self.policy = ipDiffim.makeDefaultPolicy()

        self.policy.set("kernelBasisSet", "alard-lupton")
        self.policy.set("alardNGauss", 1)
        self.policy.set("alardSigGauss", 2.5)
        self.policy.set("alardDegGauss", 2)
        self.basisList = ipDiffim.makeKernelBasisList(self.policy)


    def makeWcs(self, offset = 0):
        # taken from $AFW_DIR/tests/testMakeWcs.py
        metadata = dafBase.PropertySet()
        metadata.set("SIMPLE",                    "T") 
        metadata.set("BITPIX",                  -32) 
        metadata.set("NAXIS",                    2) 
        metadata.set("NAXIS1",                 1024) 
        metadata.set("NAXIS2",                 1153) 
        metadata.set("RADECSYS", 'FK5')
        metadata.set("EQUINOX",                2000.)
        metadata.setDouble("CRVAL1",     215.604025685476)
        metadata.setDouble("CRVAL2",     53.1595451514076)
        metadata.setDouble("CRPIX1",     1109.99981456774 + offset)
        metadata.setDouble("CRPIX2",     560.018167811613 + offset)
        metadata.set("CTYPE1", 'RA---SIN')
        metadata.set("CTYPE2", 'DEC--SIN')
        metadata.setDouble("CD1_1", 5.10808596133527E-05)
        metadata.setDouble("CD1_2", 1.85579539217196E-07)
        metadata.setDouble("CD2_2", -5.10281493481982E-05)
        metadata.setDouble("CD2_1", -8.27440751733828E-07)
        return afwImage.makeWcs(metadata)
        
    def testSubtractExposures(self, background = 100.):
        tMi, sMi, sK, kcs = diffimTools.makeFakeKernelSet(self.policy, self.basisList,
                                                          bgValue = background)
        tWcs = self.makeWcs(offset = 0)
        sWcs = self.makeWcs(offset = 1)
        tExp = afwImage.ExposureF(tMi, tWcs)
        sExp = afwImage.ExposureF(sMi, sWcs)

        self.policy.set("fitForBackground", True)
        self.policy.set("spatialKernelOrder", 1)
        self.policy.set("spatialBgOrder", 0)
        self.policy.set("spatialKernelType", "polynomial") # since that is the known function
        psfmatch = ipDiffim.ImagePsfMatch(self.policy)

        # Should fail due to registration problem
        try:
            results = psfmatch.subtractExposures(tExp, sExp, doWarping = False)
        except:
            pass
        else:
            self.fail()

        # Should work
        results = psfmatch.subtractExposures(tExp, sExp, doWarping = True)

        self.assertEqual(len(results), 4)
        self.assertEqual(type(results[0]), afwImage.ExposureF)
        self.assertEqual(type(results[1]), afwMath.LinearCombinationKernel)
        self.assertEqual(type(results[2]), afwMath.Function2D)
        self.assertEqual(type(results[3]), afwMath.SpatialCellSet)

    def testMatchExposures(self, background = 100.):
        filterPolicyFile = pexPolicy.DefaultPolicyFile("afw", "SdssFilters.paf", "tests")
        filterPolicy = pexPolicy.Policy.createPolicy(filterPolicyFile, filterPolicyFile.getRepositoryPath(), True)
        imageUtils.defineFiltersFromPolicy(filterPolicy, reset=True)

        tMi, sMi, sK, kcs = diffimTools.makeFakeKernelSet(self.policy, self.basisList,
                                                          bgValue = background)
        tWcs = self.makeWcs(offset = 0)
        sWcs = self.makeWcs(offset = 0)
        tExp = afwImage.ExposureF(tMi, tWcs)
        sExp = afwImage.ExposureF(sMi, sWcs)

        commonFilter = afwImage.Filter("r")
        tCalib = afwImage.Calib()
        tCalib.setFluxMag0(1.0e5, 1.0e3)
        tExp.setFilter(commonFilter)
        sExp.setFilter(commonFilter)
        tExp.setCalib(tCalib)

        self.policy.set("fitForBackground", True)
        self.policy.set("spatialKernelOrder", 1)
        self.policy.set("spatialBgOrder", 0)
        self.policy.set("spatialKernelType", "polynomial") # since that is the known function
        psfmatch = ipDiffim.ImagePsfMatch(self.policy)

        # Should work since already registered
        results = psfmatch.matchExposures(tExp, sExp, doWarping = False)

        # Should also work
        results = psfmatch.matchExposures(tExp, sExp, doWarping = True)
        psfMatchedExp = results[0]

        self.assertEqual(len(results), 4)
        self.assertEqual(type(results[0]), afwImage.ExposureF)
        self.assertEqual(type(results[1]), afwMath.LinearCombinationKernel)
        self.assertEqual(type(results[2]), afwMath.Function2D)
        self.assertEqual(type(results[3]), afwMath.SpatialCellSet)
        self.assertEqual(psfMatchedExp.getFilter().getName(), commonFilter.getName())
        print "Warning: testMatchExposures should test Calib object"

    def testSubtractMaskedImages(self, background = 100.):
        tMi, sMi, sK, kcs = diffimTools.makeFakeKernelSet(self.policy, self.basisList,
                                                          bgValue = background,
                                                          nCell = 9)
        self.policy.set("fitForBackground", True)
        self.policy.set("spatialKernelOrder", 1)
        self.policy.set("spatialBgOrder", 0)
        self.policy.set("spatialKernelType", "polynomial") # since that is the known function
        psfmatch = ipDiffim.ImagePsfMatch(self.policy)
        results = psfmatch.subtractMaskedImages(tMi, sMi)

        self.assertEqual(len(results), 4)
        self.assertEqual(type(tMi), type(results[0]))
        imstat = ipDiffim.ImageStatisticsF()
        imstat.apply(results[0])
        self.assertAlmostEqual(imstat.getRms(), 1.0, 1)
        self.assertEqual(type(results[1]), afwMath.LinearCombinationKernel)
        self.assertEqual(type(results[2]), afwMath.Function2D)
        self.assertEqual(type(results[3]), afwMath.SpatialCellSet)

        # Test the values of the known spatial coefficients.  Its a
        # bit tricky since noise is added and skews the results
        knownCoeffs = diffimTools.fakeCoeffs()
        fitCoeffs   = results[1].getSpatialParameters()

        for b in range(len(knownCoeffs)):
            for s in range(len(knownCoeffs[b])):
                self.assertAlmostEqual(knownCoeffs[b][s], fitCoeffs[b][s], 1)
        
    def testMatchMaskedImages(self, background = 100.):
        tMi, sMi, sK, kcs = diffimTools.makeFakeKernelSet(self.policy, self.basisList,
                                                          bgValue = background)
        self.policy.set("fitForBackground", True)
        self.policy.set("spatialKernelOrder", 1)
        self.policy.set("spatialBgOrder", 0)
        self.policy.set("spatialKernelType", "polynomial") # since that is the known function
        psfmatch = ipDiffim.ImagePsfMatch(self.policy)
        results = psfmatch.matchMaskedImages(tMi, sMi)

        self.assertEqual(len(results), 4)
        self.assertEqual(type(tMi), type(results[0]))
        self.assertAlmostEqual(afwMath.makeStatistics(results[0].getImage(), afwMath.MEDIAN).getValue(),
                               afwMath.makeStatistics(tMi.getImage(), afwMath.MEDIAN).getValue(), 1)
        self.assertEqual(type(results[1]), afwMath.LinearCombinationKernel)
        self.assertEqual(type(results[2]), afwMath.Function2D)
        self.assertEqual(type(results[3]), afwMath.SpatialCellSet)


    def testModelMatch(self):
        # Remains to be written
        psfmatch = ipDiffim.ModelPsfMatch(self.policy)
        print "Warning: ModelPsfMatch test not written"
    
    def tearDown(self):
        del self.policy
        del self.basisList

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(DiffimTestCases)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(doExit=False):
    """Run the tests"""
    tests.run(suite(), doExit)

if __name__ == "__main__":
    run(True)
