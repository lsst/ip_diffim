#!/usr/bin/env python
import unittest
import lsst.utils.tests as tests
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools

import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)

class DiffimTestCases(unittest.TestCase):
    # Some remain to be written
    def setUp(self):
        self.policy = ipDiffim.makeDefaultPolicy()

        self.policy.set("kernelBasisSet", "alard-lupton")
        self.policy.set("alardNGauss", 1)
        self.policy.set("alardSigGauss", 2.5)
        self.policy.set("alardDegGauss", 2)
        self.basisList = ipDiffim.makeKernelBasisList(self.policy)

    def testSubtractExposures(self):
        psfmatch = ipDiffim.ImagePsfMatch(self.policy)
        
    def testMatchExposures(self):
        psfmatch = ipDiffim.ImagePsfMatch(self.policy)

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
        psfmatch = ipDiffim.ModelPsfMatch(self.policy)
    
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
