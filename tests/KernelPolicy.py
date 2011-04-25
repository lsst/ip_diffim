#!/usr/bin/env python
import unittest
import eups
import os
import lsst.utils.tests as tests
import lsst.ip.diffim as ipDiffim
import lsst.pex.policy as pexPolicy

class DiffimTestCases(unittest.TestCase):
    def setUp(self):
        diffimDir         = eups.productDir('ip_diffim')
        self.policyPath   = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')
        self.p0           = pexPolicy.Policy.createPolicy(self.policyPath)

    def testNoModify(self):
        p1 = ipDiffim.makeDefaultPolicy(self.policyPath, modify=False)
        self.assertEqual(self.p0.get("kernelSize"), p1.get("kernelSize"))
        self.assertEqual(self.p0.get("fpGrowPix"), p1.get("fpGrowPix"))
        for i in range(self.p0.get("alardNGauss")):
            self.assertEqual(self.p0.getDoubleArray("alardSigGauss")[i],
                             p1.getDoubleArray("alardSigGauss")[i])

    def testModifyGreater(self, fwhm=10.):
        p1 = ipDiffim.makeDefaultPolicy(self.policyPath, fwhm=fwhm)
        self.assertTrue(self.p0.get("kernelSize") < p1.get("kernelSize"))
        self.assertTrue(self.p0.get("fpGrowPix") < p1.get("fpGrowPix"))
        for i in range(self.p0.get("alardNGauss")):
            self.assertTrue(self.p0.getDoubleArray("alardSigGauss")[i] <
                            p1.getDoubleArray("alardSigGauss")[i])

        # maxed out the sizes
        self.assertTrue(p1.get("fpGrowPix")  == self.p0.get("fpGrowMax"))
        self.assertTrue(p1.get("kernelSize") == self.p0.get("kernelSizeMax"))

    def testModifyLesser(self, fwhm=1.):
        p1 = ipDiffim.makeDefaultPolicy(self.policyPath, fwhm=fwhm)
        self.assertTrue(self.p0.get("kernelSize") > p1.get("kernelSize"))

        self.assertTrue(self.p0.get("fpGrowPix") > p1.get("fpGrowPix"))
        for i in range(self.p0.get("alardNGauss")):
            self.assertTrue(self.p0.getDoubleArray("alardSigGauss")[i] >
                            p1.getDoubleArray("alardSigGauss")[i])

        # minned out the sizes
        self.assertTrue(p1.get("fpGrowPix") == self.p0.get("fpGrowMin"))
        self.assertTrue(p1.get("kernelSize") == self.p0.get("kernelSizeMin"))

    def tearDown(self):
        del self.p0
        
#####
        
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
    
