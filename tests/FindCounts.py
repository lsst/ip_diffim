#!/usr/bin/env python
import os

import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.pex.policy as pexPolicy
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as logging

import lsst.afw.display.ds9 as ds9

Verbosity = 4
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = pexPolicy.Policy.createPolicy(diffimPolicy)
        
    def tearDown(self):
        del self.policy

    def testCountsF(self):
        mi = afwImage.MaskedImageF(20,20)
        mi.set(100, 0x0, 1)

        fc = ipDiffim.FindCountsF()
        fc.apply(mi)

        self.assertEqual(fc.getCounts(), 20*20*100)

    def testCountsMask(self):
        mi = afwImage.MaskedImageF(20,20)
        mi.set(100, 0x0, 1)
        mi.set(10, 10, (10000, 0x1, 1))

        fc = ipDiffim.FindCountsF()
        fc.apply(mi)

        self.assertEqual(fc.getCounts(), (20*20 - 1)*100)

    def testCountsI(self):
        mi = afwImage.MaskedImageI(20,20)
        mi.set(100, 0x0, 1)

        fc = ipDiffim.FindCountsI()
        fc.apply(mi)

        self.assertEqual(fc.getCounts(), 20*20*100)

#####
        
def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(DiffimTestCases)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
