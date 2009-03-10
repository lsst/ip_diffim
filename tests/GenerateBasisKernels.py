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

    def doit(self, width=10, height=10):
        ks1 = ipDiffim.generateDeltaFunctionKernelSet(width, height)
        nk  = 0
        for rowi in range(height):
            for colj in range(width):
                kernel = ks1[nk]
                kimage = afwImage.ImageD(kernel.getDimensions())
                ksum   = kernel.computeImage(kimage, False)
                assert (ksum == 1.)

                for rowk in range(height):
                    for coll in range(width):
                        if (rowi == rowk) and (colj == coll):
                            assert(kimage.get(coll, rowk) == 1.)
                        else:
                            assert(kimage.get(coll, rowk) == 0.)
                nk += 1

    def testSquare(self):
        self.doit(10, 10)

    def testNonSquare(self):
        self.doit(7, 10)
        
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
