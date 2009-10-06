#!/usr/bin/env python
import os
import numpy as num

import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.pex.policy as pexPolicy
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as logging
import lsst.afw.image.testUtils as imTestUtils

Verbosity = 4
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = pexPolicy.Policy.createPolicy(diffimPolicy)
        self.kCols  = self.policy.getInt('kernelCols')
        self.kRows  = self.policy.getInt('kernelRows')

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

    def testAlardLupton(self):
        nGauss   = self.policy.get("alardNGauss")
        sigGauss = self.policy.getDoubleArray("alardSigGauss")
        degGauss = self.policy.getIntArray("alardDegGauss")
        assert len(sigGauss) == nGauss
        assert len(degGauss) == nGauss
        assert self.kCols == self.kRows  # square
        assert self.kCols % 2 == 1       # odd sized
        kHalfWidth = int(self.kCols/2)
        self.basisList  = ipDiffim.generateAlardLuptonKernelSet(kHalfWidth, nGauss, sigGauss, degGauss)
        kim1 = afwImage.ImageD(self.basisList[0].getDimensions())
        nBasis = len(self.basisList)

        # the first one sums to 1; the rest sum to 0
        self.basisList[0].computeImage(kim1, False)
        self.assertAlmostEqual(num.sum(num.ravel(imTestUtils.arrayFromImage(kim1))), 1.0)
        for k1 in range(1, nBasis):
            self.basisList[k1].computeImage(kim1, False)
            self.assertAlmostEqual(num.sum(num.ravel(imTestUtils.arrayFromImage(kim1))), 0.0)

        # the images dotted with themselves is 1, except for the first
        for k1 in range(1, nBasis):
            self.basisList[k1].computeImage(kim1, False)
            arr1 = imTestUtils.arrayFromImage(kim1)
            self.assertAlmostEqual(num.sum(arr1*arr1), 1.0)
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
