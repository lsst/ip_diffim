#!/usr/bin/env python
import os
import numpy as num

import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as logging
import lsst.afw.image.testUtils as imTestUtils

verbosity = 4
logging.Trace_setVerbosity('lsst.ip.diffim', verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = ipDiffim.generateDefaultPolicy(diffimPolicy, modify=False)
        self.kCols  = self.policy.getInt('kernelCols')
        self.kRows  = self.policy.getInt('kernelRows')

    def tearDown(self):
        del self.policy

    def runDeltaFunction(self, width=10, height=10):
        ks1 = ipDiffim.generateDeltaFunctionBasisSet(width, height)
        nk  = 0
        for rowi in range(height):
            for colj in range(width):
                kernel = ks1[nk]
                kimage = afwImage.ImageD(kernel.getDimensions())
                ksum   = kernel.computeImage(kimage, False)
                self.assertEqual(ksum, 1.)

                for rowk in range(height):
                    for coll in range(width):
                        if (rowi == rowk) and (colj == coll):
                            self.assertEqual(kimage.get(coll, rowk), 1.)
                        else:
                            self.assertEqual(kimage.get(coll, rowk), 0.)
                nk += 1

    def testSquareDeltaFunction(self):
        self.runDeltaFunction(10, 10)

    def testNonSquareDeltaFunction(self):
        self.runDeltaFunction(7, 10)

    def testRegularization(self):

        unwrap = ipDiffim.UNWRAPPED # 0
        wrap = ipDiffim.WRAPPED     # 1
        taper = ipDiffim.TAPERED    # 2
        forw = ipDiffim.FORWARD_DIFFERENCE # 0
        cent = ipDiffim.CENTRAL_DIFFERENCE # 1
        
        # test order
        order = -1
        try:
            h = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, order, unwrap, forw)
        except Exception, e:
            pass
        else:
            self.fail()

        order = 0
        try:
            h = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, order, unwrap, forw)
        except Exception, e:
            self.fail()

        order = 3
        try:
            h = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, order, unwrap, forw)
        except Exception, e:
            pass
        else:
            self.fail()

        # test boundary_style
        boundary = -1 # too low
        try:
            h = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, 0, boundary, forw)
        except Exception, e:
            pass
        else:
            self.fail()        

        boundary = wrap
        try:
            h = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, 0, boundary, forw)
        except Exception, e:
            self.fail()        

        boundary = 3  # too high
        try:
            h = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, 0, boundary, forw)
        except Exception, e:
            pass
        else:
            self.fail()
            
        # test difference_style
        difference = -1 # too low
        try:
            h = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, 0, 0, difference)
        except Exception, e:
            pass
        else:
            self.fail()

        difference = forw
        try:
            h = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, 0, unwrap, difference)
        except Exception, e:
            self.fail()

        difference = 2 # too high
        try:
            h = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, 0, unwrap, difference)
        except Exception, e:
            pass
        else:
            self.fail()
        
            
    def testRenormalize(self):
        # inputs
        gauss1 = afwMath.GaussianFunction2D(2, 2)
        gauss2 = afwMath.GaussianFunction2D(3, 4)
        gauss3 = afwMath.GaussianFunction2D(0.2, 0.25)
        gaussKernel1 = afwMath.AnalyticKernel(self.kCols, self.kRows, gauss1)
        gaussKernel2 = afwMath.AnalyticKernel(self.kCols, self.kRows, gauss2)
        gaussKernel3 = afwMath.AnalyticKernel(self.kCols, self.kRows, gauss3)
        kimage1 = afwImage.ImageD(gaussKernel1.getDimensions())
        ksum1   = gaussKernel1.computeImage(kimage1, False)
        kimage2 = afwImage.ImageD(gaussKernel2.getDimensions())
        ksum2   = gaussKernel2.computeImage(kimage2, False)
        kimage3 = afwImage.ImageD(gaussKernel3.getDimensions())
        ksum3   = gaussKernel3.computeImage(kimage3, False)
        self.assertTrue(ksum1 != 1.)
        self.assertTrue(ksum2 != 1.)
        self.assertTrue(ksum3 != 1.)
        # no constraints on first kernels norm
        self.assertTrue(num.sum( num.ravel(ipDiffim.vectorFromImage(kimage2))**2 ) != 1.)
        self.assertTrue(num.sum( num.ravel(ipDiffim.vectorFromImage(kimage3))**2 ) != 1.)
        basisListIn = afwMath.KernelList()
        basisListIn.push_back(gaussKernel1)
        basisListIn.push_back(gaussKernel2)
        basisListIn.push_back(gaussKernel3)
        
        # outputs
        basisListOut = ipDiffim.renormalizeKernelList(basisListIn)
        gaussKernel1 = basisListOut[0]
        gaussKernel2 = basisListOut[1]
        gaussKernel3 = basisListOut[2]
        ksum1 = gaussKernel1.computeImage(kimage1, False)
        ksum2 = gaussKernel2.computeImage(kimage2, False)
        ksum3 = gaussKernel3.computeImage(kimage3, False)
        self.assertAlmostEqual(ksum1, 1.)
        self.assertAlmostEqual(ksum2, 0.)
        self.assertAlmostEqual(ksum3, 0.)
        # no constraints on first kernels norm
        self.assertAlmostEqual(num.sum( num.ravel(ipDiffim.vectorFromImage(kimage2))**2 ), 1.)
        self.assertAlmostEqual(num.sum( num.ravel(ipDiffim.vectorFromImage(kimage3))**2 ), 1.)

    def testCreateDeltaFunction(self):
        self.policy.set("kernelBasisSet", "delta-function")
        kFunctor = ipDiffim.createKernelFunctor(self.policy)
        basisList = kFunctor.getBasisList()

        nk  = 0
        for rowi in range(self.kRows):
            for colj in range(self.kCols):
                kernel = basisList[nk]
                kimage = afwImage.ImageD(kernel.getDimensions())
                ksum   = kernel.computeImage(kimage, False)
                assert (ksum == 1.)

                for rowk in range(self.kRows):
                    for coll in range(self.kCols):
                        if (rowi == rowk) and (colj == coll):
                            assert(kimage.get(coll, rowk) == 1.)
                        else:
                            assert(kimage.get(coll, rowk) == 0.)
                nk += 1

    ##### Alard Lupton

    def testAlardLupton(self):
        nGauss   = self.policy.get("alardNGauss")
        sigGauss = self.policy.getDoubleArray("alardSigGauss")
        degGauss = self.policy.getIntArray("alardDegGauss")
        assert len(sigGauss) == nGauss
        assert len(degGauss) == nGauss
        assert self.kCols == self.kRows  # square
        assert self.kCols % 2 == 1       # odd sized
        kHalfWidth = int(self.kCols/2)
        basisList  = ipDiffim.generateAlardLuptonBasisSet(kHalfWidth, nGauss, sigGauss, degGauss)
        kim1 = afwImage.ImageD(basisList[0].getDimensions())
        nBasis = len(basisList)

        # the first one sums to 1; the rest sum to 0
        basisList[0].computeImage(kim1, False)
        self.assertAlmostEqual(num.sum(num.ravel(imTestUtils.arrayFromImage(kim1))), 1.0)
        for k1 in range(1, nBasis):
            basisList[k1].computeImage(kim1, False)
            self.assertAlmostEqual(num.sum(num.ravel(imTestUtils.arrayFromImage(kim1))), 0.0)

        # the images dotted with themselves is 1, except for the first
        for k1 in range(1, nBasis):
            basisList[k1].computeImage(kim1, False)
            arr1 = imTestUtils.arrayFromImage(kim1)
            self.assertAlmostEqual(num.sum(arr1*arr1), 1.0)

    def testCreateAlardLupton(self):
        self.policy.set("kernelBasisSet", "alard-lupton")
        kFunctor = ipDiffim.createKernelFunctor(self.policy)
        basisList = kFunctor.getBasisList()

        kim1 = afwImage.ImageD(basisList[0].getDimensions())
        nBasis = len(basisList)

        # the first one sums to 1; the rest sum to 0
        basisList[0].computeImage(kim1, False)
        self.assertAlmostEqual(num.sum(num.ravel(imTestUtils.arrayFromImage(kim1))), 1.0)
        for k1 in range(1, nBasis):
            basisList[k1].computeImage(kim1, False)
            self.assertAlmostEqual(num.sum(num.ravel(imTestUtils.arrayFromImage(kim1))), 0.0)

        # the images dotted with themselves is 1, except for the first
        for k1 in range(1, nBasis):
            basisList[k1].computeImage(kim1, False)
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

def run(doExit=False):
    """Run the tests"""
    tests.run(suite(), doExit)

if __name__ == "__main__":
    run(True)
