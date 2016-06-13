#!/usr/bin/env python
import numpy as num

import unittest
import lsst.utils.tests as tests

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as logging
import lsst.pex.config as pexConfig

verbosity = 4
logging.Trace_setVerbosity("lsst.ip.diffim", verbosity)

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.configAL    = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.configAL.kernel.name = "AL"
        self.subconfigAL = self.configAL.kernel.active

        self.configDF    = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.configDF.kernel.name = "DF"
        self.subconfigDF = self.configDF.kernel.active

        self.policyAL = pexConfig.makePolicy(self.subconfigAL)
        self.policyDF = pexConfig.makePolicy(self.subconfigDF)

        self.kSize    = self.policyAL.getInt("kernelSize")

    def tearDown(self):
        del self.configAL
        del self.policyAL
        del self.configDF
        del self.policyDF

    #
    ### Delta function
    #
    
    def deltaFunctionTest(self, ks):
        # right shape
        nk  = 0
        for rowi in range(self.kSize):
            for colj in range(self.kSize):
                kernel = ks[nk]
                kimage = afwImage.ImageD(kernel.getDimensions())
                ksum   = kernel.computeImage(kimage, False)
                self.assertEqual(ksum, 1.)

                for rowk in range(self.kSize):
                    for coll in range(self.kSize):
                        if (rowi == rowk) and (colj == coll):
                            self.assertEqual(kimage.get(coll, rowk), 1.)
                        else:
                            self.assertEqual(kimage.get(coll, rowk), 0.)
                nk += 1
        

    def testMakeDeltaFunction(self):
        ks = ipDiffim.makeDeltaFunctionBasisList(self.kSize, self.kSize)

        # right size
        self.assertEqual(len(ks), self.kSize * self.kSize)

        # right shape
        self.deltaFunctionTest(ks)


    #
    ### Alard Lupton
    #

    def alardLuptonTest(self, ks):
        kim    = afwImage.ImageD(ks[0].getDimensions())
        nBasis = len(ks)

        # the first one sums to 1; the rest sum to 0
        ks[0].computeImage(kim, False)
        self.assertAlmostEqual(num.sum(num.ravel(kim.getArray())), 1.0)
        
        for k in range(1, nBasis):
            ks[k].computeImage(kim, False)
            self.assertAlmostEqual(num.sum(num.ravel(kim.getArray())), 0.0)

        # the images dotted with themselves is 1, except for the first
        for k in range(1, nBasis):
            ks[k].computeImage(kim, False)
            arr = kim.getArray()
            self.assertAlmostEqual(num.sum(arr*arr), 1.0)

    def testMakeAlardLupton(self):
        nGauss   = self.policyAL.get("alardNGauss")
        sigGauss = self.policyAL.getDoubleArray("alardSigGauss")
        degGauss = self.policyAL.getIntArray("alardDegGauss")
        self.assertEqual(len(sigGauss), nGauss)
        self.assertEqual(len(degGauss), nGauss)
        self.assertEqual(self.kSize % 2, 1)       # odd sized
        kHalfWidth = self.kSize // 2

        ks = ipDiffim.makeAlardLuptonBasisList(kHalfWidth, nGauss, sigGauss, degGauss)

        # right size
        nTot = 0
        for deg in degGauss:
            nTot += (deg + 1) * (deg + 2) / 2
        self.assertEqual(len(ks), nTot)

        # right orthogonality
        self.alardLuptonTest(ks)

    def testGenerateAlardLupton(self):
        # defaults
        ks = ipDiffim.generateAlardLuptonBasisList(self.subconfigAL)
        self.alardLuptonTest(ks)

        # send FWHM
        ks = ipDiffim.generateAlardLuptonBasisList(self.subconfigAL, targetFwhmPix=3.0,referenceFwhmPix=4.0)
        self.alardLuptonTest(ks)

        
    #
    ### Make
    #
    def testMakeKernelBasisList(self):
        ks = ipDiffim.makeKernelBasisList(self.subconfigAL)
        self.alardLuptonTest(ks)

        ks = ipDiffim.makeKernelBasisList(self.subconfigDF)
        self.deltaFunctionTest(ks)
        
    #
    ### Renormalize
    #

    def testRenormalize(self):
        # inputs
        gauss1 = afwMath.GaussianFunction2D(2, 2)
        gauss2 = afwMath.GaussianFunction2D(3, 4)
        gauss3 = afwMath.GaussianFunction2D(0.2, 0.25)
        gaussKernel1 = afwMath.AnalyticKernel(self.kSize, self.kSize, gauss1)
        gaussKernel2 = afwMath.AnalyticKernel(self.kSize, self.kSize, gauss2)
        gaussKernel3 = afwMath.AnalyticKernel(self.kSize, self.kSize, gauss3)
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
        self.assertTrue(num.sum(num.ravel(kimage2.getArray())**2 ) != 1.)
        self.assertTrue(num.sum(num.ravel(kimage3.getArray())**2 ) != 1.)
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
        self.assertAlmostEqual(num.sum(num.ravel(kimage2.getArray())**2 ), 1.)
        self.assertAlmostEqual(num.sum(num.ravel(kimage3.getArray())**2 ), 1.)

    #
    ### Regularize
    #

    def testCentralRegularization(self):
        #
        self.policyDF.set("regularizationType", "centralDifference")
        try:
            self.policyDF.set("centralRegularizationStencil", 1)
            ipDiffim.makeRegularizationMatrix(self.policyDF)
        except Exception, e: # stencil of 1 not allowed
            pass
        else:
            self.fail()
            
        self.policyDF.set("centralRegularizationStencil", 5)
        try:
            ipDiffim.makeRegularizationMatrix(self.policyDF)
        except Exception, e: # stencil of 5 allowed
            print e
            self.fail()
        else:
            pass

        self.policyDF.set("centralRegularizationStencil", 9)
        try:
            ipDiffim.makeRegularizationMatrix(self.policyDF)
        except Exception, e: # stencil of 9 allowed
            print e
            self.fail()
        else:
            pass

        self.policyDF.set("regularizationBorderPenalty", -1.0)
        try:
            ipDiffim.makeRegularizationMatrix(self.policyDF)
        except Exception, e: # border penalty > 0
            pass
        else:
            self.fail()

        self.policyDF.set("regularizationBorderPenalty", 0.0)
        try:
            ipDiffim.makeRegularizationMatrix(self.policyDF)
        except Exception, e: # border penalty > 0
            print e
            self.fail()
        else:
            pass

    def testForwardRegularization(self):
        self.policyDF.set("regularizationType", "forwardDifference")
        self.policyDF.set("forwardRegularizationOrders", 0)
        try:
            ipDiffim.makeRegularizationMatrix(self.policyDF)
        except Exception, e: # order 1..3 allowed
            pass
        else:
            self.fail()

        self.policyDF.set("forwardRegularizationOrders", 1)
        try:
            ipDiffim.makeRegularizationMatrix(self.policyDF)
        except Exception, e: # order 1..3 allowed
            print e
            self.fail()
        else:
            pass

        self.policyDF.set("forwardRegularizationOrders", 4)
        try:
            ipDiffim.makeRegularizationMatrix(self.policyDF)
        except Exception, e: # order 1..3 allowed
            pass
        else:
            self.fail()

        self.policyDF.set("forwardRegularizationOrders", 1)
        self.policyDF.add("forwardRegularizationOrders", 2)
        try:
            ipDiffim.makeRegularizationMatrix(self.policyDF)
        except Exception as e: # order 1..3 allowed
            print e
            self.fail()
        else:
            pass

    def testBadRegularization(self):
        try:
            self.policyDF.set("regularizationType", "foo")
            ipDiffim.makeRegularizationMatrix(self.policyDF)
        except Exception: # invalid option
            pass
        else:
            self.fail()

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
