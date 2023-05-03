import unittest

import numpy as num

import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.utils.logging as logUtils
import lsst.pex.config as pexConfig
import lsst.pex.exceptions

from lsst.ip.diffim import PsfMatchConfigAL, PsfMatchConfigDF

# Increase the number for more verbose messages
logUtils.trace_set_at("lsst.ip.diffim", 0)


class DiffimTestCases(unittest.TestCase):

    def setUp(self):
        self.configAL = PsfMatchConfigAL()

        self.configDF = PsfMatchConfigDF()

        self.psAL = pexConfig.makePropertySet(self.configAL)
        self.psDF = pexConfig.makePropertySet(self.configDF)

        self.kSize = self.psAL["kernelSize"]

    def tearDown(self):
        del self.configAL
        del self.psAL
        del self.configDF
        del self.psDF

    def deltaFunctionTest(self, ks):
        # right shape
        nk = 0
        for rowi in range(self.kSize):
            for colj in range(self.kSize):
                kernel = ks[nk]
                kimage = afwImage.ImageD(kernel.getDimensions())
                ksum = kernel.computeImage(kimage, False)
                self.assertEqual(ksum, 1.)

                for rowk in range(self.kSize):
                    for coll in range(self.kSize):
                        if (rowi == rowk) and (colj == coll):
                            self.assertEqual(kimage[coll, rowk, afwImage.LOCAL], 1.)
                        else:
                            self.assertEqual(kimage[coll, rowk, afwImage.LOCAL], 0.)
                nk += 1

    def testMakeDeltaFunction(self):
        ks = ipDiffim.makeDeltaFunctionBasisList(self.kSize, self.kSize)

        # right size
        self.assertEqual(len(ks), self.kSize * self.kSize)

        # right shape
        self.deltaFunctionTest(ks)

    def alardLuptonTest(self, ks):
        kim = afwImage.ImageD(ks[0].getDimensions())
        nBasis = len(ks)

        # the first one sums to 1; the rest sum to 0
        ks[0].computeImage(kim, False)
        self.assertAlmostEqual(num.sum(num.ravel(kim.array)), 1.0)

        for k in range(1, nBasis):
            ks[k].computeImage(kim, False)
            self.assertAlmostEqual(num.sum(num.ravel(kim.array)), 0.0)

        # the images dotted with themselves is 1, except for the first
        for k in range(1, nBasis):
            ks[k].computeImage(kim, False)
            arr = kim.array
            self.assertAlmostEqual(num.sum(arr*arr), 1.0)

    def testMakeAlardLupton(self):
        nGauss = self.psAL["alardNGauss"]
        sigGauss = self.psAL.getArray("alardSigGauss")
        degGauss = self.psAL.getArray("alardDegGauss")
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
        ks = ipDiffim.generateAlardLuptonBasisList(self.configAL)
        self.alardLuptonTest(ks)

        # send FWHM
        ks = ipDiffim.generateAlardLuptonBasisList(self.configAL, targetFwhmPix=3.0, referenceFwhmPix=4.0)
        self.alardLuptonTest(ks)

    def testMakeKernelBasisList(self):
        ks = ipDiffim.makeKernelBasisList(self.configAL)
        self.alardLuptonTest(ks)

        ks = ipDiffim.makeKernelBasisList(self.configDF)
        self.deltaFunctionTest(ks)

    def testRenormalize(self):
        # inputs
        gauss1 = afwMath.GaussianFunction2D(2, 2)
        gauss2 = afwMath.GaussianFunction2D(3, 4)
        gauss3 = afwMath.GaussianFunction2D(0.2, 0.25)
        gaussKernel1 = afwMath.AnalyticKernel(self.kSize, self.kSize, gauss1)
        gaussKernel2 = afwMath.AnalyticKernel(self.kSize, self.kSize, gauss2)
        gaussKernel3 = afwMath.AnalyticKernel(self.kSize, self.kSize, gauss3)
        kimage1 = afwImage.ImageD(gaussKernel1.getDimensions())
        ksum1 = gaussKernel1.computeImage(kimage1, False)
        kimage2 = afwImage.ImageD(gaussKernel2.getDimensions())
        ksum2 = gaussKernel2.computeImage(kimage2, False)
        kimage3 = afwImage.ImageD(gaussKernel3.getDimensions())
        ksum3 = gaussKernel3.computeImage(kimage3, False)
        self.assertNotEqual(ksum1, 1.)
        self.assertNotEqual(ksum2, 1.)
        self.assertNotEqual(ksum3, 1.)
        # no constraints on first kernels norm
        self.assertNotEqual(num.sum(num.ravel(kimage2.array)**2), 1.)
        self.assertNotEqual(num.sum(num.ravel(kimage3.array)**2), 1.)
        basisListIn = []
        basisListIn.append(gaussKernel1)
        basisListIn.append(gaussKernel2)
        basisListIn.append(gaussKernel3)

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
        self.assertAlmostEqual(num.sum(num.ravel(kimage2.array)**2), 1.)
        self.assertAlmostEqual(num.sum(num.ravel(kimage3.array)**2), 1.)

    def testCentralRegularization(self):
        # stencil of 1 not allowed
        self.psDF["regularizationType"] = "centralDifference"
        with self.assertRaises(lsst.pex.exceptions.Exception):
            self.psDF["centralRegularizationStencil"] = 1
            ipDiffim.makeRegularizationMatrix(self.psDF)

        # stencil of 5 allowed
        self.psDF["centralRegularizationStencil"] = 5
        try:
            ipDiffim.makeRegularizationMatrix(self.psDF)
        except lsst.pex.exceptions.Exception as e:
            self.fail("Should not raise %s: stencil of 5 is allowed."%e)

        # stencil of 9 allowed
        self.psDF["centralRegularizationStencil"] = 9
        try:
            ipDiffim.makeRegularizationMatrix(self.psDF)
        except lsst.pex.exceptions.Exception as e:
            self.fail("Should not raise %s: stencil of 9 is allowed"%e)

        # border penalty < 0
        self.psDF["regularizationBorderPenalty"] = -1.0
        with self.assertRaises(lsst.pex.exceptions.Exception):
            ipDiffim.makeRegularizationMatrix(self.psDF)

        # border penalty > 0
        self.psDF["regularizationBorderPenalty"] = 0.0
        try:
            ipDiffim.makeRegularizationMatrix(self.psDF)
        except lsst.pex.exceptions.Exception as e:
            self.fail("Should not raise %s: Border penalty > 0"%e)

    def testForwardRegularization(self):
        self.psDF["regularizationType"] = "forwardDifference"

        # order 1..3 allowed
        self.psDF["forwardRegularizationOrders"] = 0
        with self.assertRaises(lsst.pex.exceptions.Exception):
            ipDiffim.makeRegularizationMatrix(self.psDF)

        self.psDF["forwardRegularizationOrders"] = 1
        try:
            ipDiffim.makeRegularizationMatrix(self.psDF)
        except lsst.pex.exceptions.Exception as e:
            self.fail("Should not raise %s: order 1 allowed"%e)

        self.psDF["forwardRegularizationOrders"] = 4
        with self.assertRaises(lsst.pex.exceptions.Exception):
            ipDiffim.makeRegularizationMatrix(self.psDF)

        self.psDF["forwardRegularizationOrders"] = [1, 2]
        try:
            ipDiffim.makeRegularizationMatrix(self.psDF)
        except lsst.pex.exceptions.Exception as e:
            self.fail("Should not raise %s: order 1,2 allowed"%e)

    def testBadRegularization(self):
        with self.assertRaises(lsst.pex.exceptions.Exception):
            self.psDF["regularizationType"] = "foo"
            ipDiffim.makeRegularizationMatrix(self.psDF)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
