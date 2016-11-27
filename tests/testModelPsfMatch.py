#!/usr/bin/env python
import unittest
import lsst.utils.tests
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim
import lsst.log.utils as logUtils
import lsst.meas.algorithms as measAlg

logUtils.traceSetAt("ip.diffim", 4)


class PsfMatchTestCases(lsst.utils.tests.TestCase):

    def setUp(self):
        self.config = ipDiffim.ModelPsfMatchTask.ConfigClass()
        self.subconfig = self.config.kernel.active
        self.subconfig.scaleByFwhm = False
        self.subconfig.kernelSize = 9
        self.subconfig.kernelSizeMin = 9

        self.imsize = 2 * self.subconfig.sizeCellX
        self.ksize = 11
        self.sigma1 = 2.0
        self.sigma2 = 3.7
        self.exp = afwImage.ExposureF(afwGeom.Extent2I(self.imsize, self.imsize))
        self.exp.setPsf(measAlg.DoubleGaussianPsf(self.ksize, self.ksize, self.sigma1))

    def testTooBig(self):
        self.subconfig.kernelSize = self.ksize
        psf = measAlg.DoubleGaussianPsf(self.ksize, self.ksize, self.sigma2)
        psfMatch = ipDiffim.ModelPsfMatchTask(config=self.config)
        try:
            psfMatch.run(self.exp, psf)
        except Exception:
            pass
        else:
            self.fail()

    def testMatch(self):
        for order in (0, 1):
            for ksum in (0.5, 1.0, 2.7):
                self.runMatch(kOrder=order, kSumIn=ksum)

    def runMatch(self, kOrder=0, kSumIn=3.7):
        self.subconfig.spatialKernelOrder = kOrder

        psf = measAlg.DoubleGaussianPsf(self.ksize, self.ksize, self.sigma2)
        psfMatch = ipDiffim.ModelPsfMatchTask(config=self.config)
        results = psfMatch.run(self.exp, psf, kernelSum=kSumIn)

        matchingKernel = results.psfMatchingKernel

        kImage = afwImage.ImageD(matchingKernel.getDimensions())
        kSumOut = matchingKernel.computeImage(kImage, False)

        self.assertAlmostEqual(kSumIn, kSumOut)

    def testAdjustModelSize(self):
        """Test that modelPsfMatch correctly adjusts the model PSF dimensions to
        match those of the science PSF.
        """
        psfModel = measAlg.DoubleGaussianPsf(self.ksize + 2, self.ksize + 2, self.sigma2)
        psfMatch = ipDiffim.ModelPsfMatchTask(config=self.config)
        results = psfMatch.run(self.exp, psfModel)
        self.assertEqual(results.referencePsfModel.computeImage().getDimensions(),
                         self.exp.getPsf().computeImage().getDimensions())
        self.assertEqual(results.referencePsfModel.getSigma1(), self.sigma2)

    def tearDown(self):
        del self.exp
        del self.subconfig


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
