from __future__ import absolute_import, division, print_function
import unittest

import lsst.utils.tests
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim
import lsst.log.utils as logUtils
import lsst.meas.algorithms as measAlg
from lsst.ip.diffim.modelPsfMatch import nextOddInteger

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
        """Test that modelPsfMatchTask raises if kernel size is too big and
        and automatic padding disabled
        """
        self.config.doAutoPadPsf = False
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
        self.config.doAutoPadPsf = False
        psfModel = measAlg.DoubleGaussianPsf(self.ksize + 2, self.ksize + 2, self.sigma2)
        psfMatch = ipDiffim.ModelPsfMatchTask(config=self.config)
        results = psfMatch.run(self.exp, psfModel)
        self.assertEqual(results.psfMatchedExposure.getPsf().computeImage().getDimensions(),
                         self.exp.getPsf().computeImage().getDimensions())
        self.assertEqual(results.psfMatchedExposure.getPsf().getSigma1(), self.sigma2)

    def testPadPsf(self):
        """Test automatic and manual PSF padding

        Compare expected Psf size, after padding, to the reference Psf size.
        The reference Psf Size is proxy for the Sciencee Psf size here.
        """
        psfModel = measAlg.DoubleGaussianPsf(self.ksize, self.ksize, self.sigma2)

        # Test automatic padding (doAutoPadPsf is True by default)
        autoPaddedKernel = nextOddInteger(self.subconfig.kernelSize*self.config.autoPadPsfTo)
        psfMatch = ipDiffim.ModelPsfMatchTask(config=self.config)
        results = psfMatch.run(self.exp, psfModel)
        self.assertEqual(results.psfMatchedExposure.getPsf().computeImage().getWidth(), autoPaddedKernel)

        # Test manual padding
        self.config.doAutoPadPsf = False
        PAD_EVEN_VALUES = [0, 2, 4]
        for padPix in PAD_EVEN_VALUES:
            self.config.padPsfBy = padPix
            psfMatch = ipDiffim.ModelPsfMatchTask(config=self.config)
            results = psfMatch.run(self.exp, psfModel)
            self.assertEqual(results.psfMatchedExposure.getPsf().computeImage().getWidth(),
                             self.ksize + padPix)

        PAD_ODD_VALUES = [1, 3, 5]
        for padPix in PAD_ODD_VALUES:
            self.config.padPsfBy = padPix
            psfMatch = ipDiffim.ModelPsfMatchTask(config=self.config)
            with self.assertRaises(ValueError):
                results = psfMatch.run(self.exp, psfModel)

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
