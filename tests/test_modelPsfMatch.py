import unittest

import lsst.utils.tests
import lsst.daf.base as dafBase
from lsst.afw.coord import Observatory, Weather
import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.ip.diffim as ipDiffim
import lsst.utils.logging as logUtils
import lsst.meas.algorithms as measAlg
from lsst.ip.diffim.modelPsfMatch import nextOddInteger

logUtils.trace_set_at("lsst.ip.diffim", 4)


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
        self.exp = afwImage.ExposureF(geom.Extent2I(self.imsize, self.imsize))
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
        self.assertEqual(
            results.psfMatchedExposure.getPsf().computeImage(
                results.psfMatchedExposure.getPsf().getAveragePosition()
            ).getDimensions(),
            self.exp.getPsf().computeImage(
                self.exp.getPsf().getAveragePosition()
            ).getDimensions()
        )
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
        self.assertEqual(
            results.psfMatchedExposure.getPsf().computeImage(
                results.psfMatchedExposure.getPsf().getAveragePosition()
            ).getWidth(),
            autoPaddedKernel
        )

        # Test manual padding
        self.config.doAutoPadPsf = False
        PAD_EVEN_VALUES = [0, 2, 4]
        for padPix in PAD_EVEN_VALUES:
            self.config.padPsfBy = padPix
            psfMatch = ipDiffim.ModelPsfMatchTask(config=self.config)
            results = psfMatch.run(self.exp, psfModel)
            self.assertEqual(
                results.psfMatchedExposure.getPsf().computeImage(
                    results.psfMatchedExposure.getPsf().getAveragePosition()
                ).getWidth(),
                self.ksize + padPix
            )

        PAD_ODD_VALUES = [1, 3, 5]
        for padPix in PAD_ODD_VALUES:
            self.config.padPsfBy = padPix
            psfMatch = ipDiffim.ModelPsfMatchTask(config=self.config)
            with self.assertRaises(ValueError):
                results = psfMatch.run(self.exp, psfModel)

    def testPropagateVisitInfo(self):
        """Test that a PSF-matched exposure preserves the original VisitInfo.
        """
        self.exp.getInfo().setVisitInfo(makeVisitInfo())
        psfModel = measAlg.DoubleGaussianPsf(self.ksize + 2, self.ksize + 2, self.sigma2)
        psfMatch = ipDiffim.ModelPsfMatchTask(config=self.config)
        psfMatchedExposure = psfMatch.run(self.exp, psfModel).psfMatchedExposure
        self.assertEqual(psfMatchedExposure.getInfo().getVisitInfo(),
                         self.exp.getInfo().getVisitInfo())

    def tearDown(self):
        del self.exp
        del self.subconfig


def makeVisitInfo():
    """Return a non-NaN visitInfo."""
    return afwImage.VisitInfo(exposureId=10313423,
                              exposureTime=10.01,
                              darkTime=11.02,
                              date=dafBase.DateTime(65321.1, dafBase.DateTime.MJD, dafBase.DateTime.TAI),
                              ut1=12345.1,
                              era=45.1*geom.degrees,
                              boresightRaDec=geom.SpherePoint(23.1, 73.2, geom.degrees),
                              boresightAzAlt=geom.SpherePoint(134.5, 33.3, geom.degrees),
                              boresightAirmass=1.73,
                              boresightRotAngle=73.2*geom.degrees,
                              rotType=afwImage.RotType.SKY,
                              observatory=Observatory(
                                  11.1*geom.degrees, 22.2*geom.degrees, 0.333),
                              weather=Weather(1.1, 2.2, 34.5),
                              )


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
