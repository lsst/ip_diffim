import unittest
import numpy as num


import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.geom as geom
import lsst.utils.logging as logUtils
import lsst.pex.config as pexConfig

# Increase the number for more verbose messages; decrease for fewer messages
logUtils.trace_set_at("lsst.ip.diffim", 4)


class DiffimTestCases(unittest.TestCase):

    def setUp(self):
        config = ipDiffim.PsfMatchConfigAL()

        self.ps = pexConfig.makePropertySet(config)
        self.kList = ipDiffim.makeKernelBasisList(config)

        self.ksize = self.ps['kernelSize']

    def makeSpatialKernel(self, order):
        basicGaussian1 = afwMath.GaussianFunction2D(2., 2., 0.)
        basicKernel1 = afwMath.AnalyticKernel(self.ksize, self.ksize, basicGaussian1)

        basicGaussian2 = afwMath.GaussianFunction2D(5., 3., 0.5 * num.pi)
        basicKernel2 = afwMath.AnalyticKernel(self.ksize, self.ksize, basicGaussian2)

        basisList = []
        basisList.append(basicKernel1)
        basisList.append(basicKernel2)
        basisList = ipDiffim.renormalizeKernelList(basisList)

        spatialKernelFunction = afwMath.PolynomialFunction2D(order)
        spatialKernel = afwMath.LinearCombinationKernel(basisList, spatialKernelFunction)
        kCoeffs = [[0.0 for x in range(1, spatialKernelFunction.getNParameters()+1)],
                   [0.01 * x for x in range(1, spatialKernelFunction.getNParameters()+1)]]
        kCoeffs[0][0] = 1.0  # it does not vary spatially; constant across image
        spatialKernel.setSpatialParameters(kCoeffs)
        return spatialKernel

    def tearDown(self):
        del self.ps
        del self.kList

    def testGood(self):
        ti = afwImage.MaskedImageF(geom.Extent2I(100, 100))
        ti.getVariance().set(0.1)
        ti[50, 50, afwImage.LOCAL] = (1., 0x0, 1.)
        sKernel = self.makeSpatialKernel(2)
        si = afwImage.MaskedImageF(ti.getDimensions())
        convolutionControl = afwMath.ConvolutionControl()
        convolutionControl.setDoNormalize(True)
        afwMath.convolve(si, ti, sKernel, convolutionControl)

        bbox = geom.Box2I(geom.Point2I(25, 25),
                          geom.Point2I(75, 75))
        si = afwImage.MaskedImageF(si, bbox, origin=afwImage.LOCAL)
        ti = afwImage.MaskedImageF(ti, bbox, origin=afwImage.LOCAL)
        kc = ipDiffim.KernelCandidateF(50., 50., ti, si, self.ps)

        sBg = afwMath.PolynomialFunction2D(1)
        bgCoeffs = [0., 0., 0.]
        sBg.setParameters(bgCoeffs)

        # must be initialized
        bskv = ipDiffim.BuildSingleKernelVisitorF(self.kList, self.ps)
        bskv.processCandidate(kc)
        self.assertEqual(kc.isInitialized(), True)

        askv = ipDiffim.AssessSpatialKernelVisitorF(sKernel, sBg, self.ps)
        askv.processCandidate(kc)

        self.assertEqual(askv.getNProcessed(), 1)
        self.assertEqual(askv.getNRejected(), 0)
        self.assertEqual(kc.getStatus(), afwMath.SpatialCellCandidate.GOOD)

    def testBad(self):
        ti = afwImage.MaskedImageF(geom.Extent2I(100, 100))
        ti.getVariance().set(0.1)
        ti[50, 50, afwImage.LOCAL] = (1., 0x0, 1.)
        sKernel = self.makeSpatialKernel(2)
        si = afwImage.MaskedImageF(ti.getDimensions())
        convolutionControl = afwMath.ConvolutionControl()
        convolutionControl.setDoNormalize(True)
        afwMath.convolve(si, ti, sKernel, convolutionControl)

        bbox = geom.Box2I(geom.Point2I(25, 25),
                          geom.Point2I(75, 75))
        si = afwImage.MaskedImageF(si, bbox, origin=afwImage.LOCAL)
        ti = afwImage.MaskedImageF(ti, bbox, origin=afwImage.LOCAL)
        kc = ipDiffim.KernelCandidateF(50., 50., ti, si, self.ps)

        badGaussian = afwMath.GaussianFunction2D(1., 1., 0.)
        badKernel = afwMath.AnalyticKernel(self.ksize, self.ksize, badGaussian)
        basisList = []
        basisList.append(badKernel)
        badSpatialKernelFunction = afwMath.PolynomialFunction2D(0)
        badSpatialKernel = afwMath.LinearCombinationKernel(basisList, badSpatialKernelFunction)
        badSpatialKernel.setSpatialParameters([[1, ]])

        sBg = afwMath.PolynomialFunction2D(1)
        bgCoeffs = [10., 10., 10.]
        sBg.setParameters(bgCoeffs)

        # must be initialized
        bskv = ipDiffim.BuildSingleKernelVisitorF(self.kList, self.ps)
        bskv.processCandidate(kc)
        self.assertEqual(kc.isInitialized(), True)

        askv = ipDiffim.AssessSpatialKernelVisitorF(badSpatialKernel, sBg, self.ps)
        askv.processCandidate(kc)

        self.assertEqual(askv.getNProcessed(), 1)
        self.assertEqual(askv.getNRejected(), 1)
        self.assertEqual(kc.getStatus(), afwMath.SpatialCellCandidate.BAD)


#####

class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
