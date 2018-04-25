import unittest


import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.log.utils as logUtils
import lsst.pex.config as pexConfig

logUtils.traceSetAt("ip.diffim", 4)


class DiffimTestCases(lsst.utils.tests.TestCase):

    def setUp(self):
        self.config = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.config.kernel.name = "DF"
        self.subconfig = self.config.kernel.active

        self.kList = ipDiffim.makeKernelBasisList(self.subconfig)
        self.policy = pexConfig.makePolicy(self.subconfig)
        self.policy.set("useRegularization", False)

    def tearDown(self):
        del self.config
        del self.policy
        del self.kList

    def makeCandidate(self, kSum, x, y, size=51):
        mi1 = afwImage.MaskedImageF(afwGeom.Extent2I(size, size))
        mi1.getVariance().set(1.0)  # avoid NaNs
        mi1.set(size//2, size//2, (1, 0x0, 1))
        mi2 = afwImage.MaskedImageF(afwGeom.Extent2I(size, size))
        mi2.getVariance().set(1.0)  # avoid NaNs
        mi2.set(size//2, size//2, (kSum, 0x0, kSum))
        kc = ipDiffim.makeKernelCandidate(x, y, mi1, mi2, self.policy)
        return kc

    def testGaussian(self, size=51):
        gaussFunction = afwMath.GaussianFunction2D(2, 3)
        gaussKernel = afwMath.AnalyticKernel(size, size, gaussFunction)

        imagePca1 = ipDiffim.KernelPcaD()  # mean subtract
        imagePca2 = ipDiffim.KernelPcaD()  # don't mean subtract
        kpv1 = ipDiffim.KernelPcaVisitorF(imagePca1)
        kpv2 = ipDiffim.KernelPcaVisitorF(imagePca2)

        kRefIm = None

        for i in range(100):
            kImage1 = afwImage.ImageD(gaussKernel.getDimensions())
            gaussKernel.computeImage(kImage1, False)
            kImage1 *= 10000  # to get some decent peak source counts
            kImage1 += 10     # to get some sky background noise

            if kRefIm is None:
                kRefIm = kImage1

            kImage1 = diffimTools.makePoissonNoiseImage(kImage1)
            kImage2 = afwImage.ImageD(kImage1, True)

            imagePca1.addImage(kImage1, 1.0)
            imagePca2.addImage(kImage2, 1.0)

        kpv1.subtractMean()

        imagePca1.analyze()
        imagePca2.analyze()

        pcaBasisList1 = kpv1.getEigenKernels()
        pcaBasisList2 = kpv2.getEigenKernels()

        eVal1 = imagePca1.getEigenValues()
        eVal2 = imagePca2.getEigenValues()

        # First term is far more signficant without mean subtraction
        self.assertGreater(eVal2[0], eVal1[0])

        # Last term basically zero with mean subtraction
        self.assertAlmostEqual(eVal1[-1], 0.0)

        # Extra image with mean subtraction
        self.assertEqual(len(pcaBasisList1), (len(eVal1) + 1))

        # Same shape
        self.assertEqual(len(pcaBasisList2), len(eVal2))

        # Mean kernel close to kRefIm
        kImageM = afwImage.ImageD(gaussKernel.getDimensions())
        pcaBasisList1[0].computeImage(kImageM, False)
        for y in range(kRefIm.getHeight()):
            for x in range(kRefIm.getWidth()):
                self.assertLess(abs(kRefIm.get(x, y) - kImageM.get(x, y)) / kRefIm.get(x, y), 0.2)

        # First mean-unsubtracted Pca kernel close to kRefIm (normalized to peak of 1.0)
        kImage0 = afwImage.ImageD(gaussKernel.getDimensions())
        pcaBasisList2[0].computeImage(kImage0, False)
        maxVal = afwMath.makeStatistics(kRefIm, afwMath.MAX).getValue(afwMath.MAX)
        kRefIm /= maxVal
        for y in range(kRefIm.getHeight()):
            for x in range(kRefIm.getWidth()):
                self.assertLess(abs(kRefIm.get(x, y) - kImage0.get(x, y)) / kRefIm.get(x, y), 0.2)

    def testImagePca(self):
        # Test out the ImagePca behavior
        kc1 = self.makeCandidate(1, 0.0, 0.0)
        kc1.build(self.kList)
        kc2 = self.makeCandidate(2, 0.0, 0.0)
        kc2.build(self.kList)
        kc3 = self.makeCandidate(3, 0.0, 0.0)
        kc3.build(self.kList)

        imagePca = ipDiffim.KernelPcaD()
        kpv = ipDiffim.KernelPcaVisitorF(imagePca)
        kpv.processCandidate(kc1)
        kpv.processCandidate(kc2)
        kpv.processCandidate(kc3)

        imagePca.analyze()
        eigenImages = imagePca.getEigenImages()
        # NOTE : this needs to be changed once ticket #1649 is resolved
        for i in range(len(eigenImages)):
            for j in range(i, len(eigenImages)):
                print(i, j, afwImage.innerProduct(eigenImages[i], eigenImages[j]))

    def testEigenValues(self):
        kc1 = self.makeCandidate(1, 0.0, 0.0)
        kc1.build(self.kList)

        kc2 = self.makeCandidate(2, 0.0, 0.0)
        kc2.build(self.kList)

        kc3 = self.makeCandidate(3, 0.0, 0.0)
        kc3.build(self.kList)

        imagePca = ipDiffim.KernelPcaD()
        kpv = ipDiffim.KernelPcaVisitorF(imagePca)
        kpv.processCandidate(kc1)
        kpv.processCandidate(kc2)
        kpv.processCandidate(kc3)

        imagePca.analyze()
        eigenImages = imagePca.getEigenImages()
        eigenValues = imagePca.getEigenValues()

        # took in 3 images
        self.assertEqual(len(eigenImages), 3)
        self.assertEqual(len(eigenValues), 3)

        # all the same shape, only 1 eigenvalue
        self.assertAlmostEqual(eigenValues[0], 1.0)
        self.assertAlmostEqual(eigenValues[1], 0.0)
        self.assertAlmostEqual(eigenValues[2], 0.0)

    def testMeanSubtraction(self):
        kc1 = self.makeCandidate(1, 0.0, 0.0)
        kc1.build(self.kList)

        kc2 = self.makeCandidate(2, 0.0, 0.0)
        kc2.build(self.kList)

        kc3 = self.makeCandidate(3, 0.0, 0.0)
        kc3.build(self.kList)

        imagePca = ipDiffim.KernelPcaD()
        kpv = ipDiffim.KernelPcaVisitorF(imagePca)
        kpv.processCandidate(kc1)
        kpv.processCandidate(kc2)
        kpv.processCandidate(kc3)
        kpv.subtractMean()  # subtract it *from* imagePca

        imagePca.analyze()
        eigenImages = imagePca.getEigenImages()
        eigenValues = imagePca.getEigenValues()

        # took in 3 images
        self.assertEqual(len(eigenImages), 3)
        self.assertEqual(len(eigenValues), 3)

        # all the same shape, mean subtracted, so *no* eigenvalues
        self.assertAlmostEqual(eigenValues[0], 0.0)
        self.assertAlmostEqual(eigenValues[1], 0.0)
        self.assertAlmostEqual(eigenValues[2], 0.0)

        # finally, since imagePca normalizes by the sum, this should
        # have central pixel value 1.0 and the rest 0.0
        imageMean = kpv.returnMean()
        rows = imageMean.getHeight()
        cols = imageMean.getWidth()
        for y in range(rows):
            for x in range(cols):
                if x == cols // 2 and y == rows // 2:
                    self.assertAlmostEqual(imageMean.get(x, y), 1.0)
                else:
                    self.assertAlmostEqual(imageMean.get(x, y), 0.0)

    def testVisit(self, nCell=3):
        imagePca = ipDiffim.KernelPcaD()
        kpv = ipDiffim.makeKernelPcaVisitor(imagePca)

        sizeCellX = self.policy.get("sizeCellX")
        sizeCellY = self.policy.get("sizeCellY")

        kernelCellSet = afwMath.SpatialCellSet(afwGeom.Box2I(afwGeom.Point2I(0,
                                                                             0),
                                                             afwGeom.Extent2I(sizeCellX * nCell,
                                                                              sizeCellY * nCell)),
                                               sizeCellX,
                                               sizeCellY)

        for candX in range(nCell):
            for candY in range(nCell):
                if candX == nCell // 2 and candY == nCell // 2:
                    kc = self.makeCandidate(100.0,
                                            candX * sizeCellX + sizeCellX // 2,
                                            candY * sizeCellY + sizeCellY // 2)
                else:
                    kc = self.makeCandidate(1.0,
                                            candX * sizeCellX + sizeCellX // 2,
                                            candY * sizeCellY + sizeCellY // 2)
                kc.build(self.kList)
                kernelCellSet.insertCandidate(kc)

        kernelCellSet.visitCandidates(kpv, 1)
        imagePca.analyze()
        eigenImages = imagePca.getEigenImages()
        eigenValues = imagePca.getEigenValues()

        # took in 3 images
        self.assertEqual(len(eigenImages), nCell * nCell)
        self.assertEqual(len(eigenValues), nCell * nCell)

        # all the same shape, only 1 eigenvalue
        self.assertAlmostEqual(eigenValues[0], 1.0)
        self.assertAlmostEqual(eigenValues[1], 0.0)
        self.assertAlmostEqual(eigenValues[2], 0.0)

#####


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
