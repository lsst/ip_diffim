import unittest


import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.ip.diffim as ipDiffim
import lsst.utils.logging as logUtils
import lsst.pex.config as pexConfig

from lsst.ip.diffim import PsfMatchConfigAL

logUtils.trace_set_at("lsst.ip.diffim", 4)

# This tests the basics of the BuildSpatialKernelVisitor.  E.g. that
# it makes the right size solution.  For more complex behaviors such
# as reducing a delta function basis set into a Pca basis set, look at
# FitSpatialKernelFromCandidates.py


class DiffimTestCases(unittest.TestCase):

    def setUp(self):
        self.config = PsfMatchConfigAL()

        self.ps = pexConfig.makePropertySet(self.config)
        self.size = 51

    def tearDown(self):
        del self.ps

    def makeCandidate(self, kSum, x, y):
        mi1 = afwImage.MaskedImageF(geom.Extent2I(self.size, self.size))
        mi1.getVariance().set(1.0)  # avoid NaNs
        mi1[self.size//2, self.size//2, afwImage.LOCAL] = (1, 0x0, 1)
        mi2 = afwImage.MaskedImageF(geom.Extent2I(self.size, self.size))
        mi2.getVariance().set(1.0)  # avoid NaNs
        mi2[self.size//2, self.size//2, afwImage.LOCAL] = (kSum, 0x0, kSum)
        kc = ipDiffim.makeKernelCandidate(x, y, mi1, mi2, self.ps)
        return kc

    def testNoBg(self):
        self.runNoBg(0)
        self.runNoBg(1)
        self.runNoBg(2)

    def runNoBg(self, sko):
        basisList = ipDiffim.makeKernelBasisList(self.config)
        self.ps['spatialKernelOrder'] = sko
        self.ps['fitForBackground'] = False

        bbox = geom.Box2I(geom.Point2I(0, 0),
                          geom.Extent2I(self.size*10, self.size*10))

        bsikv = ipDiffim.BuildSingleKernelVisitorF(basisList, self.ps)
        bspkv = ipDiffim.BuildSpatialKernelVisitorF(basisList, bbox, self.ps)

        for x in range(1, self.size, 10):
            for y in range(1, self.size, 10):
                cand = self.makeCandidate(1.0, x, y)
                bsikv.processCandidate(cand)
                bspkv.processCandidate(cand)

        bspkv.solveLinearEquation()
        sk, sb = bspkv.getSolutionPair()

        # Kernel
        if sko == 0:
            # Specialization for speedup
            spatialKernelSolution = sk.getKernelParameters()

            # One term for each basis function
            self.assertEqual(len(spatialKernelSolution), len(basisList))

        else:
            spatialKernelSolution = sk.getSpatialParameters()

            nSpatialTerms = int(0.5 * (sko + 1) * (sko + 2))
            # One model for each basis function
            self.assertEqual(len(spatialKernelSolution), len(basisList))
            # First basis has no spatial variation
            for i in range(1, nSpatialTerms):
                self.assertEqual(spatialKernelSolution[0][i], 0.)
            # All bases have correct number of terms
            for i in range(len(spatialKernelSolution)):
                self.assertEqual(len(spatialKernelSolution[i]), nSpatialTerms)

        # Background
        spatialBgSolution = sb.getParameters()
        nBgTerms = 1
        self.assertEqual(len(spatialBgSolution), nBgTerms)

    def testModelType(self):
        bbox = geom.Box2I(geom.Point2I(10, 10),
                          geom.Extent2I(10, 10))
        basisList = ipDiffim.makeKernelBasisList(self.config)

        self.ps["spatialModelType"] = "polynomial"
        ipDiffim.BuildSpatialKernelVisitorF(basisList, bbox, self.ps)

        self.ps["spatialModelType"] = "chebyshev1"
        ipDiffim.BuildSpatialKernelVisitorF(basisList, bbox, self.ps)

        try:
            self.ps["spatialModelType"] = "foo"
            ipDiffim.BuildSpatialKernelVisitorF(basisList, bbox, self.ps)
        except Exception:
            pass
        else:
            self.fail()

    def testAlSpatialModel(self):
        self.runAlSpatialModel(0, 0)
        self.runAlSpatialModel(1, 0)
        self.runAlSpatialModel(0, 1)
        self.runAlSpatialModel(1, 1)
        self.runAlSpatialModel(2, 2)

    def runAlSpatialModel(self, sko, bgo):
        basisList = ipDiffim.makeKernelBasisList(self.config)
        self.ps['spatialKernelOrder'] = sko
        self.ps['spatialBgOrder'] = bgo
        self.ps['fitForBackground'] = True

        bbox = geom.Box2I(geom.Point2I(0, 0),
                          geom.Extent2I(self.size*10, self.size*10))

        bsikv = ipDiffim.BuildSingleKernelVisitorF(basisList, self.ps)
        bspkv = ipDiffim.BuildSpatialKernelVisitorF(basisList, bbox, self.ps)

        for x in range(1, self.size, 10):
            for y in range(1, self.size, 10):
                cand = self.makeCandidate(1.0, x, y)
                bsikv.processCandidate(cand)
                bspkv.processCandidate(cand)

        bspkv.solveLinearEquation()
        sk, sb = bspkv.getSolutionPair()

        # Kernel
        if sko == 0:
            # Specialization for speedup
            spatialKernelSolution = sk.getKernelParameters()

            # One term for each basis function
            self.assertEqual(len(spatialKernelSolution), len(basisList))

        else:
            spatialKernelSolution = sk.getSpatialParameters()

            nSpatialTerms = int(0.5 * (sko + 1) * (sko + 2))
            # One model for each basis function
            self.assertEqual(len(spatialKernelSolution), len(basisList))
            # First basis has no spatial variation
            for i in range(1, nSpatialTerms):
                self.assertEqual(spatialKernelSolution[0][i], 0.)
            # All bases have correct number of terms
            for i in range(len(spatialKernelSolution)):
                self.assertEqual(len(spatialKernelSolution[i]), nSpatialTerms)

        # Background
        spatialBgSolution = sb.getParameters()
        nBgTerms = int(0.5 * (bgo + 1) * (bgo + 2))
        self.assertEqual(len(spatialBgSolution), nBgTerms)


#####

class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
