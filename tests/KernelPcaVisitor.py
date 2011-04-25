#!/usr/bin/env python
import os, sys
import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as pexLog

pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = ipDiffim.makeDefaultPolicy()
        self.policy.set("kernelBasisSet", "delta-function")
        self.policy.set("useRegularization", False)
        self.kList = ipDiffim.makeKernelBasisList(self.policy)

    def tearDown(self):
        del self.policy
        del self.kList

    def makeCandidate(self, kSum, x, y, size = 51):
        mi1 = afwImage.MaskedImageF(afwGeom.Extent2I(size, size))
        mi1.getVariance().set(0.1) # avoid NaNs
        mi1.set(size//2, size//2, (1, 0x0, 1))
        mi2 = afwImage.MaskedImageF(afwGeom.Extent2I(size, size))
        mi2.getVariance().set(0.1) # avoid NaNs
        mi2.set(size//2, size//2, (kSum, 0x0, 1))
        kc = ipDiffim.makeKernelCandidate(x, y, mi1, mi2, self.policy)
        return kc

    def testAlardLupton(self):
        self.policy.set("kernelBasisSet", "alard-lupton")
        self.kList = ipDiffim.makeKernelBasisList(self.policy)
        nTerms = len(self.kList)
        
        

    def testEigenValues(self):
        kc1 = self.makeCandidate(1, 0.0, 0.0)
        kc1.build(self.kList)

        kc2 = self.makeCandidate(2, 0.0, 0.0)
        kc2.build(self.kList)

        kc3 = self.makeCandidate(3, 0.0, 0.0)
        kc3.build(self.kList)

        imagePca = afwImage.ImagePcaD()
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

        imagePca = afwImage.ImagePcaD()
        kpv = ipDiffim.KernelPcaVisitorF(imagePca)
        kpv.processCandidate(kc1)
        kpv.processCandidate(kc2)
        kpv.processCandidate(kc3)
        kpv.subtractMean() # subtract it *from* imagePca

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

    def xtestVisit(self, nCell = 3):
        # This currently fails since I can't get visitCandidates to
        # tell this is a pointer
        imagePca = afwImage.ImagePcaD()
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
        self.assertEqual(len(eigenImages), 3)
        self.assertEqual(len(eigenValues), 3)

        # all the same shape, only 1 eigenvalue
        self.assertAlmostEqual(eigenValues[0], 1.0)
        self.assertAlmostEqual(eigenValues[1], 0.0)
        self.assertAlmostEqual(eigenValues[2], 0.0)
        
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
