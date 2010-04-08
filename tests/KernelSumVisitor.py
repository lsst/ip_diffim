#!/usr/bin/env python
import os, sys
import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
#import lsst.ip.diffim.detail as ipDiffimDetail
import lsst.pex.logging as pexLog

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

pexLog.Trace_setVerbosity('lsst.ip.diffim', 3)

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = ipDiffim.generateDefaultPolicy(diffimPolicy)
        self.policy.set("kernelBasisSet", "delta-function")
        self.policy.set("useRegularization", False)
        self.kList = ipDiffim.makeKernelBasisList(self.policy)

    def makeCandidate(self, kSum, x, y, size = 50):
        mi1 = afwImage.MaskedImageF(size, size)
        mi1.getVariance().set(0.1) # avoid NaNs
        mi1.set(size//2, size//2, (1, 0x0, 1))
        mi2 = afwImage.MaskedImageF(size, size)
        mi2.getVariance().set(0.1) # avoid NaNs
        mi2.set(size//2, size//2, (kSum, 0x0, 1))
        # currently works, but cannot put in spatialcell
        kc = ipDiffim.KernelCandidateF(x, y, mi1, mi2, self.policy)
        # currently fails
        #kc = ipDiffim.makeKernelCandidate(x, y, mi1, mi2, self.policy)
        return kc
    
    def tearDown(self):
        del self.policy

    def testAggregate(self, kSums = [1., 1., 1., 1., 2., 3., 4.]):
        ksv = ipDiffim.KernelSumVisitorF(self.policy)
        ksv.setMode(ipDiffim.KernelSumVisitorF.AGGREGATE)

        # should fail, kernel not initialized
        kc = self.makeCandidate(1, 0.0, 0.0) 
        try:
            ksv.processCandidate(kc)
        except Exception, e:
            pass
        else:
            self.fail()
            
        for kSum in kSums:
            kc = self.makeCandidate(kSum, 0., 0.)
            kc.build(self.kList)
            self.assertAlmostEqual(kSum, kc.getKsum(ipDiffim.KernelCandidateF.RECENT))
            ksv.processCandidate(kc)

        for method in (ksv.getNRejected,
                       ksv.getkSumMean,
                       ksv.getkSumStd,
                       ksv.getdkSumMax,
                       ksv.getkSumNpts):
            self.assertEqual(method(), 0.0)

        ksv.processKsumDistribution()

        self.assertEqual(ksv.getNRejected(), 0)
        self.assertAlmostEqual(ksv.getkSumMean(),
                               afwMath.makeStatistics(kSums, afwMath.MEANCLIP).getValue(afwMath.MEANCLIP))
        self.assertAlmostEqual(ksv.getkSumStd(),
                               afwMath.makeStatistics(kSums, afwMath.STDEVCLIP).getValue(afwMath.STDEVCLIP))
        self.assertEqual(ksv.getdkSumMax(),
                         self.policy.get("maxKsumSigma") * ksv.getkSumStd())
        self.assertEqual(ksv.getkSumNpts(), len(kSums))


    def testReject(self):
        self.doReject(clipping = False)
        self.doReject(clipping = True)

    def doReject(self, clipping, kSums = [1., 1., 1., 1., 2., 3., 4., 50.]):
        self.policy.set("kernelSumClipping", clipping)
        ksv = ipDiffim.KernelSumVisitorF(self.policy)
        ksv.setMode(ipDiffim.KernelSumVisitorF.AGGREGATE)
        kcList = []

        for kSum in kSums:
            kc = self.makeCandidate(kSum, 0., 0.)
            kc.build(self.kList)
            kc.setStatus(afwMath.SpatialCellCandidate.GOOD)
            self.assertAlmostEqual(kSum, kc.getKsum(ipDiffim.KernelCandidateF.RECENT))
            ksv.processCandidate(kc)
            kcList.append(kc)

        ksv.processKsumDistribution()
        
        ksv.setMode(ipDiffim.KernelSumVisitorF.REJECT)
        for kc in kcList:
            ksv.processCandidate(kc)
            if clipping and kc == kcList[-1]:
                self.assertEqual(kc.getStatus(), afwMath.SpatialCellCandidate.BAD)
            else:
                self.assertEqual(kc.getStatus(), afwMath.SpatialCellCandidate.GOOD)

        if clipping: 
            self.assertEqual(ksv.getNRejected(), 1)
        else:
            self.assertEqual(ksv.getNRejected(), 0)


    def testVisit(self, nCell = 3):
        # instead of manually visiting, one by one, visit using spatialcell
        ksv = ipDiffim.KernelSumVisitorF(self.policy)

        sizeCellX = self.policy.get("sizeCellX")
        sizeCellY = self.policy.get("sizeCellY")
        
        kernelCellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(0, 0)),
                                               sizeCellX * nCell,
                                               sizeCellY * nCell)
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
                print kc.getXCenter(), kc.getYCenter(), kc.getKsum(ipDiffim.KernelCandidateF.RECENT)

        # NEEDS TO GET FINISHED ONCE I CAN MAKECANDIDATE
        ksv.setMode(ipDiffim.KernelSumVisitorF.AGGREGATE)
        kernelCells.visitCandidates(ksv, 1)
        ksv.processKsumDistribution()
        ksv.setMode(ipDiffim.KernelSumVisitorF.REJECT)
        kernelCells.visitCandidates(ksv, 1)

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
