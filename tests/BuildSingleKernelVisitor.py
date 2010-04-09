#!/usr/bin/env python
import os, sys
import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as pexLog

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = ipDiffim.generateDefaultPolicy(diffimPolicy)
        self.policy.set("kernelBasisSet", "delta-function")
        self.policy.set("useRegularization", False)
        self.kList = ipDiffim.makeKernelBasisList(self.policy)

    def makeCandidate(self, kSum, x, y, size = 51):
        mi1 = afwImage.MaskedImageF(size, size)
        mi1.getVariance().set(0.1) # avoid NaNs
        mi1.set(size//2, size//2, (1, 0x0, 1))
        mi2 = afwImage.MaskedImageF(size, size)
        mi2.getVariance().set(0.1) # avoid NaNs
        mi2.set(size//2, size//2, (kSum, 0x0, 1))
        kc = ipDiffim.makeKernelCandidate(x, y, mi1, mi2, self.policy)
        return kc


    def testWithOneBasis(self):
        kc1 = self.makeCandidate(1, 0.0, 0.0)
        kc2 = self.makeCandidate(2, 0.0, 0.0)
        kc3 = self.makeCandidate(3, 0.0, 0.0)

        #hMat = ipDiffim.makeRegularizationMatrix(self.policy)
        #print type(hMat)
        #bskv = ipDiffim.BuildSingleKernelVisitorF(self.kList, self.policy, hMat)
        print type(self.kList), type(self.kList[0])
        bskv = ipDiffim.BuildSingleKernelVisitorF(self.kList, self.policy)

        import pdb
        pdb.set_trace()
        bskv.processCandidate(kc1)
        bskv.processCandidate(kc2)
        bskv.processCandidate(kc3)

        # Initialized
        self.assertEqual(kc1.isInitialized(), True)
        self.assertEqual(kc2.isInitialized(), True)
        self.assertEqual(kc3.isInitialized(), True)

        # Is a solution
        try:
            kc1.getKernelSolution(ipDiffim.KernelCandidateF.RECENT)
            kc2.getKernelSolution(ipDiffim.KernelCandidateF.RECENT)
            kc3.getKernelSolution(ipDiffim.KernelCandidateF.RECENT)
        except Exception, e:
            print e
            self.fail()

        # Its not the Pca one
        try:
            kc1.getKernelSolution(ipDiffim.KernelCandidateF.PCA)
            kc2.getKernelSolution(ipDiffim.KernelCandidateF.PCA)
            kc3.getKernelSolution(ipDiffim.KernelCandidateF.PCA)
        except Exception, e:
            pass
        else:
            self.fail()

        # Processed all of them
        self.assertEqual(bskv.getNProcessed(), 3)

        # Rejected none
        self.assertEqual(bskv.getNRejected(), 0)

        # Skips built candidates
        bskv.reset()
        bskv.setSkipBuilt(True)
        bskv.processCandidate(kc1)
        bskv.processCandidate(kc2)
        bskv.processCandidate(kc3)
        # Processed none of them
        self.assertEqual(bskv.getNProcessed(), 0)
        
        
    def testWithThreeBases(self):
        kc1 = self.makeCandidate(1, 0.0, 0.0)
        kc2 = self.makeCandidate(2, 0.0, 0.0)
        kc3 = self.makeCandidate(3, 0.0, 0.0)
        bskv1 = ipDiffim.BuildSingleKernelVisitorF(self.kList, self.policy)
        bskv1.processCandidate(kc1)
        bskv1.processCandidate(kc2)
        bskv1.processCandidate(kc3)
        return

        # do pca basis; visit manually since visitCandidates is still broken
        imagePca = afwImage.ImagePcaD()
        kpv = ipDiffim.KernelPcaVisitorF(imagePca)
        kpv.processCandidate(kc1)
        kpv.processCandidate(kc2)
        kpv.processCandidate(kc3)
        imagePca.analyze()
        eigenImages = imagePca.getEigenImages()

        # do twice to mimic a Pca loop
        bskv2 = ipDiffim.BuildSingleKernelVisitorF(eigenImages, self.policy)
        bskv2.processCandidate(kc1)
        bskv2.processCandidate(kc2)
        bskv2.processCandidate(kc3)

        # do twice to mimic a Pca loop
        bskv3 = ipDiffim.BuildSingleKernelVisitorF(eigenImages, self.policy)
        bskv3.processCandidate(kc1)
        bskv3.processCandidate(kc2)
        bskv3.processCandidate(kc3)
        
        

#    def testRejection
        

#        bskv.setCandidateKernel(True)
        
#        bskv.setCandidateKernel(False)

#    def testRegularization(self):
        

#    def testVisit(self):
#        # test that visitCandidates calls reset() first
        
    def tearDown(self):
        del self.policy

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
