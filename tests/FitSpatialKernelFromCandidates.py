#!/usr/bin/env python
import os, sys
import numpy
import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as pexLog

pexLog.Trace_setVerbosity('lsst.ip.diffim', 9)
rdm = afwMath.Random(afwMath.Random.MT19937, 10101)

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = ipDiffim.makeDefaultPolicy()
        self.policy.set("checkConditionNumber", False) # these images have been hand-constructed
        self.policy.set("fitForBackground", True)      # with background testing
        self.size   = 50
        self.policy.set("sizeCellX", 1)
        self.policy.set("sizeCellY", 1)

        # Don't let the random noise we add cause any problems
        self.policy.set("singleKernelClipping", False)        
        self.policy.set("kernelSumClipping", False)        
        self.policy.set("spatialKernelClipping", False)        
        
        self.kernelCellSet = afwMath.SpatialCellSet(afwGeom.Box2I(afwGeom.Point2I(0,0),
                                                                  afwGeom.Extent2I(self.size, self.size)),
                                                    self.policy.getInt("sizeCellX"),
                                                    self.policy.getInt("sizeCellY"))

    def tearDown(self):
        del self.policy
        del self.kernelCellSet
        
    def addNoise(self, mi):
        img       = mi.getImage()
        rdmImage  = img.Factory(img.getDimensions())
        afwMath.randomGaussianImage(rdmImage, rdm)
        img      += rdmImage

    def makeCandidate(self, kSum, x, y, addNoise = True):
        mi1 = afwImage.MaskedImageF(afwGeom.Extent2I(self.size, self.size))
        mi1.getVariance().set(1.0) # level of addNoise
        mi1.set(self.size//2, self.size//2, (5.0, 0x0, 5.0))
        mi2 = afwImage.MaskedImageF(afwGeom.Extent2I(self.size, self.size))
        mi2.getVariance().set(1.0) # level of addNoise
        mi2.set(self.size//2, self.size//2, (kSum, 0x0, kSum))
        if addNoise:
            self.addNoise(mi1)
            self.addNoise(mi2)
        kc = ipDiffim.makeKernelCandidate(x, y, mi1, mi2, self.policy)
        return kc

    def testAlardLupton(self):
        self.runAlardLupton(0, 0)
        self.runAlardLupton(1, 1)
        self.runAlardLuptonPca(False)
        self.runAlardLuptonPca(True)

    def runAlardLupton(self, sko, bgo):
        self.policy.set('kernelBasisSet', 'alard-lupton')
        self.policy.set('spatialKernelOrder', sko)
        self.policy.set('spatialBgOrder', bgo)
        self.policy.set('usePcaForSpatialKernel', False)
        basisList = ipDiffim.makeKernelBasisList(self.policy)
        
        for x in range(1, self.size, 10):
            for y in range(1, self.size, 10):
                cand = self.makeCandidate(10.0, x, y)
                self.kernelCellSet.insertCandidate(cand)

        sk, sb = ipDiffim.fitSpatialKernelFromCandidates(self.kernelCellSet, self.policy)

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

    def runAlardLuptonPca(self, subtractMean, sko = 1, bgo = 1):
        # Just test the Pca part
        self.policy.set('kernelBasisSet', 'alard-lupton')
        self.policy.set('spatialKernelOrder', sko)
        self.policy.set('spatialBgOrder', bgo)
        self.policy.set('usePcaForSpatialKernel', True)
        self.policy.set('subtractMeanForPca', subtractMean)

        count = 0
        for x in range(1, self.size, 10):
            for y in range(1, self.size, 10):
                cand = self.makeCandidate(10.0, x, y)
                self.kernelCellSet.insertCandidate(cand)
                count += 1
        if subtractMean:
            # all the components are the same!  just keep mean and first component
            self.policy.set('numPrincipalComponents', 2)
        else:
            self.policy.set('numPrincipalComponents', count)
            
        sk, sb = ipDiffim.fitSpatialKernelFromCandidates(self.kernelCellSet, self.policy)

        spatialKernelSolution = sk.getSpatialParameters()

        nPca = self.policy.get('numPrincipalComponents')
        self.assertEqual(len(spatialKernelSolution), nPca)
        
        nSpatialTerms = int(0.5 * (sko + 1) * (sko + 2))
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

    def testDfTooSmall(self):
        # Test too small of stamp size for kernel size
        self.size = self.policy.get('kernelSize')
        try:
            self.runDfPca(False, False)
        except Exception, e:
            pass
        else:
            self.fail()
        
    def testDf(self):
        self.runDfPca(False, False)
        self.runDfPca(False, True)
        self.runDfPca(True, False)
        self.runDfPca(True, True)


    def runDfPca(self, useRegularization, subtractMean, sko = 1, bgo = 1):
        # Just test the Pca part
        self.policy.set('kernelBasisSet', 'delta-function')
        self.policy.set('spatialKernelOrder', sko)
        self.policy.set('spatialBgOrder', bgo)
        self.policy.set('usePcaForSpatialKernel', True)
        basisList = ipDiffim.makeKernelBasisList(self.policy)
        
        # Switched on and off
        self.policy.set('subtractMeanForPca', subtractMean)
        self.policy.set('useRegularization', useRegularization)

        # We need to add noise to create some variation of the kernels for Pca
        count = 0
        for x in range(1, self.size, 10):
            for y in range(1, self.size, 10):
                cand = self.makeCandidate(10.0, x, y)
                self.kernelCellSet.insertCandidate(cand)
                count += 1

        if subtractMean:
            # all the components are the same!  just keep mean and first component
            self.policy.set('numPrincipalComponents', 2)
        else:
            self.policy.set('numPrincipalComponents', count)

        sk, sb = ipDiffim.fitSpatialKernelFromCandidates(self.kernelCellSet, self.policy)

        spatialKernelSolution = sk.getSpatialParameters()

        nPca = self.policy.get('numPrincipalComponents')
        self.assertEqual(len(spatialKernelSolution), nPca)
        
        nSpatialTerms = int(0.5 * (sko + 1) * (sko + 2))
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
