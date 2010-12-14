#!/usr/bin/env python
import os, sys
import numpy
import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as pexLog

pexLog.Trace_setVerbosity('lsst.ip.diffim', 7)

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = ipDiffim.createDefaultPolicy()
        self.policy.set("checkConditionNumber", False) # these images have been hand-constructed
        self.size   = 30
        self.policy.set("sizeCellX", self.size//3)
        self.policy.set("sizeCellY", self.size//3)
        self.kernelCellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(0,0),
                                                                  self.size,
                                                                  self.size),
                                                    self.policy.getInt("sizeCellX"),
                                                    self.policy.getInt("sizeCellY"))

    def tearDown(self):
        del self.policy
        del self.kernelCellSet
        
    def addNoise(self, mi):
        img       = mi.getImage()
        seed      = int(afwMath.makeStatistics(mi.getVariance(), afwMath.MAX).getValue())
        rdm       = afwMath.Random(afwMath.Random.MT19937, seed)
        rdmImage  = img.Factory(img.getDimensions())
        afwMath.randomGaussianImage(rdmImage, rdm)
        img      += rdmImage

    def makeCandidate(self, kSum, x, y, addNoise = True):
        mi1 = afwImage.MaskedImageF(self.size, self.size)
        mi1.getVariance().set(0.1) # avoid NaNs
        mi1.set(self.size//2, self.size//2, (1, 0x0, 1))
        mi2 = afwImage.MaskedImageF(self.size, self.size)
        mi2.getVariance().set(0.1) # avoid NaNs
        mi2.set(self.size//2, self.size//2, (kSum, 0x0, 1))
        if addNoise:
            self.addNoise(mi1)
            self.addNoise(mi2)
        kc = ipDiffim.makeKernelCandidate(x, y, mi1, mi2, self.policy)
        return kc

    def testAlardLupton(self):
        #self.runAlardLupton(0, 0)
        #self.runAlardLupton(1, 1)
        #self.runAlardLuptonPca(False)
        self.runAlardLuptonPca(True)

    def runAlardLupton(self, sko, bgo):
        self.policy.set('kernelBasisSet', 'alard-lupton')
        self.policy.set('spatialKernelOrder', sko)
        self.policy.set('spatialBgOrder', bgo)
        self.policy.set('usePcaForSpatialKernel', False)
        basisList = ipDiffim.makeKernelBasisList(self.policy)
        
        for x in numpy.arange(1, self.size, 10):
            for y in numpy.arange(1, self.size, 10):
                cand = self.makeCandidate(1.0, x, y)
                self.kernelCellSet.insertCandidate(cand)

        result = ipDiffim.fitSpatialKernelFromCandidates(self.kernelCellSet, self.policy)
        sk = result.first
        sb = result.second

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
        for x in numpy.arange(1, self.size, 10):
            for y in numpy.arange(1, self.size, 10):
                cand = self.makeCandidate(1.0, x, y)
                self.kernelCellSet.insertCandidate(cand)
                count += 1
        #if subtractMean:
        #    # all the components are the same!  just keep mean and first component
        #    self.policy.set('numPrincipalComponents', 2)
        #else:
        #    self.policy.set('numPrincipalComponents', count)
            
        self.policy.set('numPrincipalComponents', count)

        result = ipDiffim.fitSpatialKernelFromCandidates(self.kernelCellSet, self.policy)
        sk = result.first
        sb = result.second

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

    def xtestDf(self):
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
        
        self.policy.set('subtractMeanForPca', subtractMean)
        self.policy.set('useRegularization', useRegularization)

        # We need to add noise to create some variation of the kernels for Pca
        kSum = 1.
        for x in numpy.arange(1, self.size, 10):
            for y in numpy.arange(1, self.size, 10):
                cand = self.makeCandidate(kSum, x, y, True)
                self.kernelCellSet.insertCandidate(cand)
                kSum += 1

        result = ipDiffim.fitSpatialKernelFromCandidates(self.kernelCellSet, self.policy)
        sk = result.first
        sb = result.second

        spatialKernelSolution = sk.getSpatialParameters()

        nPca = self.policy.get('numPrincipalComponents')
        if subtractMean:
            nPca += 1
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
