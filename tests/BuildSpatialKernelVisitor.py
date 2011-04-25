#!/usr/bin/env python
import os, sys
import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.pex.logging as pexLog

import numpy
import pylab

pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)

# This tests the basics of the BuildSpatialKernelVisitor.  E.g. that
# it makes the right size solution.  For more complex behaviors such
# as reducing a delta function basis set into a Pca basis set, look at
# FitSpatialKernelFromCandidates.py

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = ipDiffim.makeDefaultPolicy()
        self.size   = 51

    def tearDown(self):
        del self.policy

    def makeCandidate(self, kSum, x, y):
        mi1 = afwImage.MaskedImageF(afwGeom.Extent2I(self.size, self.size))
        mi1.getVariance().set(0.1) # avoid NaNs
        mi1.set(self.size//2, self.size//2, (1, 0x0, 1))
        mi2 = afwImage.MaskedImageF(afwGeom.Extent2I(self.size, self.size))
        mi2.getVariance().set(0.1) # avoid NaNs
        mi2.set(self.size//2, self.size//2, (kSum, 0x0, 1))
        kc = ipDiffim.makeKernelCandidate(x, y, mi1, mi2, self.policy)
        return kc

    def testNoBg(self):
        self.policy.set('kernelBasisSet', 'alard-lupton')
        self.policy.set('spatialKernelOrder', 1)
        self.policy.set('fitForBackground', False)
        basisList = ipDiffim.makeKernelBasisList(self.policy)

        bsikv = ipDiffim.BuildSingleKernelVisitorF(basisList, self.policy)
        bspkv = ipDiffim.BuildSpatialKernelVisitorF(basisList, self.policy)
        
        for x in numpy.arange(1, self.size, 10):
            for y in numpy.arange(1, self.size, 10):
                cand = self.makeCandidate(1.0, x, y)
                bsikv.processCandidate(cand)
                bspkv.processCandidate(cand)
                
        bspkv.solveLinearEquation()
        sk, sb = bspkv.getSolutionPair()

    def testAlSpatialModel(self):
        self.runAlSpatialModel(0, 0)
        self.runAlSpatialModel(1, 0)
        self.runAlSpatialModel(0, 1)
        self.runAlSpatialModel(1, 1)
        self.runAlSpatialModel(2, 2)
        
    def runAlSpatialModel(self, sko, bgo):
        self.policy.set('kernelBasisSet', 'alard-lupton')
        self.policy.set('spatialKernelOrder', sko)
        self.policy.set('spatialBgOrder', bgo)
        basisList = ipDiffim.makeKernelBasisList(self.policy)

        bsikv = ipDiffim.BuildSingleKernelVisitorF(basisList, self.policy)
        bspkv = ipDiffim.BuildSpatialKernelVisitorF(basisList, self.policy)
        
        for x in numpy.arange(1, self.size, 10):
            for y in numpy.arange(1, self.size, 10):
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
