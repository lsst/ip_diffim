#!/usr/bin/env python
import os

import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.policy as pexPolicy
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as logging

import lsst.afw.display.ds9 as ds9
import numpy as num

Verbosity = 4
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = pexPolicy.Policy.createPolicy(diffimPolicy)
        
    def tearDown(self):
        del self.policy

    def testImageStatisticsZero(self):
        numArray = num.zeros((20,20))
        mi       = afwImage.MaskedImageF(20,20)
        for j in range(mi.getHeight()):
            for i in range(mi.getWidth()):
                mi.set( i, j, (numArray[j][i], 0x0, 1) )

        imstat = ipDiffim.ImageStatisticsF()
        imstat.apply(mi)

        self.assertEqual(imstat.getMean(), 0)
        self.assertEqual(imstat.getRms(), 0)
        self.assertEqual(imstat.getNpix(), 20*20)

    def testImageStatisticsOne(self):
        numArray = num.ones((20,20))
        mi       = afwImage.MaskedImageF(20,20)
        for j in range(mi.getHeight()):
            for i in range(mi.getWidth()):
                mi.set( i, j, (numArray[j][i], 0x0, 1) )

        imstat = ipDiffim.ImageStatisticsF()
        imstat.apply(mi)

        self.assertEqual(imstat.getMean(), 1)
        self.assertEqual(imstat.getRms(), 0)
        self.assertEqual(imstat.getNpix(), 20*20)

    def testImageStatisticsOne(self):
        numArray = num.ones((20,20))
        mi       = afwImage.MaskedImageF(20,20)
        for j in range(mi.getHeight()):
            for i in range(mi.getWidth()):
                mi.set( i, j, (numArray[j][i], 0x0, 1) )

        imstat = ipDiffim.ImageStatisticsF()
        imstat.apply(mi)

        self.assertEqual(imstat.getMean(), 1)
        self.assertEqual(imstat.getRms(), 0)
        self.assertEqual(imstat.getNpix(), 20*20)
        
    def testImageStatisticsGeneral(self):
        numArray = num.ones((20,20))
        mi       = afwImage.MaskedImageF(20,20)
        for j in range(mi.getHeight()):
            for i in range(mi.getWidth()):
                val = i + 2.3 * j
                mi.set( i, j, (val, 0x0, 1) )
                numArray[j][i] = val

        imstat = ipDiffim.ImageStatisticsF()
        imstat.apply(mi)

        self.assertAlmostEqual(imstat.getMean(), numArray.mean())
        # note that these don't agree exactly...
        self.assertAlmostEqual(imstat.getRms(), numArray.std(), 1)
        self.assertEqual(imstat.getNpix(), 20*20)

        afwStat = afwMath.makeStatistics(mi.getImage(), afwMath.MEAN | afwMath.STDEV)
        self.assertAlmostEqual(imstat.getMean(), afwStat.getValue(afwMath.MEAN))
        # even though these do
        self.assertAlmostEqual(imstat.getRms(), afwStat.getValue(afwMath.STDEV))

    def testImageStatisticsMask(self):
        numArray = num.ones((20,19))
        mi       = afwImage.MaskedImageF(20,20)
        for j in range(mi.getHeight()):
            for i in range(mi.getWidth()):
                val = i + 2.3 * j
                
                if i == 19:
                    mi.set( i, j, (val, 0x1, 1) )
                else:
                    mi.set( i, j, (val, 0x0, 1) )
                    numArray[j][i] = val

        imstat = ipDiffim.ImageStatisticsF()
        imstat.apply(mi)

        self.assertAlmostEqual(imstat.getMean(), numArray.mean())
        # note that these don't agree exactly...
        self.assertAlmostEqual(imstat.getRms(), numArray.std(), 1)
        self.assertEqual(imstat.getNpix(), 20*(20-1))


#####
        
def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(DiffimTestCases)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
