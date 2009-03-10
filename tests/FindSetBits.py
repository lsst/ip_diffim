#!/usr/bin/env python
import os

import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.pex.policy as pexPolicy
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as logging

import lsst.afw.display.ds9 as ds9

Verbosity = 4
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = pexPolicy.Policy.createPolicy(diffimPolicy)
        
    def tearDown(self):
        del self.policy

    def testNoMask(self):
        mask = afwImage.MaskU(20,20)
        mask.set(0)
        fsb  = ipDiffim.FindSetBitsU(mask)

        bbox     = afwImage.BBox(afwImage.PointI(0,10),
                                 afwImage.PointI(9,12))
        fp       = afwDetection.Footprint(bbox)
        fsb.apply(fp)

        self.assertEqual(fsb.getBits(), 0)

    def testOneMask(self):
        mask = afwImage.MaskU(20,20)
        mask.set(0)
        bitmaskBad = mask.getPlaneBitMask('BAD')
        fsb = ipDiffim.FindSetBitsU(mask)

        bbox     = afwImage.BBox(afwImage.PointI(9,10),
                                 afwImage.PointI(11,12))
        submask  = afwImage.MaskU(mask, bbox)
        submask |= bitmaskBad

        bbox2    = afwImage.BBox(afwImage.PointI(8,8),
                                 afwImage.PointI(19,19))
        fp       = afwDetection.Footprint(bbox2)
        fsb.apply(fp)

        self.assertEqual(fsb.getBits(), bitmaskBad)

    def testManyMask(self):
        mask = afwImage.MaskU(20,20)
        mask.set(0)
        bitmaskBad = mask.getPlaneBitMask('BAD')
        bitmaskSat = mask.getPlaneBitMask('SAT')
        fsb = ipDiffim.FindSetBitsU(mask)

        bbox      = afwImage.BBox(afwImage.PointI(9,10),
                                  afwImage.PointI(11,12))
        submask   = afwImage.MaskU(mask, bbox)
        submask  |= bitmaskBad

        bbox2     = afwImage.BBox(afwImage.PointI(8,8),
                                  afwImage.PointI(19,19))
        submask2  = afwImage.MaskU(mask, bbox2)
        submask2 |= bitmaskSat

        bbox3     = afwImage.BBox(afwImage.PointI(0,0),
                                  afwImage.PointI(19,19))
        fp       = afwDetection.Footprint(bbox3)
        fsb.apply(fp)

        self.assertEqual(fsb.getBits(), bitmaskBad | bitmaskSat)

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
