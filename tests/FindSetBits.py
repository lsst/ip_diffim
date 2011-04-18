#!/usr/bin/env python
import os

import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as logging

verbosity = 1
logging.Trace_setVerbosity('lsst.ip.diffim', verbosity)

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = ipDiffim.createDefaultPolicy()
        
    def tearDown(self):
        del self.policy

    def testNoMask(self):
        mask = afwImage.MaskU(afwGeom.Extent2I(20, 20))
        mask.set(0)
        fsb  = ipDiffim.FindSetBitsU()

        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 10),
                             afwGeom.Point2I(9, 12))
        fsb.apply(afwImage.MaskU(mask, bbox, afwImage.LOCAL))

        self.assertEqual(fsb.getBits(), 0)

    def testOneMask(self):
        mask = afwImage.MaskU(afwGeom.Extent2I(20, 20))
        mask.set(0)
        bitmaskBad = mask.getPlaneBitMask('BAD')
        fsb = ipDiffim.FindSetBitsU()

        bbox     = afwGeom.Box2I(afwGeom.Point2I(9, 10),
                                 afwGeom.Point2I(11, 12))
        submask  = afwImage.MaskU(mask, bbox, afwImage.LOCAL)
        submask |= bitmaskBad

        bbox2    = afwGeom.Box2I(afwGeom.Point2I(8, 8),
                                 afwGeom.Point2I(19, 19))
        fsb.apply(afwImage.MaskU(mask, bbox2, afwImage.LOCAL))

        self.assertEqual(fsb.getBits(), bitmaskBad)

    def testManyMask(self):
        mask = afwImage.MaskU(afwGeom.Extent2I(20, 20))
        mask.set(0)
        bitmaskBad = mask.getPlaneBitMask('BAD')
        bitmaskSat = mask.getPlaneBitMask('SAT')
        fsb = ipDiffim.FindSetBitsU()

        bbox      = afwGeom.Box2I(afwGeom.Point2I(9, 10),
                                  afwGeom.Point2I(11, 12))
        submask   = afwImage.MaskU(mask, bbox, afwImage.LOCAL)
        submask  |= bitmaskBad

        bbox2     = afwGeom.Box2I(afwGeom.Point2I(8, 8),
                                  afwGeom.Point2I(19, 19))
        submask2  = afwImage.MaskU(mask, bbox2, afwImage.LOCAL)
        submask2 |= bitmaskSat

        bbox3     = afwGeom.Box2I(afwGeom.Point2I(0, 0),
                                  afwGeom.Point2I(19, 19))
        fsb.apply(afwImage.MaskU(mask, bbox3, afwImage.LOCAL))

        self.assertEqual(fsb.getBits(), bitmaskBad | bitmaskSat)

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
