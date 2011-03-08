#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import os

import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as logging

verbosity = 1
logging.Trace_setVerbosity('lsst.ip.diffim', verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'policy', 'ImageSubtractStageDictionary.paf')

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = ipDiffim.generateDefaultPolicy(diffimPolicy)
        
    def tearDown(self):
        del self.policy

    def testNoMask(self):
        mask = afwImage.MaskU(20, 20)
        mask.set(0)
        fsb  = ipDiffim.FindSetBitsU()

        bbox     = afwImage.BBox(afwImage.PointI(0, 10),
                                 afwImage.PointI(9, 12))
        fsb.apply(afwImage.MaskU(mask, bbox))

        self.assertEqual(fsb.getBits(), 0)

    def testOneMask(self):
        mask = afwImage.MaskU(20, 20)
        mask.set(0)
        bitmaskBad = mask.getPlaneBitMask('BAD')
        fsb = ipDiffim.FindSetBitsU()

        bbox     = afwImage.BBox(afwImage.PointI(9, 10),
                                 afwImage.PointI(11, 12))
        submask  = afwImage.MaskU(mask, bbox)
        submask |= bitmaskBad

        bbox2    = afwImage.BBox(afwImage.PointI(8, 8),
                                 afwImage.PointI(19, 19))
        fsb.apply(afwImage.MaskU(mask, bbox2))

        self.assertEqual(fsb.getBits(), bitmaskBad)

    def testManyMask(self):
        mask = afwImage.MaskU(20, 20)
        mask.set(0)
        bitmaskBad = mask.getPlaneBitMask('BAD')
        bitmaskSat = mask.getPlaneBitMask('SAT')
        fsb = ipDiffim.FindSetBitsU()

        bbox      = afwImage.BBox(afwImage.PointI(9, 10),
                                  afwImage.PointI(11, 12))
        submask   = afwImage.MaskU(mask, bbox) 
        submask  |= bitmaskBad

        bbox2     = afwImage.BBox(afwImage.PointI(8, 8),
                                  afwImage.PointI(19, 19))
        submask2  = afwImage.MaskU(mask, bbox2)
        submask2 |= bitmaskSat

        bbox3     = afwImage.BBox(afwImage.PointI(0, 0),
                                  afwImage.PointI(19, 19))
        fsb.apply(afwImage.MaskU(mask, bbox3))

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
