#
# LSST Data Management System
# Copyright 2008-2016 LSST Corporation.
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
import unittest

import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.ip.diffim as ipDiffim
import lsst.log.utils as logUtils

verbosity = 0
logUtils.traceSetAt("lsst.ip.diffim", verbosity)


class DiffimTestCases(unittest.TestCase):

    def testNoMask(self):
        mask = afwImage.Mask(geom.Extent2I(20, 20))
        mask.set(0)
        fsb = ipDiffim.FindSetBitsU()

        bbox = geom.Box2I(geom.Point2I(0, 10),
                          geom.Point2I(9, 12))
        fsb.apply(afwImage.Mask(mask, bbox, afwImage.LOCAL))

        self.assertEqual(fsb.getBits(), 0)

    def testOneMask(self):
        mask = afwImage.Mask(geom.Extent2I(20, 20))
        mask.set(0)
        bitmaskBad = mask.getPlaneBitMask('BAD')
        fsb = ipDiffim.FindSetBitsU()

        bbox = geom.Box2I(geom.Point2I(9, 10),
                          geom.Point2I(11, 12))
        submask = afwImage.Mask(mask, bbox, afwImage.LOCAL)
        submask |= bitmaskBad

        bbox2 = geom.Box2I(geom.Point2I(8, 8),
                           geom.Point2I(19, 19))
        fsb.apply(afwImage.Mask(mask, bbox2, afwImage.LOCAL))

        self.assertEqual(fsb.getBits(), bitmaskBad)

    def testManyMask(self):
        mask = afwImage.Mask(geom.Extent2I(20, 20))
        mask.set(0)
        bitmaskBad = mask.getPlaneBitMask('BAD')
        bitmaskSat = mask.getPlaneBitMask('SAT')
        fsb = ipDiffim.FindSetBitsU()

        bbox = geom.Box2I(geom.Point2I(9, 10),
                          geom.Point2I(11, 12))
        submask = afwImage.Mask(mask, bbox, afwImage.LOCAL)
        submask |= bitmaskBad

        bbox2 = geom.Box2I(geom.Point2I(8, 8),
                           geom.Point2I(19, 19))
        submask2 = afwImage.Mask(mask, bbox2, afwImage.LOCAL)
        submask2 |= bitmaskSat

        bbox3 = geom.Box2I(geom.Point2I(0, 0),
                           geom.Point2I(19, 19))
        fsb.apply(afwImage.Mask(mask, bbox3, afwImage.LOCAL))

        self.assertEqual(fsb.getBits(), bitmaskBad | bitmaskSat)

#####


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
