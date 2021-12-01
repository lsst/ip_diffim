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
import lsst.afw.math as afwMath
import lsst.geom as geom
import lsst.ip.diffim as ipDiffim
import lsst.utils.logging as logUtils
import lsst.pex.config as pexConfig
import numpy as num

verbosity = 0
logUtils.trace_set_at("lsst.ip.diffim", verbosity)


class DiffimTestCases(unittest.TestCase):

    def setUp(self):
        self.config = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.subconfig = self.config.kernel["DF"]
        self.ps = pexConfig.makePropertySet(self.subconfig)

    def tearDown(self):
        del self.ps

    def testImageStatisticsNan(self, core=3):
        numArray = num.zeros((20, 20))
        mi = afwImage.MaskedImageF(geom.Extent2I(20, 20))
        for j in range(mi.getHeight()):
            for i in range(mi.getWidth()):
                mi[i, j, afwImage.LOCAL] = (numArray[j][i], 0x0, 0)

        # inverse variance weight of 0 is NaN
        imstat = ipDiffim.ImageStatisticsF(self.ps)
        imstat.apply(mi)
        self.assertEqual(imstat.getNpix(), 0)

        imstat = ipDiffim.ImageStatisticsF(self.ps)
        imstat.apply(mi, core)
        self.assertEqual(imstat.getNpix(), 0)

    def testImageStatisticsZero(self):
        numArray = num.zeros((20, 20))
        mi = afwImage.MaskedImageF(geom.Extent2I(20, 20))
        for j in range(mi.getHeight()):
            for i in range(mi.getWidth()):
                mi[i, j, afwImage.LOCAL] = (numArray[j][i], 0x0, 1)

        imstat = ipDiffim.ImageStatisticsF(self.ps)
        imstat.apply(mi)

        self.assertEqual(imstat.getMean(), 0)
        self.assertEqual(imstat.getRms(), 0)
        self.assertEqual(imstat.getNpix(), 20*20)

    def testImageStatisticsOne(self):
        numArray = num.ones((20, 20))
        mi = afwImage.MaskedImageF(geom.Extent2I(20, 20))
        for j in range(mi.getHeight()):
            for i in range(mi.getWidth()):
                mi[i, j, afwImage.LOCAL] = (numArray[j][i], 0x0, 1)

        imstat = ipDiffim.ImageStatisticsF(self.ps)
        imstat.apply(mi)

        self.assertEqual(imstat.getMean(), 1)
        self.assertEqual(imstat.getRms(), 0)
        self.assertEqual(imstat.getNpix(), 20*20)

    def testImageStatisticsCore(self, core=3):
        numArray = num.ones((20, 20))
        mi = afwImage.MaskedImageF(geom.Extent2I(20, 20))
        for j in range(mi.getHeight()):
            for i in range(mi.getWidth()):
                mi[i, j, afwImage.LOCAL] = (numArray[j][i], 0x0, 1)

        imstat = ipDiffim.ImageStatisticsF(self.ps)
        imstat.apply(mi, core)

        self.assertEqual(imstat.getMean(), 1)
        self.assertEqual(imstat.getRms(), 0)
        self.assertEqual(imstat.getNpix(), (2*core+1)**2)

    def testImageStatisticsGeneral(self):
        numArray = num.ones((20, 20))
        mi = afwImage.MaskedImageF(geom.Extent2I(20, 20))
        for j in range(mi.getHeight()):
            for i in range(mi.getWidth()):
                val = i + 2.3 * j
                mi[i, j, afwImage.LOCAL] = (val, 0x0, 1)
                numArray[j][i] = val

        imstat = ipDiffim.ImageStatisticsF(self.ps)
        imstat.apply(mi)

        self.assertAlmostEqual(imstat.getMean(), numArray.mean())
        # note that these don't agree exactly...
        self.assertAlmostEqual(imstat.getRms(), numArray.std(), 1)
        self.assertEqual(imstat.getNpix(), 20 * 20)

        afwStat = afwMath.makeStatistics(mi.getImage(), afwMath.MEAN | afwMath.STDEV)
        self.assertAlmostEqual(imstat.getMean(), afwStat.getValue(afwMath.MEAN))
        # even though these do
        self.assertAlmostEqual(imstat.getRms(), afwStat.getValue(afwMath.STDEV))

    def testImageStatisticsMask1(self):
        # Mask value that gets ignored
        maskPlane = self.ps.getArray("badMaskPlanes")[0]
        maskVal = afwImage.Mask.getPlaneBitMask(maskPlane)
        numArray = num.ones((20, 19))
        mi = afwImage.MaskedImageF(geom.Extent2I(20, 20))
        for j in range(mi.getHeight()):
            for i in range(mi.getWidth()):
                val = i + 2.3 * j

                if i == 19:
                    mi[i, j, afwImage.LOCAL] = (val, maskVal, 1)
                else:
                    mi[i, j, afwImage.LOCAL] = (val, 0x0, 1)
                    numArray[j][i] = val

        imstat = ipDiffim.ImageStatisticsF(self.ps)
        imstat.apply(mi)

        self.assertAlmostEqual(imstat.getMean(), numArray.mean())
        # note that these don't agree exactly...
        self.assertAlmostEqual(imstat.getRms(), numArray.std(), 1)
        self.assertEqual(imstat.getNpix(), 20 * (20 - 1))

    def testImageStatisticsMask2(self):
        # Mask value that does not get ignored
        maskPlanes = self.ps.getArray("badMaskPlanes")
        for maskPlane in ("BAD", "EDGE", "CR", "SAT", "INTRP"):
            if maskPlane not in maskPlanes:
                maskVal = afwImage.Mask.getPlaneBitMask(maskPlane)
                break
        self.assertGreater(maskVal, 0)

        numArray = num.ones((20, 20))
        mi = afwImage.MaskedImageF(geom.Extent2I(20, 20))
        for j in range(mi.getHeight()):
            for i in range(mi.getWidth()):
                val = i + 2.3 * j

                if i == 19:
                    mi[i, j, afwImage.LOCAL] = (val, maskVal, 1)
                    numArray[j][i] = val
                else:
                    mi[i, j, afwImage.LOCAL] = (val, 0x0, 1)
                    numArray[j][i] = val

        imstat = ipDiffim.ImageStatisticsF(self.ps)
        imstat.apply(mi)

        self.assertAlmostEqual(imstat.getMean(), numArray.mean())
        # note that these don't agree exactly...
        self.assertAlmostEqual(imstat.getRms(), numArray.std(), 1)
        self.assertEqual(imstat.getNpix(), 20 * 20)


#####

class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
