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
import lsst.ip.diffim as ipDiffim


class DiffimTestCases(lsst.utils.tests.TestCase):

    def testDiaSourceAnalystConfig(self):
        config = ipDiffim.DiaSourceAnalystConfig()
        config.validate()

    def testImagePsfMatchConfig(self):
        config = ipDiffim.ImagePsfMatchConfig()
        config.validate()

    def testSnapPsfMatchConfig(self):
        config = ipDiffim.SnapPsfMatchConfig()
        config.validate()

    def testSnapPsfMatchConfigDF(self):
        config = ipDiffim.SnapPsfMatchConfigDF()
        config.validate()

    def testSnapPsfMatchConfigAL(self):
        config = ipDiffim.SnapPsfMatchConfigAL()
        config.validate()

    def testModelPsfMatchConfig(self):
        ipDiffim.ModelPsfMatchConfig()

    def testDetectionConfig(self):
        config = ipDiffim.DetectionConfig()
        config.fpGrowPix = 10
        config.validate()

    #

    def testPsfMatchConfig(self):
        config = ipDiffim.PsfMatchConfig()
        config.validate()

    def testPsfMatchConfigAL(self):
        config = ipDiffim.PsfMatchConfigAL()
        config.validate()

    def testPsfMatchConfigDF(self):
        config = ipDiffim.PsfMatchConfigDF()
        config.validate()

    #

#####


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
