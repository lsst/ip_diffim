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

import unittest
import eups
import os
import lsst.utils.tests as tests
import lsst.ip.diffim as ipDiffim

class DiffimTestCases(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testDiaSourceAnalystConfig(self):
        config = ipDiffim.DiaSourceAnalystConfig()
        config.validate()

    def testImagePsfMatchConfig(self):
        config = ipDiffim.ImagePsfMatchConfig()
        config.validate()

    def testSnapPsfMatchConfigDF(self):
        config = ipDiffim.SnapPsfMatchConfigDF()
        config.validate()

    def testSnapPsfMatchConfigAL(self):
        config = ipDiffim.SnapPsfMatchConfigAL()
        config.validate()

    def testModelPsfMatchConfig(self):
        config = ipDiffim.ModelPsfMatchConfig()

    def testDetectionConfig(self):
        config = ipDiffim.DetectionConfig()

        # Try invalid entry
        config.fpGrowPix = 1
        config.fpGrowMin = 2
        self.assertEqual(config.validate(), False)

    def testAfwBackgroundConfig(self):
        config = ipDiffim.AfwBackgroundConfig()
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

    def testSnapPsfMatchConfig(self):
        config = ipDiffim.SnapPsfMatchConfig()
        config.validate()

    def testSnapPsfMatchConfigAL(self):
        config = ipDiffim.SnapPsfMatchConfigAL()
        config.validate()

    def testSnapPsfMatchConfigDF(self):
        config = ipDiffim.SnapPsfMatchConfigDF()
        config.validate()

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
    
