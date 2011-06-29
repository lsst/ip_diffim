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
import lsst.pex.policy as pexPolicy

class DiffimTestCases(unittest.TestCase):
    def setUp(self):
        self.p0      = pexPolicy.Policy()
        diffimDir    = eups.productDir("ip_diffim")
        policyPath   = os.path.join(diffimDir, "policy", "PsfMatchingDictionary.paf")
        defPolicy    = pexPolicy.Policy.createPolicy(policyPath)
        self.p0.mergeDefaults(defPolicy.getDictionary())

    def tearDown(self):
        del self.p0
        
    def testMerge(self):
        # Invalid; parameter does not exist
        mergePolicy = pexPolicy.Policy()
        mergePolicy.set("testMerge", True)
        try:
            p1 = ipDiffim.makeDefaultPolicy(mergePolicy = mergePolicy)
        except:
            pass
        else:
            self.fail()

        # Change parameter by sending policy
        mergePolicy = pexPolicy.Policy()
        mergePolicy.set("fitForBackground", not self.p0.get("fitForBackground"))
        p2 = ipDiffim.makeDefaultPolicy(mergePolicy = mergePolicy)
        self.assertTrue(p2.get("fitForBackground") == (not self.p0.get("fitForBackground")))

        # Change parameter by sending path
        mergePolicyPath = os.path.join(os.getenv("IP_DIFFIM_DIR"), "policy", "DeconvolutionPolicy.paf")
        p3 = ipDiffim.makeDefaultPolicy(mergePolicy = mergePolicyPath)
        self.assertTrue(p3.get("modifyForDeconvolution") == True)


    def testModifyImagePsfMatch(self):
        p1 = ipDiffim.modifyForImagePsfMatch(self.p0, 3, 4)
        self.assertTrue(p1.get("modifyForImagePsfMatch") == True)

        p2 = ipDiffim.modifyForImagePsfMatch(self.p0, 4, 3)
        self.assertTrue(p2.get("modifyForImagePsfMatch") == True)
        self.assertTrue(p2.get("modifyForDeconvolution") == True)

    def testModifyDeconv(self):
        p1 = ipDiffim.modifyForDeconvolution(self.p0)
        self.assertTrue(p1.get("modifyForDeconvolution") == True)
        
    def testModifyModelPsfMatch(self):
        p1 = ipDiffim.modifyForModelPsfMatch(self.p0)
        self.assertTrue(p1.get("modifyForModelPsfMatch") == True)

    def testModifySnap(self):
        p1 = ipDiffim.modifyForSnapSubtraction(self.p0)
        self.assertTrue(p1.get("modifyForSnapSubtraction") == True)
        
    def tearDown(self):
        del self.p0
        
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
    
