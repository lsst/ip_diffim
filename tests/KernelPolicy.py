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
        diffimDir    = eups.productDir('ip_diffim')
        policyPath   = os.path.join(diffimDir, 'policy', 'PsfMatchingDictionary.paf')
        defPolicy    = pexPolicy.Policy.createPolicy(policyPath)
        self.p0.mergeDefaults(defPolicy.getDictionary())

    def tearDown(self):
        del self.p0
        
    def testNoModify(self):
        p1 = ipDiffim.makeDefaultPolicy(modify=False)
        self.assertEqual(self.p0.get("kernelSize"), p1.get("kernelSize"))
        self.assertEqual(self.p0.getPolicy("detectionPolicy").get("fpGrowPix"),
                         p1.getPolicy("detectionPolicy").get("fpGrowPix"))
        
        for i in range(self.p0.get("alardNGauss")):
            self.assertEqual(self.p0.getDoubleArray("alardSigGauss")[i],
                             p1.getDoubleArray("alardSigGauss")[i])

    def testModifyGreater(self, fwhm=20.):
        p1 = ipDiffim.makeDefaultPolicy(fwhm = fwhm, modify = True)
        
        self.assertTrue(self.p0.get("kernelSize") < p1.get("kernelSize"))
        self.assertTrue(self.p0.getPolicy("detectionPolicy").get("fpGrowPix") <
                        p1.getPolicy("detectionPolicy").get("fpGrowPix"))

        for i in range(self.p0.get("alardNGauss")):
            self.assertTrue(self.p0.getDoubleArray("alardSigGauss")[i] <
                            p1.getDoubleArray("alardSigGauss")[i])

        # maxed out the sizes
        self.assertTrue(p1.getPolicy("detectionPolicy").get("fpGrowPix") ==
                        p1.getPolicy("detectionPolicy").get("fpGrowMax"))
        self.assertTrue(p1.get("kernelSize") == p1.get("kernelSizeMax"))

    def testModifyLesser(self, fwhm=1.):
        p1 = ipDiffim.makeDefaultPolicy(fwhm = fwhm, modify = True)

        self.assertTrue(self.p0.get("kernelSize") > p1.get("kernelSize"))
        self.assertTrue(self.p0.getPolicy("detectionPolicy").get("fpGrowPix") >
                        p1.getPolicy("detectionPolicy").get("fpGrowPix"))
        
        for i in range(self.p0.get("alardNGauss")):
            self.assertTrue(self.p0.getDoubleArray("alardSigGauss")[i] >
                            p1.getDoubleArray("alardSigGauss")[i])

        # minned out the sizes
        self.assertTrue(p1.getPolicy("detectionPolicy").get("fpGrowPix") ==
                        p1.getPolicy("detectionPolicy").get("fpGrowMin"))
        self.assertTrue(p1.get("kernelSize") == p1.get("kernelSizeMin"))

    def testModify1(self):
        p1 = ipDiffim.modifyForDeconvolution(self.p0)
        self.assertTrue(self.p0.getIntArray('alardDegGauss') != p1.getIntArray('alardDegGauss'))
        
    def testModify2(self):
        p1 = ipDiffim.modifyForModelPsfMatch(self.p0)
        self.assertTrue(self.p0.get('singleKernelClipping') != p1.get('singleKernelClipping'))
        
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
    
