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
        diffimDir         = eups.productDir('ip_diffim')
        self.policyPath   = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')
        self.p0           = pexPolicy.Policy.createPolicy(self.policyPath)

    def tearDown(self):
        del self.p0
        
    def testNoModify(self):
        p1 = ipDiffim.generateDefaultPolicy(self.policyPath, modify=False)
        self.assertEqual(self.p0.get("kernelRows"), p1.get("kernelRows"))
        self.assertEqual(self.p0.get("kernelCols"), p1.get("kernelCols"))
        self.assertEqual(self.p0.get("fpGrowPix"), p1.get("fpGrowPix"))
        for i in range(self.p0.get("alardNGauss")):
            self.assertEqual(self.p0.getDoubleArray("alardSigGauss")[i],
                             p1.getDoubleArray("alardSigGauss")[i])

    def testModifyGreater(self, fwhm=10.):
        p1 = ipDiffim.generateDefaultPolicy(self.policyPath, fwhm=fwhm)
        self.assertTrue(self.p0.get("kernelRows") < p1.get("kernelRows"))
        self.assertTrue(self.p0.get("kernelCols") < p1.get("kernelCols"))
        self.assertTrue(self.p0.get("fpGrowPix") < p1.get("fpGrowPix"))
        for i in range(self.p0.get("alardNGauss")):
            self.assertTrue(self.p0.getDoubleArray("alardSigGauss")[i] <
                            p1.getDoubleArray("alardSigGauss")[i])

        # maxed out the sizes
        self.assertTrue(p1.get("fpGrowPix") == self.p0.get("fpGrowMax"))
        self.assertTrue(p1.get("kernelRows")//2 == self.p0.get("kernelRadiusMax"))
        self.assertTrue(p1.get("kernelCols")//2 == self.p0.get("kernelRadiusMax"))

    def testModifyLesser(self, fwhm=1.):
        p1 = ipDiffim.generateDefaultPolicy(self.policyPath, fwhm=fwhm)
        self.assertTrue(self.p0.get("kernelRows") > p1.get("kernelRows"))
        self.assertTrue(self.p0.get("kernelCols") > p1.get("kernelCols"))

        self.assertTrue(self.p0.get("fpGrowPix") > p1.get("fpGrowPix"))
        for i in range(self.p0.get("alardNGauss")):
            self.assertTrue(self.p0.getDoubleArray("alardSigGauss")[i] >
                            p1.getDoubleArray("alardSigGauss")[i])

        # minned out the sizes
        self.assertTrue(p1.get("fpGrowPix") == self.p0.get("fpGrowMin"))
        self.assertTrue(p1.get("kernelRows")//2 == self.p0.get("kernelRadiusMin"))
        self.assertTrue(p1.get("kernelCols")//2 == self.p0.get("kernelRadiusMin"))
            
        
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
    
