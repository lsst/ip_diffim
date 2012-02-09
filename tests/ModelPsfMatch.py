#!/usr/bin/env python
import unittest
import lsst.utils.tests as tests
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
import lsst.ip.diffim as ipDiffim

import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)

class PsfMatchTestCases(unittest.TestCase):

    def setUp(self):
        self.config = ipDiffim.ModelPsfMatchConfig()
        self.size   = 100
        self.exp    = afwImage.ExposureF(afwGeom.Extent2I(self.size, self.size))
        self.exp.setPsf(afwDet.createPsf("DoubleGaussian", 11, 11, 1.0))

    def testMatch(self):
        psf = afwDet.createPsf("DoubleGaussian", 11, 11, 3.0)
        psfMatch = ipDiffim.ModelPsfMatch(self.config)
        psfMatch.matchExposure(self.exp, psf)

    def testKsum(self):
        pass

    def tearDown(self):
        del self.exp
        del self.config

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(PsfMatchTestCases)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(doExit=False):
    """Run the tests"""
    tests.run(suite(), doExit)

if __name__ == "__main__":
    run(True)
