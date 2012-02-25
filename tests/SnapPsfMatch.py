#!/usr/bin/env python
import unittest
import lsst.utils.tests as tests
import lsst.afw.image as afwImage
import lsst.afw.image.utils as imageUtils
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.daf.base as dafBase
import lsst.pex.policy as pexPolicy

import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)

class PsfMatchTestCases(unittest.TestCase):

    def setUp(self):
        self.configAL    = ipDiffim.SnapPsfMatchTask.ConfigClass()
        self.configAL.kernel.name = "AL"
        self.subconfigAL = self.configAL.kernel.active

        self.configDF    = ipDiffim.SnapPsfMatchTask.ConfigClass()
        self.configDF.kernel.name = "DF"
        self.subconfigDF = self.configDF.kernel.active

        self.configDFr    = ipDiffim.SnapPsfMatchTask.ConfigClass()
        self.configDFr.kernel.name = "DF"
        self.subconfigDFr = self.configDFr.kernel.active

        self.subconfigDF.useRegularization = False
        self.subconfigDFr.useRegularization = True

        # variance is a hack
        self.subconfigAL.singleKernelClipping   = False   
        self.subconfigAL.spatialKernelClipping  = False  
        self.subconfigDF.singleKernelClipping   = False   
        self.subconfigDF.spatialKernelClipping  = False  
        self.subconfigDFr.singleKernelClipping  = False   
        self.subconfigDFr.spatialKernelClipping = False  

        # Send fake kernel a differential background
        self.bgValue = 100.
        self.subconfigAL.fitForBackground = True
        self.subconfigDF.fitForBackground = True
        self.subconfigDFr.fitForBackground = True

    def testSnap(self):
        tMi, sMi, sK, kcs, confake = diffimTools.makeFakeKernelSet(bgValue = self.bgValue)

        psfMatchAL   = ipDiffim.SnapPsfMatchTask(self.subconfigAL)
        psfMatchDF   = ipDiffim.SnapPsfMatchTask(self.subconfigDF)
        psfMatchDFr  = ipDiffim.SnapPsfMatchTask(self.subconfigDFr)
        resultsAL    = psfMatchAL.run(tMi, sMi, "subtractMaskedImages")
        resultsDF    = psfMatchDF.run(tMi, sMi, "subtractMaskedImages")
        resultsDFr   = psfMatchDFr.run(tMi, sMi, "subtractMaskedImages")

    def tearDown(self):
        del self.configAL 
        del self.configDF 
        del self.configDFr

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
