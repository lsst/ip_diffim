#!/usr/bin/env python
import unittest
import lsst.utils.tests as tests
import lsst.afw.image as afwImage
import lsst.afw.image.utils as imageUtils
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDet
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

        # Make ideal PSF
        self.ksize  = 21
        self.sigma = 2.0
        self.psf = afwDet.createPsf("DoubleGaussian", self.ksize, self.ksize, self.sigma)

    def makeWcs(self, offset = 0):
        # taken from $AFW_DIR/tests/testMakeWcs.py
        metadata = dafBase.PropertySet()
        metadata.set("SIMPLE",                    "T") 
        metadata.set("BITPIX",                  -32) 
        metadata.set("NAXIS",                    2) 
        metadata.set("NAXIS1",                 1024) 
        metadata.set("NAXIS2",                 1153) 
        metadata.set("RADECSYS", 'FK5')
        metadata.set("EQUINOX",                2000.)
        metadata.setDouble("CRVAL1",     215.604025685476)
        metadata.setDouble("CRVAL2",     53.1595451514076)
        metadata.setDouble("CRPIX1",     1109.99981456774 + offset)
        metadata.setDouble("CRPIX2",     560.018167811613 + offset)
        metadata.set("CTYPE1", 'RA---SIN')
        metadata.set("CTYPE2", 'DEC--SIN')
        metadata.setDouble("CD1_1", 5.10808596133527E-05)
        metadata.setDouble("CD1_2", 1.85579539217196E-07)
        metadata.setDouble("CD2_2", -5.10281493481982E-05)
        metadata.setDouble("CD2_1", -8.27440751733828E-07)
        return afwImage.makeWcs(metadata)

    def testSnap(self):
        tMi, sMi, sK, kcs, confake = diffimTools.makeFakeKernelSet(bgValue = self.bgValue)

        tWcs = self.makeWcs(offset = 0)
        sWcs = self.makeWcs(offset = 0)
        tExp = afwImage.ExposureF(tMi, tWcs)
        sExp = afwImage.ExposureF(sMi, sWcs)
	sExp.setPsf(self.psf)
        psfMatchAL   = ipDiffim.SnapPsfMatchTask(config=self.configAL)
        psfMatchDF   = ipDiffim.SnapPsfMatchTask(config=self.configDF)
        psfMatchDFr  = ipDiffim.SnapPsfMatchTask(config=self.configDFr)
        candlist = psfMatchAL.makeCandidateList(tExp, sExp)
        print len(candlist)
        resultsAL    = psfMatchAL.subtractMaskedImages(tMi, sMi, psfMatchAL.makeCandidateList(tExp, sExp))
        resultsDF    = psfMatchDF.subtractMaskedImages(tMi, sMi, psfMatchDF.makeCandidateList(tExp, sExp))
        resultsDFr   = psfMatchDFr.subtractMaskedImages(tMi, sMi, psfMatchDFr.makeCandidateList(tExp, sExp))

    def tearDown(self):
        del self.configAL 
        del self.configDF 
        del self.configDFr
	del self.psf

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
