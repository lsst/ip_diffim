#!/usr/bin/env python
import os
import pdb
import sys
import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as logging
import lsst.ip.diffim.diffimTools as diffimTools

verbosity = 5
logging.Trace_setVerbosity('lsst.ip.diffim', verbosity)

class DiffimTestCases(unittest.TestCase):
    def setUp(self):
        self.policy       = ipDiffim.makeDefaultPolicy()

        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:

            defTemplatePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                           "v5-e0-c011-a00.sci")
            defSciencePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                          "v26-e0-c011-a00.sci")
            
            self.scienceImage   = afwImage.ExposureF(defSciencePath)
            self.templateImage  = afwImage.ExposureF(defTemplatePath)
            
            diffimTools.backgroundSubtract(self.policy.getPolicy("afwBackgroundPolicy"),
                                           [self.templateImage.getMaskedImage(),
                                            self.scienceImage.getMaskedImage()])
   
    def tearDown(self):
        del self.policy
        if self.defDataDir:
            del self.scienceImage
            del self.templateImage

 
    def testRatings(self):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up"
            return

        psfmatch = ipDiffim.ImagePsfMatch(self.policy)
        results  = psfmatch.subtractExposures(self.templateImage, self.scienceImage,
                                              doWarping = True)
        differenceExposure, spatialKernel, backgroundModel, kernelCellSet = results
        ratings = ipDiffim.makeRatingVector(kernelCellSet, spatialKernel, backgroundModel)

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
