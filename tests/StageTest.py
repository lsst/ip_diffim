#!/usr/bin/env python
"""
Run with:
   python DiffimStageTest.py
"""

import sys, os, math

import eups
import pdb
import unittest

import lsst.utils.tests as utilsTests
import lsst.pex.harness.Queue as pexQueue
import lsst.pex.harness.Clipboard as pexClipboard
import lsst.pex.policy as pexPolicy
import lsst.pex.logging as pexLog
import lsst.ip.diffim.diffimStages as diffimStages
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Verbosity = 4
pexLog.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

class DiffimStageTestCase(unittest.TestCase):

    def setUp(self):
        # processing policy
        self.policy = pexPolicy.Policy()

        # PIPELINE INPUTS
        self.policy.add('scienceExposureKey', 'scienceExposure0')
        self.policy.add('templateExposureKey', 'templateExposure0')
        # ISR PROCESSING
        path = os.path.join(eups.productDir("ip_diffim"),
                            "pipeline", 
                            "ImageSubtractStageDictionary.paf")
        diffimPolicyFile = pexPolicy.PolicyFile(path)
        
        diffimPolicy = pexPolicy.Policy(diffimPolicyFile)
        self.policy.add('diffimPolicy', diffimPolicy)
         
        # OUTPUTS
        self.policy.add('differenceExposureKey', 'differenceExposure0')
        self.policy.add('sdqaRatingSetKey',      'sdqaRatingSet0')
              
        # create clipboard and fill 'er up!
        defDataDir = eups.productDir('afwdata')
        defSciencePath = os.path.join(defDataDir, "CFHT", "D4", 
                                      "cal-53535-i-797722_1")
        defTemplatePath = defSciencePath + "_tmpl"
        
        scienceExposure = afwImage.ExposureF(defSciencePath)
        templateExposure = afwImage.ExposureF(defTemplatePath)
        
        clipboard = pexClipboard.Clipboard()
        clipboard.put(self.policy.get('scienceExposureKey'), scienceExposure)
        clipboard.put(self.policy.get('templateExposureKey'), templateExposure)

        inQueue = pexQueue.Queue()
        inQueue.addDataset(clipboard)
        self.outQueue = pexQueue.Queue()
       
        self.stage = diffimStages.DiffimStage(0, self.policy)
        self.stage.initialize(self.outQueue, inQueue)
        self.stage.setUniverseSize(1)
        self.stage.setRun('SingleExposureTest')


    def tearDown(self):
        del self.stage
        del self.outQueue

    def testSingleInputExposure(self):
        self.stage.process()
        clipboard = self.outQueue.getNextDataset()
        assert(clipboard.contains(self.policy.getString('differenceExposureKey')))
        assert(clipboard.contains(self.policy.getString('sdqaRatingSetKey')))

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(DiffimStageTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
