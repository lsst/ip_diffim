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
import lsst.afw.math as afwMath
import lsst.daf.base as dafBase

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Verbosity = 1
pexLog.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

class DiffimStageTestCase(unittest.TestCase):

    def setUp(self):
        # processing policy
        self.policy = pexPolicy.Policy()

        # PIPELINE INPUTS
        self.policy.add('scienceExposureKey', 'scienceExposure0')
        self.policy.add('templateExposureKey', 'templateExposure0')
        # ISR PROCESSING
        path = os.path.join(os.environ["IP_DIFFIM_DIR"],
                            "pipeline", 
                            "ImageSubtractStageDictionary.paf")
        diffimPolicyFile = pexPolicy.PolicyFile(path)
        
        diffimPolicy = pexPolicy.Policy(diffimPolicyFile)
        self.policy.add('diffimPolicy', diffimPolicy)
         
        # OUTPUTS
        self.policy.add('differenceExposureKey', 'differenceExposure0')
        self.policy.add('sdqaRatingSetKey',      'sdqaRatingSet0')

        clipboard = pexClipboard.Clipboard()
              
        # create clipboard and fill 'er up!
        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:
            defSciencePath = os.path.join(self.defDataDir, "CFHT", "D4", 
                                          "cal-53535-i-797722_1")
            defTemplatePath = defSciencePath + "_tmpl"
            
            self.policy.set("diffimPolicy.spatialKernelOrder", 1)
            self.policy.set("diffimPolicy.sizeCellX", 128)
            self.policy.set("diffimPolicy.sizeCellY", 128)
            bbox = afwImage.BBox(afwImage.PointI(32,32), 512, 512)
            scienceExposure = afwImage.ExposureF(defSciencePath, 0, bbox)
            templateExposure = afwImage.ExposureF(defTemplatePath)

            # NOTE - you need to subtract off background from the image
            # you run detection on.  Here it is the template.
            algorithm = self.policy.get("diffimPolicy.backgroundPolicy.algorithm")
            binsize   = self.policy.get("diffimPolicy.backgroundPolicy.binsize")
            bctrl     = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE)
            bctrl.setNxSample(int(templateExposure.getWidth()//binsize) + 1)
            bctrl.setNySample(int(templateExposure.getHeight()//binsize) + 1)

            image   = templateExposure.getMaskedImage().getImage() 
            backobj = afwMath.makeBackground(image, bctrl)
            image  -= backobj.getImageF()

            image   = scienceExposure.getMaskedImage().getImage() 
            backobj = afwMath.makeBackground(image, bctrl)
            image  -= backobj.getImageF()
            
            del image; del backobj
        
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
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up; not running StageTest.py"
            return
        
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
