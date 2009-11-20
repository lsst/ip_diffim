#!/usr/bin/env python
import unittest
import lsst.utils.tests as tests
import lsst.ip.diffim as ipDiffim
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import eups
import os

class DiffimTestCases(unittest.TestCase):
    def setUp(self):
        self.diffimDir    = eups.productDir('ip_diffim')
        self.diffimPolicy = os.path.join(self.diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')
        self.policy       = ipDiffim.generateDefaultPolicy(self.diffimPolicy)
        
        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:

            defTemplatePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                           "v5-e0-c011-a00.sci")
            defSciencePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                          "v26-e0-c011-a00.sci")
            
            self.scienceImage   = afwImage.ExposureF(defSciencePath)
            self.templateImage  = afwImage.ExposureF(defTemplatePath)
 
    def tearDown(self):
        del self.scienceImage
        del self.templateImage

    def testWarp(self):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata not set up; not running WarpTemplateExposure.py"
            return

        # image 1 gets remapped to match up with image 2
        remappedImage = ipDiffim.warpTemplateExposure(self.templateImage,
                                                      self.scienceImage,
                                                      self.policy)


        # sizes in pixels
        self.assertEqual(remappedImage.getHeight(),
                         self.scienceImage.getHeight())
        self.assertEqual(remappedImage.getWidth(),
                         self.scienceImage.getWidth())

        # sizes on the sky
        wcs1 = remappedImage.getWcs()
        wcs2 = self.scienceImage.getWcs()

        self.assertEqual(wcs1.xyToRaDec(0,0)[0],
                         wcs2.xyToRaDec(0,0)[0])
        self.assertEqual(wcs1.xyToRaDec(0,0)[1],
                         wcs2.xyToRaDec(0,0)[1])


        self.assertEqual(wcs1.xyToRaDec(remappedImage.getWidth(),
                                        remappedImage.getHeight())[0],
                         wcs2.xyToRaDec(remappedImage.getWidth(),
                                        remappedImage.getHeight())[0])
        self.assertEqual(wcs1.xyToRaDec(remappedImage.getWidth(),
                                        remappedImage.getHeight())[1],
                         wcs2.xyToRaDec(remappedImage.getWidth(),
                                        remappedImage.getHeight())[1])
                         
                         
#####
        
def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(DiffimTestCases)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
    
        


     
