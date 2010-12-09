#!/usr/bin/env python
import unittest
import lsst.utils.tests as tests
import lsst.ip.diffim as ipDiffim
import lsst.afw.image as afwImage
import eups
import os
import sys
#import lsst.afw.display.ds9 as ds9

class DiffimTestCases(unittest.TestCase):
    def setUp(self):
        self.policy       = ipDiffim.createDefaultPolicy()
        
        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:

            defTemplatePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                           "v5-e0-c011-a00.sci")
            defSciencePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                          "v26-e0-c011-a00.sci")
            
            self.scienceImage   = afwImage.ExposureF(defSciencePath)
            self.templateImage  = afwImage.ExposureF(defTemplatePath)
 
    def tearDown(self):
        del self.policy
        if self.defDataDir:
            del self.scienceImage
            del self.templateImage

    def testWarp(self):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata not set up; not running WarpTemplateExposure.py"
            return

        # image 1 gets remapped to match up with image 2
        remappedImage = ipDiffim.warpTemplateExposure(self.templateImage,
                                                      self.scienceImage,
                                                      self.policy.getPolicy("warpingPolicy"))


        # sizes in pixels
        self.assertEqual(remappedImage.getHeight(),
                         self.scienceImage.getHeight())
        self.assertEqual(remappedImage.getWidth(),
                         self.scienceImage.getWidth())

        # sizes on the sky
        wcs1 = remappedImage.getWcs()
        wcs2 = self.scienceImage.getWcs()

        self.assertEqual(wcs1.pixelToSky(0, 0)[0],
                         wcs2.pixelToSky(0, 0)[0])
        self.assertEqual(wcs1.pixelToSky(0, 0)[1],
                         wcs2.pixelToSky(0, 0)[1])


        self.assertEqual(wcs1.pixelToSky(remappedImage.getWidth(),
                                         remappedImage.getHeight())[0],
                         wcs2.pixelToSky(remappedImage.getWidth(),
                                         remappedImage.getHeight())[0])
        self.assertEqual(wcs1.pixelToSky(remappedImage.getWidth(),
                                         remappedImage.getHeight())[1],
                         wcs2.pixelToSky(remappedImage.getWidth(),
                                         remappedImage.getHeight())[1])
                         

    def testXY0(self):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata not set up; not running WarpTemplateExposure.py"
            return

        bbox     = afwImage.BBox(afwImage.PointI(7, 900),
                                 afwImage.PointI(102, 1000))
        templateSubImage = afwImage.ExposureF(self.templateImage, bbox)
        scienceSubImage  = afwImage.ExposureF(self.scienceImage, bbox)

        # image 1 gets remapped to match up with image 2
        remappedImage = ipDiffim.warpTemplateExposure(templateSubImage,
                                                      scienceSubImage,
                                                      self.policy)


        # sizes in pixels
        self.assertEqual(remappedImage.getHeight(),
                         scienceSubImage.getHeight())
        self.assertEqual(remappedImage.getWidth(),
                         scienceSubImage.getWidth())

        # sizes on the sky
        wcs1 = remappedImage.getWcs()
        wcs2 = scienceSubImage.getWcs()

        self.assertEqual(wcs1.pixelToSky(0, 0)[0],
                         wcs2.pixelToSky(0, 0)[0])
        self.assertEqual(wcs1.pixelToSky(0, 0)[1],
                         wcs2.pixelToSky(0, 0)[1])


        self.assertEqual(wcs1.pixelToSky(remappedImage.getWidth(),
                                         remappedImage.getHeight())[0],
                         wcs2.pixelToSky(remappedImage.getWidth(),
                                         remappedImage.getHeight())[0])
        self.assertEqual(wcs1.pixelToSky(remappedImage.getWidth(),
                                         remappedImage.getHeight())[1],
                         wcs2.pixelToSky(remappedImage.getWidth(),
                                         remappedImage.getHeight())[1])
        
       
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
    
        


     
