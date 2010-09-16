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
import lsst.utils.tests as tests
import lsst.ip.diffim as ipDiffim
import lsst.afw.image as afwImage
import eups
import os
import sys
#import lsst.afw.display.ds9 as ds9

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
                                                      self.policy)


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

        bbox     = afwImage.BBox(afwImage.PointI(2, 900),
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

        #ds9.mtv(scienceSubImage, frame=1)
        #ds9.mtv(remappedImage, frame=2)
        
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
    
        


     
