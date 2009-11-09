#!/usr/bin/env python
import os
import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.policy as pexPolicy
import lsst.pex.logging as logging

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy      = ipDiffim.generateDefaultPolicy(diffimPolicy)
        self.kCols       = self.policy.getInt('kernelCols')
        self.kRows       = self.policy.getInt('kernelRows')

        # gaussian reference kernel
        self.gSize         = self.kCols
        self.gaussFunction = afwMath.GaussianFunction2D(2, 3)
        self.gaussKernel   = afwMath.AnalyticKernel(self.gSize, self.gSize, self.gaussFunction)

        # known input images
        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:
            defImagePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                        "v5-e0-c011-a00.sci")
            self.templateImage  = afwImage.MaskedImageF(defImagePath)
            self.scienceImage   = self.templateImage.Factory( self.templateImage.getDimensions() )
            
            afwMath.convolve(self.scienceImage, self.templateImage, self.gaussKernel, False)

    def tearDown(self):
        del self.policy

    def testGetCollection(self):

        # NOTE - you need to subtract off background from the image
        # you run detection on.  Here it is the template.
        algorithm   = self.policy.get("backgroundPolicy.algorithm")
        binsize     = self.policy.get("backgroundPolicy.binsize")
        undersample = self.policy.get("backgroundPolicy.undersample")
        bctrl       = afwMath.BackgroundControl(algorithm)
        bctrl.setNxSample(self.templateImage.getWidth()//binsize + 1)
        bctrl.setNySample(self.templateImage.getHeight()//binsize + 1)
        bctrl.setUndersampleStyle(undersample)

        image   = self.templateImage.getImage() 
        backobj = afwMath.makeBackground(image, bctrl)
        image  -= backobj.getImageF()
        del image; del backobj
        
        fpList1 = ipDiffim.getCollectionOfFootprintsForPsfMatching(self.templateImage,
                                                                   self.scienceImage,
                                                                   self.policy)
        self.assertTrue(len(fpList1) != 0)

        for fp in fpList1:
            bbox = fp.getBBox()
            tmi  = afwImage.MaskedImageF(self.templateImage, bbox)
            smi  = afwImage.MaskedImageF(self.scienceImage, bbox)
            tmask = tmi.getMask()
            smask = smi.getMask()

            for j in range(tmask.getHeight()):
                for i in range(tmask.getWidth()):
                    # No masked pixels in either image
                    self.assertEqual(tmask.get(i, j), 0)
                    self.assertEqual(smask.get(i, j), 0)

        # add a masked pixel to the template image and make sure you don't get it
        afwImage.MaskedImageF(self.templateImage, fpList1[0].getBBox()).getMask().set(tmask.getWidth()//2,
                                                                                      tmask.getHeight()//2,
                                                                                      0x1)
        fpList2 = ipDiffim.getCollectionOfFootprintsForPsfMatching(self.templateImage,
                                                                   self.scienceImage,
                                                                   self.policy)
        self.assertTrue(len(fpList2) == (len(fpList1)-1))

        # add a masked pixel to the science image and make sure you don't get it
        afwImage.MaskedImageF(self.scienceImage, fpList1[1].getBBox()).getMask().set(smask.getWidth()//2,
                                                                                     smask.getHeight()//2,
                                                                                     0x1)
        afwImage.MaskedImageF(self.scienceImage, fpList1[2].getBBox()).getMask().set(smask.getWidth()//2,
                                                                                     smask.getHeight()//2,
                                                                                     0x1)
        fpList3 = ipDiffim.getCollectionOfFootprintsForPsfMatching(self.templateImage,
                                                                   self.scienceImage,
                                                                   self.policy)
        self.assertTrue(len(fpList3) == (len(fpList1)-3))

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
