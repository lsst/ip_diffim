#!/usr/bin/env python
import os
import pdb
import sys
import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as logging
import lsst.ip.diffim.diffimTools as diffimTools

verbosity = 3
logging.Trace_setVerbosity('lsst.ip.diffim', verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

# This one tests convolve and subtract of subimages / XY0

display = False

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
        
    def setUp(self):
        self.diffimDir    = eups.productDir('ip_diffim')
        self.diffimPolicy = os.path.join(self.diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')
        self.policy       = ipDiffim.createDefaultPolicy(self.diffimPolicy)
        
        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:

            defTemplatePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                           "v5-e0-c011-a00.sci")
            defSciencePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                          "v26-e0-c011-a00.sci")
            
            self.scienceImage   = afwImage.ExposureF(defSciencePath)
            self.templateImage  = afwImage.ExposureF(defTemplatePath)
            
            diffimTools.backgroundSubtract(self.policy, [self.templateImage.getMaskedImage(),
                                                         self.scienceImage.getMaskedImage()])
            self.offset   = 1500
            self.bbox     = afwImage.BBox(afwImage.PointI(0, self.offset),
                                          afwImage.PointI(511, 2046))

    def tearDown(self):
        del self.policy
        if self.defDataDir:
            del self.scienceImage
            del self.templateImage

    def testAL(self):
        self.policy.set('kernelBasisSet', 'alard-lupton')
        self.policy.set('spatialKernelOrder', 1)
        self.policy.set('spatialBgOrder', 0) # already bg-subtracted
        self.policy.set('usePcaForSpatialKernel', False)
        self.runXY0()

    def testWarping(self):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up"
            return

        templateSubImage = afwImage.ExposureF(self.templateImage, self.bbox)
        scienceSubImage  = afwImage.ExposureF(self.scienceImage, self.bbox)
        try:
            ipDiffim.subtractExposures(templateSubImage, scienceSubImage, self.policy, doWarping = False)
        except Exception, e:
            pass
        else:
            self.fail()


    def runXY0(self):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up"
            return

        templateSubImage = afwImage.ExposureF(self.templateImage, self.bbox)
        scienceSubImage  = afwImage.ExposureF(self.scienceImage, self.bbox)

        # Have an XY0
        results1 = ipDiffim.subtractExposures(templateSubImage, scienceSubImage, self.policy,
                                              display=display, frame=0)
        differenceExposure1, spatialKernel1, backgroundModel1, kernelCellSet1 = results1

        # And then take away XY0
        templateSubImage.setXY0(0, 0) # do it to the exposure so the Wcs gets modified too
        scienceSubImage.setXY0(0, 0)
        results2 = ipDiffim.subtractExposures(templateSubImage, scienceSubImage, self.policy,
                                              display=display, frame=3)
        differenceExposure2, spatialKernel2, backgroundModel2, kernelCellSet2 = results2

        # need to count up the candidates first, since its a running tally
        count = 0
        for cell in kernelCellSet1.getCellList():
            for cand1 in cell.begin(False): 
                count += 1

        for cell in kernelCellSet1.getCellList():
            for cand1 in cell.begin(False): 
                cand1 = ipDiffim.cast_KernelCandidateF(cand1)
                cand2 = ipDiffim.cast_KernelCandidateF(kernelCellSet2.getCandidateById(cand1.getId() + count))

                # positions are the same
                self.assertEqual(cand1.getXCenter(), cand2.getXCenter())
                self.assertEqual(cand1.getYCenter(), cand2.getYCenter() + self.offset)

                # kernels are the same
                im1   = cand1.getKernelImage(ipDiffim.KernelCandidateF.RECENT)
                im2   = cand2.getKernelImage(ipDiffim.KernelCandidateF.RECENT)
                for y in range(im1.getHeight()):
                    for x in range(im1.getWidth()):
                        self.assertAlmostEqual(im1.get(x, y), im2.get(x, y))

        # Spatial fits are the same
        skp1 = spatialKernel1.getSpatialParameters()
        skp2 = spatialKernel2.getSpatialParameters()
        bgp1 = backgroundModel1.getParameters()
        bgp2 = backgroundModel2.getParameters()

        # first term = kernel sum, 0, 0
        self.assertAlmostEqual(skp1[0][0], skp2[0][0])
        # on other terms, the spatial terms are the same, the zpt terms are different
        for nk in range(1, len(skp1)):
            self.assertNotEqual(skp1[nk][0], skp2[nk][0])
            for np in range(1, len(skp1[nk])):
                self.assertAlmostEqual(skp1[nk][np], skp2[nk][np])

        for np in range(len(bgp1)):
            self.assertAlmostEqual(bgp1[np], bgp2[np])

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
    
        


     
