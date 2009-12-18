#!/usr/bin/env python
import os
import sys
import pdb

import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.policy as pexPolicy
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as logging

import lsst.afw.display.ds9 as ds9

Verbosity = 1
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

# This tests core functionality of PsfMatchingFunctor

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
        
    def setUp(self):
        self.policy    = ipDiffim.generateDefaultPolicy(diffimPolicy, modify=False)
        self.kCols     = self.policy.getInt('kernelCols')
        self.kRows     = self.policy.getInt('kernelRows')
        self.basisList = ipDiffim.generateDeltaFunctionBasisSet(self.kCols, self.kRows)

        # gaussian reference kernel
        self.gSize         = self.kCols
        self.gaussFunction = afwMath.GaussianFunction2D(2, 3)
        self.gaussKernel   = afwMath.AnalyticKernel(self.gSize, self.gSize, self.gaussFunction)
        self.kImageIn      = afwImage.ImageD(self.gSize, self.gSize)
        self.gaussKernel.computeImage(self.kImageIn, False)

        # difference imaging functor
        self.kFunctor      = ipDiffim.PsfMatchingFunctorF(self.basisList)

        # known input images
        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:
            defSciencePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                          "v26-e0-c011-a00.sci")
            defTemplatePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                           "v5-e0-c011-a00.sci")
            
            self.scienceImage   = afwImage.ExposureF(defSciencePath)
            self.templateImage  = afwImage.ExposureF(defTemplatePath)
            self.templateImage  = ipDiffim.warpTemplateExposure(self.templateImage,
                                                                self.scienceImage,
                                                                self.policy)

    def tearDown(self):
        del self.policy

    def doNormalize(self, xloc=397, yloc=580, imscale=2.):
        imsize = int(3 * self.kCols)
        bbox = afwImage.BBox( afwImage.PointI(xloc - imsize/2,
                                              yloc - imsize/2),
                              afwImage.PointI(xloc + imsize/2,
                                              yloc + imsize/2) )
        
        smi  = afwImage.MaskedImageF(self.scienceImage.getMaskedImage(),  bbox)
        tmi  = afwImage.MaskedImageF(self.templateImage.getMaskedImage(), bbox)

        # make sure its a valid subregion!
        mask     = smi.getMask()
        for j in range(mask.getHeight()):
            for i in range(mask.getWidth()):
                self.assertEqual(smi.getMask().get(i, j), 0)
                self.assertEqual(tmi.getMask().get(i, j), 0)

        smi *= imscale
        
        # estimate of the variance
        var  = afwImage.MaskedImageF(smi, True)
        var -= tmi

        # accepts : imageToConvolve, imageToNotConvolve
        self.kFunctor.apply(tmi.getImage(), smi.getImage(), var.getVariance(), self.policy)
        kb1       = self.kFunctor.getSolution()
        kSoln1    = kb1.first
        bgSoln1   = kb1.second
        kImageOut = afwImage.ImageD(self.kCols, self.kRows)
        kSum1     = kSoln1.computeImage(kImageOut, False)
        self.assertNotEqual(kSum1, 1.)

        # test normalization of the coefficient vectors
        self.kFunctor.normalizeKernel()

        kb2       = self.kFunctor.getSolution()
        kSoln2    = kb2.first
        bgSoln2   = kb2.second
        kSum2     = kSoln2.computeImage(kImageOut, False)
        self.assertAlmostEqual(kSum2, 1.)

        # the background should not change!
        self.assertEqual(bgSoln1, bgSoln2)

    def testNormalize(self):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up"
            return
        
        self.doNormalize()
        
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
