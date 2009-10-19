#!/usr/bin/env python
import os, pdb, sys
import numpy as num
import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.policy as pexPolicy
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.pex.logging as logging

import lsst.afw.display.ds9 as ds9

Verbosity = 4
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

# This one tests convolve and subtract

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
        
    def setUp(self):
        self.policy      = pexPolicy.Policy.createPolicy(diffimPolicy)
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

    def runConvolveAndSubtract1(self, bgVal = 0, xloc = 408, yloc = 580):
        imsize = int(5 * self.kCols)

        p0 = afwImage.PointI(xloc - imsize/2, yloc - imsize/2)
        p1 = afwImage.PointI(xloc + imsize/2, yloc + imsize/2)
        bbox = afwImage.BBox(p0, p1)

        tmi     = afwImage.MaskedImageF(self.templateImage, bbox)
        smi     = afwImage.MaskedImageF(self.scienceImage, bbox)
        diffIm  = ipDiffim.convolveAndSubtract(tmi, smi, self.gaussKernel, bgVal)

        p0 = afwImage.PointI(self.gaussKernel.getCtrX(),
                             self.gaussKernel.getCtrY())
        p1 = afwImage.PointI(imsize - (self.gaussKernel.getWidth()  - self.gaussKernel.getCtrX()),
                             imsize - (self.gaussKernel.getHeight() - self.gaussKernel.getCtrY()))
        bbox    = afwImage.BBox(p0, p1)
        diffIm2 = afwImage.MaskedImageF(diffIm, bbox)

        # image is empty (or the additional background you subtracted off)
        for j in range(diffIm2.getHeight()):
            for i in range(diffIm2.getWidth()):
                self.assertAlmostEqual(diffIm2.getImage().get(i, j), -1.*bgVal, 4)

    def runConvolveAndSubtract2(self, bgOrder=0, xloc = 408, yloc = 580):
        imsize = int(5 * self.kCols)

        p0 = afwImage.PointI(xloc - imsize/2, yloc - imsize/2)
        p1 = afwImage.PointI(xloc + imsize/2, yloc + imsize/2)
        bbox = afwImage.BBox(p0, p1)

        tmi     = afwImage.MaskedImageF(self.templateImage, bbox)
        smi     = afwImage.MaskedImageF(self.scienceImage, bbox)
        bgFunc  = afwMath.PolynomialFunction2D(bgOrder)  # coeffs are 0. by default
        diffIm  = ipDiffim.convolveAndSubtract(tmi, smi, self.gaussKernel, bgFunc)

        p0 = afwImage.PointI(self.gaussKernel.getCtrX(),
                             self.gaussKernel.getCtrY())
        p1 = afwImage.PointI(imsize - (self.gaussKernel.getWidth()  - self.gaussKernel.getCtrX()),
                             imsize - (self.gaussKernel.getHeight() - self.gaussKernel.getCtrY()))
        bbox    = afwImage.BBox(p0, p1)
        diffIm2 = afwImage.MaskedImageF(diffIm, bbox)
        for j in range(diffIm2.getHeight()):
            for i in range(diffIm2.getWidth()):
                self.assertAlmostEqual(diffIm2.getImage().get(i, j), 0., 4)


    def testConvolveAndSubtract(self):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up"
            return

        self.runConvolveAndSubtract1(bgVal=0)
        self.runConvolveAndSubtract1(bgVal=10)
        # this one uses a function
        self.runConvolveAndSubtract2(bgOrder=0)
        self.runConvolveAndSubtract2(bgOrder=2)

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
