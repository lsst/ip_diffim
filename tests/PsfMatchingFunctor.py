#!/usr/bin/env python
import os

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

######### INCOMPLETE

Verbosity = 4
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = pexPolicy.Policy.createPolicy(diffimPolicy)
        self.kCols  = self.policy.getInt('kernelCols')
        self.kRows  = self.policy.getInt('kernelRows')
        self.basisList = ipDiffim.generateDeltaFunctionKernelSet(self.kCols, self.kRows)

        self.gSize         = 10
        self.gaussFunction = afwMath.GaussianFunction2D(2, 3)
        self.gaussKernel   = afwMath.AnalyticKernel(self.gSize, self.gSize, self.gaussFunction)
        
    def tearDown(self):
        del self.policy

    def doDeltaFunction(self, pixX, pixY, scaling, background, imsize=50):
        # image to convolve with kernel
        mi = afwImage.MaskedImageF(imsize, imsize)
        mi.set(1, 0x0, 1)

        # hot pixel 
        mi.set(pixX, pixY, (scaling, 0x0, scaling))

        # convolve with gaussian
        cmi = afwImage.MaskedImageF(imsize, imsize)
        afwMath.convolve(cmi, mi, self.gaussKernel, False)

        # grab only the non-masked subregion
        bbox     = afwImage.BBox(afwImage.PointI(self.gSize/2, self.gSize/2),
                                 afwImage.PointI(imsize-self.gSize/2,
                                                 imsize-self.gSize/2))
        cmi1     = afwImage.MaskedImageF(cmi, bbox)
        mask     = cmi1.getMask()
        for j in range(mask.getHeight()):
            for i in range(mask.getWidth()):
                self.assertEqual(mask.get(i, j), 0)
        # should be 41 x 41 now for imsize=50
        self.assertEqual(mask.getHeight(), 41)
        self.assertEqual(mask.getWidth(),  41)

        # deep copy
        cmi2   = afwImage.MaskedImageF(cmi1, True)
        cmi2  *= scaling
        cmi2  += background

        # variance estimate
        cmi3   = afwImage.MaskedImageF(cmi1, True)
        cmi3  -= cmi2

        # D = I - (K.x.T + bg)
        # 
        # accepts : imageToConvolve, imageToNotConvolve
        kFunctor = ipDiffim.PsfMatchingFunctorGslF(self.basisList)
        kFunctor.apply(cmi2, cmi1, cmi3.getVariance(), self.policy)

        kImage   = afwImage.ImageD(self.kCols, self.kRows)
        kernel   = kFunctor.getKernel()
        kSum     = kernel.computeImage(kImage, False)

        # only 3 sig digits
        self.assertAlmostEqual(kFunctor.getBackground(), background)
        self.assertAlmostEqual(kSum, scaling)

        for j in range(kImage.getHeight()):
            for i in range(kImage.getWidth()):
                if i == pixX and j == pixY:
                    self.assertAlmostEqual(kImage.get(i, j), scaling)
                else:
                    self.assertAlmostEqual(kImage.get(i, j), 0.)
                    
            
    def testDeltaFunction1(self):
        # hot central pixel
        self.doDeltaFunction(25, 25, 100, 0)

    def testDeltaFunction2(self):
        # hot central pixel + bg
        self.doDeltaFunction(25, 25, 100, 7)

    def testDeltaFunction3(self):
        # hot central pixel + bg
        self.doDeltaFunction(25, 25, 0.1, 7)

#    def testDeltaFunction2(self):
#        # hot off center pixel
#        self.doDeltaFunction(15, 15, 100, 0)


        
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
