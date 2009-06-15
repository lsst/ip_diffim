#!/usr/bin/env python
import os, pdb

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

Verbosity = 4
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
        
    def setUp(self):
        self.policy    = pexPolicy.Policy.createPolicy(diffimPolicy)
        self.kCols     = self.policy.getInt('kernelCols')
        self.kRows     = self.policy.getInt('kernelRows')
        self.basisList = ipDiffim.generateDeltaFunctionKernelSet(self.kCols, self.kRows)

        # gaussian reference kernel
        self.gSize         = self.kCols
        self.gaussFunction = afwMath.GaussianFunction2D(2, 3)
        self.gaussKernel   = afwMath.AnalyticKernel(self.gSize, self.gSize, self.gaussFunction)
        self.kImageIn      = afwImage.ImageD(self.gSize, self.gSize)
        self.gaussKernel.computeImage(self.kImageIn, False)

        # difference imaging functor
        self.kFunctor      = ipDiffim.PsfMatchingFunctorF(self.basisList)

        # known input images
        defDataDir = eups.productDir('afwdata')
        defSciencePath = os.path.join(defDataDir, "CFHT", "D4", 
                                      "cal-53535-i-797722_1")
        self.scienceImage  = afwImage.MaskedImageF(defSciencePath)
        
        
    def tearDown(self):
        del self.policy

    def doGaussian(self, imsize=50, xloc=1118, yloc=2483):
        # NOTE : the size of these images have to be bigger
        #        size you lose pixels due to the convolution with the gaussian
        #        so adjust the size a bit to compensate 
        imsize += self.gSize

        # chop out a region around a known object
        bbox = afwImage.BBox( afwImage.PointI(xloc - imsize/2,
                                              yloc - imsize/2),
                              afwImage.PointI(xloc + imsize/2,
                                              yloc + imsize/2) )
        tmi  = afwImage.MaskedImageF(self.scienceImage, bbox)

        # now convolve it with a gaussian to make a science image
        cmi = afwImage.MaskedImageF(imsize, imsize)
        afwMath.convolve(cmi, tmi, self.gaussKernel, False)

        # grab only the non-masked subregion
        bbox     = afwImage.BBox(afwImage.PointI(self.gaussKernel.getCtrX(),
                                                 self.gaussKernel.getCtrY()) ,
                                 afwImage.PointI(imsize+1 - (self.gaussKernel.getWidth()  - self.gaussKernel.getCtrX()),
                                                 imsize+1 - (self.gaussKernel.getHeight() - self.gaussKernel.getCtrY())))
                                 
        tmi2     = afwImage.MaskedImageF(tmi, bbox)
        cmi2     = afwImage.MaskedImageF(cmi, bbox)

        # OUTPUT
        self.kImageIn.writeFits('k1.fits')
        tmi2.writeFits('t')
        cmi2.writeFits('c')

        # make sure its a valid subregion!
        mask     = cmi2.getMask()
        for j in range(mask.getHeight()):
            for i in range(mask.getWidth()):
                self.assertEqual(mask.get(i, j), 0)
                
        # estimate of the variance
        var  = afwImage.MaskedImageF(cmi2, True)
        var -= tmi2

        # accepts : imageToConvolve, imageToNotConvolve
        self.kFunctor.apply(tmi2.getImage(), cmi2.getImage(), var.getVariance(), self.policy)
        kernel    = self.kFunctor.getKernel()
        kImageOut = afwImage.ImageD(self.kCols, self.kRows)
        kSum      = kernel.computeImage(kImageOut, False)
        diffIm    = ipDiffim.convolveAndSubtract(tmi2, cmi2, kernel, self.kFunctor.getBackground())

        # OUTPUT
        kImageOut.writeFits('k2a.fits')
        diffIm.writeFits('da')

        # Second iteration
        self.kFunctor.apply(tmi2.getImage(), cmi2.getImage(), diffIm.getVariance(), self.policy)
        kernel    = self.kFunctor.getKernel()
        kSum      = kernel.computeImage(kImageOut, False)
        diffIm    = ipDiffim.convolveAndSubtract(tmi2, cmi2, kernel, self.kFunctor.getBackground())
        kImageOut.writeFits('k2b.fits')
        diffIm.writeFits('db')

        # OUTPUT
        tmi2.getImage().writeFits('t.fits')
        cmi2.getImage().writeFits('c.fits')

        # make sure the derived kernel looks like the input kernel
        for j in range(kImageOut.getHeight()):
            for i in range(kImageOut.getWidth()):
                if self.kImageIn.get(i,j) > 1e-3:
                    self.assertAlmostEqual(kImageOut.get(i, j)/self.kImageIn.get(i,j), 1.0, 1)


            
    def testGaussian(self):
        self.doGaussian()

        
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
