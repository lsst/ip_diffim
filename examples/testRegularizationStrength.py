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

import lsst.afw.display.ds9 as ds9

verbosity = 5
logging.Trace_setVerbosity('lsst.ip.diffim', verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

display = True
writefits = False

# This one just creates example convolution and deconvolution kernels

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)

    def diffimQuality(self, kFunctor, tmi, smi, var, foffset=0):
        kFunctor.apply(tmi.getImage(), smi.getImage(), var.getVariance(), self.policy)
        kb         = kFunctor.getSolution()
        kernel     = kb.first
        background = kb.second
        kImageOut  = afwImage.ImageD(self.kCols, self.kRows)

        kSum      = kernel.computeImage(kImageOut, False)
        diffIm    = ipDiffim.convolveAndSubtract(tmi, smi, kernel, background)
        bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                  kernel.getCtrY()) ,
                                  afwImage.PointI(diffIm.getWidth() - (kernel.getWidth()  -
                                                                       kernel.getCtrX()),
                                                  diffIm.getHeight() - (kernel.getHeight() -
                                                                        kernel.getCtrY())))
        diffIm2   = afwImage.MaskedImageF(diffIm, bbox)
        self.dStats.apply( diffIm2 )

        if display:
            ds9.mtv(kImageOut, frame=foffset)
            ds9.mtv(diffIm2, frame=foffset+1)
        return kSum
        
        
    def setUp(self):
        self.policy      = ipDiffim.generateDefaultPolicy(diffimPolicy)
        self.kCols       = self.policy.getInt('kernelCols')
        self.kRows       = self.policy.getInt('kernelRows')
        self.basisList   = ipDiffim.generateDeltaFunctionBasisSet(self.kCols, self.kRows)

        # Regularization terms
        forward, central     = 0, 1
        noWrap, wrap, taper  = 0, 1, 2
        boundaryStyle, diffStyle = taper, forward
        self.h0 = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, 0,
                                                                  boundaryStyle, diffStyle)
        self.h1 = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, 1,
                                                                  boundaryStyle, diffStyle)
        self.h2 = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, 2,
                                                                  boundaryStyle, diffStyle)
        
        # difference imaging functor
        self.kFunctor      = ipDiffim.PsfMatchingFunctorF(self.basisList)
        self.kFunctor0     = ipDiffim.PsfMatchingFunctorF(self.basisList, self.h0)
        self.kFunctor1     = ipDiffim.PsfMatchingFunctorF(self.basisList, self.h1)
        self.kFunctor2     = ipDiffim.PsfMatchingFunctorF(self.basisList, self.h2)

        # known input images
        defDataDir = eups.productDir('afwdata')
        defSciencePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                      "v26-e0-c011-a00.sci")
        defTemplatePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                       "v5-e0-c011-a00.sci")
        self.scienceImage   = afwImage.ExposureF(defSciencePath)
        self.templateImage  = afwImage.ExposureF(defTemplatePath)

        # Remap the template to the image; replace self.templateImage with warped image
        wKernel = afwMath.makeWarpingKernel('lanczos4')
        self.remappedImage = self.templateImage.Factory(
            self.scienceImage.getWidth(), 
            self.scienceImage.getHeight(),
            self.scienceImage.getWcs())
        self.remappedImage.getMaskedImage().setXY0( self.scienceImage.getMaskedImage().getXY0() )
        afwMath.warpExposure(self.remappedImage, 
                             self.templateImage, 
                             wKernel)
        self.templateImage = self.remappedImage

        # image statistics
        self.dStats  = ipDiffim.ImageStatisticsF()
        
    def tearDown(self):
        del self.policy

    def applyFunctor(self, imscale=4, invert=False, foffset=0, xloc=397, yloc=580):
        imsize = int(imscale * self.kCols)

        # chop out a region around a known object
        bbox = afwImage.BBox( afwImage.PointI(xloc - imsize/2,
                                              yloc - imsize/2),
                              afwImage.PointI(xloc + imsize/2,
                                              yloc + imsize/2) )

        if invert:
            tmi  = afwImage.MaskedImageF(self.scienceImage.getMaskedImage(),  bbox)
            smi  = afwImage.MaskedImageF(self.templateImage.getMaskedImage(), bbox)
        else:
            smi  = afwImage.MaskedImageF(self.scienceImage.getMaskedImage(),  bbox)
            tmi  = afwImage.MaskedImageF(self.templateImage.getMaskedImage(), bbox)

        # estimate of the variance
        var  = afwImage.MaskedImageF(smi, True)
        var -= tmi

        # accepts : imageToConvolve, imageToNotConvolve
        kSum = self.diffimQuality(self.kFunctor, tmi, smi, var, foffset=foffset+0)
        print 'DF : %.2f +/- %.2f (%.3f)' % (self.dStats.getMean(), self.dStats.getRms(), kSum)

        kSum = self.diffimQuality(self.kFunctor0, tmi, smi, var, foffset=foffset+2)
        print 'DFr0 : %.2f +/- %.2f (%.3f)' % (self.dStats.getMean(), self.dStats.getRms(), kSum)

        kSum = self.diffimQuality(self.kFunctor1, tmi, smi, var, foffset=foffset+4)
        print 'DFr1 : %.2f +/- %.2f (%.3f)' % (self.dStats.getMean(), self.dStats.getRms(), kSum)

        kSum = self.diffimQuality(self.kFunctor2, tmi, smi, var, foffset=foffset+6)
        print 'DFr2 : %.2f +/- %.2f (%.3f)' % (self.dStats.getMean(), self.dStats.getRms(), kSum)

    def testFunctor(self):
        lams = range(-8, 1, 2)
        #lams = [-6]
        for i in range(len(lams)):
            lam = lams[i]
            self.policy.set('regularizationScaling', 1.0 * 10**lam)
            self.applyFunctor(foffset=i*8)
                
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
    if '-d' in sys.argv:
        display = True
    if '-w' in sys.argv:
        writefits = True
        
    run(True)
