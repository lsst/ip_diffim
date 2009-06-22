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
import lsst.pex.logging as logging

import lsst.afw.display.ds9 as ds9

Verbosity = 4
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

display = False
writefits = False

# Recover a known smoothing kernel applied to real data

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

        # edge bit
        self.edgeBit = afwImage.MaskU().getMaskPlane('EDGE')
        
    def tearDown(self):
        del self.policy

    def addNoise(self, mi):
        # use median of variance for seed
        # also the sqrt of this is used to scale the noise image
        img       = mi.getImage()

        seed      = int(afwMath.makeStatistics(mi.getVariance(), afwMath.MEDIAN).getValue())
        rdm       = afwMath.Random(afwMath.Random.MT19937, seed)
        rdmImage  = img.Factory(img.getDimensions())
        afwMath.randomGaussianImage(rdmImage, rdm)
        rdmImage *= num.sqrt(seed)
        img      += rdmImage

    def doGaussian(self, kNorm, imsize=60, xloc=1118, yloc=2483, addNoise=False):
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
        smi = afwImage.MaskedImageF(imsize, imsize)
        afwMath.convolve(smi, tmi, self.gaussKernel, kNorm, self.edgeBit)

        if addNoise:
            self.addNoise(smi)
            
        # grab only the non-masked subregion
        bbox     = afwImage.BBox(afwImage.PointI(self.gaussKernel.getCtrX(),
                                                 self.gaussKernel.getCtrY()) ,
                                 afwImage.PointI(imsize - (self.gaussKernel.getWidth()  - self.gaussKernel.getCtrX()),
                                                 imsize - (self.gaussKernel.getHeight() - self.gaussKernel.getCtrY())))

        # For testing only
        invert = False
        if invert:
            tmi2     = afwImage.MaskedImageF(smi, bbox)
            smi2     = afwImage.MaskedImageF(tmi, bbox)
        else:
            tmi2     = afwImage.MaskedImageF(tmi, bbox)
            smi2     = afwImage.MaskedImageF(smi, bbox)

        # OUTPUT
        if display:
            ds9.mtv(tmi2, frame=1)
            ds9.mtv(smi2, frame=2)
            ds9.mtv(self.kImageIn, frame=3)
        if writefits:
            self.kImageIn.writeFits('kIn.fits')
            tmi2.writeFits('t2')
            smi2.writeFits('s2')
            
        # make sure its a valid subregion!
        mask     = smi2.getMask()
        for j in range(mask.getHeight()):
            for i in range(mask.getWidth()):
                self.assertEqual(mask.get(i, j), 0)
                
        # estimate of the variance
        var  = afwImage.MaskedImageF(smi2, True)
        var -= tmi2

        # accepts : imageToConvolve, imageToNotConvolve
        self.kFunctor.apply(tmi2.getImage(), smi2.getImage(), var.getVariance(), self.policy)
        kernel    = self.kFunctor.getKernel()
        kImageOut = afwImage.ImageD(self.kCols, self.kRows)
        kSum      = kernel.computeImage(kImageOut, False)
        diffIm    = ipDiffim.convolveAndSubtract(tmi2, smi2, kernel, self.kFunctor.getBackground())
        bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                  kernel.getCtrY()) ,
                                  afwImage.PointI(diffIm.getWidth() - (kernel.getWidth()  - kernel.getCtrX()),
                                                  diffIm.getHeight() - (kernel.getHeight() - kernel.getCtrY())))
        diffIm2   = afwImage.MaskedImageF(diffIm, bbox)

        # OUTPUT
        if display:
            ds9.mtv(kImageOut, frame=4)
            ds9.mtv(diffIm2, frame=5)
        if writefits:
            kImageOut.writeFits('k1.fits')
            diffIm2.writeFits('dA2')

        # Second iteration
        self.kFunctor.apply(tmi2.getImage(), smi2.getImage(), diffIm.getVariance(), self.policy)
        kernel    = self.kFunctor.getKernel()
        kSum      = kernel.computeImage(kImageOut, False)
        diffIm    = ipDiffim.convolveAndSubtract(tmi2, smi2, kernel, self.kFunctor.getBackground())
        bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                  kernel.getCtrY()) ,
                                  afwImage.PointI(diffIm.getWidth() - (kernel.getWidth()  - kernel.getCtrX()),
                                                  diffIm.getHeight() - (kernel.getHeight() - kernel.getCtrY())))
        diffIm2   = afwImage.MaskedImageF(diffIm, bbox)
        
        # OUTPUT
        if display:
            ds9.mtv(kImageOut, frame=6)
            ds9.mtv(diffIm2, frame=7)
        if writefits:
            kImageOut.writeFits('k2.fits')
            diffIm2.writeFits('dB2')

        # kernel sum should be 1.0 if kNorm
        if kNorm:
            if addNoise:
                self.assertAlmostEqual(kSum, 1.0, 1)
            else:
                self.assertAlmostEqual(kSum, 1.0, 5)
        
        # make sure the derived kernel looks like the input kernel
        # only if you haven't normalized the kernel sum to be 1.0 during the initial convolution
        for j in range(kImageOut.getHeight()):
            for i in range(kImageOut.getWidth()):
                if not kNorm:
                    if addNoise:
                        self.assertAlmostEqual(kImageOut.get(i, j), self.kImageIn.get(i,j), 1)
                    else:
                        self.assertAlmostEqual(kImageOut.get(i, j), self.kImageIn.get(i,j), 4)


    def testGaussian(self):
        self.doGaussian(kNorm=True, addNoise=False)
        self.doGaussian(kNorm=False, addNoise=False)
        
        self.doGaussian(kNorm=True, addNoise=True)
        self.doGaussian(kNorm=False, addNoise=True)

        
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
    if '-d' in sys.argv:
        display = True
    if '-w' in sys.argv:
        writefits = True
        
    run(True)
