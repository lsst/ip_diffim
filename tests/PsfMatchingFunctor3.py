#!/usr/bin/env python
import os, pdb, sys

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

# This one just creates example convolution and deconvolution kernels

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
        
    def setUp(self):
        self.policy      = pexPolicy.Policy.createPolicy(diffimPolicy)
        self.kCols       = self.policy.getInt('kernelCols')
        self.kRows       = self.policy.getInt('kernelRows')
        self.fpGrowKsize = self.policy.getDouble('fpGrowKsize')
        self.basisList   = ipDiffim.generateDeltaFunctionKernelSet(self.kCols, self.kRows)
        self.H           = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, 0)
        
        self.gSize         = self.kCols
        self.gaussFunction = afwMath.GaussianFunction2D(2, 2)
        self.gaussKernel   = afwMath.AnalyticKernel(self.gSize, self.gSize, self.gaussFunction)
        self.kImageIn      = afwImage.ImageD(self.gSize, self.gSize)
        self.gaussKernel.computeImage(self.kImageIn, False)

        # difference imaging functor
        self.kFunctor      = ipDiffim.PsfMatchingFunctorF(self.basisList, self.H)

        # known input images
        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:
            defSciencePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                          "v26-e0-c011-a00.sci")
            defTemplatePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                           "v5-e0-c011-a00.sci")
            self.scienceImage   = afwImage.ExposureF(defSciencePath)
            self.templateImage  = afwImage.ExposureF(defTemplatePath)
            self.templateImage  = ipDiffim.warpTemplateExposure(self.templateImage, self.scienceImage, self.policy)
            
            # image statistics
            self.dStats  = ipDiffim.ImageStatisticsF()
        
    def tearDown(self):
        del self.policy

    def applyFunctor(self, invert=False, foffset=0, xloc=397, yloc=580):
        # NOTE : the size of these images have to be bigger
        #        size you lose pixels due to the convolution with the gaussian
        #        so adjust the size a bit to compensate 
        imsize = int(3.5 * self.kCols)

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

        # convolve science image with a gaussian for testing...
        #cmi = smi.Factory(smi.getDimensions())
        #afwMath.convolve(cmi, smi, self.gaussKernel, False)
        #smi = cmi
        
        # OUTPUT
        if display:
            ds9.mtv(tmi, frame=1+foffset)
            ds9.mtv(smi, frame=2+foffset)
        if writefits:
            tmi.writeFits('t')
            smi.writeFits('s')
            
        # estimate of the variance
        var  = afwImage.MaskedImageF(smi, True)
        var -= tmi

        #print 'Estimated variance   : %.2f %.2f -> %.2f' % (afwMath.makeStatistics(tmi.getVariance(), afwMath.MEDIAN).getValue(),
        #                                                    afwMath.makeStatistics(smi.getVariance(), afwMath.MEDIAN).getValue(),
        #                                                    afwMath.makeStatistics(var.getVariance(), afwMath.MEDIAN).getValue())
        
        # accepts : imageToConvolve, imageToNotConvolve
        self.kFunctor.apply(tmi.getImage(), smi.getImage(), var.getVariance(), self.policy)
        kernel    = self.kFunctor.getKernel()
        kImageOut = afwImage.ImageD(self.kCols, self.kRows)
        kSum      = kernel.computeImage(kImageOut, False)
        diffIm    = ipDiffim.convolveAndSubtract(tmi, smi, kernel, self.kFunctor.getBackground())
        bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                  kernel.getCtrY()) ,
                                  afwImage.PointI(diffIm.getWidth() - (kernel.getWidth()  - kernel.getCtrX()),
                                                  diffIm.getHeight() - (kernel.getHeight() - kernel.getCtrY())))
        diffIm2   = afwImage.MaskedImageF(diffIm, bbox)
        self.dStats.apply( diffIm2 )
        print 'Diffim residuals1 : %.2f +/- %.2f (%.3f)' % (self.dStats.getMean(), self.dStats.getRms(), kSum)
        #print 'Diffim variance   : %.2f %.2f -> %.2f' % (afwMath.makeStatistics(tmi.getVariance(), afwMath.MEDIAN).getValue(),
        #                                                 afwMath.makeStatistics(smi.getVariance(), afwMath.MEDIAN).getValue(),
        #                                                 afwMath.makeStatistics(diffIm2.getVariance(), afwMath.MEDIAN).getValue())
                                                       
        
        # OUTPUT
#        if display:
#            ds9.mtv(kImageOut, frame=3)
#            ds9.mtv(diffIm, frame=4)
#            ds9.mtv(diffIm2, frame=5)
#        if writefits:
#            kImageOut.writeFits('k1.fits')
#            diffIm.writeFits('dA1')
#            diffIm2.writeFits('dA2')

        # Second iteration
        self.kFunctor.apply(tmi.getImage(), smi.getImage(), diffIm.getVariance(), self.policy)
        kernel    = self.kFunctor.getKernel()
        kSum      = kernel.computeImage(kImageOut, False)
        diffIm    = ipDiffim.convolveAndSubtract(tmi, smi, kernel, self.kFunctor.getBackground())
        bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                  kernel.getCtrY()) ,
                                  afwImage.PointI(diffIm.getWidth() - (kernel.getWidth()  - kernel.getCtrX()),
                                                  diffIm.getHeight() - (kernel.getHeight() - kernel.getCtrY())))
        diffIm2   = afwImage.MaskedImageF(diffIm, bbox)
        self.dStats.apply( diffIm2 )
        print 'Diffim residuals2 : %.2f +/- %.2f (%.2f)' % (self.dStats.getMean(), self.dStats.getRms(), kSum)
        #print 'Diffim variance   : %.2f %.2f -> %.2f' % (afwMath.makeStatistics(tmi.getVariance(), afwMath.MEDIAN).getValue(),
        #                                                 afwMath.makeStatistics(smi.getVariance(), afwMath.MEDIAN).getValue(),
        #                                                 afwMath.makeStatistics(diffIm2.getVariance(), afwMath.MEDIAN).getValue())
        
        # OUTPUT
        if display:
            ds9.mtv(kImageOut, frame=3+foffset)
            ds9.mtv(diffIm2, frame=4+foffset)
        if writefits:
            kImageOut.writeFits('k2.fits')
            diffIm2.writeFits('dB2')

    def testFunctor(self):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up"
            return
        
        #self.applyFunctor(invert=False, foffset=0)
        #self.applyFunctor(invert=True, foffset=4)
        
        #self.applyFunctor(invert=False, foffset=0, xloc=460, yloc=1656)
        #self.applyFunctor(invert=True, foffset=4, xloc=460, yloc=1656)
        
        #self.applyFunctor(invert=False, foffset=0, xloc=288, yloc=952)
        #self.applyFunctor(invert=True, foffset=4, xloc=288, yloc=952)

        self.applyFunctor(invert=False, foffset=0, xloc=460, yloc=1656)
        self.applyFunctor(invert=True, foffset=4, xloc=460, yloc=1656)
        
        #self.applyFunctor(invert=False, foffset=0, xloc=251, yloc=1950)
        #self.applyFunctor(invert=True, foffset=4, xloc=251, yloc=1950)
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
