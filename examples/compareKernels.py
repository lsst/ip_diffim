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

Verbosity = 5
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

display = True
writefits = False
iterate = False

# This one compares DeltaFunction and AlardLupton kernels

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
        
    def setUp(self):
        self.policy      = pexPolicy.Policy.createPolicy(diffimPolicy)
        self.kCols       = self.policy.getInt('kernelCols')
        self.kRows       = self.policy.getInt('kernelRows')
        self.fpGrowKsize = self.policy.getDouble('fpGrowKsize')

        # Delta function basis set
        self.basisList1  = ipDiffim.generateDeltaFunctionKernelSet(self.kCols, self.kRows)
        self.kFunctor1   = ipDiffim.PsfMatchingFunctorF(self.basisList1)

        # Alard-Lupton basis set
        nGauss   = self.policy.get("alardNGauss")
        sigGauss = self.policy.getDoubleArray("alardSigGauss")
        degGauss = self.policy.getIntArray("alardDegGauss")
        assert len(sigGauss) == nGauss
        assert len(degGauss) == nGauss
        assert self.kCols == self.kRows  # square
        assert self.kCols % 2 == 1  # odd sized
        kHalfWidth = int(self.kCols/2)
        self.basisList2  = ipDiffim.generateAlardLuptonKernelSet(kHalfWidth, nGauss, sigGauss, degGauss)
        self.kFunctor2   = ipDiffim.PsfMatchingFunctorF(self.basisList2)

        # Regularized delta function basis set
        self.H = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, 0, 2, 0)
        self.kFunctor3   = ipDiffim.PsfMatchingFunctorF(self.basisList1, self.H)

        # known input images
        defDataDir = eups.productDir('afwdata')
        defSciencePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                      "v26-e0-c011-a00.sci")
        defTemplatePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                       "v5-e0-c011-a00.sci")
        self.scienceImage   = afwImage.ExposureF(defSciencePath)
        self.templateImage  = afwImage.ExposureF(defTemplatePath)
        self.templateImage  = ipDiffim.warpTemplateExposure(self.templateImage, self.scienceImage, self.policy)

        # image statistics
        self.dStats  = ipDiffim.ImageStatisticsF()

        # footprints
        self.detSet     = afwDetection.makeDetectionSet(self.scienceImage.getMaskedImage(), afwDetection.Threshold(575))
        self.footprints = self.detSet.getFootprints()
        # BLOCKED BY TICKET 911
        #self.footprints  = ipDiffim.getCollectionOfFootprintsForPsfMatching(self.templateImage.getMaskedImage(),
        #                                                                    self.scienceImage.getMaskedImage(),
        #                                                                    self.policy)
        
    def tearDown(self):
        del self.policy

    def applyFunctor(self, invert=False, xloc=397, yloc=580):
        print '# %.2f %.2f' % (xloc, yloc)
        
        imsize = int(3 * self.kCols)

        # chop out a region around a known object
        bbox = afwImage.BBox( afwImage.PointI(xloc - imsize/2,
                                              yloc - imsize/2),
                              afwImage.PointI(xloc + imsize/2,
                                              yloc + imsize/2) )

        # sometimes the box goes off the image; no big deal...
        try:
            if invert:
                tmi  = afwImage.MaskedImageF(self.scienceImage.getMaskedImage(),  bbox)
                smi  = afwImage.MaskedImageF(self.templateImage.getMaskedImage(), bbox)
            else:
                smi  = afwImage.MaskedImageF(self.scienceImage.getMaskedImage(),  bbox)
                tmi  = afwImage.MaskedImageF(self.templateImage.getMaskedImage(), bbox)
        except:
            return None

        # estimate of the variance
        var  = afwImage.MaskedImageF(smi, True)
        var -= tmi



        # delta function kernel
        for func in (self.kFunctor1.apply,):
            func(tmi.getImage(), smi.getImage(), var.getVariance(), self.policy)
            kernel    = self.kFunctor1.getKernel()
            kImageOut = afwImage.ImageD(self.kCols, self.kRows)
            kSum      = kernel.computeImage(kImageOut, False)
            diffIm    = ipDiffim.convolveAndSubtract(tmi, smi, kernel, self.kFunctor1.getBackground())
            bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                      kernel.getCtrY()) ,
                                      afwImage.PointI(diffIm.getWidth() - (kernel.getWidth()  - kernel.getCtrX()),
                                                      diffIm.getHeight() - (kernel.getHeight() - kernel.getCtrY())))
            diffIm2   = afwImage.MaskedImageF(diffIm, bbox)
            self.dStats.apply( diffIm2 )
    
            dmean1 = afwMath.makeStatistics(diffIm2.getImage(),    afwMath.MEAN).getValue()
            dstd1  = afwMath.makeStatistics(diffIm2.getImage(),    afwMath.STDEV).getValue()
            vmean1 = afwMath.makeStatistics(diffIm2.getVariance(), afwMath.MEAN).getValue()
            
            print 'DF Diffim residuals : %.2f +/- %.2f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(), self.dStats.getRms(),
                                                                                        kSum, self.kFunctor1.getBackground(),
                                                                                        dmean1, dstd1, vmean1)
        # outputs
        if display:
            ds9.mtv(tmi, frame=0)
            ds9.mtv(smi, frame=1)
            ds9.mtv(kImageOut, frame=2)
            ds9.mtv(diffIm2, frame=3)
        if writefits:
            tmi.writeFits('t')
            smi.writeFits('s')
            diffIm2.writeFits('d1')
            kImageOut.writeFits('k1.fits')

        # alard-lupton kernel
        for func in (self.kFunctor2.apply,):
            func(tmi.getImage(), smi.getImage(), var.getVariance(), self.policy)
            kernel    = self.kFunctor2.getKernel()
            kImageOut = afwImage.ImageD(self.kCols, self.kRows)
            kSum      = kernel.computeImage(kImageOut, False)
            diffIm    = ipDiffim.convolveAndSubtract(tmi, smi, kernel, self.kFunctor2.getBackground())
            bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                      kernel.getCtrY()) ,
                                      afwImage.PointI(diffIm.getWidth() - (kernel.getWidth()  - kernel.getCtrX()),
                                                      diffIm.getHeight() - (kernel.getHeight() - kernel.getCtrY())))
            diffIm2   = afwImage.MaskedImageF(diffIm, bbox)
            self.dStats.apply( diffIm2 )
    
            dmean2 = afwMath.makeStatistics(diffIm2.getImage(),    afwMath.MEAN).getValue()
            dstd2  = afwMath.makeStatistics(diffIm2.getImage(),    afwMath.STDEV).getValue()
            vmean2 = afwMath.makeStatistics(diffIm2.getVariance(), afwMath.MEAN).getValue()
            
            print 'AL Diffim residuals : %.2f +/- %.2f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(), self.dStats.getRms(),
                                                                                        kSum, self.kFunctor2.getBackground(),
                                                                                        dmean2, dstd2, vmean2)
            diffImAL = diffIm2
            
        # outputs
        if display:
            ds9.mtv(tmi, frame=4)
            ds9.mtv(smi, frame=5)
            ds9.mtv(kImageOut, frame=6)
            ds9.mtv(diffIm2, frame=7)
        if writefits:
            diffIm2.writeFits('d2')
            kImageOut.writeFits('k2.fits')

        # regularized delta function kernel
        for func in (self.kFunctor3.apply,):
            func(tmi.getImage(), smi.getImage(), var.getVariance(), self.policy)
            kernel    = self.kFunctor3.getKernel()
            kImageOut = afwImage.ImageD(self.kCols, self.kRows)
            kSum      = kernel.computeImage(kImageOut, False)
            diffIm    = ipDiffim.convolveAndSubtract(tmi, smi, kernel, self.kFunctor3.getBackground())
            bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                      kernel.getCtrY()) ,
                                      afwImage.PointI(diffIm.getWidth() - (kernel.getWidth()  - kernel.getCtrX()),
                                                      diffIm.getHeight() - (kernel.getHeight() - kernel.getCtrY())))
            diffIm2   = afwImage.MaskedImageF(diffIm, bbox)
            self.dStats.apply( diffIm2 )
    
            dmean1 = afwMath.makeStatistics(diffIm2.getImage(),    afwMath.MEAN).getValue()
            dstd1  = afwMath.makeStatistics(diffIm2.getImage(),    afwMath.STDEV).getValue()
            vmean1 = afwMath.makeStatistics(diffIm2.getVariance(), afwMath.MEAN).getValue()
            
            print 'DFr Diffim residuals : %.2f +/- %.2f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(), self.dStats.getRms(),
                                                                                         kSum, self.kFunctor3.getBackground(),
                                                                                         dmean1, dstd1, vmean1)
            diffImDRr = diffIm2

        diffImAL -= diffImDRr
        self.dStats.apply(diffImAL)
        dmeanD = afwMath.makeStatistics(diffImAL.getImage(),    afwMath.MEAN).getValue()
        dstdD  = afwMath.makeStatistics(diffImAL.getImage(),    afwMath.STDEV).getValue()
        print 'AL-DFr Diffim residuals : %.2f +/- %.2f' % (self.dStats.getMean(), self.dStats.getRms())
        print
        
        # outputs
        if display:
            ds9.mtv(tmi, frame=8)
            ds9.mtv(smi, frame=9)
            ds9.mtv(kImageOut, frame=10)
            ds9.mtv(diffIm2, frame=11)
        if writefits:
            diffIm2.writeFits('d3')
            kImageOut.writeFits('k3.fits')


    def testFunctor(self):
        for object in self.footprints:
            # note this returns the kernel images
            self.applyFunctor(invert=False, 
                              xloc= int(0.5 * ( object.getBBox().getX0() + object.getBBox().getX1() )),
                              yloc= int(0.5 * ( object.getBBox().getY0() + object.getBBox().getY1() )))
       
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
    if '-i' in sys.argv:
        iterate = True
        
    run(True)
