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

# This one compares DeltaFunction and AlardLupton kernels

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
    def setUp(self, CFHT=True):
        self.policy      = ipDiffim.generateDefaultPolicy(diffimPolicy)
        self.kCols       = self.policy.getInt('kernelCols')
        self.kRows       = self.policy.getInt('kernelRows')

        # Delta function basis set
        self.basisList1  = ipDiffim.generateDeltaFunctionBasisSet(self.kCols, self.kRows)
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
        self.basisList2  = ipDiffim.generateAlardLuptonBasisSet(kHalfWidth, nGauss, sigGauss, degGauss)
        self.kFunctor2   = ipDiffim.PsfMatchingFunctorF(self.basisList2)

        # Regularized delta function basis set
        self.h = ipDiffim.generateFiniteDifferenceRegularization(self.kCols, self.kRows, 2, 1, 0)
        self.kFunctor3   = ipDiffim.PsfMatchingFunctorF(self.basisList1, self.h)

        # known input images
        defDataDir = eups.productDir('afwdata')
        if CFHT:
            defSciencePath  = os.path.join(defDataDir, 'CFHT', 'D4', 'cal-53535-i-797722_1')
            defTemplatePath = os.path.join(defDataDir, 'CFHT', 'D4', 'cal-53535-i-797722_1_tmpl')

            # no need to remap
            self.scienceImage   = afwImage.ExposureF(defSciencePath)
            self.templateImage  = afwImage.ExposureF(defTemplatePath)
        else:
            defSciencePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                          "v26-e0-c011-a00.sci")
            defTemplatePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                           "v5-e0-c011-a00.sci")
            
            self.scienceImage   = afwImage.ExposureF(defSciencePath)
            self.templateImage  = afwImage.ExposureF(defTemplatePath)
            self.templateImage  = ipDiffim.warpTemplateExposure(self.templateImage,
                                                                self.scienceImage,
                                                                self.policy)


        # image statistics
        self.dStats  = ipDiffim.ImageStatisticsF()

        #
        tmi = self.templateImage.getMaskedImage()
        smi = self.scienceImage.getMaskedImage()
        self.footprints = ipDiffim.getCollectionOfFootprintsForPsfMatching(tmi, smi, self.policy)
        
    def tearDown(self):
        del self.policy

    def apply(self, functor, tmi, smi, var):
        functor.apply(tmi.getImage(), smi.getImage(), var.getVariance(), self.policy)
        pair      = functor.getSolution()
        kernel    = pair.first
        bg        = pair.second
        kImageOut = afwImage.ImageD(self.kCols, self.kRows)
        kSum      = kernel.computeImage(kImageOut, False)
        diffIm    = ipDiffim.convolveAndSubtract(tmi, smi, kernel, bg)
        bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                  kernel.getCtrY()) ,
                                  afwImage.PointI(diffIm.getWidth()-(kernel.getWidth()-kernel.getCtrX()),
                                                  diffIm.getHeight()-(kernel.getHeight()-kernel.getCtrY())))
        diffIm2   = afwImage.MaskedImageF(diffIm, bbox)
        self.dStats.apply( diffIm2 )
        
        dmean = afwMath.makeStatistics(diffIm2.getImage(),    afwMath.MEAN).getValue()
        dstd  = afwMath.makeStatistics(diffIm2.getImage(),    afwMath.STDEV).getValue()
        vmean = afwMath.makeStatistics(diffIm2.getVariance(), afwMath.MEAN).getValue()
        return kSum, bg, dmean, dstd, vmean, kImageOut, diffIm2
        
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
        except Exception, e:
            return None

        # estimate of the variance
        var  = afwImage.MaskedImageF(smi, True)
        var -= tmi

        # delta function kernel
        kSum1, bg1, dmean1, dstd1, vmean1, kImageOut1, diffIm1 = self.apply(self.kFunctor1, tmi, smi, var)
        print 'DF Diffim residuals : %.2f +/- %.2f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(),
                                                                                    self.dStats.getRms(),
                                                                                    kSum1, bg1,
                                                                                    dmean1, dstd1, vmean1)
        if display:
            ds9.mtv(tmi, frame=1) # ds9 switches frame 0 and 1 for some reason
            ds9.mtv(smi, frame=0)
            ds9.mtv(kImageOut1, frame=2)
            ds9.mtv(diffIm1, frame=3)
        if writefits:
            tmi.writeFits('t')
            smi.writeFits('s')
            kImageOut1.writeFits('k1.fits')
            diffIm1.writeFits('d1')

        # alard-lupton kernel
        kSum2, bg2, dmean2, dstd2, vmean2, kImageOut2, diffIm2 = self.apply(self.kFunctor2, tmi, smi, var)
        print 'AL Diffim residuals : %.2f +/- %.2f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(),
                                                                                    self.dStats.getRms(),
                                                                                    kSum2, bg2,
                                                                                    dmean2, dstd2, vmean2)
        if display:
            ds9.mtv(tmi, frame=4)
            ds9.mtv(smi, frame=5)
            ds9.mtv(kImageOut2, frame=6)
            ds9.mtv(diffIm2, frame=7)
        if writefits:
            kImageOut2.writeFits('k2.fits')
            diffIm2.writeFits('d2')

        # regularized delta function kernel
        kSum3, bg3, dmean3, dstd3, vmean3, kImageOut3, diffIm3 = self.apply(self.kFunctor3, tmi, smi, var)
        print 'DFr Diffim residuals : %.2f +/- %.2f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(),
                                                                                     self.dStats.getRms(),
                                                                                     kSum3, bg3,
                                                                                     dmean3, dstd3, vmean3)
        # outputs
        if display:
            ds9.mtv(tmi, frame=8)
            ds9.mtv(smi, frame=9)
            ds9.mtv(kImageOut3, frame=10)
            ds9.mtv(diffIm3, frame=11)
        if writefits:
            kImageOut3.writeFits('k3.fits')
            diffIm3.writeFits('d3')

        raw_input('Next: ')

    def testFunctor(self):
        for fp in self.footprints:
            # note this returns the kernel images
            self.applyFunctor(invert=False, 
                              xloc= int(0.5 * ( fp.getBBox().getX0() + fp.getBBox().getX1() )),
                              yloc= int(0.5 * ( fp.getBBox().getY0() + fp.getBBox().getY1() )))
       
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
