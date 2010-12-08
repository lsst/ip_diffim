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
        self.policy      = ipDiffim.createDefaultPolicy(diffimPolicy)
        self.kSize       = self.policy.getInt('kernelSize')

        # Delta function basis set
        self.basisList1  = ipDiffim.makeDeltaFunctionBasisSet(self.kSize, self.kSize)
        self.kFunctor1   = ipDiffim.PsfMatchingFunctorF(self.basisList1)

        # Alard-Lupton basis set
        nGauss   = self.policy.get("alardNGauss")
        sigGauss = self.policy.getDoubleArray("alardSigGauss")
        degGauss = self.policy.getIntArray("alardDegGauss")
        assert len(sigGauss) == nGauss
        assert len(degGauss) == nGauss
        assert self.kSize % 2 == 1  # odd sized
        kHalfWidth = iself.kSize // 2
        self.basisList2  = ipDiffim.makeAlardLuptonBasisSet(kHalfWidth, nGauss, sigGauss, degGauss)
        self.kFunctor2   = ipDiffim.PsfMatchingFunctorF(self.basisList2)

        # Regularized delta function basis set using default forward diff
        self.policy.set("regularizationType", "forwardDifference")
        h3 = ipDiffim.makeRegularizationMatrix(self.policy)
        self.kFunctor3   = ipDiffim.PsfMatchingFunctorF(self.basisList1, h3)

        # Regularized delta function basis set using default central diff
        self.policy.set("regularizationType", "centralDifference")
        h4 = ipDiffim.makeRegularizationMatrix(self.policy)
        self.kFunctor4   = ipDiffim.PsfMatchingFunctorF(self.basisList1, h4)

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

        self.policy.set("detThreshold", 100.)
        kcDetect = ipDiffim.KernelCandidateDetectionF(self.policy)
        kcDetect.apply(tmi, smi)
        self.footprints = kcDetect.getFootprints()
        
    def tearDown(self):
        del self.policy

    def apply(self, functor, tmi, smi, var):
        functor.apply(tmi.getImage(), smi.getImage(), var.getVariance(), self.policy)
        pair      = functor.getSolution()
        kernel    = pair.first
        bg        = pair.second
        kImageOut = afwImage.ImageD(self.kSize, self.kSize)
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
        if invert:
            frame0 = 16
        else:
            frame0 = 0
            
        print '# %.2f %.2f %s' % (xloc, yloc, invert)
        
        imsize = int(3 * self.kSize)

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
            ds9.mtv(tmi, frame=frame0+0)
            ds9.mtv(smi, frame=frame0+1)
            ds9.mtv(kImageOut1, frame=frame0+2)
            ds9.mtv(diffIm1, frame=frame0+3)
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
            ds9.mtv(tmi, frame=frame0+4)
            ds9.mtv(smi, frame=frame0+5)
            ds9.mtv(kImageOut2, frame=frame0+6)
            ds9.mtv(diffIm2, frame=frame0+7)
        if writefits:
            kImageOut2.writeFits('k2.fits')
            diffIm2.writeFits('d2')

        # regularized delta function kernel 1
        kSum3, bg3, dmean3, dstd3, vmean3, kImageOut3, diffIm3 = self.apply(self.kFunctor3, tmi, smi, var)
        print 'DFr Diffim residuals : %.2f +/- %.2f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(),
                                                                                     self.dStats.getRms(),
                                                                                     kSum3, bg3,
                                                                                     dmean3, dstd3, vmean3)
        # outputs
        if display:
            ds9.mtv(tmi, frame=frame0+8)
            ds9.mtv(smi, frame=frame0+9)
            ds9.mtv(kImageOut3, frame=frame0+10)
            ds9.mtv(diffIm3, frame=frame0+11)
        if writefits:
            kImageOut3.writeFits('k3.fits')
            diffIm3.writeFits('d3')

        # regularized delta function kernel 2
        kSum4, bg4, dmean4, dstd4, vmean4, kImageOut4, diffIm4 = self.apply(self.kFunctor4, tmi, smi, var)
        print 'DFr Diffim residuals : %.2f +/- %.2f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(),
                                                                                     self.dStats.getRms(),
                                                                                     kSum4, bg4,
                                                                                     dmean4, dstd4, vmean4)
        # outputs
        if display:
            ds9.mtv(tmi, frame=frame0+12)
            ds9.mtv(smi, frame=frame0+13)
            ds9.mtv(kImageOut4, frame=frame0+14)
            ds9.mtv(diffIm4, frame=frame0+15)
        if writefits:
            kImageOut4.writeFits('k4.fits')
            diffIm4.writeFits('d4')

    def testFunctor(self):
        for fp in self.footprints:
            # note this returns the kernel images
            self.applyFunctor(invert=False, 
                              xloc= int(0.5 * ( fp.getBBox().getX0() + fp.getBBox().getX1() )),
                              yloc= int(0.5 * ( fp.getBBox().getY0() + fp.getBBox().getY1() )))
            self.applyFunctor(invert=True, 
                              xloc= int(0.5 * ( fp.getBBox().getX0() + fp.getBBox().getX1() )),
                              yloc= int(0.5 * ( fp.getBBox().getY0() + fp.getBBox().getY1() )))
            raw_input('Next: ')

       
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
