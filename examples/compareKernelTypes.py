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
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.pex.logging as logging

import lsst.afw.display.ds9 as ds9

verbosity = 5
logging.Trace_setVerbosity('lsst.ip.diffim', verbosity)

display = True
writefits = False

# This one compares DeltaFunction and AlardLupton kernels

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
    def setUp(self, CFHT=False):
        self.policy1     = ipDiffim.createDefaultPolicy()
        self.policy2     = ipDiffim.createDefaultPolicy()
        self.policy3     = ipDiffim.createDefaultPolicy()

        self.policy1.set("kernelBasisSet", "delta-function")
        self.policy1.set("useRegularization", False)
        self.policy1.set("maxConditionNumber", 5.0e6)
        self.policy1.set("checkConditionNumber", False)
        self.policy1.set('fitForBackground', False)
        self.policy1.set('constantVarianceWeighting', True)
        self.kList1 = ipDiffim.makeKernelBasisList(self.policy1)
        self.bskv1  = ipDiffim.BuildSingleKernelVisitorF(self.kList1, self.policy1)
        
        self.policy2.set("kernelBasisSet", "delta-function")
        self.policy2.set("useRegularization", True)
        self.policy2.set("maxConditionNumber", 5.0e6)
        self.policy2.set("checkConditionNumber", False)
        self.policy2.set('fitForBackground', False)
        self.policy2.set('lambdaValue', 1000.0)
        self.policy2.set('constantVarianceWeighting', True)
        self.kList2 = ipDiffim.makeKernelBasisList(self.policy2)
        self.hMat2  = ipDiffim.makeRegularizationMatrix(self.policy2)
        self.bskv2  = ipDiffim.BuildSingleKernelVisitorF(self.kList2, self.policy2, self.hMat2)
        
        self.policy3.set("kernelBasisSet", "alard-lupton")
        self.policy3.set("maxConditionNumber", 5.0e7)
        self.policy3.set("checkConditionNumber", False)
        self.policy3.set('fitForBackground', False)
        self.policy3.set('constantVarianceWeighting', True)
        self.kList3 = ipDiffim.makeKernelBasisList(self.policy3)
        self.bskv3  = ipDiffim.BuildSingleKernelVisitorF(self.kList3, self.policy3)

        # known input images
        defDataDir = eups.productDir('afwdata')
        if CFHT:
            defSciencePath  = os.path.join(defDataDir, 'CFHT', 'D4', 'cal-53535-i-797722_1')
            defTemplatePath = os.path.join(defDataDir, 'CFHT', 'D4', 'cal-53535-i-797722_1_tmpl')

            # no need to remap
            self.scienceImage   = afwImage.ExposureF(defSciencePath)
            self.templateImage  = afwImage.ExposureF(defTemplatePath)
        elif True:
            #self.scienceImage  = afwImage.ExposureF('s2L006430_0106g4TANSIPwInv.fits')
            #self.templateImage = afwImage.ExposureF('s2L200006_00770078g4.fits')
            self.scienceImage  = afwImage.ExposureF('s2L007173_0100g4TANSIPwInv.fits')
            self.templateImage = afwImage.ExposureF('oneTemplate100006_0072g4.fits')
            diffimTools.backgroundSubtract(self.policy1.getPolicy("afwBackgroundPolicy"),
                                           [self.scienceImage.getMaskedImage(),])
            self.templateImage = ipDiffim.warpTemplateExposure(self.templateImage,
                                                               self.scienceImage,
                                                               self.policy1.getPolicy("warpingPolicy"))
            ### reverse order!
            #foo = self.templateImage
            #self.templateImage = self.scienceImage
            #self.scienceImage = foo
            
        else:
            defSciencePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                          "v26-e0-c011-a00.sci")
            defTemplatePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                           "v5-e0-c011-a00.sci")
            
            self.scienceImage   = afwImage.ExposureF(defSciencePath)
            self.templateImage  = afwImage.ExposureF(defTemplatePath)
            self.templateImage  = ipDiffim.warpTemplateExposure(self.templateImage,
                                                                self.scienceImage,
                                                                self.policy1.getPolicy("warpingPolicy"))


        # image statistics
        self.dStats  = ipDiffim.ImageStatisticsF()

        #
        tmi = self.templateImage.getMaskedImage()
        smi = self.scienceImage.getMaskedImage()

        detPolicy = self.policy1.getPolicy("detectionPolicy")
        detPolicy.set("detThreshold", 100.)
        detPolicy.set("detOnTemplate", False)
        kcDetect = ipDiffim.KernelCandidateDetectionF(detPolicy)
        kcDetect.apply(tmi, smi)
        self.footprints = kcDetect.getFootprints()

        
    def tearDown(self):
        del self.policy1
        del self.policy2
        del self.policy3
        del self.kList1
        del self.kList2
        del self.kList3
        del self.hMat2
        del self.bskv1
        del self.bskv2
        del self.bskv3
        del self.scienceImage
        del self.templateImage

    def apply(self, policy, visitor, xloc, yloc, tmi, smi):
        kc     = ipDiffim.makeKernelCandidate(xloc, yloc, tmi, smi, policy)
        visitor.processCandidate(kc)
        kim    = kc.getKernelImage(ipDiffim.KernelCandidateF.RECENT)
        diffIm = kc.getDifferenceImage(ipDiffim.KernelCandidateF.RECENT)
        kSum   = kc.getKsum(ipDiffim.KernelCandidateF.RECENT)
        bg     = kc.getBackground(ipDiffim.KernelCandidateF.RECENT)

        p0, p1 = diffimTools.getConvolvedImageLimits(kc.getKernel(ipDiffim.KernelCandidateF.RECENT), diffIm)
        bbox   = afwImage.BBox(p0, p1)
        diffIm = afwImage.MaskedImageF(diffIm, bbox)
        self.dStats.apply(diffIm)
        
        dmean = afwMath.makeStatistics(diffIm.getImage(),    afwMath.MEAN).getValue()
        dstd  = afwMath.makeStatistics(diffIm.getImage(),    afwMath.STDEV).getValue()
        vmean = afwMath.makeStatistics(diffIm.getVariance(), afwMath.MEAN).getValue()
        return kSum, bg, dmean, dstd, vmean, kim, diffIm, kc
        
    def applyVisitor(self, invert=False, xloc=397, yloc=580):
        print '# %.2f %.2f' % (xloc, yloc)

        #if (int(xloc) != 1149) and (int(yloc) != 179):
        #    return
            
        imsize = int(3 * self.policy1.get("kernelSize"))

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

        # delta function kernel
        results1 = self.apply(self.policy1, self.bskv1, xloc, yloc, tmi, smi)
        kSum1, bg1, dmean1, dstd1, vmean1, kImageOut1, diffIm1, kc1 = results1
        kc1.getKernelSolution(ipDiffim.KernelCandidateF.RECENT).getConditionNumber(ipDiffim.KernelSolution.EIGENVALUE)
        kc1.getKernelSolution(ipDiffim.KernelCandidateF.RECENT).getConditionNumber(ipDiffim.KernelSolution.SVD)
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

        # regularized delta function kernel
        results2 = self.apply(self.policy2, self.bskv2, xloc, yloc, tmi, smi)
        kSum2, bg2, dmean2, dstd2, vmean2, kImageOut2, diffIm2, kc2 = results2
        kc2.getKernelSolution(ipDiffim.KernelCandidateF.RECENT).getConditionNumber(ipDiffim.KernelSolution.EIGENVALUE)
        kc2.getKernelSolution(ipDiffim.KernelCandidateF.RECENT).getConditionNumber(ipDiffim.KernelSolution.SVD)
        print 'DFr Diffim residuals : %.2f +/- %.2f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(),
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

        # alard-lupton kernel
        results3 = self.apply(self.policy3, self.bskv3, xloc, yloc, tmi, smi)
        kSum3, bg3, dmean3, dstd3, vmean3, kImageOut3, diffIm3, kc3 = results3
        kc3.getKernelSolution(ipDiffim.KernelCandidateF.RECENT).getConditionNumber(ipDiffim.KernelSolution.EIGENVALUE)
        kc3.getKernelSolution(ipDiffim.KernelCandidateF.RECENT).getConditionNumber(ipDiffim.KernelSolution.SVD)
        print 'AL Diffim residuals : %.2f +/- %.2f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(),
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
            self.applyVisitor(invert=False, 
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
