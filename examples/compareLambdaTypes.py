#!/usr/bin/env python
import os
import pdb
import sys
import unittest
import lsst.utils.tests as tests
import numpy as num
import eups
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.pex.logging as logging

import lsst.afw.display.ds9 as ds9

verbosity = 7
logging.Trace_setVerbosity('lsst.ip.diffim', verbosity)

display = False
writefits = False

# This one compares DeltaFunction and AlardLupton kernels

CFHTTORUN = 'cal-53535-i-797722_1'

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
    def setUp(self, CFHT=True):
        self.policy1     = ipDiffim.createDefaultPolicy()
        self.policy2     = ipDiffim.createDefaultPolicy()
        self.policy3     = ipDiffim.createDefaultPolicy()
        self.policy4     = ipDiffim.createDefaultPolicy()

        self.policy1.set("kernelBasisSet", "delta-function")
        self.policy1.set("useRegularization", False)
        self.policy1.set("fitForBackground", False)
        self.kList1 = ipDiffim.makeKernelBasisList(self.policy1)
        self.bskv1  = ipDiffim.BuildSingleKernelVisitorF(self.kList1, self.policy1)
        
        self.policy2.set("kernelBasisSet", "delta-function")
        self.policy2.set("useRegularization", True)
        self.policy2.set("lambdaType", "absolute")
        self.policy2.set("lambdaValue", 1.0)
        self.policy2.set("regularizationType", "centralDifference")
        self.policy2.set("centralRegularizationStencil", 9)
        self.policy2.set("fitForBackground", False)
        self.kList2 = ipDiffim.makeKernelBasisList(self.policy2)
        self.hMat2  = ipDiffim.makeRegularizationMatrix(self.policy2)
        self.bskv2  = ipDiffim.BuildSingleKernelVisitorF(self.kList2, self.policy2, self.hMat2)

        self.policy3.set("kernelBasisSet", "delta-function")
        self.policy3.set("useRegularization", True)
        self.policy3.set("lambdaType", "minimizeBiasedRisk")
        self.policy3.set("regularizationConditionTolerance", 1.e5)
        self.policy3.set("regularizationType", "centralDifference")
        self.policy3.set("centralRegularizationStencil", 9)
        self.policy3.set("fitForBackground", False)
        self.kList3 = ipDiffim.makeKernelBasisList(self.policy3)
        self.hMat3  = ipDiffim.makeRegularizationMatrix(self.policy3)
        self.bskv3  = ipDiffim.BuildSingleKernelVisitorF(self.kList3, self.policy3, self.hMat3)

        self.policy4.set("kernelBasisSet", "delta-function")
        self.policy4.set("useRegularization", True)
        self.policy4.set("lambdaType", "minimizeUnbiasedRisk")
        self.policy4.set("regularizationType", "centralDifference")
        self.policy4.set("centralRegularizationStencil", 9)
        self.policy4.set("fitForBackground", False)
        self.kList4 = ipDiffim.makeKernelBasisList(self.policy4)
        self.hMat4  = ipDiffim.makeRegularizationMatrix(self.policy4)
        self.bskv4  = ipDiffim.BuildSingleKernelVisitorF(self.kList4, self.policy4, self.hMat4)

        # known input images
        defDataDir = eups.productDir('afwdata')
        if CFHT:
            defSciencePath  = os.path.join(defDataDir, 'CFHT', 'D4', CFHTTORUN)
            defTemplatePath = os.path.join(defDataDir, 'CFHT', 'D4', CFHTTORUN+'_tmpl')

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
                                                                self.policy1.getPolicy("warpingPolicy"))

        diffimTools.backgroundSubtract(self.policy1.getPolicy("afwBackgroundPolicy"),
                                       [self.scienceImage.getMaskedImage(),
                                        self.templateImage.getMaskedImage()])

        # image statistics
        self.dStats  = ipDiffim.ImageStatisticsF()

        #
        tmi = self.templateImage.getMaskedImage()
        smi = self.scienceImage.getMaskedImage()

        self.policy1.set("detThreshold", 100.)
        kcDetect = ipDiffim.KernelCandidateDetectionF(self.policy1.getPolicy("detectionPolicy"))
        kcDetect.apply(tmi, smi)
        self.footprints = kcDetect.getFootprints()

        
    def tearDown(self):
        del self.policy1
        del self.policy2
        del self.policy3
        del self.policy4
        del self.kList1
        del self.kList2
        del self.kList3
        del self.kList4
        del self.hMat2
        del self.bskv1
        del self.bskv2
        del self.bskv3
        del self.bskv4
        del self.scienceImage
        del self.templateImage

    def apply(self, policy, visitor, xloc, yloc, tmi, smi):
        kc     = ipDiffim.makeKernelCandidate(xloc, yloc, tmi, smi, policy)
        visitor.processCandidate(kc)
        kim    = kc.getKernelImage(ipDiffim.KernelCandidateF.RECENT)
        diffIm = kc.getDifferenceImage(ipDiffim.KernelCandidateF.RECENT)
        kSum   = kc.getKsum(ipDiffim.KernelCandidateF.RECENT)
        bg     = kc.getBackground(ipDiffim.KernelCandidateF.RECENT)

        bbox = kc.getKernel(ipDiffim.KernelCandidateF.RECENT).shrinkBBox(diffim.getBBox(afwImage.LOCAL))
        diffIm = afwImage.MaskedImageF(diffIm, bbox, afwImage.LOCAL)
        self.dStats.apply(diffIm)
        
        dmean = afwMath.makeStatistics(diffIm.getImage(),    afwMath.MEAN).getValue()
        dstd  = afwMath.makeStatistics(diffIm.getImage(),    afwMath.STDEV).getValue()
        vmean = afwMath.makeStatistics(diffIm.getVariance(), afwMath.MEAN).getValue()
        return kSum, bg, dmean, dstd, vmean, kim, diffIm, kc
        
    def applyVisitor(self, invert=False, xloc=397, yloc=580):
        print '# %.2f %.2f' % (xloc, yloc)
        #if xloc != 1312 and yloc != 160:
        #    return
        
        imsize = int(3 * self.policy1.get("kernelSize"))

        # chop out a region around a known object
        bbox = afwGeom.Box2I(afwGeom.Point2I(xloc - imsize/2,
                                             yloc - imsize/2),
                             afwGeom.Point2I(xloc + imsize/2,
                                             yloc + imsize/2) )

        # sometimes the box goes off the image; no big deal...
        try:
            if invert:
                tmi  = afwImage.MaskedImageF(self.scienceImage.getMaskedImage(), bbox, afwImage.LOCAL)
                smi  = afwImage.MaskedImageF(self.templateImage.getMaskedImage(), bbox, afwImage.LOCAL)
            else:
                smi  = afwImage.MaskedImageF(self.scienceImage.getMaskedImage(), bbox, afwImage.LOCAL)
                tmi  = afwImage.MaskedImageF(self.templateImage.getMaskedImage(), bbox, afwImage.LOCAL)
        except Exception, e:
            return None

        # delta function kernel
        logging.Trace("lsst.ip.diffim.compareLambdaTypes", 1, 'DF run')
        results1 = self.apply(self.policy1, self.bskv1, xloc, yloc, tmi, smi)
        kSum1, bg1, dmean1, dstd1, vmean1, kImageOut1, diffIm1, kc1 = results1
        kc1.getKernelSolution(ipDiffim.KernelCandidateF.RECENT).getConditionNumber(ipDiffim.KernelSolution.EIGENVALUE)
        kc1.getKernelSolution(ipDiffim.KernelCandidateF.RECENT).getConditionNumber(ipDiffim.KernelSolution.SVD)
        res = 'DF residuals : %.3f +/- %.3f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(),
                                                                             self.dStats.getRms(),
                                                                             kSum1, bg1,
                                                                             dmean1, dstd1, vmean1)
        logging.Trace("lsst.ip.diffim.compareLambdaTypes", 1, res)

        if display:
            ds9.mtv(tmi, frame=1) # ds9 switches frame 0 and 1 for some reason
            ds9.mtv(smi, frame=0)
            ds9.mtv(kImageOut1, frame=2)
            ds9.mtv(diffIm1, frame=3)
        if writefits:
            tmi.writeFits('t.fits')
            smi.writeFits('s.fits')
            kImageOut1.writeFits('k1.fits')
            diffIm1.writeFits('d1.fits')

        # regularized delta function kernel
        logging.Trace("lsst.ip.diffim.compareLambdaTypes", 1, 'DFr1 run')
        results2 = self.apply(self.policy2, self.bskv2, xloc, yloc, tmi, smi)
        kSum2, bg2, dmean2, dstd2, vmean2, kImageOut2, diffIm2, kc2 = results2
        res = 'DFr1 residuals : %.3f +/- %.3f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(),
                                                                               self.dStats.getRms(),
                                                                               kSum2, bg2,
                                                                               dmean2, dstd2, vmean2)
        logging.Trace("lsst.ip.diffim.compareLambdaTypes", 1, res)
        
        if display:
            ds9.mtv(tmi, frame=4)
            ds9.mtv(smi, frame=5)
            ds9.mtv(kImageOut2, frame=6)
            ds9.mtv(diffIm2, frame=7)
        if writefits:
            kImageOut2.writeFits('k2.fits')
            diffIm2.writeFits('d2')


        # minimize Biased Risk
        for logTol in num.arange(4, 6, 0.1):
            self.policy3.set("regularizationConditionTolerance", 10**logTol)
            logging.Trace("lsst.ip.diffim.compareLambdaTypes", 1,
                         'DFrBr run %.3e' % (self.policy3.get("regularizationConditionTolerance")))
            results3 = self.apply(self.policy3, self.bskv3, xloc, yloc, tmi, smi)
            kSum3, bg3, dmean3, dstd3, vmean3, kImageOut3, diffIm3, kc3 = results3
            res = 'DFrBr %.3e : %.3f +/- %.3f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.policy3.get("regularizationConditionTolerance"),
                                                                               self.dStats.getMean(),
                                                                               self.dStats.getRms(),
                                                                               kSum3, bg3,
                                                                               dmean3, dstd3, vmean3)
            logging.Trace("lsst.ip.diffim.compareLambdaTypes", 1, res)
 
        # outputs
        if display:
            ds9.mtv(tmi, frame=8)
            ds9.mtv(smi, frame=9)
            ds9.mtv(kImageOut3, frame=10)
            ds9.mtv(diffIm3, frame=11)
        if writefits:
            kImageOut3.writeFits('k3.fits')
            diffIm3.writeFits('d3')


        # Minimize Unbised Rick
        logging.Trace("lsst.ip.diffim.compareLambdaTypes", 1, 'DFrUr run')
        results4 = self.apply(self.policy4, self.bskv4, xloc, yloc, tmi, smi)
        kSum4, bg4, dmean4, dstd4, vmean4, kImageOut4, diffIm4, kc4 = results4
        res = 'DFrUr residuals : %.3f +/- %.3f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(),
                                                                                self.dStats.getRms(),
                                                                                kSum4, bg4,
                                                                                dmean4, dstd4, vmean4)
        logging.Trace("lsst.ip.diffim.compareLambdaTypes", 1, res)

        # outputs
        if display:
            ds9.mtv(tmi, frame=12)
            ds9.mtv(smi, frame=13)
            ds9.mtv(kImageOut4, frame=14)
            ds9.mtv(diffIm4, frame=15)
        if writefits:
            kImageOut4.writeFits('k4.fits')
            diffIm4.writeFits('d4')

        #raw_input('Next: ')

    def testFunctor(self):
        for fp in self.footprints:
            # note this returns the kernel images
            self.applyVisitor(invert=False, 
                              xloc= int(0.5 * ( fp.getBBox().getMinX() + fp.getBBox().getMaxX() )),
                              yloc= int(0.5 * ( fp.getBBox().getMinY() + fp.getBBox().getMaxY() )))
       
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

    if len(sys.argv) > 1:
        CFHTTORUN = sys.argv[1]

    run(True)
