#!/usr/bin/env python
from __future__ import print_function
from builtins import input
import os
import sys
import unittest

import lsst.utils.tests as tests
import lsst.utils
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.pex.logging as logging
import lsst.pex.config as pexConfig

import lsst.afw.display.ds9 as ds9

verbosity = 7
logging.Trace_setVerbosity('lsst.ip.diffim', verbosity)

display = True
writefits = False

# This one compares DeltaFunction kernels of different types; iterate lambdaVal for different strengths

CFHTTORUN = 'cal-53535-i-797722_1'


class DiffimTestCases(unittest.TestCase):

    # D = I - (K.x.T + bg)
    def setUp(self, CFHT=True):
        lambdaValue = 1.0

        self.config1 = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.config1.kernel.name = "DF"
        self.subconfig1 = self.config1.kernel.active

        self.config2 = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.config2.kernel.name = "DF"
        self.subconfig2 = self.config2.kernel.active

        self.config3 = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.config3.kernel.name = "DF"
        self.subconfig3 = self.config3.kernel.active

        self.config4 = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.config4.kernel.name = "DF"
        self.subconfig4 = self.config4.kernel.active

        self.subconfig1.useRegularization = False

        self.subconfig2.useRegularization = True
        self.subconfig2.lambdaType = "absolute"
        self.subconfig2.lambdaValue = lambdaValue
        self.subconfig2.regularizationType = "centralDifference"
        self.subconfig2.centralRegularizationStencil = 5

        self.subconfig3.useRegularization = True
        self.subconfig3.lambdaType = "absolute"
        self.subconfig3.lambdaValue = lambdaValue
        self.subconfig3.regularizationType = "centralDifference"
        self.subconfig3.centralRegularizationStencil = 9

        self.subconfig4.useRegularization = True
        self.subconfig4.lambdaType = "absolute"
        self.subconfig4.lambdaValue = lambdaValue
        self.subconfig4.regularizationType = "forwardDifference"
        self.subconfig4.forwardRegularizationOrders = [1, 2]

        self.kList1 = ipDiffim.makeKernelBasisList(self.subconfig1)
        self.bskv1 = ipDiffim.BuildSingleKernelVisitorF(self.kList1, pexConfig.makePolicy(self.subconfig1))

        self.kList2 = ipDiffim.makeKernelBasisList(self.subconfig2)
        self.hMat2 = ipDiffim.makeRegularizationMatrix(pexConfig.makePolicy(self.subconfig2))
        self.bskv2 = ipDiffim.BuildSingleKernelVisitorF(self.kList2, 
                                                        pexConfig.makePolicy(self.subconfig2), self.hMat2)

        self.kList3 = ipDiffim.makeKernelBasisList(self.subconfig3)
        self.hMat3 = ipDiffim.makeRegularizationMatrix(pexConfig.makePolicy(self.subconfig3))
        self.bskv3 = ipDiffim.BuildSingleKernelVisitorF(self.kList3, 
                                                        pexConfig.makePolicy(self.subconfig3), self.hMat3)

        self.kList4 = ipDiffim.makeKernelBasisList(self.subconfig4)
        self.hMat4 = ipDiffim.makeRegularizationMatrix(pexConfig.makePolicy(self.subconfig4))
        self.bskv4 = ipDiffim.BuildSingleKernelVisitorF(self.kList4, 
                                                        pexConfig.makePolicy(self.subconfig4), self.hMat4)

        # known input images
        defDataDir = lsst.utils.getPackageDir('afwdata')
        if CFHT:
            defSciencePath = os.path.join(defDataDir, 'CFHT', 'D4', CFHTTORUN+'.fits')
            defTemplatePath = os.path.join(defDataDir, 'CFHT', 'D4', CFHTTORUN+'_tmpl.fits')

            # no need to remap
            self.scienceExposure = afwImage.ExposureF(defSciencePath)
            self.templateExposure = afwImage.ExposureF(defTemplatePath)
        else:
            defSciencePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                          "v26-e0-c011-a00.sci")
            defTemplatePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                           "v5-e0-c011-a00.sci")

            self.scienceExposure = afwImage.ExposureF(defSciencePath)
            self.templateExposure = afwImage.ExposureF(defTemplatePath)
            warper = afwMath.Warper.fromConfig(self.subconfig1.warpingConfig)
            self.templateExposure = warper.warpExposure(self.scienceExposure.getWcs(), self.templateExposure,
                                                        destBBox=self.scienceExposure.getBBox())

        diffimTools.backgroundSubtract(self.subconfig1.afwBackgroundConfig,
                                       [self.scienceExposure.getMaskedImage(),
                                        self.templateExposure.getMaskedImage()])

        # image statistics
        self.dStats = ipDiffim.ImageStatisticsF()

        #
        tmi = self.templateExposure.getMaskedImage()
        smi = self.scienceExposure.getMaskedImage()

        detConfig = self.subconfig1.detectionConfig
        detPolicy = pexConfig.makePolicy(detConfig)
        detPolicy.set("detThreshold", 50.)
        detPolicy.set("detOnTemplate", False)
        kcDetect = ipDiffim.KernelCandidateDetectionF(detPolicy)
        kcDetect.apply(tmi, smi)
        self.footprints = kcDetect.getFootprints()

    def tearDown(self):
        del self.subconfig1
        del self.subconfig2
        del self.subconfig3
        del self.subconfig4
        del self.kList1
        del self.kList2
        del self.kList3
        del self.kList4
        del self.hMat2
        del self.bskv1
        del self.bskv2
        del self.bskv3
        del self.bskv4
        del self.scienceExposure
        del self.templateExposure

    def apply(self, policy, visitor, xloc, yloc, tmi, smi):
        kc = ipDiffim.makeKernelCandidate(xloc, yloc, tmi, smi, policy)
        visitor.processCandidate(kc)
        kim = kc.getKernelImage(ipDiffim.KernelCandidateF.RECENT)
        diffIm = kc.getDifferenceImage(ipDiffim.KernelCandidateF.RECENT)
        kSum = kc.getKsum(ipDiffim.KernelCandidateF.RECENT)
        bg = kc.getBackground(ipDiffim.KernelCandidateF.RECENT)

        bbox = kc.getKernel(ipDiffim.KernelCandidateF.RECENT).shrinkBBox(diffIm.getBBox(afwImage.LOCAL))
        diffIm = afwImage.MaskedImageF(diffIm, bbox, afwImage.LOCAL)
        self.dStats.apply(diffIm)

        dmean = afwMath.makeStatistics(diffIm.getImage(), afwMath.MEAN).getValue()
        dstd = afwMath.makeStatistics(diffIm.getImage(), afwMath.STDEV).getValue()
        vmean = afwMath.makeStatistics(diffIm.getVariance(), afwMath.MEAN).getValue()
        return kSum, bg, dmean, dstd, vmean, kim, diffIm, kc

    def applyVisitor(self, invert=False, xloc=397, yloc=580):
        print('# %.2f %.2f' % (xloc, yloc))
        #if xloc != 1312 and yloc != 160:
        #    return

        imsize = int(3 * self.subconfig1.kernelSize)

        # chop out a region around a known object
        bbox = afwGeom.Box2I(afwGeom.Point2I(xloc - imsize//2,
                                             yloc - imsize//2),
                             afwGeom.Point2I(xloc + imsize//2,
                                             yloc + imsize//2))

        # sometimes the box goes off the image; no big deal...
        try:
            if invert:
                tmi = afwImage.MaskedImageF(self.scienceExposure.getMaskedImage(), bbox, afwImage.LOCAL)
                smi = afwImage.MaskedImageF(self.templateExposure.getMaskedImage(), bbox, afwImage.LOCAL)
            else:
                smi = afwImage.MaskedImageF(self.scienceExposure.getMaskedImage(), bbox, afwImage.LOCAL)
                tmi = afwImage.MaskedImageF(self.templateExposure.getMaskedImage(), bbox, afwImage.LOCAL)
        except Exception:
            return None

        # delta function kernel
        logging.Trace("lsst.ip.diffim.compareLambdaTypes", 1, 'DF run')
        results1 = self.apply(pexConfig.makePolicy(self.subconfig1), self.bskv1, xloc, yloc, tmi, smi)
        kSum1, bg1, dmean1, dstd1, vmean1, kImageOut1, diffIm1, kc1 = results1
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
        logging.Trace("lsst.ip.diffim.compareLambdaTypes", 1, 'DFrC5 run')
        results2 = self.apply(pexConfig.makePolicy(self.subconfig2), self.bskv2, xloc, yloc, tmi, smi)
        kSum2, bg2, dmean2, dstd2, vmean2, kImageOut2, diffIm2, kc2 = results2
        res = 'DFrC5 residuals : %.3f +/- %.3f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(),
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

        # regularized delta function kernel
        logging.Trace("lsst.ip.diffim.compareLambdaTypes", 1, 'DFrC9 run')
        results3 = self.apply(pexConfig.makePolicy(self.subconfig3), self.bskv3, xloc, yloc, tmi, smi)
        kSum3, bg3, dmean3, dstd3, vmean3, kImageOut3, diffIm3, kc3 = results3
        res = 'DFrC9 residuals : %.3f +/- %.3f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(),
                                                                                self.dStats.getRms(),
                                                                                kSum3, bg3,
                                                                                dmean3, dstd3, vmean3)
        logging.Trace("lsst.ip.diffim.compareLambdaTypes", 1, res)
        if display:
            ds9.mtv(tmi, frame=8)
            ds9.mtv(smi, frame=9)
            ds9.mtv(kImageOut3, frame=10)
            ds9.mtv(diffIm3, frame=11)
        if writefits:
            kImageOut2.writeFits('k3.fits')
            diffIm2.writeFits('d3')

        # regularized delta function kernel
        logging.Trace("lsst.ip.diffim.compareLambdaTypes", 1, 'DFrF12 run')
        results4 = self.apply(pexConfig.makePolicy(self.subconfig4), self.bskv4, xloc, yloc, tmi, smi)
        kSum4, bg4, dmean4, dstd4, vmean4, kImageOut4, diffIm4, kc4 = results4
        res = 'DFrF12 residuals : %.3f +/- %.3f; %.2f, %.2f; %.2f %.2f, %.2f' % (self.dStats.getMean(),
                                                                                 self.dStats.getRms(),
                                                                                 kSum4, bg4,
                                                                                 dmean4, dstd4, vmean4)
        logging.Trace("lsst.ip.diffim.compareLambdaTypes", 1, res)
        if display:
            ds9.mtv(tmi, frame=12)
            ds9.mtv(smi, frame=13)
            ds9.mtv(kImageOut4, frame=14)
            ds9.mtv(diffIm4, frame=15)
        if writefits:
            kImageOut2.writeFits('k4.fits')
            diffIm2.writeFits('d4')

        input('Next: ')

    def testFunctor(self):
        for fp in self.footprints:
            # note this returns the kernel images
            self.applyVisitor(invert=False,
                              xloc=int(0.5 * (fp.getBBox().getMinX() + fp.getBBox().getMaxX())),
                              yloc=int(0.5 * (fp.getBBox().getMinY() + fp.getBBox().getMaxY())))

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
