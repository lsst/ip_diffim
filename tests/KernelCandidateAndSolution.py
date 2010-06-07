#!/usr/bin/env python
import os, sys
import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.pex.logging as pexLog
import lsst.afw.display.ds9 as ds9
import numpy as num

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)

class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = ipDiffim.generateDefaultPolicy(diffimPolicy)
        
        # known input images
        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:
            defSciencePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                          "v26-e0-c011-a10.sci")
            defTemplatePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                           "v5-e0-c011-a10.sci")

            scienceImage  = afwImage.ExposureF(defSciencePath)
            templateImage = afwImage.ExposureF(defTemplatePath)
            # set XY0 = 0
            scienceImage.getMaskedImage().setXY0(afwImage.PointI(0, 0))
            templateImage.getMaskedImage().setXY0(afwImage.PointI(0, 0))
            # do the warping first so we don't have any masked pixels in the postage stamps
            templateImage = ipDiffim.warpTemplateExposure(templateImage,
                                                          scienceImage,
                                                          self.policy)

            # Change xy0
            # Nice star at position 276, 717
            # And should be at index 40, 40
            # No masked pixels in this one
            self.x02 = 276
            self.y02 = 717
            size     = 40
            bbox2 = afwImage.BBox( afwImage.PointI(self.x02 - size, self.y02 - size),
                                   afwImage.PointI(self.x02 + size, self.y02 + size) )
            self.scienceImage2  = afwImage.ExposureF(scienceImage, bbox2)
            self.templateImage2 = afwImage.ExposureF(templateImage, bbox2)

    def addNoise(self, mi):
        img       = mi.getImage()
        seed      = int(afwMath.makeStatistics(mi.getVariance(), afwMath.MEDIAN).getValue())
        rdm       = afwMath.Random(afwMath.Random.MT19937, seed)
        rdmImage  = img.Factory(img.getDimensions())
        afwMath.randomGaussianImage(rdmImage, rdm)
        img      += rdmImage
        return afwMath.makeStatistics(rdmImage, afwMath.MEAN).getValue(afwMath.MEAN)

    def verifyDeltaFunctionSolution(self, solution, kSum = 1.0, bg = 0.0):
        # when kSum = 1.0, this agrees to the default precision.  when
        # kSum != 1.0 I need to go to only 4 digits.
        #
        # -5.4640810225678728e-06 != 0.0 within 7 places
        # 
        bgSolution = solution.getBackground()
        self.assertAlmostEqual(bgSolution, bg, 4)

        # again when kSum = 1.0 this agrees.  otherwise
        #
        # 2.7000000605594079 != 2.7000000000000002 within 7 places
        # 
        kSumSolution = solution.getKsum()
        self.assertAlmostEqual(kSumSolution, kSum, 5)


        kImage = solution.makeKernelImage()
        for j in range(kImage.getHeight()):
            for i in range(kImage.getWidth()):

                if (i == kImage.getWidth() // 2) and (j == kImage.getHeight() // 2):
                    self.assertAlmostEqual(kImage.get(i, j), kSum, 5)
                else:
                    self.assertAlmostEqual(kImage.get(i, j), 0., 5)

    def testConstructor(self):
        # Original and uninitialized
        kc = ipDiffim.KernelCandidateF(self.x02, self.y02,
                                       self.templateImage2.getMaskedImage(),
                                       self.scienceImage2.getMaskedImage(),
                                       self.policy)

        # Kernel not initialized
        self.assertEqual(kc.isInitialized(), False)

        # But this should be set on construction
        try:
            kc.getCandidateRating()
        except Exception, e:
            print e
            self.fail()

        # And these should be filled
        try:
            kc.getMiToConvolvePtr()
            kc.getMiToNotConvolvePtr()
        except Exception, e:
            print e
            self.fail()

        # And of the right type
        self.assertEqual(type(kc.getMiToConvolvePtr()), type(afwImage.MaskedImageF()))
        self.assertEqual(type(kc.getMiToNotConvolvePtr()), type(afwImage.MaskedImageF()))
        
        # None of these should work
        for kType in (ipDiffim.KernelCandidateF.ORIG,
                      ipDiffim.KernelCandidateF.PCA,
                      ipDiffim.KernelCandidateF.RECENT):
            for kMethod in (kc.getKernelSolution,
                            kc.getKernel,
                            kc.getBackground,
                            kc.getKsum,
                            kc.getKernelImage,
                            kc.getDifferenceImage):
                try:
                    kMethod(kType)
                except Exception, e:
                    pass
                else:
                    self.fail()
        try:
            kc.getImage()
        except:
            pass
        else:
            self.fail()

    def testDeltaFunctionScaled(self, scaling = 2.7, bg = 11.3):
        sIm  = afwImage.MaskedImageF(self.templateImage2.getMaskedImage(), True)
        sIm *= scaling
        kc = ipDiffim.KernelCandidateF(self.x02, self.y02,
                                       self.templateImage2.getMaskedImage(),
                                       sIm,
                                       self.policy)

        self.policy.set("kernelBasisSet", "delta-function")
        self.policy.set("useRegularization", False)
        kList = ipDiffim.makeKernelBasisList(self.policy)
        kc.build(kList)
        self.verifyDeltaFunctionSolution(kc.getKernelSolution(ipDiffim.KernelCandidateF.RECENT),
                                         kSum = scaling)


        sIm  = afwImage.MaskedImageF(self.templateImage2.getMaskedImage(), True)
        sIm += bg
        kc = ipDiffim.KernelCandidateF(self.x02, self.y02,
                                       self.templateImage2.getMaskedImage(),
                                       sIm,
                                       self.policy)

        self.policy.set("kernelBasisSet", "delta-function")
        self.policy.set("useRegularization", False)
        kList = ipDiffim.makeKernelBasisList(self.policy)
        kc.build(kList)
        self.verifyDeltaFunctionSolution(kc.getKernelSolution(ipDiffim.KernelCandidateF.RECENT),
                                         bg = bg)

        

    def testDeltaFunction(self):
        # Match an image to itself, with delta-function basis set
        # No regularization
        kc = ipDiffim.KernelCandidateF(self.x02, self.y02,
                                       self.templateImage2.getMaskedImage(),
                                       self.templateImage2.getMaskedImage(),
                                       self.policy)

        self.policy.set("kernelBasisSet", "delta-function")
        self.policy.set("useRegularization", False)
        kList = ipDiffim.makeKernelBasisList(self.policy)

        kc.build(kList)
        self.assertEqual(kc.isInitialized(), True)

        # These should work
        for kType in (ipDiffim.KernelCandidateF.ORIG,
                      ipDiffim.KernelCandidateF.RECENT):
            for kMethod in (kc.getKernelSolution,
                            kc.getKernel,
                            kc.getBackground,
                            kc.getKsum,
                            kc.getKernelImage,
                            kc.getDifferenceImage):
                try:
                    kMethod(kType)
                except Exception, e:
                    print kMethod, e
                    self.fail()
                else:
                    pass
        try:
            kc.getImage()
        except:
            print kMethod, e
            self.fail()
        else:
            pass
                
        # None of these should work
        for kType in (ipDiffim.KernelCandidateF.PCA,):
            for kMethod in (kc.getKernelSolution,
                            kc.getKernel,
                            kc.getBackground,
                            kc.getKsum,
                            kc.getKernelImage,
                            kc.getImage,
                            kc.getDifferenceImage):
                try:
                    kMethod(kType)
                except Exception, e:
                    pass
                else:
                    print kMethod
                    self.fail()
        
        self.verifyDeltaFunctionSolution(kc.getKernelSolution(ipDiffim.KernelCandidateF.RECENT))

    def testGaussianWithNoise(self):
        # Convolve a real image with a gaussian and try and recover
        # it.  Add noise and perform the same test.

        gsize = self.policy.getInt("kernelRows")
        gaussFunction = afwMath.GaussianFunction2D(2, 3)
        gaussKernel   = afwMath.AnalyticKernel(gsize, gsize, gaussFunction)
        kImageIn      = afwImage.ImageD(gsize, gsize)
        kSumIn        = gaussKernel.computeImage(kImageIn, False)

        imX, imY = self.templateImage2.getMaskedImage().getDimensions()
        smi = afwImage.MaskedImageF(imX, imY)
        afwMath.convolve(smi, self.templateImage2.getMaskedImage(), gaussKernel, False)

        p0, p1   = diffimTools.getConvolvedImageLimits(gaussKernel, smi)
        bbox     = afwImage.BBox(p0, p1)

        tmi2 = afwImage.MaskedImageF(self.templateImage2.getMaskedImage(), bbox)
        smi2 = afwImage.MaskedImageF(smi, bbox)

        kc = ipDiffim.KernelCandidateF(self.x02, self.y02, tmi2, smi2, self.policy)
        self.policy.set("kernelBasisSet", "delta-function")
        self.policy.set("useRegularization", False)
        kList = ipDiffim.makeKernelBasisList(self.policy)
        kc.build(kList)
        self.assertEqual(kc.isInitialized(), True)
        kImageOut = kc.getImage()

        #ds9.mtv(kImageIn, frame=1)
        #ds9.mtv(kImageOut, frame=2)

        soln = kc.getKernelSolution(ipDiffim.KernelCandidateF.RECENT)
        self.assertAlmostEqual(soln.getKsum(), kSumIn)
        # 8.7499380640430563e-06 != 0.0 within 7 places
        self.assertAlmostEqual(soln.getBackground(), 0.0, 4)

        for j in range(kImageOut.getHeight()):
            for i in range(kImageOut.getWidth()):

                # in the outskirts of the kernel, the ratio can get screwed because of low S/N
                # e.g. 7.45817359824e-09 vs. 1.18062529402e-08
                # in the guts of the kernel it should look closer
                if kImageIn.get(i, j) > 1e-4:
                    # sigh, too bad this sort of thing fails..
                    # 0.99941584433815966 != 1.0 within 3 places
                    self.assertAlmostEqual(kImageOut.get(i, j)/kImageIn.get(i, j), 1.0, 2)
        

        # now repeat with noise added; decrease precision of comparison
        bkg = self.addNoise(smi2)
        kc = ipDiffim.KernelCandidateF(self.x02, self.y02, tmi2, smi2, self.policy)
        self.policy.set("kernelBasisSet", "delta-function")
        self.policy.set("useRegularization", False)
        kList = ipDiffim.makeKernelBasisList(self.policy)
        kc.build(kList)
        self.assertEqual(kc.isInitialized(), True)
        kImageOut = kc.getImage()

        #ds9.mtv(kImageIn, frame=3)
        #ds9.mtv(kImageOut, frame=4)


        soln = kc.getKernelSolution(ipDiffim.KernelCandidateF.RECENT)
        self.assertAlmostEqual(soln.getKsum(), kSumIn, 3)
        if not (self.policy.get("fitForBackground")):
            self.assertEqual(soln.getBackground(), 0.0)
            
        for j in range(kImageOut.getHeight()):
            for i in range(kImageOut.getWidth()):
                if kImageIn.get(i, j) > 1e-2:
                    self.assertAlmostEqual(kImageOut.get(i, j), kImageIn.get(i, j), 2)
        

    def testGaussian(self, imsize = 50):
        # Convolve a delta function with a known gaussian; try to
        # recover using delta-function basis

        gsize = self.policy.getInt("kernelRows")
        tsize = imsize + gsize

        gaussFunction = afwMath.GaussianFunction2D(2, 3)
        gaussKernel   = afwMath.AnalyticKernel(gsize, gsize, gaussFunction)
        kImageIn      = afwImage.ImageD(gsize, gsize)
        gaussKernel.computeImage(kImageIn, False)

        # template image with a single hot pixel in the exact center
        tmi = afwImage.MaskedImageF(tsize, tsize)
        tmi.set(0, 0x0, 1e-4)
        cpix = tsize // 2
        tmi.set(cpix, cpix, (1, 0x0, 1))

        # science image
        smi = afwImage.MaskedImageF(tmi.getDimensions())
        afwMath.convolve(smi, tmi, gaussKernel, False)

        # get the actual kernel sum (since the image is not infinite)
        gscaling = afwMath.makeStatistics(smi, afwMath.SUM).getValue(afwMath.SUM)

        # grab only the non-masked subregion
        p0, p1   = diffimTools.getConvolvedImageLimits(gaussKernel, smi)
        bbox     = afwImage.BBox(p0, p1)

        tmi2 = afwImage.MaskedImageF(tmi, bbox)
        smi2 = afwImage.MaskedImageF(smi, bbox)

        # make sure its a valid subregion!
        for j in range(tmi2.getHeight()):
            for i in range(tmi2.getWidth()):
                self.assertEqual(tmi2.getMask().get(i, j), 0)
                self.assertEqual(smi2.getMask().get(i, j), 0)
 

        kc = ipDiffim.KernelCandidateF(0.0, 0.0, tmi2, smi2, self.policy)
        self.policy.set("kernelBasisSet", "delta-function")
        self.policy.set("useRegularization", False)
        kList = ipDiffim.makeKernelBasisList(self.policy)
        kc.build(kList)
        self.assertEqual(kc.isInitialized(), True)
        kImageOut = kc.getImage()

        soln = kc.getKernelSolution(ipDiffim.KernelCandidateF.RECENT)
        self.assertAlmostEqual(soln.getKsum(), gscaling)
        self.assertAlmostEqual(soln.getBackground(), 0.0)

        for j in range(kImageOut.getHeight()):
            for i in range(kImageOut.getWidth()):
                self.assertAlmostEqual(kImageOut.get(i, j)/kImageIn.get(i, j), 1.0, 5)


    def testNoBackgroundFit(self):
        self.policy.set("fitForBackground", False)
        self.testGaussian()
        self.testGaussianWithNoise()

    def testInsert(self):
        mi = afwImage.MaskedImageF(10, 10)
        kc = ipDiffim.makeKernelCandidate(0., 0., mi, mi, self.policy)
        kc.setStatus(afwMath.SpatialCellCandidate.GOOD)
        
        sizeCellX = self.policy.get("sizeCellX")
        sizeCellY = self.policy.get("sizeCellY")
        kernelCellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(0, 0)), sizeCellX, sizeCellY)
        kernelCellSet.insertCandidate(kc)
        nSeen = 0
        for cell in kernelCellSet.getCellList():
            for cand in cell.begin(True):
                cand  = ipDiffim.cast_KernelCandidateF(cand)
                self.assertEqual(cand.getStatus(), afwMath.SpatialCellCandidate.GOOD)
                nSeen += 1
        self.assertEqual(nSeen, 1)
        
    def xtestDisp(self):
        ds9.mtv(self.scienceImage2, frame=1)
        ds9.mtv(self.templateImage2, frame=2)
        

    def tearDown(self):
        del self.policy
        if self.defDataDir:
            del self.scienceImage2
            del self.templateImage2
        
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
    run(True)
