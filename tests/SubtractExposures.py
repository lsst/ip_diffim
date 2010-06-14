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
import lsst.ip.diffim.diffimTools as diffimTools

verbosity = 8
logging.Trace_setVerbosity('lsst.ip.diffim', verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

# This one tests convolve and subtract of subimages / XY0

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
        
    def setUp(self):
        self.diffimDir    = eups.productDir('ip_diffim')
        self.diffimPolicy = os.path.join(self.diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')
        self.policy       = ipDiffim.generateDefaultPolicy(self.diffimPolicy)
        
        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:

            defTemplatePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                           "v5-e0-c011-a00.sci")
            defSciencePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                          "v26-e0-c011-a00.sci")
            
            self.scienceImage   = afwImage.ExposureF(defSciencePath)
            self.templateImage  = afwImage.ExposureF(defTemplatePath)
            
            diffimTools.backgroundSubtract(self.policy, [self.templateImage.getMaskedImage(),
                                                         self.scienceImage.getMaskedImage()])

            self.bbox     = afwImage.BBox(afwImage.PointI(0, 1500),
                                          afwImage.PointI(511, 2046))

    def tearDown(self):
        del self.policy
        if self.defDataDir:
            del self.scienceImage
            del self.templateImage

    def testAL(self):
        self.policy.set('kernelBasisSet', 'alard-lupton')
        self.policy.set('spatialKernelOrder', 1)
        self.policy.set('spatialBgOrder', 0) # already bg-subtracted
        self.policy.set('usePcaForSpatialKernel', False)
        self.runXY0()

    def xtestDFr(self):
        self.policy.set('kernelBasisSet', 'delta-function')
        self.policy.set('spatialKernelOrder', 1)
        self.policy.set('spatialBgOrder', 0) # already bg-subtracted
        self.policy.set('usePcaForSpatialKernel', True)
        self.policy.set('useRegularization', True)
        self.runXY0()

    def runXY0(self):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up"
            return

        templateSubImage = afwImage.ExposureF(self.templateImage, self.bbox)
        scienceSubImage  = afwImage.ExposureF(self.scienceImage, self.bbox)

        # Have an XY0
        results1 = ipDiffim.subtractExposures(templateSubImage, scienceSubImage, self.policy,
                                              display=True, frame=0)
        differenceExposure1, spatialKernel1, backgroundModel1, kernelCellSet1 = results1

        # And then take away XY0
        # Crap, the warping makes things slighly offset w.r.t. the previous run
        # So you don't get *exact* same coords
        templateSubImage.getMaskedImage().setXY0(0, 0)
        scienceSubImage.getMaskedImage().setXY0(0, 0)
        results2 = ipDiffim.subtractExposures(templateSubImage, scienceSubImage, self.policy,
                                              display=True, frame=3)
        differenceExposure2, spatialKernel2, backgroundModel2, kernelCellSet2 = results2

        skp1 = spatialKernel1.getSpatialParameters()
        skp2 = spatialKernel2.getSpatialParameters()
        # not all these will be the same since the coords are
        # different; at least the kernel sum (first coeff) should be
        # the same
        self.assertAlmostEqual(skp1[0][0], skp2[0][0])


        sys.exit(1)
        # DIES here 



        # and compare candidate quality
        kImage1 = afwImage.ImageD(spatialKernel1.getDimensions())
        kImage2 = afwImage.ImageD(spatialKernel2.getDimensions())
        imstats = ipDiffim.ImageStatisticsF()
        # need to count up the candidates first, since its a running tally
        count = 0
        for cell in kernelCellSet1.getCellList():
            for cand1 in cell.begin(False): 
                count += 1
                    
        # Id in the second set is id+count
        for cell in kernelCellSet1.getCellList():
            for cand1 in cell.begin(True): 
                cand1 = ipDiffim.cast_KernelCandidateF(cand1)
                
                if cand1.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                    cand2 = kernelCellSet2.getCandidateById(cand1.getId() + count)
                    cand2 = ipDiffim.cast_KernelCandidateF(cand2)
                    
                    # evaluate kernel and background at position of candidate 1
                    xCand1 = int(cand1.getXCenter())
                    yCand1 = int(cand1.getYCenter())
                    kSum1  = spatialKernel1.computeImage(kImage1, False,
                                                         afwImage.indexToPosition(xCand1),
                                                         afwImage.indexToPosition(yCand1))
                    kernel1 = afwMath.FixedKernel(kImage1)
                    background1 = backgroundModel1(afwImage.indexToPosition(xCand1),
                                                   afwImage.indexToPosition(yCand1))
                    diffIm1 = cand1.returnDifferenceImage(kernel1, background1)
                    imstats.apply(diffIm1)
                    candMean1   = imstats.getMean()
                    candRms1    = imstats.getRms()

                    # evaluate kernel and background at position of candidate 2
                    xCand2 = int(cand2.getXCenter())
                    yCand2 = int(cand2.getYCenter())
                    kSum2  = spatialKernel2.computeImage(kImage2, False,
                                                         afwImage.indexToPosition(xCand2),
                                                         afwImage.indexToPosition(yCand2))
                    kernel2 = afwMath.FixedKernel(kImage2)
                    background2 = backgroundModel2(afwImage.indexToPosition(xCand2),
                                                   afwImage.indexToPosition(yCand2))
                    diffIm2 = cand2.returnDifferenceImage(kernel2, background2)
                    imstats.apply(diffIm2)
                    candMean2   = imstats.getMean()
                    candRms2    = imstats.getRms()

                    self.assertAlmostEqual(candMean1, candMean2)
                    self.assertAlmostEqual(candRms1, candRms2)
                    self.assertAlmostEqual(kSum1, kSum2)
                    self.assertAlmostEqual(background1, background2)
        
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
    
        


     
