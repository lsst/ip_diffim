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

Verbosity = 1
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

display = False
writefits = False

# This one looks for the PCA of a known Gaussian smearing kernel with
# noise added to the smeared image

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
        
    def setUp(self):
        self.policy      = ipDiffim.generateDefaultPolicy(diffimPolicy, modify=False)
        self.kCols       = self.policy.getInt('kernelCols')
        self.kRows       = self.policy.getInt('kernelRows')
        
        # gaussian reference kernel
        self.gSize         = self.kCols
        self.gaussFunction = afwMath.GaussianFunction2D(2, 3)
        self.gaussKernel   = afwMath.AnalyticKernel(self.gSize, self.gSize, self.gaussFunction)
        self.kImageIn      = afwImage.ImageD(self.gSize, self.gSize)
        self.kSumIn        = self.gaussKernel.computeImage(self.kImageIn, False)

        # known input images
        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:
            defSciencePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                          "v5-e0-c011-a00.sci")
            self.templateImage = afwImage.MaskedImageF(defSciencePath)
            
            # image statistics
            self.dStats  = ipDiffim.ImageStatisticsF()

            self.scienceImage = self.templateImage.Factory( self.templateImage.getDimensions() )
            afwMath.convolve(self.scienceImage, self.templateImage, self.gaussKernel, False)
            self.addNoise(self.scienceImage)

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

    def testPca(self):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up"
            return

        self.policy.set("kernelBasisSet", "delta-function")
        self.policy.set("useRegularization", False)
        self.policy.set("spatialKernelOrder", 0)
        self.policy.set("spatialBgOrder", 0)

        # NOTE - you need to subtract off background from the image
        # you run detection on.  Here it is the template.
        algorithm = self.policy.get("backgroundPolicy.algorithm")
        binsize   = self.policy.get("backgroundPolicy.binsize")

        # hack to get around ticket #987
        # hack to get around ticket #987
        # hack to get around ticket #987
        binsize   = binsize//2
        # hack to get around ticket #987
        # hack to get around ticket #987
        # hack to get around ticket #987
        
        bctrl     = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE)
        bctrl.setNxSample(int(self.templateImage.getWidth()//binsize) + 1)
        bctrl.setNySample(int(self.templateImage.getHeight()//binsize) + 1)
        image   = self.templateImage.getImage() 
        backobj = afwMath.makeBackground(image, bctrl)
        image  -= backobj.getImageF()
        del image; del backobj

        result = ipDiffim.createPsfMatchingKernel(self.templateImage,
                                                  self.scienceImage,
                                                  self.policy)
        
        spatialKernel, spatialBg, kernelCellSet = result
        basisList = spatialKernel.getKernelList()
        kernel0 = basisList[0]
        im0     = afwImage.ImageD(spatialKernel.getDimensions())
        ksum0   = kernel0.computeImage(im0, False)    

        if display:
            frame=0
            ds9.mtv(self.kImageIn, frame=0)
            frame += 1
            
            for idx in range(min(5, len(basisList))):
                kernel = basisList[idx]
                im     = afwImage.ImageD(spatialKernel.getDimensions())
                ksum   = kernel.computeImage(im, False)    
                ds9.mtv(im, frame=frame)
                frame += 1
                
        # mean kernel is the known applied kernel
        self.assertAlmostEqual(ksum0, self.kSumIn, 2)
        for j in range(self.kImageIn.getHeight()):
            for i in range(self.kImageIn.getWidth()):
                self.assertAlmostEqual(self.kImageIn.get(i,j), im0.get(i, j), 2)

        # and the eigen kernel sums are close to 0
        # tough to test each pixel since there are fluctuations
        # so this test just says, they should carry no overall power
        for idx in range(1, len(basisList)):
            kernel = basisList[idx]
            im     = afwImage.ImageD(spatialKernel.getDimensions())
            ksum   = kernel.computeImage(im, False)    

            self.assertAlmostEqual(ksum, 0, 2)
                    
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
