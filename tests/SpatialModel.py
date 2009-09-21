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
import lsst.ip.diffim.spatialKernelFit as spatialKernelFit 
import lsst.pex.logging as logging

import lsst.afw.display.ds9 as ds9

Verbosity = 1
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

display = False
writefits = False
iterate = False

# This one looks for the spatial model of a known Gaussian smearing
# kernel with noise added to the images.

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
        
    def setUp(self):
        self.policy      = pexPolicy.Policy.createPolicy(diffimPolicy)
        self.kCols       = self.policy.getInt('kernelCols')
        self.kRows       = self.policy.getInt('kernelRows')
        self.basisList   = ipDiffim.generateDeltaFunctionKernelSet(self.kCols, self.kRows)

        # difference imaging functor
        self.kFunctor      = ipDiffim.PsfMatchingFunctorF(self.basisList)

        # gaussian reference kernel
        self.gSize         = self.kCols
        self.gaussFunction = afwMath.GaussianFunction2D(2, 3)
        self.gaussKernel   = afwMath.AnalyticKernel(self.gSize, self.gSize, self.gaussFunction)
        self.kImageIn      = afwImage.ImageD(self.gSize, self.gSize)
        self.gaussKernel.computeImage(self.kImageIn, False)

        # edge bit
        self.edgeBit = afwImage.MaskU().getMaskPlane('EDGE')

        # known input images
        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:
            defImagePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                        "v5-e0-c011-a00.sci")
            self.templateImage  = afwImage.MaskedImageF(defImagePath)
            self.scienceImage   = self.templateImage.Factory( self.templateImage.getDimensions() )
            afwMath.convolve(self.scienceImage, self.templateImage, self.gaussKernel, False, self.edgeBit)
            self.addNoise(self.scienceImage) # entire image!
            
            # image statistics
            self.dStats  = ipDiffim.ImageStatisticsF()
            
            self.policy.set('minCleanFp', 10)
            self.policy.set('detThresholdType', 'value')
            self.policy.set('detThreshold', 750.)
            
            self.footprints = ipDiffim.getCollectionOfFootprintsForPsfMatching(
                self.templateImage,
                self.scienceImage,
                self.policy)

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

    def testSpatialModel(self):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up"
            return
        
        spatialCells = diffimTools.createSpatialModelKernelCells(
            self.templateImage,
            self.scienceImage,
            self.footprints,
            self.kFunctor,
            self.policy
            )
        
        # no kernels should have failed
        for scID, scPtr in enumerate(spatialCells):
            if scPtr.getNModels() == 0:
                continue
            self.assertEqual(scPtr.isUsable(), True)
        
        # no kernels should be rejected based on kSum
        diffimTools.rejectKernelSumOutliers(spatialCells, self.policy)
        for scID, scPtr in enumerate(spatialCells):
            if scPtr.getNModels() == 0:
                continue
            self.assertEqual(scPtr.isUsable(), True)

        if writefits:
            self.policy.set('debugIO', True)
            
        # spatial fit by Pca
        mKernel, eKernelVector, eVal, eCoeff = \
                 spatialKernelFit.spatialModelKernelPca(spatialCells, self.policy)

        # mean kernel is the known applied kernel
        mImage = afwImage.ImageD(self.kCols, self.kRows)
        mKernel.computeImage(mImage, False)
        for j in range(self.kImageIn.getHeight()):
            for i in range(self.kImageIn.getWidth()):
                self.assertAlmostEqual(mImage.get(i, j), self.kImageIn.get(i,j), 2)

        # do lots of spatial fitting
        for order in range(0,2):
            self.policy.set('kernelSpatialOrder', order)
            self.policy.set('bgSpatialOrder', order)
            
            for nEig in range(len(self.footprints)-1):
                
                sKernel, bgFunction = spatialKernelFit.spatialModelByPca(
                    spatialCells,
                    mKernel,
                    eKernelVector,
                    eCoeff,
                    nEig,
                    self.policy)
                
                kPar = sKernel.getSpatialParameters()

                # mean kernel has literally *no* spatial variation
                self.assertEqual(kPar[0][0], 1.)
                for i in range( 1, sKernel.getNSpatialParameters() ):
                    self.assertEqual(kPar[0][i], 0.)

                # eigen kernels should have small mean and smaller spatial variation
                for i in range( 1, nEig+1 ):
                    self.assertAlmostEqual(kPar[i][0], 0., 1)
                    for j in range( 1, sKernel.getNSpatialParameters() ):
                        self.assertAlmostEqual(kPar[i][j], 0., 4)

                # none should be rejected
                nRejected, sdqaList = spatialKernelFit.evaluateModelByPca(
                    spatialCells,
                    bgFunction, 
                    sKernel,
                    self.policy, 
                    reject=True)
                self.assertEqual(nRejected, 0)

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
