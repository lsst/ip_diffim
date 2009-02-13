import unittest

import os
import random
import numpy as num
import pdb

import eups
import lsst.pex.policy
import lsst.utils.tests as tests
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.ip.diffim as ipDiffim

ipDiffimDir = eups.productDir("ip_diffim")
if not ipDiffimDir:
    raise RuntimeError("Could not get path to ip_diffim")

policyPath = os.path.join(ipDiffimDir, "pipeline", "ImageSubtractStageDictionary.paf")
policy = lsst.pex.policy.Policy.createPolicy(policyPath)

class FindSetBitsUT(unittest.TestCase):
    """Test case for FindSetBits functor class"""
    def testFunctor(self):
        mask = afwImage.MaskU(10, 10)
        fp   = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(0,0),
                                                    mask.getWidth(),
                                                    mask.getHeight()))
        
        count = ipDiffim.FindSetBits(mask)
        count.apply(fp)
        assert(count.getBits() == 0)
        
        mask.set(0, 0, 0x1)
        count = ipDiffim.FindSetBits(mask)
        count.apply(fp)
        assert(count.getBits() == 1)
        
        mask.set(1, 1, 0x2)
        mask.set(2, 2, 0x5)
        count = ipDiffim.FindSetBits(mask)
        count.apply(fp)
        assert(count.getBits() == 8)

class FindCountsUT(unittest.TestCase):
    """Test case for FindCounts functor class"""
    def testFunctorF(self):
        image = afwImage.MaskedImageF(10, 10)
        fp    = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(0,0),
                                                     image.getWidth(),
                                                     image.getHeight()))

        image.getImage().set(100)
        image.getVariance().set(0)
        image.getMask().set(0)

        count = ipDiffim.FindCounts(image)
        count.apply(fp)
        assert(count.getCounts() == 10*10*100)
        
        image.set(0, 0, (0,0,0))
        count = ipDiffim.FindCounts(image)
        count.apply(fp)
        assert(count.getCounts() == 10*9*100)
        
        image.set(1, 1, (7,0,0))
        image.set(2, 2, (7,1,2))
        count = ipDiffim.FindCounts(image)
        count.apply(fp)
        assert(count.getCounts() == (10*7*100 + 14))

        rd = random.random()
        image.set(3, 3, (rd,0,0))
        count = ipDiffim.FindCounts(image)
        count.apply(fp)
        assert(count.getCounts() == (10*6*100 + 14 + rd))
        
class ImageStatisticsUT(unittest.TestCase):
    """Test case for ImageStatistics functor class"""
    def testFunctorF(self):
        image = afwImage.MaskedImageF(10, 10)
        fp    = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(0,0),
                                                     image.getWidth(),
                                                     image.getHeight()))

        imageData    = num.zeros((10,10))
        varianceData = num.zeros((10,10))
        maskData     = num.zeros((10,10))

        for i in range(10):
            for j in range(10):
                rd1 = random.random()
                rd2 = abs(random.random())

                imageData[i][j]    = rd1
                varianceData[i][j] = rd2

                image.set(i, j, (rd1, rd2, 0))

        count = ipDiffim.ImageStatistics(image)
        count.apply(fp)
        assert(count.getMean()     == (imageData/varianceData).mean())
        assert(count.getVariance() == (imageData/varianceData).std())

        image.getMask().set(3, 3, 1)
        maskData[3][3] = 1
        
        idx = num.where(maskData == 0)

        count = ipDiffim.ImageStatistics(image)
        count.apply(fp)
        assert(count.getMean()     == (imageData[idx]/varianceData[idx]).mean())
        assert(count.getVariance() == (imageData[idx]/varianceData[idx]).std())
        
class generateDeltaFunctionKernelSetUT(unittest.TestCase):
    """Test case for generateDeltaFunctionKernelSet"""
    def doit(self, width=10, height=10):
        ks1 = ipDiffim.generateDeltaFunctionKernelSet(width, height)
        nk  = 0
        for rowi in range(height):
            for colj in range(width):
                kernel = ks1[nk]
                kimage = afwImage.ImageD(kernel.getDimensions())
                ksum   = kernel.computeImage(kimage, False)
                assert (ksum == 1.)

                for rowk in range(height):
                    for coll in range(width):
                        if (rowi == rowk) and (colj == coll):
                            assert(kimage.get(coll, rowk) == 1.)
                        else:
                            assert(kimage.get(coll, rowk) == 0.)
                nk += 1

    def testSquare(self):
        self.doit(10, 10)

    def testNonSquare(self):
        self.doit(7, 10)
        
class testgenerateAlardLuptonKernelSetUT(unittest.TestCase):
    """Not implemented in code so an empty unit test"""
    def testDefault(self):
        pass

class getCollectionOfFootprintsForPsfMatchingUT(unittest.TestCase):
    def testDefault(self):
        pass

class convolveAndSubtractUT(unittest.TestCase):
    def testDefault(self):
        pass

class computePsfMatchingKernelForFootprintUT(unittest.TestCase):
    def testDefault(self):
        pass

#
###
#

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []

    # functors
    suites += unittest.makeSuite(FindSetBitsUT)
    suites += unittest.makeSuite(FindCountsUT)
    suites += unittest.makeSuite(ImageStatisticsUT)

    # classes
    # this might end up being deprecated in favor of ImageStatistics functors
    # suites += unittest.makeSuite(DifferenceImageStatisticsUT)

    # subroutines
    suites += unittest.makeSuite(generateDeltaFunctionKernelSetUT)
    suites += unittest.makeSuite(testgenerateAlardLuptonKernelSetUT)
    suites += unittest.makeSuite(getCollectionOfFootprintsForPsfMatchingUT)
    suites += unittest.makeSuite(convolveAndSubtractUT)
    suites += unittest.makeSuite(computePsfMatchingKernelForFootprintUT)
        
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    """Tests each unit of C-code"""
    run(True)
