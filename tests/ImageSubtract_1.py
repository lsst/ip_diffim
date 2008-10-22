import os
import math
import pdb
import unittest
import eups

import lsst.pex.policy
import lsst.utils.tests as tests
import lsst.pex.logging as logging
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.detection.detectionLib as detection
import lsst.ip.diffim as ipDiff
import lsst.afw.image.testUtils as imageTest

import lsst.ip.diffim.diffimTools as ipDiffimTools

try:
    type(verbosity)
except NameError:
    verbosity = 0                       # increase to see trace
logging.Trace_setVerbosity("lsst.afw", verbosity)

try:
    type(debugIO)
except NameError:
    debugIO = 0

import lsst.afw.display.ds9 as ds9
try:
    type(display)
except NameError:
    display = False

dataDir = eups.productDir("afwdata")
if not dataDir:
    raise RuntimeError("Must set up afwdata to run these tests")
imageProcDir = eups.productDir("ip_diffim")
if not imageProcDir:
    raise RuntimeError("Could not get path to ip_diffim")
policyPath = os.path.join(imageProcDir, "pipeline", "ImageSubtractStageDictionary.paf")
policy = lsst.pex.policy.Policy.createPolicy(policyPath)

InputMaskedImagePath = os.path.join(dataDir, "CFHT", "D4", "cal-53535-i-797722_1")
TemplateMaskedImagePath = os.path.join(dataDir, "CFHT", "D4", "cal-53535-i-797722_1_tmpl")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def initializeTestCases():
    MaskedImage = afwImage.MaskedImageD       # the desired type of MaskedImage

    templateMaskedImage2 = MaskedImage()
    templateMaskedImage2.readFits(TemplateMaskedImagePath)

    templateMaskedImage = MaskedImage()
    templateMaskedImage.readFits(InputMaskedImagePath)
    scienceMaskedImage = MaskedImage()
    scienceMaskedImage.readFits(InputMaskedImagePath)
    
    kernelCols = policy.get('kernelCols')
    kernelRows = policy.get('kernelRows')
    kernelSpatialOrder = policy.get('kernelSpatialOrder')
    backgroundSpatialOrder = policy.get('backgroundSpatialOrder')

    # create basis vectors
    kernelBasisList = ipDiff.generateDeltaFunctionKernelSet(kernelCols, kernelRows)
    
    # create output kernel pointer
    kernelPtr = afwMath.LinearCombinationKernelPtr(afwMath.LinearCombinationKernel())
    
    # and its function for spatial variation
    kernelFunction = afwMath.PolynomialFunction2D(kernelSpatialOrder)
    
    # and background function
    backgroundFunction = afwMath.PolynomialFunction2D(backgroundSpatialOrder)
    
    # make single good footprint at known object position in cal-53535-i-797722_1
    size = 40
    footprintList = detection.FootprintContainerT()
    footprint = detection.FootprintPtrT(detection.Footprint( afwImage.BBox2i(366 - size/2,
                                                                             364 - size/2,
                                                                             size,
                                                                             size) )
                                        )
    footprintList.push_back(footprint)

    return templateMaskedImage2, templateMaskedImage, scienceMaskedImage, kernelBasisList, kernelPtr, kernelFunction, backgroundFunction, footprintList

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class ConvolveTestCase(unittest.TestCase):
    """Test case for deriving Delta Function as best kernel to match two of the same images"""
    def setUp(self):
        testObjects = initializeTestCases()
        self.templateMaskedImage2  = testObjects[0]
        self.templateMaskedImage   = testObjects[1]
        self.scienceMaskedImage    = testObjects[2]
        self.kernelBasisList       = testObjects[3]
        self.kernelPtr             = testObjects[4]
        self.kernelFunction        = testObjects[5]
        self.backgroundFunction    = testObjects[6]
        self.footprintList         = testObjects[7]
               
    def tearDown(self):
        del self.templateMaskedImage2    
        del self.templateMaskedImage    
        del self.scienceMaskedImage    
        del self.kernelBasisList        
        del self.kernelPtr             
        del self.kernelFunction     
        del self.backgroundFunction 
        del self.footprintList         

    def testConvolve(self, sigmaX=2, sigmaY=3):
        """Make sure that you recover a known convolution kernel"""
        kernelCols = policy.get('kernelCols')
        kernelRows = policy.get('kernelRows')
        gaussFunction = afwMath.GaussianFunction2D(sigmaX,sigmaY)
        gaussKernel = afwMath.AnalyticKernel(gaussFunction, kernelCols, kernelRows)
        convolvedScienceMaskedImage = afwMath.convolveNew(self.scienceMaskedImage, gaussKernel, 0, False)

        kImageIn  = afwImage.ImageD(kernelCols, kernelRows)
        kSumIn    = gaussKernel.computeImage(kImageIn, 0.0, 0.0, False)
        if debugIO:
            kImageIn.writeFits('kiFits.fits')

        kImageOut = afwImage.ImageD(kernelCols, kernelRows)
        
        for footprintID, iFootprintPtr in enumerate(self.footprintList):
            footprintBBox              = iFootprintPtr.getBBox()
            imageToConvolveStampPtr    = self.templateMaskedImage.getSubImage(footprintBBox)
            imageToNotConvolveStampPtr = convolvedScienceMaskedImage.getSubImage(footprintBBox)
            
            kernelCoeffList, background = ipDiff.computePsfMatchingKernelForFootprint(
                imageToConvolveStampPtr.get(),
                imageToNotConvolveStampPtr.get(),
                self.kernelBasisList,
                policy
                )

            footprintKernelPtr = afwMath.LinearCombinationKernelPtr(
                afwMath.LinearCombinationKernel(self.kernelBasisList, kernelCoeffList)
                )

            kSumOut = footprintKernelPtr.computeImage(kImageOut, 0.0, 0.0, False)

            if debugIO:
                imageToConvolveStampPtr.writeFits('tFits_%d' % (footprintID,))
                imageToNotConvolveStampPtr.writeFits('sFits_%d' % (footprintID,))
                kImageOut.writeFits('koFits_%d.fits' % (footprintID,))
            if display:
                ds9.mtv(kImageIn, frame=0)
                ds9.mtv(kImageOut, frame=1)

            # make sure it matches the known kernel
            # N.b. background test fails with floating point images 
            for i in range(kImageOut.getCols()):
                for j in range(kImageOut.getRows()):
                    if False:
                        tol = 2e-4
                        self.assertEqual(math.fabs(kImageIn.getVal(i,j) - kImageOut.getVal(i,j)) < tol, True,
                                         "K(%d,%d): |%g - %g| < %g" %
                                         (i, j, kImageIn.getVal(i,j), kImageOut.getVal(i,j), tol))
                    else:
                        self.assertAlmostEqual(kImageIn.getVal(i,j), kImageOut.getVal(i,j))                        

            # make sure that the background is zero
            self.assertAlmostEqual(background, 0.0)

            # make sure that the kSum is scaling
            self.assertAlmostEqual(kSumIn, kSumOut)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class DeltaFunctionTestCase(unittest.TestCase):
    """Test case for deriving Delta Function as best kernel to match two of the same images"""
    def setUp(self):
        testObjects = initializeTestCases()
        self.templateMaskedImage2  = testObjects[0]
        self.templateMaskedImage   = testObjects[1]
        self.scienceMaskedImage    = testObjects[2]
        self.kernelBasisList       = testObjects[3]
        self.kernelPtr             = testObjects[4]
        self.kernelFunction        = testObjects[5]
        self.backgroundFunction    = testObjects[6]
        self.footprintList         = testObjects[7]
               
    def tearDown(self):
        del self.templateMaskedImage2    
        del self.templateMaskedImage    
        del self.scienceMaskedImage    
        del self.kernelBasisList        
        del self.kernelPtr             
        del self.kernelFunction     
        del self.backgroundFunction 
        del self.footprintList         

    def testDeltaFunction(self, bg=0.0, scaling=1.0):
        """Make sure that the output kernels are delta functions"""

        kernelCols = policy.get('kernelCols')
        kernelRows = policy.get('kernelRows')

        kImage = afwImage.ImageD(kernelCols, kernelRows)
        
        for footprintID, iFootprintPtr in enumerate(self.footprintList):
            footprintBBox              = iFootprintPtr.getBBox()
            imageToConvolveStampPtr    = self.templateMaskedImage.getSubImage(footprintBBox)

            imageToNotConvolveStampPtr = self.scienceMaskedImage.getSubImage(footprintBBox)
            # this is a bit of a hack to deal with -= and *= problems with swigged masked images
            imageArray, varianceArray, maskArray = imageTest.arraysFromMaskedImage(imageToNotConvolveStampPtr.get())
            imageArray += bg
            imageArray *= scaling
            varianceArray *= scaling**2
            imageToNotConvolveStamp = imageTest.maskedImageFromArrays( (imageArray, varianceArray, maskArray) )

            if debugIO:
                imageToConvolveStampPtr.writeFits('tFits_%d' % (footprintID,))
                imageToNotConvolveStamp.writeFits('sFits_%d' % (footprintID,))

            kernelCoeffList, background = ipDiff.computePsfMatchingKernelForFootprint(
                imageToConvolveStampPtr.get(),
                imageToNotConvolveStamp,
                self.kernelBasisList,
                policy
                )

            footprintKernelPtr = afwMath.LinearCombinationKernelPtr(
                afwMath.LinearCombinationKernel(self.kernelBasisList, kernelCoeffList)
                )

            kSum = footprintKernelPtr.computeImage(kImage, 0.0, 0.0, False)
            if debugIO:
                kImage.writeFits('kFits_%d.fits' % (footprintID,))

            # make sure its a delta function
            for i in range(kImage.getCols()):
                for j in range(kImage.getRows()):
                    if i==j and i==kImage.getCols()/2:
                        self.assertAlmostEqual(kImage.getVal(i,j), scaling)
                    else:
                        self.assertAlmostEqual(kImage.getVal(i,j), 0.0)

            # make sure that the background is zero
            self.assertAlmostEqual(background, bg, places=6)

            # make sure that the kSum is scaling
            self.assertAlmostEqual(kSum, scaling)
            
    def testBackground(self, bg=17.5):
        """Make sure that the background is correctly determined"""
        self.testDeltaFunction(bg=bg)

    def testScaling(self, scaling=1.75):
        """Make sure that the output kernel is scaled correctly"""
        self.testDeltaFunction(scaling=scaling)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class DeconvolveTestCase(unittest.TestCase):
    """Make sure that the deconvolution kernel convolved with convolution kernel is delta function"""
    def setUp(self):
        testObjects = initializeTestCases()
        self.templateMaskedImage2  = testObjects[0]
        self.templateMaskedImage   = testObjects[1]
        self.scienceMaskedImage    = testObjects[2]
        self.kernelBasisList       = testObjects[3]
        self.kernelPtr             = testObjects[4]
        self.kernelFunction        = testObjects[5]
        self.backgroundFunction    = testObjects[6]
        self.footprintList         = testObjects[7]
               
    def tearDown(self):
        del self.templateMaskedImage2    
        del self.templateMaskedImage    
        del self.scienceMaskedImage    
        del self.kernelBasisList        
        del self.kernelPtr             
        del self.kernelFunction     
        del self.backgroundFunction 
        del self.footprintList         

    def testDeconvolve(self, sigma=3.0, scale=4):
        kernelCols = policy.get('kernelCols')
        kernelRows = policy.get('kernelRows')

        gaussFunction1 = afwMath.GaussianFunction2D(sigma,sigma)
        gaussKernel1   = afwMath.AnalyticKernel(gaussFunction1, scale*kernelCols+1, scale*kernelRows+1)
        gaussFunction2 = afwMath.GaussianFunction2D(sigma+2,sigma*2)
        gaussKernel2   = afwMath.AnalyticKernel(gaussFunction2, scale*kernelCols+1, scale*kernelRows+1)

        dfMaskedImage  = afwImage.MaskedImageD(2*scale*kernelCols+1, 2*scale*kernelRows+1)
        dfMaskedImage.getImage().set(scale*kernelCols+1, scale*kernelRows+1, 1)
        dfMaskedImage.getVariance().set(scale*kernelCols+1, scale*kernelRows+1, 1)
        
        maskedImage1   = afwMath.convolveNew(dfMaskedImage, gaussKernel1, 0, False)
        maskedImage2   = afwMath.convolveNew(dfMaskedImage, gaussKernel2, 0, False)

        # give it some sky so that there is no zero-valued variance
        img1  = maskedImage1.getImage()
        img1 += 1.e-4
        var1  = maskedImage1.getVariance()
        var1 += 1.e-4
        img2  = maskedImage2.getImage()
        img2 += 1.e-4
        var2  = maskedImage2.getVariance()
        var2 += 1.e-4
        # give the pixels realistic values
        maskedImage1 *= 1.e4
        maskedImage2 *= 1.e4 

        if debugIO:
            maskedImage1.writeFits('MI1a')
            maskedImage2.writeFits('MI2a')

        goodData = afwImage.BBox2i(scale*kernelCols/2+1, scale*kernelRows/2+1, scale*kernelCols, scale*kernelRows)
        maskedSubImage1Ptr = maskedImage1.getSubImage(goodData)
        maskedSubImage2Ptr = maskedImage2.getSubImage(goodData)
        kMaskPtr           = afwImage.MaskUPtr( afwImage.MaskU(kernelCols, kernelRows) )
        
        emptyStamp  = afwImage.MaskedImageD(maskedSubImage1Ptr.getCols(), maskedSubImage1Ptr.getRows())
        emptyStamp += maskedSubImage1Ptr.get()
        emptyStamp -= maskedSubImage2Ptr.get()

        if debugIO:
            maskedSubImage1Ptr.writeFits('MI1b')
            maskedSubImage2Ptr.writeFits('MI2b')

        # convolve one way
        vectorPair1 = ipDiff.computePsfMatchingKernelForFootprint2(
            maskedSubImage1Ptr.get(), maskedSubImage2Ptr.get(),
            emptyStamp, self.kernelBasisList, policy
            )
        kernelVector1, kernelErrorVector1, background1, backgroundError1 = ipDiffimTools.vectorPairToVectors(vectorPair1)
        kernelPtr1 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(self.kernelBasisList, kernelVector1)
            )
        diffIm1 = ipDiff.convolveAndSubtract(maskedSubImage1Ptr.get(),
                                             maskedSubImage2Ptr.get(),
                                             kernelPtr1,
                                             background1)
        kImage1Ptr    = afwImage.ImageDPtr( afwImage.ImageD(kernelCols, kernelRows) )
        kSum1 = kernelPtr1.computeImage(kImage1Ptr.get(), 0.0, 0.0, False)
        kMaskedImage1 = afwImage.MaskedImageD(kImage1Ptr, kMaskPtr)
        if debugIO:
            kImage1Ptr.writeFits('kFits1.fits')
            kMaskedImage1.writeFits('kFits1_Mi')
                                    

        # convolve the other way
        vectorPair2 = ipDiff.computePsfMatchingKernelForFootprint2(
            maskedSubImage2Ptr.get(), maskedSubImage1Ptr.get(),
            emptyStamp, self.kernelBasisList, policy
            )
        kernelVector2, kernelErrorVector2, background2, backgroundError2 = ipDiffimTools.vectorPairToVectors(vectorPair2)
        kernelPtr2 = afwMath.KernelPtr(
            afwMath.LinearCombinationKernel(self.kernelBasisList, kernelVector2)
            )
        diffIm2 = ipDiff.convolveAndSubtract(maskedSubImage2Ptr.get(),
                                             maskedSubImage1Ptr.get(),
                                             kernelPtr2,
                                             background2)
        kImage2Ptr    = afwImage.ImageDPtr( afwImage.ImageD(kernelCols, kernelRows) )
        kSum2 = kernelPtr2.computeImage(kImage2Ptr.get(), 0.0, 0.0, False)
        kMaskedImage2 = afwImage.MaskedImageD(kImage2Ptr, kMaskPtr)
        if debugIO:
            kImage2Ptr.writeFits('kFits2.fits')
            kMaskedImage2.writeFits('kFits2_Mi')

        # check difference images
        stats1  = ipDiff.DifferenceImageStatisticsD(diffIm1)
        self.assertAlmostEqual(stats1.getResidualMean(), 0.0)
        stats2  = ipDiff.DifferenceImageStatisticsD(diffIm2)
        self.assertAlmostEqual(stats2.getResidualMean(), 0.0)
        
        if debugIO:
            diffIm1.writeFits('DI1')
            diffIm2.writeFits('DI2')
        
        # check that you get a delta function
        testConv12  = afwMath.convolveNew(kMaskedImage1, kernelPtr2.get(), 0, False)
        testConv21  = afwMath.convolveNew(kMaskedImage2, kernelPtr1.get(), 0, False)
        testImage12 = testConv12.getImage()
        testImage21 = testConv21.getImage()
        # normalize to sum = 1.0
        sum12 = 0.0
        sum21 = 0.0
        for i in range(testImage12.getCols()):
            for j in range(testImage12.getRows()):
                sum12 += testImage12.getVal(i,j)
                sum21 += testImage21.getVal(i,j)
                
        testConv12 /= sum12
        testConv21 /= sum21

        if debugIO:
            testConv12.writeFits('deltaFunc12')
            testConv21.writeFits('deltaFunc21')
        
        testImage12 = testConv12.getImage()
        testImage21 = testConv21.getImage()
        # In practice these are close but not exact due to noise
        for i in range(testImage12.getCols()):
            for j in range(testImage12.getRows()):
                if i==j and i==testImage12.getCols()/2:
                    self.assertAlmostEqual(testImage12.getVal(i,j), 1.0, places=2)
                    self.assertAlmostEqual(testImage21.getVal(i,j), 1.0, places=2)
                else:
                    self.assertAlmostEqual(testImage12.getVal(i,j), 0.0, places=2)
                    self.assertAlmostEqual(testImage21.getVal(i,j), 0.0, places=2)
                
            
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(ConvolveTestCase)
    suites += unittest.makeSuite(DeltaFunctionTestCase)
    suites += unittest.makeSuite(DeconvolveTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
