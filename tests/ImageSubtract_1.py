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
            self.assertAlmostEqual(background, bg)

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
    """Make sure that the deconvolution kernel convolved with convolution kernel is delta function; NOT FINISHED"""
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

    def testDeconvolve(self):

        kernelCols = policy.get('kernelCols')
        kernelRows = policy.get('kernelRows')

        kImageCPtr    = afwImage.ImageDPtr( afwImage.ImageD(kernelCols, kernelRows) )
        kImageDPtr    = afwImage.ImageDPtr( afwImage.ImageD(kernelCols, kernelRows) )

        # hack to turn kernel image into maskedimage for convolution testing
        kMaskPtr      = afwImage.MaskUPtr( afwImage.MaskU(kernelCols, kernelRows) )
        
        for footprintID, iFootprintPtr in enumerate(self.footprintList):
            footprintBBox              = iFootprintPtr.getBBox()
            imageToConvolveStampPtr    = self.templateMaskedImage2.getSubImage(footprintBBox)
            imageToNotConvolveStampPtr = self.scienceMaskedImage.getSubImage(footprintBBox)

            if debugIO:
                imageToConvolveStampPtr.writeFits('tFits_%d' % (footprintID,))
                imageToNotConvolveStampPtr.writeFits('sFits_%d' % (footprintID,))

            # convolve
            kernelCoeffListC, backgroundC = ipDiff.computePsfMatchingKernelForFootprint(
                imageToConvolveStampPtr.get(),
                imageToNotConvolveStampPtr.get(),
                self.kernelBasisList,
                policy
                )
            footprintKernelPtrC = afwMath.LinearCombinationKernelPtr(
                afwMath.LinearCombinationKernel(self.kernelBasisList, kernelCoeffListC)
                )
            kSumC = footprintKernelPtrC.computeImage(kImageCPtr.get(), 0.0, 0.0, False)
            kMaskedImageC = afwImage.MaskedImageD(kImageCPtr, kMaskPtr)
            if debugIO:
                kImageCPtr.writeFits('kFitsC_%d.fits' % (footprintID,))
                kMaskedImageC.writeFits('kFitsMiC_%d' % (footprintID,))
            
            # deconvolve
            kernelCoeffListD, backgroundD = ipDiff.computePsfMatchingKernelForFootprint(
                imageToNotConvolveStampPtr.get(),
                imageToConvolveStampPtr.get(),
                self.kernelBasisList,
                policy
                )
            footprintKernelPtrD = afwMath.LinearCombinationKernelPtr(
                afwMath.LinearCombinationKernel(self.kernelBasisList, kernelCoeffListD)
                )
            kSumD = footprintKernelPtrD.computeImage(kImageDPtr.get(), 0.0, 0.0, False)
            if debugIO:
                kImageDPtr.writeFits('kFitsD_%d.fits' % (footprintID,))

            # check that you get a delta function
            testMaskedImage = afwMath.convolveNew(kMaskedImageC, footprintKernelPtrD.get(), 0, False)
            print testMaskedImage
            if debugIO:
                testMaskedImage.writeFits('deltaFunc_%d' % (footprintID,))
            testImage = testMaskedImage.getImage()
                
            # make sure its a delta function
            for i in range(testImage.getCols()):
                for j in range(testImage.getRows()):
                    print i, j, testImage.getVal(i,j)
                    #if i==j and i==testImage.getCols()/2:
                    #    self.assertAlmostEqual(testImage.getVal(i,j), 1.0)
                    #else:
                    #    self.assertAlmostEqual(testImage.getVal(i,j), 0.0)

            # make sure that the background is the same
            print backgroundC, backgroundD
            #self.assertAlmostEqual(backgroundC, -1 * backgroundD)

            # make sure that the kSum is scaling
            print kSumC, kSumD
            #self.assertAlmostEqual(kSumC, 1./kSumD)
            
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(ConvolveTestCase)
    suites += unittest.makeSuite(DeltaFunctionTestCase)
    #suites += unittest.makeSuite(DeconvolveTestCase) # TBD
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
