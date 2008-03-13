import os
import pdb
import unittest
import eups

import lsst.mwi.policy
import lsst.mwi.tests as tests
import lsst.mwi.utils as mwiu
import lsst.fw.Core.fwLib as fw
import lsst.detection.detectionLib as detection
import lsst.imageproc.imageprocLib as imageproc
import lsst.fw.Core.imageTestUtils as imTestUtils

verbosity = 0 # increase to see trace
mwiu.Trace_setVerbosity("lsst.fw", verbosity)

debugIO = 1

dataDir = os.environ.get("FWDATA_DIR", "")
imageProcDir = eups.productDir("imageproc", "setup")
policyPath = os.path.join(imageProcDir, "pipeline", "ImageSubtractStageDictionary.paf")
policy = lsst.mwi.policy.Policy.createPolicy(policyPath)

if not dataDir:
    raise RuntimeError("Must set up fwData to run these tests")

InputMaskedImagePath = os.path.join(dataDir, "CFHT", "D4", "cal-53535-i-797722_1")
TemplateMaskedImagePath = os.path.join(dataDir, "CFHT", "D4", "cal-53535-i-797722_1_tmpl")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def initializeTestCases():
    templateMaskedImage2 = fw.MaskedImageD()
    templateMaskedImage2.readFits(TemplateMaskedImagePath)

    templateMaskedImage = fw.MaskedImageD()
    templateMaskedImage.readFits(InputMaskedImagePath)
    scienceMaskedImage = fw.MaskedImageD()
    scienceMaskedImage.readFits(InputMaskedImagePath)
    
    kernelCols = policy.get('kernelCols')
    kernelRows = policy.get('kernelRows')
    kernelSpatialOrder = policy.get('kernelSpatialOrder')
    backgroundSpatialOrder = policy.get('backgroundSpatialOrder')

    # create basis vectors
    kernelBasisList = imageproc.generateDeltaFunctionKernelSetD(kernelCols, kernelRows)
    
    # create output kernel pointer
    kernelPtr = fw.LinearCombinationKernelDPtr(fw.LinearCombinationKernelD())
    
    # and its function for spatial variation
    kernelFunctionPtr = fw.Function2DPtr(fw.PolynomialFunction2D(kernelSpatialOrder))
    
    # and background function
    backgroundFunctionPtr = fw.Function2DPtr(fw.PolynomialFunction2D(backgroundSpatialOrder))
    
    # make single good footprint at known object position in cal-53535-i-797722_1
    size = 40
    footprintList = detection.FootprintContainerT()
    footprint = detection.FootprintPtrT(detection.Footprint( fw.BBox2i(366 - size/2,
                                                                       364 - size/2,
                                                                       size,
                                                                       size) )
                                        )
    footprintList.push_back(footprint)

    return templateMaskedImage2, templateMaskedImage, scienceMaskedImage, kernelBasisList, kernelPtr, kernelFunctionPtr, backgroundFunctionPtr, footprintList

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
        self.kernelFunctionPtr     = testObjects[5]
        self.backgroundFunctionPtr = testObjects[6]
        self.footprintList         = testObjects[7]
               
    def tearDown(self):
        del self.templateMaskedImage2    
        del self.templateMaskedImage    
        del self.scienceMaskedImage    
        del self.kernelBasisList        
        del self.kernelPtr             
        del self.kernelFunctionPtr     
        del self.backgroundFunctionPtr 
        del self.footprintList         

    def testConvolve(self, sigmaX=2, sigmaY=3):
        """Make sure that you recover a known convolution kernel"""
        kernelCols = policy.get('kernelCols')
        kernelRows = policy.get('kernelRows')
        gaussFunctionPtr = fw.Function2DPtr(fw.GaussianFunction2D(sigmaX,sigmaY))
        gaussKernel = fw.AnalyticKernelD(gaussFunctionPtr, kernelCols, kernelRows)
        convolvedScienceMaskedImage = fw.convolve(self.scienceMaskedImage, gaussKernel, 0, False)

        kImageIn  = fw.ImageD(kernelCols, kernelRows)
        kSumIn    = gaussKernel.computeImage(kImageIn, 0.0, 0.0, False)
        if debugIO:
            kImageIn.writeFits('kiFits.fits')

        kImageOut = fw.ImageD(kernelCols, kernelRows)
        
        for footprintID, iFootprintPtr in enumerate(self.footprintList):
            footprintBBox              = iFootprintPtr.getBBox()
            imageToConvolveStampPtr    = self.templateMaskedImage.getSubImage(footprintBBox)
            imageToNotConvolveStampPtr = convolvedScienceMaskedImage.getSubImage(footprintBBox)
            
            kernelCoeffList, background = imageproc.computePsfMatchingKernelForPostageStamp(
                imageToConvolveStampPtr.get(),
                imageToNotConvolveStampPtr.get(),
                self.kernelBasisList,
                policy
                )

            footprintKernelPtr = fw.LinearCombinationKernelDPtr(
                fw.LinearCombinationKernelD(self.kernelBasisList, kernelCoeffList)
                )

            kSumOut = footprintKernelPtr.computeImage(kImageOut, 0.0, 0.0, False)

            if debugIO:
                imageToConvolveStampPtr.writeFits('tFits_%d' % (footprintID,))
                imageToNotConvolveStampPtr.writeFits('sFits_%d' % (footprintID,))
                kImageOut.writeFits('koFits_%d.fits' % (footprintID,))
                
            # make sure it matches the known kernel
            for i in range(kImageOut.getCols()):
                for j in range(kImageOut.getRows()):
                    self.assertAlmostEqual(kImageIn.getPtr(i,j), kImageOut.getPtr(i,j))

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
        self.kernelFunctionPtr     = testObjects[5]
        self.backgroundFunctionPtr = testObjects[6]
        self.footprintList         = testObjects[7]
               
    def tearDown(self):
        del self.templateMaskedImage2    
        del self.templateMaskedImage    
        del self.scienceMaskedImage    
        del self.kernelBasisList        
        del self.kernelPtr             
        del self.kernelFunctionPtr     
        del self.backgroundFunctionPtr 
        del self.footprintList         

    def testDeltaFunction(self, bg=0.0, scaling=1.0):
        """Make sure that the output kernels are delta functions"""

        kernelCols = policy.get('kernelCols')
        kernelRows = policy.get('kernelRows')

        kImage = fw.ImageD(kernelCols, kernelRows)
        
        for footprintID, iFootprintPtr in enumerate(self.footprintList):
            footprintBBox              = iFootprintPtr.getBBox()
            imageToConvolveStampPtr    = self.templateMaskedImage.getSubImage(footprintBBox)

            imageToNotConvolveStampPtr = self.scienceMaskedImage.getSubImage(footprintBBox)
            # this is a bit of a hack to deal with -= and *= problems with swigged masked images
            imageArray, varianceArray, maskArray = imTestUtils.arraysFromMaskedImage(imageToNotConvolveStampPtr.get())
            imageArray += bg
            imageArray *= scaling
            varianceArray *= scaling**2
            imageToNotConvolveStamp = imTestUtils.maskedImageFromArrays( (imageArray, varianceArray, maskArray) )

            if debugIO:
                imageToConvolveStampPtr.writeFits('tFits_%d' % (footprintID,))
                imageToNotConvolveStamp.writeFits('sFits_%d' % (footprintID,))

            kernelCoeffList, background = imageproc.computePsfMatchingKernelForPostageStamp(
                imageToConvolveStampPtr.get(),
                imageToNotConvolveStamp,
                self.kernelBasisList,
                policy
                )

            footprintKernelPtr = fw.LinearCombinationKernelDPtr(
                fw.LinearCombinationKernelD(self.kernelBasisList, kernelCoeffList)
                )

            kSum = footprintKernelPtr.computeImage(kImage, 0.0, 0.0, False)
            if debugIO:
                kImage.writeFits('kFits_%d.fits' % (footprintID,))

            # make sure its a delta function
            for i in range(kImage.getCols()):
                for j in range(kImage.getRows()):
                    if i==j and i==kImage.getCols()/2:
                        self.assertAlmostEqual(kImage.getPtr(i,j), scaling)
                    else:
                        self.assertAlmostEqual(kImage.getPtr(i,j), 0.0)

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
        self.kernelFunctionPtr     = testObjects[5]
        self.backgroundFunctionPtr = testObjects[6]
        self.footprintList         = testObjects[7]
               
    def tearDown(self):
        del self.templateMaskedImage2    
        del self.templateMaskedImage    
        del self.scienceMaskedImage    
        del self.kernelBasisList        
        del self.kernelPtr             
        del self.kernelFunctionPtr     
        del self.backgroundFunctionPtr 
        del self.footprintList         

    def testDeconvolve(self):

        kernelCols = policy.get('kernelCols')
        kernelRows = policy.get('kernelRows')

        kImageCPtr    = fw.ImageDPtr( fw.ImageD(kernelCols, kernelRows) )
        kImageDPtr    = fw.ImageDPtr( fw.ImageD(kernelCols, kernelRows) )

        # hack to turn kernel image into maskedimage for convolution testing
        kMaskPtr      = fw.MaskUPtr( fw.MaskU(kernelCols, kernelRows) )
        
        for footprintID, iFootprintPtr in enumerate(self.footprintList):
            footprintBBox              = iFootprintPtr.getBBox()
            imageToConvolveStampPtr    = self.templateMaskedImage2.getSubImage(footprintBBox)
            imageToNotConvolveStampPtr = self.scienceMaskedImage.getSubImage(footprintBBox)

            if debugIO:
                imageToConvolveStampPtr.writeFits('tFits_%d' % (footprintID,))
                imageToNotConvolveStampPtr.writeFits('sFits_%d' % (footprintID,))

            # convolve
            kernelCoeffListC, backgroundC = imageproc.computePsfMatchingKernelForPostageStamp(
                imageToConvolveStampPtr.get(),
                imageToNotConvolveStampPtr.get(),
                self.kernelBasisList,
                policy
                )
            footprintKernelPtrC = fw.LinearCombinationKernelDPtr(
                fw.LinearCombinationKernelD(self.kernelBasisList, kernelCoeffListC)
                )
            kSumC = footprintKernelPtrC.computeImage(kImageCPtr.get(), 0.0, 0.0, False)
            kMaskedImageC = fw.MaskedImageD(kImageCPtr, kMaskPtr)
            if debugIO:
                kImageCPtr.writeFits('kFitsC_%d.fits' % (footprintID,))
                kMaskedImageC.writeFits('kFitsMiC_%d' % (footprintID,))
            
            # deconvolve
            kernelCoeffListD, backgroundD = imageproc.computePsfMatchingKernelForPostageStamp(
                imageToNotConvolveStampPtr.get(),
                imageToConvolveStampPtr.get(),
                self.kernelBasisList,
                policy
                )
            footprintKernelPtrD = fw.LinearCombinationKernelDPtr(
                fw.LinearCombinationKernelD(self.kernelBasisList, kernelCoeffListD)
                )
            kSumD = footprintKernelPtrD.computeImage(kImageDPtr.get(), 0.0, 0.0, False)
            if debugIO:
                kImageDPtr.writeFits('kFitsD_%d.fits' % (footprintID,))

            # check that you get a delta function
            testMaskedImage = fw.convolve(kMaskedImageC, footprintKernelPtrD.get(), 0, False)
            print testMaskedImage
            if debugIO:
                testMaskedImage.writeFits('deltaFunc_%d' % (footprintID,))
            testImage = testMaskedImage.getImage()
                
            # make sure its a delta function
            for i in range(testImage.getCols()):
                for j in range(testImage.getRows()):
                    print i, j, testImage.getPtr(i,j)
                    #if i==j and i==testImage.getCols()/2:
                    #    self.assertAlmostEqual(testImage.getPtr(i,j), 1.0)
                    #else:
                    #    self.assertAlmostEqual(testImage.getPtr(i,j), 0.0)

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
