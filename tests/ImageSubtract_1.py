import os
import pdb
import unittest
import lsst.mwi.tests as tests
import lsst.mwi.utils as mwiu
import lsst.fw.Core.fwLib as fw
import lsst.detection.detectionLib as detection

try:
    type(verbose)
except NameError:
    verbose = 0
    mwiu.Trace_setVerbosity("imageproc.ImageSubtract", verbose)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def initializeTestCases(kernelRows=7, kernelCols=7, kernelSpatialOrder=0, backgroundSpatialOrder=0):
    scienceMaskedImage  = fw.MaskedImageD()
    templateMaskedImage = fw.MaskedImageD()

    scienceMaskedImage.readFits(os.path.join(os.environ['FWDATA_DIR'], '871034p_1_MI'))
    templateMaskedImage.readFits(os.path.join(os.environ['FWDATA_DIR'], '871034p_1_MI'))
    
    # create basis vectors
    kernelBasisVec = fw.vectorKernelPtrD()
    imageproc.generateDeltaFunctionKernelSet_D(kernelRows, kernelCols, kernelBasisVec)
    
    # create output kernel pointer
    kernelPtr = imageproc.LinearCombinationKernelPtrTypeD(fw.LinearCombinationKernelD())
    
    # and its function for spatial variation
    kernelFunctionPtr = fw.Function2PtrTypeD(fw.PolynomialFunction2D(kernelSpatialOrder))
    
    # and background function
    backgroundFunctionPtr = fw.Function2PtrTypeD(fw.PolynomialFunction2D(backgroundSpatialOrder))
    
    # get good footprints
    footprintList = detection.FootprintContainerT()
    imageproc.getCollectionOfMaskedImagesForPsfMatching(footprintList)

    return scienceMaskedImage, templateMaskedImage, kernelBasisVec, kernelPtr, kernelFunctionPtr, backgroundFunctionPtr, footprintList

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class DeltaFunctionTestCase(unittest.TestCase):
    """Test case for deriving Delta Function as best kernel to match two of the same images"""
    def setUp(self):
        testObjects = initializeTestCases()
        self.scienceMaskedImage    = testObjects[0]
        self.templateMaskedImage   = testObjects[1]
        self.kernelBasisVec        = testObjects[2]
        self.kernelPtr             = testObjects[3]
        self.kernelFunctionPtr     = testObjects[4]
        self.backgroundFunctionPtr = testObjects[5]
        self.footprintList         = testObjects[6]

        imageproc.computePsfMatchingKernelForMaskedImage_FU8DD(self.templateMaskedImage,
                                                               self.scienceMaskedImage,
                                                               self.kernelBasisVec,
                                                               self.footprintList,
                                                               self.kernelPtr,
                                                               self.kernelFunctionPtr,
                                                               self.backgroundFunctionPtr,
                                                               self.policy)
        
    def tearDown(self):
        del self.scienceMaskedImage    
        del self.templateMaskedImage   
        del self.kernelBasisVec        
        del self.kernelPtr             
        del self.kernelFunctionPtr     
        del self.backgroundFunctionPtr 
        del self.footprintList         

    def testDeltaFunction(self):
        """Make sure that the mean output kernel is a delta function"""
        meanKernel = kernelPtr[0].computeNewImage()
        for i in range(meanKernel.getCols()):
            for j in range(meanKernel.getRows()):
                if i==j and i==meanKernel.getCols()/2:
                    self.assertAlmostEqual(meanKernel[i][j], 1.0)
                else:
                    self.assertAlmostEqual(meanKernel[i][j], 0.0)

    def testKernelModel(self):
        """Make sure that the basis functions have the correct power"""
        parameters = kernelPtr.getKernelFunction().getParameters()
        self.assertAlmostEqual(parameters[0], 1.0)
        for i in range(1, len(parameters)):
            self.assertAlmostEqual(parameters[0], 0.0)

    def testBackground(self):
        """Make sure that the background is zero"""
        parameters = backgroundFunctionPtr.getParameters()
        self.assertAlmostEqual(parameters[0], 0.0)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class BackgroundTestCase(unittest.TestCase):
    """Test case for determining a DC offset (background) between 2 images"""
    def setUp(self):
        testObjects = initializeTestCases()
        self.scienceMaskedImage    = testObjects[0]
        self.templateMaskedImage   = testObjects[1]
        self.kernelBasisVec        = testObjects[2]
        self.kernelPtr             = testObjects[3]
        self.kernelFunctionPtr     = testObjects[4]
        self.backgroundFunctionPtr = testObjects[5]
        self.footprintList         = testObjects[6]

        self.backgroundOffset      = 100
        self.scienceMaskedImage   += self.backgroundOffset 

        imageproc.computePsfMatchingKernelForMaskedImage_FU8DD(self.templateMaskedImage,
                                                               self.scienceMaskedImage,
                                                               self.kernelBasisVec,
                                                               self.footprintList,
                                                               self.kernelPtr,
                                                               self.kernelFunctionPtr,
                                                               self.backgroundFunctionPtr,
                                                               self.policy)


    def tearDown(self):
        del self.scienceMaskedImage    
        del self.templateMaskedImage   
        del self.kernelBasisVec        
        del self.kernelPtr             
        del self.kernelFunctionPtr     
        del self.backgroundFunctionPtr 
        del self.footprintList         

    def testDeltaFunction(self):
        """Make sure that the mean output kernel is a delta function"""
        meanKernel = kernelPtr[0].computeNewImage()
        for i in range(meanKernel.getCols()):
            for j in range(meanKernel.getRows()):
                if i==j and i==meanKernel.getCols()/2:
                    self.assertAlmostEqual(meanKernel[i][j], 1.0)
                else:
                    self.assertAlmostEqual(meanKernel[i][j], 0.0)

    def testKernelModel(self):
        """Make sure that the basis functions have the correct power"""
        parameters = kernelPtr.getKernelFunction().getParameters()
        self.assertAlmostEqual(parameters[0], 1.0)
        for i in range(1, len(parameters)):
            self.assertAlmostEqual(parameters[0], 0.0)

    def testBackground(self):
        """Make sure that the background is correct"""
        parameters = backgroundFunctionPtr.getParameters()
        self.assertAlmostEqual(parameters[0], self.backgroundOffset)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class GaussianTestCase(unittest.TestCase):
    """Test case for recovering the function you use to convolve an input image"""
    def setUp(self):
        testObjects = initializeTestCases()
        self.scienceMaskedImage    = testObjects[0]
        self.templateMaskedImage   = testObjects[1]
        self.kernelBasisVec        = testObjects[2]
        self.kernelPtr             = testObjects[3]
        self.kernelFunctionPtr     = testObjects[4]
        self.backgroundFunctionPtr = testObjects[5]
        self.footprintList         = testObjects[6]

        sigmaX = 1.0
        sigmaY = 2.0
        gaussFunctionPtr = fw.Function2PtrTypeD(fw.GaussianFunction2D(sigmaX,sigmaY))
        #gaussKernel      = fw.KernelPtrTypeD()

        kernelCols = 7
        kernelRows = 7
        self.gaussKernel = fw.AnalyticKernelD(gaussFunctionPtr, kernelCols, kernelRows)
        
        fw.convolve(self.scienceMaskedImage, self.gaussKernel, 0.0, 1)

        imageproc.computePsfMatchingKernelForMaskedImage_FU8DD(self.templateMaskedImage,
                                                               self.scienceMaskedImage,
                                                               self.kernelBasisVec,
                                                               self.footprintList,
                                                               self.kernelPtr,
                                                               self.kernelFunctionPtr,
                                                               self.backgroundFunctionPtr,
                                                               self.policy)

    def tearDown(self):
        del self.scienceMaskedImage    
        del self.templateMaskedImage   
        del self.kernelBasisVec        
        del self.kernelPtr             
        del self.kernelFunctionPtr     
        del self.backgroundFunctionPtr 
        del self.footprintList
        del self.gaussKernel

    def testDeltaFunction(self):
        """Make sure that the mean output kernel is a Gaussian"""
        meanKernel = kernelPtr[0].computeNewImage()
        convKernel = self.gaussKernel.computeNewImage()
        for i in range(meanKernel.getCols()):
            for j in range(meanKernel.getRows()):
                self.assertAlmostEqual(meanKernel[i][j], convKernel[i][j])

    def testKernelModel(self):
        """Make sure that the basis functions have the correct power"""
        parameters = kernelPtr.getKernelFunction().getParameters()
        self.assertAlmostEqual(parameters[0], 1.0)
        for i in range(1, len(parameters)):
            self.assertAlmostEqual(parameters[0], 0.0)

    def testBackground(self):
        """Make sure that the background is correct"""
        parameters = backgroundFunctionPtr.getParameters()
        self.assertAlmostEqual(parameters[0], 0.0)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(DeltaFunctionTestCase)
    suites += unittest.makeSuite(BackgroundTestCase)
    suites += unittest.makeSuite(GaussianTestCase)
    return unittest.TestSuite(suites)

if __name__ == "__main__":
    tests.run(suite())
