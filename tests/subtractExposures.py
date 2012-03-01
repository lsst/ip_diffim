#!/usr/bin/env python
import os
import pdb
import sys
import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as logging
import lsst.pex.config as pexConfig
import lsst.ip.diffim.diffimTools as diffimTools

verbosity = 5
logging.Trace_setVerbosity('lsst.ip.diffim', verbosity)
logging.Trace_setVerbosity('ImagePsfMatchTask', verbosity)

display = False

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
        
    def setUp(self):
        self.config    = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.config.kernel.name = "AL"
        self.subconfig = self.config.kernel.active

        # Some of the tests are sensitive to the centroids returned by
        # "stdev" vs "pixel_stdev"
        self.subconfig.detectionConfig.detThresholdType = "stdev"

        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:

            defTemplatePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                           "v5-e0-c011-a00.sci")
            defSciencePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                          "v26-e0-c011-a00.sci")
            
            self.scienceImage   = afwImage.ExposureF(defSciencePath)
            self.templateImage  = afwImage.ExposureF(defTemplatePath)

            bgConfig = self.subconfig.afwBackgroundConfig
            diffimTools.backgroundSubtract(bgConfig,
                                           [self.templateImage.getMaskedImage(),
                                            self.scienceImage.getMaskedImage()])

            self.offset   = 1500
            self.bbox     = afwGeom.Box2I(afwGeom.Point2I(0, self.offset),
                                          afwGeom.Point2I(511, 2046))
            self.subconfig.spatialKernelOrder = 1
            self.subconfig.spatialBgOrder = 0

    def tearDown(self):
        del self.config
        if self.defDataDir:
            del self.scienceImage
            del self.templateImage

    def testModelType(self):
        self.runModelType(fitForBackground = True)
        self.runModelType(fitForBackground = False)
        
    def runModelType(self, fitForBackground):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up"
            return

        self.subconfig.fitForBackground = fitForBackground
            
        templateSubImage = afwImage.ExposureF(self.templateImage, self.bbox, afwImage.LOCAL)
        scienceSubImage  = afwImage.ExposureF(self.scienceImage, self.bbox, afwImage.LOCAL)

        self.subconfig.spatialModelType = 'chebyshev1'
        psfmatch1 = ipDiffim.ImagePsfMatchTask(config=self.config)
        results1 = psfmatch1.run(templateSubImage, scienceSubImage, "subtractExposures", 
                                 doWarping = True)
        differenceExposure1 = results1.subtractedImage
        spatialKernel1      = results1.psfMatchingKernel
        backgroundModel1    = results1.backgroundModel
        kernelCellSet1      = results1.kernelCellSet

        self.subconfig.spatialModelType = 'polynomial'
        psfmatch2 = ipDiffim.ImagePsfMatchTask(config=self.config)
        results2 = psfmatch2.run(templateSubImage, scienceSubImage, "subtractExposures",
                                 doWarping = True)
        differenceExposure2 = results2.subtractedImage
        spatialKernel2      = results2.psfMatchingKernel
        backgroundModel2    = results2.backgroundModel
        kernelCellSet2      = results2.kernelCellSet

        # Got the types right?
        self.assertTrue(
            spatialKernel1.getSpatialFunctionList()[0].toString().startswith('Chebyshev1Function2')
            )
        self.assertTrue(
            spatialKernel2.getSpatialFunctionList()[0].toString().startswith('PolynomialFunction2')
            )

        # First order term has zero spatial variation and sum = kernel sum
        kp1par0 = spatialKernel1.getSpatialFunctionList()[0].getParameters()
        kp2par0 = spatialKernel2.getSpatialFunctionList()[0].getParameters()
        self.assertAlmostEqual(kp1par0[0], kp2par0[0], 5)
        for i in range(1, len(kp1par0)):
            self.assertAlmostEqual(kp1par0[i], 0.0, 5)
            self.assertAlmostEqual(kp1par0[i], kp2par0[i], 5)
            
        if fitForBackground:
            # Nterms (zeroth order model)
            self.assertEqual(backgroundModel1.getNParameters(), 1)
            self.assertEqual(backgroundModel2.getNParameters(), 1)

            # Same value in function
            self.assertAlmostEqual(backgroundModel1.getParameters()[0], 
                                   backgroundModel2.getParameters()[0], 4)

            # Functions evaluates to each other
            self.assertAlmostEqual(backgroundModel1(0, 0), backgroundModel2(0, 0), 4)
            
            # Spatially...
            self.assertAlmostEqual(backgroundModel1(10, 10), backgroundModel2(10, 10), 4)

        else:
            # Kernel weights (end up different when fitting for background)
            nPar1 = len(spatialKernel1.getKernelParameters())    
            nPar2 = len(spatialKernel2.getKernelParameters())    
            kp1   = afwMath.vectorD(nPar1)
            kp2   = afwMath.vectorD(nPar2)
            spatialKernel1.computeKernelParametersFromSpatialModel(kp1, 0.0, 0.0)
            spatialKernel2.computeKernelParametersFromSpatialModel(kp2, 0.0, 0.0)
            for i in range(len(kp1)):
                self.assertAlmostEqual(kp1[i], kp2[i], 2)

            # Nterms (zeroth order)
            self.assertEqual(backgroundModel1.getNParameters(), 1)
            self.assertEqual(backgroundModel2.getNParameters(), 1)

            # Zero value in function
            self.assertAlmostEqual(backgroundModel1.getParameters()[0], 0.0)
            self.assertAlmostEqual(backgroundModel2.getParameters()[0], 0.0)

            # Function evaluates to zero
            self.assertAlmostEqual(backgroundModel1(0, 0), 0.0)
            self.assertAlmostEqual(backgroundModel2(0, 0), 0.0)
            
            # Spatially...
            self.assertAlmostEqual(backgroundModel1(10, 10), 0.0)
            self.assertAlmostEqual(backgroundModel2(10, 10), 0.0)
        
    def testWarping(self):
        # Should fail since images are not aligned
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up"
            return

        templateSubImage = afwImage.ExposureF(self.templateImage, self.bbox, afwImage.LOCAL)
        scienceSubImage  = afwImage.ExposureF(self.scienceImage, self.bbox, afwImage.LOCAL)
        psfmatch = ipDiffim.ImagePsfMatchTask(config=self.config)
        try:
            psfmatch.run(templateSubImage, scienceSubImage, "subtractExposures", doWarping = False)
        except Exception, e:
            pass
        else:
            self.fail()

    def testXY0(self):
        self.runXY0('polynomial')
        self.runXY0('chebyshev1')

    def runXY0(self, poly, fitForBackground = False):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up"
            return

        self.subconfig.spatialModelType = poly
        self.subconfig.fitForBackground = fitForBackground

        templateSubImage = afwImage.ExposureF(self.templateImage, self.bbox, afwImage.PARENT)
        scienceSubImage  = afwImage.ExposureF(self.scienceImage, self.bbox, afwImage.PARENT)

        # Have an XY0
        psfmatch  = ipDiffim.ImagePsfMatchTask(config=self.config)
        results1  = psfmatch.run(templateSubImage, scienceSubImage, "subtractExposures",
                                 doWarping = True)
        differenceExposure1 = results1.subtractedImage
        spatialKernel1      = results1.psfMatchingKernel
        backgroundModel1    = results1.backgroundModel
        kernelCellSet1      = results1.kernelCellSet

        # And then take away XY0
        templateSubImage.setXY0(afwGeom.Point2I(0, 0)) 
        scienceSubImage.setXY0(afwGeom.Point2I(0, 0))
        results2  = psfmatch.run(templateSubImage, scienceSubImage, "subtractExposures",
                                 doWarping = True)
        differenceExposure2 = results2.subtractedImage
        spatialKernel2      = results2.psfMatchingKernel
        backgroundModel2    = results2.backgroundModel
        kernelCellSet2      = results2.kernelCellSet

        # need to count up the candidates first, since its a running tally
        count = 0
        for cell in kernelCellSet1.getCellList():
            for cand1 in cell.begin(False):
                count += 1

        for cell in kernelCellSet1.getCellList():
            for cand1 in cell.begin(False):
                if cand1.getStatus() == afwMath.SpatialCellCandidate.UNKNOWN:
                    continue
                if cand1.getStatus() == afwMath.SpatialCellCandidate.BAD:
                    continue
                
                cand1 = ipDiffim.cast_KernelCandidateF(cand1)
                cand2 = ipDiffim.cast_KernelCandidateF(kernelCellSet2.getCandidateById(cand1.getId()+count))

                # positions are the same
                self.assertEqual(cand1.getXCenter(), cand2.getXCenter())
                self.assertEqual(cand1.getYCenter(), cand2.getYCenter() + self.offset)

                # kernels are the same
                im1   = cand1.getKernelImage(ipDiffim.KernelCandidateF.RECENT)
                im2   = cand2.getKernelImage(ipDiffim.KernelCandidateF.RECENT)
                for y in range(im1.getHeight()):
                    for x in range(im1.getWidth()):
                        self.assertAlmostEqual(im1.get(x, y), im2.get(x, y))

        # Spatial fits are the same
        skp1 = spatialKernel1.getSpatialParameters()
        skp2 = spatialKernel2.getSpatialParameters()
        bgp1 = backgroundModel1.getParameters()
        bgp2 = backgroundModel2.getParameters()

        # first term = kernel sum, 0, 0
        self.assertAlmostEqual(skp1[0][0], skp2[0][0])

        # On other terms, the spatial terms are the same, the zpt terms are different
        for nk in range(1, len(skp1)):

            # Zeropoint
            if poly == 'polynomial':
                self.assertNotEqual(skp1[nk][0], skp2[nk][0])
            elif poly == 'chebyshev1':
                # Cheby remaps coords, so model should be the same!
                self.assertAlmostEqual(skp1[nk][0], skp2[nk][0])
            else:
                self.fail()

            # Spatial terms
            for np in range(1, len(skp1[nk])):
                self.assertAlmostEqual(skp1[nk][np], skp2[nk][np], 4)

        for np in range(len(bgp1)):
            self.assertAlmostEqual(bgp1[np], bgp2[np], 4)

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
    
        


     
