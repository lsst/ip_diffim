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

Verbosity = 4
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

display = False
writefits = False
iterate = False

# This one looks for the PCA of the convolution and devoncolution kernels

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
        
    def setUp(self):
        self.policy      = pexPolicy.Policy.createPolicy(diffimPolicy)
        self.kCols       = self.policy.getInt('kernelCols')
        self.kRows       = self.policy.getInt('kernelRows')
        self.fpGrowKsize = self.policy.getDouble('fpGrowKsize')
        self.basisList   = ipDiffim.generateDeltaFunctionKernelSet(self.kCols, self.kRows)

        # difference imaging functor
        self.kFunctor      = ipDiffim.PsfMatchingFunctorF(self.basisList)

        # known input images
        defDataDir = eups.productDir('afwdata')
        defSciencePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                      "v26-e0-c011-a00.sci")
        defTemplatePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                       "v5-e0-c011-a00.sci")
        self.scienceImage   = afwImage.ExposureF(defSciencePath)
        self.templateImage  = afwImage.ExposureF(defTemplatePath)
        self.templateImage  = ipDiffim.warpTemplateExposure(self.templateImage, self.scienceImage, self.policy)

        # image statistics
        self.dStats  = ipDiffim.ImageStatisticsF()

        # footprints
        self.detSet     = afwDetection.makeDetectionSet(self.scienceImage.getMaskedImage(), afwDetection.Threshold(575))
        self.footprints = self.detSet.getFootprints()
        
    def tearDown(self):
        del self.policy

    def applyFunctor(self, invert=False, foffset=0, xloc=397, yloc=580):
        imsize = int(3 * self.kCols)

        # chop out a region around a known object
        bbox = afwImage.BBox( afwImage.PointI(xloc - imsize/2,
                                              yloc - imsize/2),
                              afwImage.PointI(xloc + imsize/2,
                                              yloc + imsize/2) )

        # sometimes the box goes off the image; no big deal...
        try:
            if invert:
                tmi  = afwImage.MaskedImageF(self.scienceImage.getMaskedImage(),  bbox)
                smi  = afwImage.MaskedImageF(self.templateImage.getMaskedImage(), bbox)
            else:
                smi  = afwImage.MaskedImageF(self.scienceImage.getMaskedImage(),  bbox)
                tmi  = afwImage.MaskedImageF(self.templateImage.getMaskedImage(), bbox)
        except:
            return None

        # convolve science image with a gaussian for testing...
        #cmi = smi.Factory(smi.getDimensions())
        #afwMath.convolve(cmi, smi, self.gaussKernel, False)
        #smi = cmi
        
        # estimate of the variance
        var  = afwImage.MaskedImageF(smi, True)
        var -= tmi

        # accepts : imageToConvolve, imageToNotConvolve
        self.kFunctor.apply(tmi.getImage(), smi.getImage(), var.getVariance(), self.policy)
        kernel    = self.kFunctor.getKernel()
        kImageOut = afwImage.ImageD(self.kCols, self.kRows)
        kSum      = kernel.computeImage(kImageOut, False)
        diffIm    = ipDiffim.convolveAndSubtract(tmi, smi, kernel, self.kFunctor.getBackground())
        bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                  kernel.getCtrY()) ,
                                  afwImage.PointI(diffIm.getWidth() - (kernel.getWidth()  - kernel.getCtrX()),
                                                  diffIm.getHeight() - (kernel.getHeight() - kernel.getCtrY())))
        diffIm2   = afwImage.MaskedImageF(diffIm, bbox)
        self.dStats.apply( diffIm2 )
        print 'Diffim residuals1 : %.2f +/- %.2f; %.2f, %.2f' % (self.dStats.getMean(), self.dStats.getRms(),
                                                                 kSum, self.kFunctor.getBackground())
        # OUTPUT
        if display:
            ds9.mtv(tmi, frame=1+foffset)
            ds9.mtv(smi, frame=2+foffset)
        if writefits:
            tmi.writeFits('t')
            smi.writeFits('s')
            diffIm2.writeFits('d')

        if not iterate:
            return kImageOut

        # OUTPUT
#        if display:
#            ds9.mtv(kImageOut, frame=3)
#            ds9.mtv(diffIm, frame=4)
#            ds9.mtv(diffIm2, frame=5)
#        if writefits:
#            kImageOut.writeFits('k1.fits')
#            diffIm.writeFits('dA1')
#            diffIm2.writeFits('dA2')

        # Second iteration
        self.kFunctor.apply(tmi.getImage(), smi.getImage(), diffIm.getVariance(), self.policy)
        kernel    = self.kFunctor.getKernel()
        kSum      = kernel.computeImage(kImageOut, False)
        diffIm    = ipDiffim.convolveAndSubtract(tmi, smi, kernel, self.kFunctor.getBackground())
        bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                  kernel.getCtrY()) ,
                                  afwImage.PointI(diffIm.getWidth() - (kernel.getWidth()  - kernel.getCtrX()),
                                                  diffIm.getHeight() - (kernel.getHeight() - kernel.getCtrY())))
        diffIm2   = afwImage.MaskedImageF(diffIm, bbox)
        self.dStats.apply( diffIm2 )
        print 'Diffim residuals2 : %.2f +/- %.2f; %.2f, %.2f' % (self.dStats.getMean(), self.dStats.getRms(),
                                                                 kSum, self.kFunctor.getBackground())
        
        # OUTPUT
        if display:
            ds9.mtv(kImageOut, frame=3+foffset)
            ds9.mtv(diffIm2, frame=4+foffset)
        if writefits:
            kImageOut.writeFits('k2.fits')
            diffIm2.writeFits('dB2')

        return kImageOut
    
    def testFunctor(self):
        cKernels = []
        dKernels = []
        for object in self.footprints:
            # note this returns the kernel images
            cKernel = self.applyFunctor(invert=False, foffset=0,
                                        xloc= int(0.5 * ( object.getBBox().getX0() + object.getBBox().getX1() )),
                                        yloc= int(0.5 * ( object.getBBox().getY0() + object.getBBox().getY1() )))
            dKernel = self.applyFunctor(invert=True,  foffset=4, 
                                        xloc= int(0.5 * ( object.getBBox().getX0() + object.getBBox().getX1() )),
                                        yloc= int(0.5 * ( object.getBBox().getY0() + object.getBBox().getY1() )))
            if cKernel and dKernel:
                cKernels.append(cKernel)
                dKernels.append(dKernel)

        nKernels = len(cKernels)
        cM       = num.zeros((self.kCols*self.kRows, nKernels))
        dM       = num.zeros((self.kCols*self.kRows, nKernels))
        kImage   = afwImage.ImageD(self.kCols, self.kRows)
        for idx in range(nKernels):
            cM[:,idx] = diffimTools.vectorFromImage(cKernels[idx])
            dM[:,idx] = diffimTools.vectorFromImage(dKernels[idx])
            
        cMean, cU, cVal, cCoeff = diffimTools.runPca(cM, None)
        dMean, dU, dVal, dCoeff = diffimTools.runPca(dM, None)
        
        ceKernels = []
        deKernels = []
        # mean image
        ceKernels.append ( diffimTools.imageFromVector(cMean, self.kCols, self.kRows, retType=afwImage.ImageD) )
        deKernels.append ( diffimTools.imageFromVector(dMean, self.kCols, self.kRows, retType=afwImage.ImageD) )
        # eigen images
        for i in range(cU.shape[1]):
            ceKernels.append ( diffimTools.imageFromVector(cU[:,i], self.kCols, self.kRows, retType=afwImage.ImageD) )
            deKernels.append ( diffimTools.imageFromVector(dU[:,i], self.kCols, self.kRows, retType=afwImage.ImageD) )

        if writefits:
            ceKernels[0].writeFits('Mc.fits')
            deKernels[0].writeFits('Md.fits')

        if display:
            # display them
            neKernels = len(ceKernels)
            frame = 1
            for i in range(neKernels):
                ds9.mtv(ceKernels[i], frame=frame)
                frame += 1
            for i in range(neKernels):
                ds9.mtv(deKernels[i], frame=frame)
                frame += 1
            
            #pylab.plot(range(len(cVal)), cVal, 'r-')
            #pylab.plot(range(len(dVal)), dVal, 'b--')
            #pylab.show()

        cSum  = num.cumsum(cVal)
        cSum /= cSum[-1]
        dSum  = num.cumsum(dVal)
        dSum /= dSum[-1]
        print '#C', cSum
        print '#D', dSum
        
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
