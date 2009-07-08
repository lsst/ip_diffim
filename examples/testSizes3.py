#!/usr/bin/env python
import os, pdb, sys

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
import numpy as num
import pylab
import lsst.afw.display.ds9 as ds9

Verbosity = 4
logging.Trace_setVerbosity('lsst.ip.diffim', Verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

display = False
writefits = False

# This one just creates example convolution and deconvolution kernels

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
        
    def setUp(self):
        self.policy      = pexPolicy.Policy.createPolicy(diffimPolicy)
        self.kCols       = self.policy.getInt('kernelCols')
        self.kRows       = self.policy.getInt('kernelRows')
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

        # Remap the template to the image; replace self.templateImage with warped image
        wKernel = afwMath.makeWarpingKernel('lanczos4')
        self.remappedImage = self.templateImage.Factory(
            self.scienceImage.getWidth(), 
            self.scienceImage.getHeight(),
            self.scienceImage.getWcs())
        self.remappedImage.getMaskedImage().setXY0( self.scienceImage.getMaskedImage().getXY0() )
        afwMath.warpExposure(self.remappedImage, 
                             self.templateImage, 
                             wKernel)
        self.templateImage = self.remappedImage

        # edge bit
        self.edgeBit = afwImage.MaskU().getMaskPlane('EDGE')

        # image statistics
        self.dStats  = ipDiffim.ImageStatisticsF()

        # footprints
        self.detSet     = afwDetection.makeDetectionSet(self.scienceImage.getMaskedImage(), afwDetection.Threshold(575))
        self.footprints = self.detSet.getFootprints()
        print '#', len(self.footprints)
        
    def tearDown(self):
        del self.policy

    def applyFunctor(self, imscale=4, invert=False, foffset=0, xloc=397, yloc=580):
        imsize = int(imscale * self.kCols)

        # chop out a region around a known object
        bbox = afwImage.BBox( afwImage.PointI(xloc - imsize/2,
                                              yloc - imsize/2),
                              afwImage.PointI(xloc + imsize/2,
                                              yloc + imsize/2) )

        if invert:
            try:
                tmi  = afwImage.MaskedImageF(self.scienceImage.getMaskedImage(),  bbox)
                smi  = afwImage.MaskedImageF(self.templateImage.getMaskedImage(), bbox)
            except:
                return None
        else:
            try:
                smi  = afwImage.MaskedImageF(self.scienceImage.getMaskedImage(),  bbox)
                tmi  = afwImage.MaskedImageF(self.templateImage.getMaskedImage(), bbox)
            except:
                return None

        # OUTPUT
        if display:
            ds9.mtv(tmi, frame=0+foffset)
            ds9.mtv(smi, frame=1+foffset)
        if writefits:
            tmi.writeFits('t')
            smi.writeFits('s')
            
        # estimate of the variance
        var  = afwImage.MaskedImageF(smi, True)
        var -= tmi

        # accepts : imageToConvolve, imageToNotConvolve
        try:
            self.kFunctor.apply(tmi.getImage(), smi.getImage(), var.getVariance(), self.policy)
        except:
            return None
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
        print 'Diffim residuals1 : %.2f +/- %.2f (%.3f)' % (self.dStats.getMean(), self.dStats.getRms(), kSum)
                                                       
        # Second iteration
        try:
            self.kFunctor.apply(tmi.getImage(), smi.getImage(), diffIm.getVariance(), self.policy)
        except:
            return kImageOut
        kernel    = self.kFunctor.getKernel()
        kSum      = kernel.computeImage(kImageOut, False)
        diffIm    = ipDiffim.convolveAndSubtract(tmi, smi, kernel, self.kFunctor.getBackground())
        bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                  kernel.getCtrY()) ,
                                  afwImage.PointI(diffIm.getWidth() - (kernel.getWidth()  - kernel.getCtrX()),
                                                  diffIm.getHeight() - (kernel.getHeight() - kernel.getCtrY())))
        diffIm2   = afwImage.MaskedImageF(diffIm, bbox)
        self.dStats.apply( diffIm2 )
        print 'Diffim residuals2 : %.2f +/- %.2f (%.2f)' % (self.dStats.getMean(), self.dStats.getRms(), kSum)
        
        # OUTPUT
        if display:
            ds9.mtv(kImageOut, frame=2+foffset)
            ds9.mtv(diffIm2, frame=3+foffset)
        if writefits:
            kImageOut.writeFits('k2.fits')
            diffIm2.writeFits('dB2')

        return kImageOut

    def testFunctor(self):
        frame = 1
        
        for i in range(5, 24, 2):
            self.kCols       = i
            self.kRows       = i
            self.policy.set('kernelCols', i)
            self.policy.set('kernelRows', i)

            self.basisList   = ipDiffim.generateDeltaFunctionKernelSet(self.kCols, self.kRows)
            self.kFunctor    = ipDiffim.PsfMatchingFunctorF(self.basisList)
            
            for j in range(2, 11, 2):
                # make a mean kernel
                kernels = []
                for object in self.footprints:
                    # note this returns the kernel images
                    kernel = self.applyFunctor(xloc= int(0.5 * ( object.getBBox().getX0() + object.getBBox().getX1() )),
                                               yloc= int(0.5 * ( object.getBBox().getY0() + object.getBBox().getY1() )),
                                               imscale=j)
        
                    if kernel:
                        kernels.append(kernel)
        
                nKernels = len(kernels)
                kPca     = num.zeros((self.kCols*self.kRows, nKernels))
                for idx in range(nKernels):
                    kPca[:,idx] = diffimTools.vectorFromImage(kernels[idx])
                    
                kMean, kU, kVal, kCoeff = diffimTools.runPca(kPca, None)

                # new
                kSum  = num.cumsum(kVal)
                kSum /= kSum[-1]
                pylab.clf()
                pylab.plot(num.arange(len(kSum))+1, kSum, 'r-')
                pylab.xlim(1,30)
                pylab.ylim(0,1)
                pylab.savefig('%d_%d.png' % (j, i))
                
                kMeanIm = diffimTools.imageFromVector(kMean, self.kCols, self.kRows, retType=afwImage.ImageD) 
                ds9.mtv(kMeanIm, frame=frame, title='%d %d' % (j, i))
                frame += 1
                
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
