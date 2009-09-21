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

# This one looks for the PCA of a known Gaussian smearing kernel with
# noise added to the images.

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

        # gaussian reference kernel
        self.gSize         = self.kCols
        self.gaussFunction = afwMath.GaussianFunction2D(2, 3)
        self.gaussKernel   = afwMath.AnalyticKernel(self.gSize, self.gSize, self.gaussFunction)
        self.kImageIn      = afwImage.ImageD(self.gSize, self.gSize)
        self.gaussKernel.computeImage(self.kImageIn, False)

        # known input images
        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:
            defSciencePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                          "v5-e0-c011-a00.sci")
            self.scienceImage  = afwImage.MaskedImageF(defSciencePath)
            
            # image statistics
            self.dStats  = ipDiffim.ImageStatisticsF()
            
            # footprints
            self.detSet     = afwDetection.makeDetectionSet(self.scienceImage, afwDetection.Threshold(800))
            self.footprints = self.detSet.getFootprints()

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

    def applyFunctor(self, xloc=397, yloc=580):
        imsize  = int(3 * self.kCols)
        imsize += self.gSize   # since we are convolving

        # chop out a region around a known object
        bbox = afwImage.BBox( afwImage.PointI(xloc - imsize/2,
                                              yloc - imsize/2),
                              afwImage.PointI(xloc + imsize/2,
                                              yloc + imsize/2) )

        # sometimes the box goes off the image; no big deal...
        try:
            tmi  = afwImage.MaskedImageF(self.scienceImage, bbox)
        except:
            return None

        # now convolve it with a gaussian to make a science image
        smi = tmi.Factory( tmi.getDimensions() )
        afwMath.convolve(smi, tmi, self.gaussKernel, False)
        # and add noise
        self.addNoise(smi)

        # grab only the non-masked subregion
        bbox     = afwImage.BBox(afwImage.PointI(self.gaussKernel.getCtrX(),
                                                 self.gaussKernel.getCtrY()) ,
                                 afwImage.PointI(imsize - (self.gaussKernel.getWidth()  - self.gaussKernel.getCtrX()),
                                                 imsize - (self.gaussKernel.getHeight() - self.gaussKernel.getCtrY())))
        tmi2     = afwImage.MaskedImageF(tmi, bbox)
        smi2     = afwImage.MaskedImageF(smi, bbox)
        

        # OUTPUT
        if display:
            ds9.mtv(tmi2, frame=1)
            ds9.mtv(smi2, frame=2)
            ds9.mtv(self.kImageIn, frame=3)
        if writefits:
            self.kImageIn.writeFits('kIn.fits')
            tmi2.writeFits('t2_%d' % (xloc))
            smi2.writeFits('s2_%d' % (xloc))
            
        # make sure its a valid subregion!
        mask     = smi2.getMask()
        for j in range(mask.getHeight()):
            for i in range(mask.getWidth()):
                self.assertEqual(mask.get(i, j), 0)
                
        # estimate of the variance
        var  = afwImage.MaskedImageF(smi2, True)
        var -= tmi2

        # accepts : imageToConvolve, imageToNotConvolve
        self.kFunctor.apply(tmi2.getImage(), smi2.getImage(), var.getVariance(), self.policy)
        kernel    = self.kFunctor.getKernel()
        kImageOut = afwImage.ImageD(self.kCols, self.kRows)
        kSum      = kernel.computeImage(kImageOut, False)
        diffIm    = ipDiffim.convolveAndSubtract(tmi2, smi2, kernel, self.kFunctor.getBackground())
        bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                  kernel.getCtrY()) ,
                                  afwImage.PointI(diffIm.getWidth() - (kernel.getWidth()  - kernel.getCtrX()),
                                                  diffIm.getHeight() - (kernel.getHeight() - kernel.getCtrY())))
        diffIm2   = afwImage.MaskedImageF(diffIm, bbox)


        # OUTPUT
        if display:
            ds9.mtv(kImageOut, frame=4)
            ds9.mtv(diffIm2, frame=5)
        if writefits:
            kImageOut.writeFits('k1_%d.fits' % (xloc))
            diffIm2.writeFits('dA2_%d' % (xloc))

        if not iterate:
            return kImageOut

        # Second iteration
        self.kFunctor.apply(tmi2.getImage(), smi2.getImage(), diffIm.getVariance(), self.policy)
        kernel    = self.kFunctor.getKernel()
        kSum      = kernel.computeImage(kImageOut, False)
        diffIm    = ipDiffim.convolveAndSubtract(tmi2, smi2, kernel, self.kFunctor.getBackground())
        bbox      = afwImage.BBox(afwImage.PointI(kernel.getCtrX(),
                                                  kernel.getCtrY()) ,
                                  afwImage.PointI(diffIm.getWidth() - (kernel.getWidth()  - kernel.getCtrX()),
                                                  diffIm.getHeight() - (kernel.getHeight() - kernel.getCtrY())))
        diffIm2   = afwImage.MaskedImageF(diffIm, bbox)
        
        # OUTPUT
        if display:
            ds9.mtv(kImageOut, frame=6)
            ds9.mtv(diffIm2, frame=7)
        if writefits:
            kImageOut.writeFits('k2_%d.fits' % (xloc))
            diffIm2.writeFits('dB2_%d' % (xloc))

        return kImageOut
    
    def testFunctor(self):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up"
            return
        
        kernels = []
        for object in self.footprints:
            # note this returns the kernel images
            kernel = self.applyFunctor(xloc= int(0.5 * ( object.getBBox().getX0() + object.getBBox().getX1() )),
                                       yloc= int(0.5 * ( object.getBBox().getY0() + object.getBBox().getY1() )))

            if kernel:
                kernels.append(kernel)

        nKernels = len(kernels)
        kPca     = num.zeros((self.kCols*self.kRows, nKernels))
        for idx in range(nKernels):
            kPca[:,idx] = diffimTools.vectorFromImage(kernels[idx])
            
            
        kMean, kU, kVal, kCoeff = ipDiffim.runPca(kPca, None)
        
        eKernels = []
        # mean image
        eKernels.append ( diffimTools.imageFromVector(kMean, self.kCols, self.kRows, retType=afwImage.ImageD) )
        # eigen images
        for i in range(kU.shape[1]):
            eKernels.append ( diffimTools.imageFromVector(kU[:,i], self.kCols, self.kRows, retType=afwImage.ImageD) )

        neKernels = len(eKernels)
        
        if display:
            # display them
            frame = 0
            for i in range(neKernels):
                ds9.mtv(eKernels[i], frame=frame)
                frame += 1
            
        kSum  = num.cumsum(kVal)
        kSum /= kSum[-1]
        #print '#eval', kSum

        if writefits:
            eKernels[0].writeFits('mK.fits')
            for i in range(1,neKernels):
                eKernels[i].writeFits('eK%d.fits' % (i))

        # mean kernel is the known applied kernel
        for j in range(self.kImageIn.getHeight()):
            for i in range(self.kImageIn.getWidth()):
                self.assertAlmostEqual(eKernels[0].get(i, j), self.kImageIn.get(i,j), 2)

        
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
