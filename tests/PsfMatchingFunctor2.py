#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import os
import pdb
import sys
import numpy as num
import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as logging

import lsst.afw.display.ds9 as ds9

verbosity = 4
logging.Trace_setVerbosity('lsst.ip.diffim', verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'policy', 'ImageSubtractStageDictionary.paf')

display = False
writefits = False

# Recover a known smoothing kernel applied to real CFHT data with
# additional noise added

class DiffimTestCases(unittest.TestCase):
    
    # D = I - (K.x.T + bg)
        
    def setUp(self):
        self.policy    = ipDiffim.generateDefaultPolicy(diffimPolicy, modify=False)
        self.kCols     = self.policy.getInt('kernelCols')
        self.kRows     = self.policy.getInt('kernelRows')
        self.basisList = ipDiffim.generateDeltaFunctionBasisSet(self.kCols, self.kRows)

        # gaussian reference kernel
        self.gSize         = self.kCols
        self.gaussFunction = afwMath.GaussianFunction2D(2, 3)
        self.gaussKernel   = afwMath.AnalyticKernel(self.gSize, self.gSize, self.gaussFunction)
        self.kImageIn      = afwImage.ImageD(self.gSize, self.gSize)
        self.gaussKernel.computeImage(self.kImageIn, False)

        # difference imaging functor
        self.kFunctor      = ipDiffim.PsfMatchingFunctorF(self.basisList)

        # known input images
        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:
            defSciencePath = os.path.join(self.defDataDir, "CFHT", "D4", 
                                          "cal-53535-i-797722_1")
            self.scienceImage  = afwImage.MaskedImageF(defSciencePath)
            
    def tearDown(self):
        del self.policy
        del self.basisList
        del self.gaussFunction
        del self.gaussKernel
        del self.kImageIn
        del self.kFunctor
        if self.defDataDir:
            del self.scienceImage

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

    def doGaussian(self, kNorm, imsize=60, xloc=1118, yloc=2483, addNoise=False):
        # NOTE : the size of these images have to be bigger
        #        size you lose pixels due to the convolution with the gaussian
        #        so adjust the size a bit to compensate 
        imsize += self.gSize

        # chop out a region around a known object
        bbox = afwImage.BBox( afwImage.PointI(xloc - imsize/2,
                                              yloc - imsize/2),
                              afwImage.PointI(xloc + imsize/2,
                                              yloc + imsize/2) )
        tmi  = afwImage.MaskedImageF(self.scienceImage, bbox)

        # now convolve it with a gaussian to make a science image
        smi = afwImage.MaskedImageF(imsize, imsize)
        afwMath.convolve(smi, tmi, self.gaussKernel, kNorm)

        if addNoise:
            self.addNoise(smi)
            
        # grab only the non-masked subregion
        bbox     = afwImage.BBox(afwImage.PointI(self.gaussKernel.getCtrX(),
                                                 self.gaussKernel.getCtrY()) ,
                                 afwImage.PointI(imsize - (self.gaussKernel.getWidth()  -
                                                           self.gaussKernel.getCtrX()),
                                                 imsize - (self.gaussKernel.getHeight() -
                                                           self.gaussKernel.getCtrY())))

        # For testing only
        invert = False
        if invert:
            tmi2     = afwImage.MaskedImageF(smi, bbox)
            smi2     = afwImage.MaskedImageF(tmi, bbox)
        else:
            tmi2     = afwImage.MaskedImageF(tmi, bbox)
            smi2     = afwImage.MaskedImageF(smi, bbox)

        # OUTPUT
        if display:
            ds9.mtv(tmi2, frame=1)
            ds9.mtv(smi2, frame=2)
            ds9.mtv(self.kImageIn, frame=3)
        if writefits:
            self.kImageIn.writeFits('kIn.fits')
            tmi2.writeFits('t2')
            smi2.writeFits('s2')
            
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
        kb     = self.kFunctor.getSolution()
        kSoln  = kb.first
        bgSoln = kb.second
        kImageOut = afwImage.ImageD(self.kCols, self.kRows)
        kSum      = kSoln.computeImage(kImageOut, False)
        diffIm    = ipDiffim.convolveAndSubtract(tmi2, smi2, kSoln, bgSoln)
        bbox      = afwImage.BBox(afwImage.PointI(kSoln.getCtrX(),
                                                  kSoln.getCtrY()) ,
                                  afwImage.PointI(diffIm.getWidth() - (kSoln.getWidth()  -
                                                                       kSoln.getCtrX()),
                                                  diffIm.getHeight() - (kSoln.getHeight() -
                                                                        kSoln.getCtrY())))
        diffIm2   = afwImage.MaskedImageF(diffIm, bbox)

        # OUTPUT
        if display:
            ds9.mtv(kImageOut, frame=4)
            ds9.mtv(diffIm2, frame=5)
        if writefits:
            kImageOut.writeFits('k1.fits')
            diffIm2.writeFits('dA2')

        # kernel sum should be 1.0 if kNorm
        if kNorm:
            if addNoise:
                self.assertAlmostEqual(kSum, 1.0, 1)
            else:
                self.assertAlmostEqual(kSum, 1.0, 5)
        
        # make sure the derived kernel looks like the input kernel
        # only if you haven't normalized the kernel sum to be 1.0 during the initial convolution
        for j in range(kImageOut.getHeight()):
            for i in range(kImageOut.getWidth()):
                if not kNorm:
                    if addNoise:
                        self.assertAlmostEqual(kImageOut.get(i, j), self.kImageIn.get(i, j), 1)
                    else:
                        self.assertAlmostEqual(kImageOut.get(i, j), self.kImageIn.get(i, j), 4)

        # finally, stats on the diffim
        imstat = ipDiffim.ImageStatisticsF()
        imstat.apply(diffIm2)
        self.assertAlmostEqual(imstat.getMean(), 0, 4)
        if not addNoise:
            self.assertAlmostEqual(imstat.getRms(), 0, 4)
            
        
    def testGaussian(self):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up"
            return
        
        self.doGaussian(kNorm=True, addNoise=False)
        self.doGaussian(kNorm=False, addNoise=False)
        
        self.doGaussian(kNorm=True, addNoise=True)
        self.doGaussian(kNorm=False, addNoise=True)

        
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
    if '-d' in sys.argv:
        display = True
    if '-w' in sys.argv:
        writefits = True
        
    run(True)
