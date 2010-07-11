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

import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as logging

verbosity = 1
logging.Trace_setVerbosity('lsst.ip.diffim', verbosity)

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

# This tests delta function kernels with delta-function images, and
# known gaussian smoothing kernels

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

    def tearDown(self):
        del self.policy
        del self.basisList
        del self.gaussFunction
        del self.gaussKernel
        del self.kImageIn
        del self.kFunctor

    def doDeltaFunction(self, dX, dY, scaling, background, imsize=50):
        # template image with a single hot pixel in the exact center
        tmi = afwImage.MaskedImageF(imsize, imsize)
        tmi.set(0, 0x0, 1e-4)
        # central pixel
        cpix = int(imsize/2)
        # where is the hot pixel?
        tmi.set(cpix, cpix, (1, 0x0, 1))

        # science image with a potentially offset and scaled hot pixel
        smi = afwImage.MaskedImageF(imsize, imsize)
        smi.set(0, 0x0, 1e-4)
        xpix = cpix + dX
        ypix = cpix + dY
        # where is the hot pixel?
        smi.set(xpix, ypix, (scaling, 0x0, scaling))
        # any additional background
        smi += background

        # estimate of the variance
        var  = afwImage.MaskedImageF(smi, True)
        var -= tmi

        # accepts : imageToConvolve, imageToNotConvolve
        self.kFunctor.apply(tmi.getImage(), smi.getImage(), var.getVariance(), self.policy)

        kb     = self.kFunctor.getSolution()
        kSoln  = kb.first
        bgSoln = kb.second
        kImage = afwImage.ImageD(self.kCols, self.kRows)
        kSum   = kSoln.computeImage(kImage, False)

        # make sure the correct background and scaling have been determined
        self.assertAlmostEqual(bgSoln, background, 5)
        self.assertAlmostEqual(kSum, scaling, 5)

        # make sure the delta function is in the right place
        xpix = kSoln.getCtrX() - dX
        ypix = kSoln.getCtrY() - dY
        for j in range(kImage.getHeight()):
            for i in range(kImage.getWidth()):

                if i == xpix and j == ypix:
                    self.assertAlmostEqual(kImage.get(i, j), scaling, 5)
                else:
                    self.assertAlmostEqual(kImage.get(i, j), 0., 5)
        
        
        
    def doGaussian(self, scaling, background, imsize=50):
        # NOTE : the size of these images have to be bigger
        #        size you lose pixels due to the convolution with the gaussian
        #        so adjust the size a bit to compensate 
        imsize += self.gSize
        
        # template image with a single hot pixel in the exact center
        tmi = afwImage.MaskedImageF(imsize, imsize)
        tmi.set(0, 0x0, 1e-4)
        # central pixel
        cpix = int(imsize/2)
        # where is the hot pixel?
        tmi.set(cpix, cpix, (1, 0x0, 1))
        
        # science image with a potentially offset and scaled hot pixel
        smi = afwImage.MaskedImageF(imsize, imsize)
        smi.set(0, 0x0, 1e-4)
        xpix = cpix
        ypix = cpix
        # where is the hot pixel?
        smi.set(xpix, ypix, (scaling, 0x0, scaling))
        # convolve with gaussian
        cmi = afwImage.MaskedImageF(imsize, imsize)
        afwMath.convolve(cmi, smi, self.gaussKernel, False)
        # this will adjust the kernel sum a bit
        # lose some at the outskirts of the kernel
        stats = afwMath.makeStatistics(cmi, afwMath.SUM)
        cscaling = stats.getValue(afwMath.SUM)
        # any additional background
        cmi += background

        # grab only the non-masked subregion
        bbox     = afwImage.BBox(afwImage.PointI(self.gaussKernel.getCtrX(),
                                                 self.gaussKernel.getCtrY()) ,
                                 afwImage.PointI(imsize - (self.gaussKernel.getWidth()  -
                                                           self.gaussKernel.getCtrX()),
                                                 imsize - (self.gaussKernel.getHeight() -
                                                           self.gaussKernel.getCtrY())))
                                 
        tmi2     = afwImage.MaskedImageF(tmi, bbox)
        cmi2     = afwImage.MaskedImageF(cmi, bbox)

        # make sure its a valid subregion!
        mask     = cmi2.getMask()
        for j in range(mask.getHeight()):
            for i in range(mask.getWidth()):
                self.assertEqual(mask.get(i, j), 0)
                
        # estimate of the variance
        var  = afwImage.MaskedImageF(cmi2, True)
        var -= tmi2

        # accepts : imageToConvolve, imageToNotConvolve
        self.kFunctor.apply(tmi2.getImage(), cmi2.getImage(), var.getVariance(), self.policy)

        kb        = self.kFunctor.getSolution()
        kSoln     = kb.first
        bgSoln    = kb.second
        kImageOut = afwImage.ImageD(self.kCols, self.kRows)
        kSum      = kSoln.computeImage(kImageOut, False)

        # make sure the correct background and scaling have been determined
        self.assertAlmostEqual(bgSoln, background, 4)
        self.assertAlmostEqual(kSum, cscaling, 4)

        # make sure the derived kernel looks like the input kernel
        for j in range(kImageOut.getHeight()):
            for i in range(kImageOut.getWidth()):

                # once we start to add in a background, the outer
                # portions of the kernel start to get a bit noisy.
                #
                # print i, j, self.kImageIn.get(i,j), kImageOut.get(i, j),
                # kImageOut.get(i, j)/self.kImageIn.get(i,j)
                #
                # however, where the power is, the results are the
                # same

                if self.kImageIn.get(i, j) > 0.01:
                    self.assertAlmostEqual(kImageOut.get(i, j)/self.kImageIn.get(i, j), scaling, 4)
                elif self.kImageIn.get(i, j) > 0.001:
                    self.assertAlmostEqual(kImageOut.get(i, j)/self.kImageIn.get(i, j), scaling, 3)


            
    def testDeltaFunction(self):
        # hot central pixel 
        self.doDeltaFunction(0, 0, 1, 0)
        self.doDeltaFunction(0, 0, 7, 0)
        self.doDeltaFunction(0, 0, 0.375, 0)

        # hot central pixel with background
        self.doDeltaFunction(0, 0, 1, 100)
        self.doDeltaFunction(0, 0, 1, 0.391)
        self.doDeltaFunction(0, 0, 1, -17.9)

        # mixture of hot central pixel, scaling, and background
        self.doDeltaFunction(0, 0, 0.735, 14.5)
        self.doDeltaFunction(0, 0, 12.20, 14.6)
        self.doDeltaFunction(0, 0, 0.735, -12.1)
        self.doDeltaFunction(0, 0, 12.20, -12.1)

        # offset delta function
        self.doDeltaFunction(-3, 2, 1, 0)
        self.doDeltaFunction(1, -2, 1, 0)

        # offset, scaling, and background
        self.doDeltaFunction(-3, 2, 10,  0)
        self.doDeltaFunction(-3, 2, 0.1, 0)
        self.doDeltaFunction(-4, 1, 1, +10)
        self.doDeltaFunction(-4, 1, 1, -10)
        self.doDeltaFunction(-4, 1, 0.1, +10)
        self.doDeltaFunction(-4, 1, 10,  -10)
        
        
    def testGaussian(self):
        # different scalings
        self.doGaussian(1, 0)
        self.doGaussian(7, 0)
        self.doGaussian(0.375, 0)

        # different backgrounds
        self.doGaussian(1, 3.3)
        self.doGaussian(1, 0.18)
        self.doGaussian(1, -12.1)

        # scaling and background
        self.doGaussian(0.735, 14.5)
        self.doGaussian(12.20, 14.6)
        self.doGaussian(0.735, -12.1)
        self.doGaussian(12.20, -12.1)

        
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
