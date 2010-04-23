#!/usr/bin/env python
import os, sys
import unittest
import lsst.utils.tests as tests

import eups
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.pex.logging as pexLog
import lsst.afw.display.ds9 as ds9

diffimDir    = eups.productDir('ip_diffim')
diffimPolicy = os.path.join(diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')

pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)
display = True
class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = ipDiffim.generateDefaultPolicy(diffimPolicy)
        self.kSize = 11
        self.policy.set('kernelRows', self.kSize)
        self.policy.set('kernelCols', self.kSize)
        self.nCell = 5
        self.sizeCell = 64
        self.policy.set("sizeCellX", self.sizeCell)
        self.policy.set("sizeCellY", self.sizeCell)

        # we want an image that is nCell by nCell wide.
        # after convolution by the first gaussian
        # and after convolution by the spatial model

        # this sets the final extent of each convolved delta function
        self.gaussKernelWidth   = self.sizeCell//2
        # this sets how many pixels are correlated in tmi and smi
        #self.spatialKernelWidth = 2 * self.kSize + 1
        self.spatialKernelWidth = self.kSize 
        
        border = (self.gaussKernelWidth + self.spatialKernelWidth)//2
        
        totalSize = self.nCell * self.sizeCell + 2 * border
        tim = afwImage.ImageF(totalSize, totalSize)

        # matrix of delta functions
        for x in range(self.nCell):
            for y in range(self.nCell):
                # we need enough points such that a convolution with a
                # basis yields signal in all pixels; make a 3x3 grid.
                # maybe try a narrow gaussian as well?
                tim.set(x * self.sizeCell + self.sizeCell // 2 + border - 1,
                        y * self.sizeCell + self.sizeCell // 2 + border - 1,
                        1000)
        gaussFunction = afwMath.GaussianFunction2D(1.0, 1.0)
        gaussKernel   = afwMath.AnalyticKernel(self.gaussKernelWidth,
                                               self.gaussKernelWidth,
                                               gaussFunction)
        cim = afwImage.ImageF(tim.getDimensions())
        afwMath.convolve(cim, tim, gaussKernel, True)
        p0, p1   = diffimTools.getConvolvedImageLimits(gaussKernel, cim)
        bbox     = afwImage.BBox(p0, p1)

        # lets set the variance equal to the science image values
        # you have to make sure there are no negative values!
        vmin = afwMath.makeStatistics(cim, afwMath.MIN).getValue()
        self.assertTrue(vmin >= 0.0)
        tim  = afwImage.ImageF(cim, bbox)
        mask = afwImage.MaskU(tim.getDimensions())
        self.ti = afwImage.MaskedImageF(tim, mask, tim)
                
    def makeKernelCellSet(self, ti, si):
        kernelCellSet = afwMath.SpatialCellSet(afwImage.BBox(afwImage.PointI(0, 0),
                                                             self.sizeCell * self.nCell,
                                                             self.sizeCell * self.nCell),
                                               self.sizeCell,
                                               self.sizeCell)
        stampHalfWidth = 1 * self.kSize
        for x in range(self.nCell):
            for y in range(self.nCell):
                xCoord = x * self.sizeCell + self.sizeCell // 2
                yCoord = y * self.sizeCell + self.sizeCell // 2
                p0 = afwImage.PointI(xCoord - stampHalfWidth,
                                     yCoord - stampHalfWidth)
                p1 = afwImage.PointI(xCoord + stampHalfWidth,
                                     yCoord + stampHalfWidth)
                bbox = afwImage.BBox(p0, p1)
                tsi = afwImage.MaskedImageF(ti, bbox)
                ssi = afwImage.MaskedImageF(si, bbox)

                kc = ipDiffim.makeKernelCandidate(xCoord, yCoord, tsi, ssi, self.policy)
                kernelCellSet.insertCandidate(kc)

                if display:
                    ds9.mtv(tsi, frame=4)
                    ds9.mtv(ssi, frame=5)
        
        return kernelCellSet

    def tearDown(self):
        del self.policy
        del self.ti

    def addNoise(self, mi):
        img       = mi.getImage()
        seed      = int(10. * afwMath.makeStatistics(mi.getImage(), afwMath.MAX).getValue() + 1)
        print seed
        rdm       = afwMath.Random(afwMath.Random.MT19937, seed)
        rdmImage  = img.Factory(img.getDimensions())
        afwMath.randomGaussianImage(rdmImage, rdm)
        img      += rdmImage
        return afwMath.makeStatistics(rdmImage, afwMath.MEAN).getValue(afwMath.MEAN)

    def xtestSolve(self):
        pass

    def testGaussianField(self):
        #self.runGaussianField(0)
        self.runGaussianField(1)
        #self.runGaussianField(2)
        
    def xrunGaussianFieldNan(self, order, nGauss = 3, wGauss = [1.5, 2.5, 3.5]):
        self.policy.set('kernelBasisSet', 'alard-lupton')
        self.policy.set('spatialKernelOrder', order)
        self.policy.set('usePcaForSpatialKernel', False)
        
        self.policy.set('alardNGauss', nGauss)
        self.policy.set('alardSigGauss', wGauss[0])
        self.policy.set('alardDegGauss', 0)
        for i in range(1, len(wGauss)):
            self.policy.add('alardSigGauss', wGauss[i])
            self.policy.add('alardDegGauss', 0)

        # set up kernel coefficients
        polyFunc = afwMath.PolynomialFunction2D(order)
        nParams  = len(polyFunc.getParameters())
        kCoeffs  = []
        for i in range(nGauss):
            kCoeffs.append([])
            for j in range(nParams):
                kCoeffs[i].append(0.01 * (1.5 * (-1)**j + i)) # this number should be different for each basis
        
        # set up basis list
        # make kernel a little bigger than our fit so we avoid edge effects
        basisList = afwMath.KernelList()
        for i in range(nGauss):
            gaussFunction = afwMath.GaussianFunction2D(wGauss[i], wGauss[i])
            gaussKernel   = afwMath.AnalyticKernel(self.spatialKernelWidth,
                                                   self.spatialKernelWidth,
                                                   gaussFunction)
            basisList.append(gaussKernel)

        # should we normalize list so all the powers in the first one,
        # but it doesn't vary spatially?

        #
        kl1 = afwMath.KernelList()
        kl1.append(basisList[0])
        k1 = afwMath.LinearCombinationKernel(kl1, polyFunc)
        k1.setSpatialParameters([kCoeffs[0],])
        
        kl2 = afwMath.KernelList()
        kl2.append(basisList[1])
        k2 = afwMath.LinearCombinationKernel(kl2, polyFunc)
        k2.setSpatialParameters([kCoeffs[1],])

        kl3 = afwMath.KernelList()
        kl3.append(basisList[2])
        k3 = afwMath.LinearCombinationKernel(kl3, polyFunc)
        k3.setSpatialParameters([kCoeffs[2],])

        c1 = afwImage.MaskedImageF(self.ti.getDimensions())
        afwMath.convolve(c1, self.ti, k1, True)
        c2 = afwImage.MaskedImageF(self.ti.getDimensions())
        afwMath.convolve(c2, self.ti, k2, True)
        c3 = afwImage.MaskedImageF(self.ti.getDimensions())
        afwMath.convolve(c3, self.ti, k3, True)

        p0, p1 = diffimTools.getConvolvedImageLimits(k2, c2)
        bbox = afwImage.BBox(p0, p1)
        foo = afwImage.MaskedImageF(c2, bbox)
        for x in range(foo.getWidth()):
            for y in range(foo.getHeight()):
                if str(foo.getImage().get(x,y)) == 'nan':
                    print x, y, foo.getImage().get(x,y), k2.getSpatialFunction(0)(x,y)

#        for x in range(294-self.spatialKernelWidth//2, 294+self.spatialKernelWidth//2 + 1):
#            for y in range(29-self.spatialKernelWidth//2, 29+self.spatialKernelWidth//2 + 1):
#                ki = afwImage.ImageD(k3.getDimensions())
#                k3.computeImage(ki, x, y, True)
#                print x, y, c3.getImage().get(x,y),
#                print afwMath.makeStatistics(ki, afwMath.MEAN).getValue(afwMath.MEAN),
#                print k3.getSpatialFunction(0)(x,y)


        c1.writeFits('c1.fits')
        c2.writeFits('c2.fits')
        c3.writeFits('c3.fits')

        # make the full spatial kernel
        kernelOrig = afwMath.LinearCombinationKernel(basisList, polyFunc)
        kernelOrig.setSpatialParameters(kCoeffs)
        si = afwImage.MaskedImageF(self.ti.getDimensions())
        afwMath.convolve(si, self.ti, kernelOrig, True)  # True = Constant sum across the image

        # get the good subregion
        p0, p1   = diffimTools.getConvolvedImageLimits(kernelOrig, si)
        bbox     = afwImage.BBox(p0, p1)
        si       = afwImage.MaskedImageF(si, bbox)
        ti       = afwImage.MaskedImageF(self.ti, bbox, True)  # copy pixels since i modify them

        si.setXY0(afwImage.PointI(0,0))
        ti.setXY0(afwImage.PointI(0,0))

        print kCoeffs
        #diffimTools.displayKernelList(basisList)
        diffimTools.displaySpatialKernelMosaic(kernelOrig, si.getWidth(), si.getHeight(),
                                               frame = len(basisList),
                                               doNorm = True)

        # add noise?
        self.addNoise(ti)
        self.addNoise(si)

        if display:
            #ds9.mtv(self.ti, frame = 1)
            ds9.mtv(ti, frame = 4)
            ds9.mtv(si, frame = 5)
        sys.exit(1)
        # make candidates
        kernelCellSet = self.makeKernelCellSet(ti, si)

        # single kernel visitor
        basisListToFit = ipDiffim.makeKernelBasisList(self.policy)
        diffimTools.displayKernelList(basisListToFit, 4)
        
        bskv = ipDiffim.BuildSingleKernelVisitorF(basisListToFit, self.policy)

        # visit candidates by hand
        frame = 0 
        for cell in kernelCellSet.getCellList():
            for cand in cell.begin(False): # False = include bad candidates
                cand = ipDiffim.cast_KernelCandidateF(cand)
                bskv.processCandidate(cand)

                ds9.mtv(cand.getMiToConvolvePtr(), frame = frame); frame += 1
                ds9.mtv(cand.getMiToNotConvolvePtr(), frame = frame); frame += 1
                ds9.mtv(cand.getKernelImage(ipDiffim.KernelCandidateF.ORIG), frame = frame); frame += 1
                #ds9.mtv(cand.getDifferenceImage(ipDiffim.KernelCandidateF.ORIG), frame = frame); frame += 1

                ki = afwImage.ImageD(kernelOrig.getDimensions())
                kernelOrig.computeImage(ki, cand.getXCenter(), cand.getYCenter(), False)
                ds9.mtv(ki, frame = frame); frame += 1






 
    def runGaussianField(self, order, nGauss = 3, wGauss = [1.5, 2.5, 3.5]):
        self.policy.set('kernelBasisSet', 'alard-lupton')
        self.policy.set('spatialKernelOrder', order)
        self.policy.set('usePcaForSpatialKernel', False)
        
        self.policy.set('alardNGauss', nGauss)
        self.policy.set('alardSigGauss', wGauss[0])
        self.policy.set('alardDegGauss', 0)
        for i in range(1, len(wGauss)):
            self.policy.add('alardSigGauss', wGauss[i])
            self.policy.add('alardDegGauss', 0)

        
        # set up basis list
        # make kernel a little bigger than our fit so we avoid edge effects
        basisList = afwMath.KernelList()
        for i in range(nGauss):
            gaussFunction = afwMath.GaussianFunction2D(wGauss[i], wGauss[i])
            gaussKernel   = afwMath.AnalyticKernel(self.spatialKernelWidth,
                                                   self.spatialKernelWidth,
                                                   gaussFunction)
            basisList.append(gaussKernel)

        # first kernel has no spatial variation; the others have no kernel sum
        basisList = afwMath.KernelList(ipDiffim.renormalizeKernelList(basisList))

        # set up kernel coefficients
        polyFunc = afwMath.PolynomialFunction2D(order)
        nParams  = len(polyFunc.getParameters())
        kCoeffs  = []
        # first one does not vary spatially
        kCoeffs.append([])
        kCoeffs[0].append(1.)
        for i in range(1, nParams):
            kCoeffs[0].append(0.)
            
        for i in range(1, nGauss):
            kCoeffs.append([])
            for j in range(nParams):
                kCoeffs[i].append(0.001 * (1.5 * (-1)**j + i))

        # estimate of the variance comes from main gaussian
        var = afwImage.ImageF(self.ti.getDimensions())
        afwMath.convolve(var, self.ti.getImage(), basisList[0], False)

        # make the full spatial kernel
        kernelOrig = afwMath.LinearCombinationKernel(basisList, polyFunc)
        kernelOrig.setSpatialParameters(kCoeffs)
        im = afwImage.ImageF(self.ti.getDimensions())
        afwMath.convolve(im, self.ti.getImage(), kernelOrig, False)

        # get the good subregion
        p0, p1   = diffimTools.getConvolvedImageLimits(kernelOrig, im)
        bbox     = afwImage.BBox(p0, p1)
        sim      = afwImage.ImageF(im, bbox)
        var      = afwImage.ImageF(var, bbox)
        mask     = afwImage.MaskU(sim.getDimensions())
        si       = afwImage.MaskedImageF(sim, mask, var)
        ti       = afwImage.MaskedImageF(self.ti, bbox, True)  # copy pixels since i modify them

        si.setXY0(afwImage.PointI(0,0))
        ti.setXY0(afwImage.PointI(0,0))

        diffimTools.displaySpatialKernelMosaic(kernelOrig, si.getWidth(), si.getHeight(),
                                               frame = len(basisList),
                                               doNorm = True)

        # add noise
        #self.addNoise(ti)
        #self.addNoise(si)
        var = ti.getVariance(); var += 1.0
        var = si.getVariance(); var += 1.0

        if display:
            #ds9.mtv(self.ti, frame = 1)
            ds9.mtv(ti, frame = 4)
            ds9.mtv(si, frame = 5)

        # make candidates
        kernelCellSet = self.makeKernelCellSet(ti, si)

        # single kernel visitor
        basisListToFit = ipDiffim.makeKernelBasisList(self.policy)
        diffimTools.displayKernelList(basisListToFit, 4)
        
        bskv = ipDiffim.BuildSingleKernelVisitorF(basisListToFit, self.policy)

        # visit candidates by hand
        frame = 0 
        for cell in kernelCellSet.getCellList():
            for cand in cell.begin(False): # False = include bad candidates
                cand = ipDiffim.cast_KernelCandidateF(cand)
                bskv.processCandidate(cand)

                #ds9.mtv(cand.getMiToConvolvePtr(), frame = frame); frame += 1
                #ds9.mtv(cand.getMiToNotConvolvePtr(), frame = frame); frame += 1
                #ds9.mtv(cand.getKernelImage(ipDiffim.KernelCandidateF.ORIG), frame = frame); frame += 1
                #ds9.mtv(cand.getDifferenceImage(ipDiffim.KernelCandidateF.ORIG), frame = frame); frame += 1
                kernelComp = cand.getKernel(ipDiffim.KernelCandidateF.ORIG)

                print kernelComp.getKernelParameters()
                print kernelOrig.getSpatialFunction(0)(cand.getXCenter(), cand.getYCenter()),
                print kernelOrig.getSpatialFunction(1)(cand.getXCenter(), cand.getYCenter()),
                print kernelOrig.getSpatialFunction(2)(cand.getXCenter(), cand.getYCenter())
                
                ki = afwImage.ImageD(kernelOrig.getDimensions())
                kernelOrig.computeImage(ki, cand.getXCenter(), cand.getYCenter(), False)
                #ds9.mtv(ki, frame = frame); frame += 1
                cand.getMiToConvolvePtr().writeFits('t%d.fits' % (frame))
                cand.getMiToNotConvolvePtr().writeFits('s%d.fits' % (frame))
                cand.getDifferenceImage(ipDiffim.KernelCandidateF.ORIG).writeFits('d%d.fits' % (frame))
                frame += 1
 
    def xtestDeltaFunction(self):
        pass
    
    def xtestVisit(self):
        pass

        
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
