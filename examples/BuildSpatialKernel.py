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

import numpy
import pylab

pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)
display = False
class DiffimTestCases(unittest.TestCase):
    
    def setUp(self):
        self.policy = ipDiffim.createDefaultPolicy()
        self.kSize = 19
        self.policy.set("kernelSize", self.kSize)
        self.sizeCell = 64
        self.policy.set("sizeCellX", self.sizeCell)
        self.policy.set("sizeCellY", self.sizeCell)
        self.policy.set('fitForBackground', False)

    def tearDown(self):
        del self.policy

    def xtestSolve(self):
        pass

    def xtestDeltaFunction(self):
        self.policy.set('kernelBasisSet', 'delta-function')
        self.policy.set('useRegularization', False)
        kset = self.runGaussianField(0, title = 'DF')
        diffimTools.displayKernelMosaic(kset, 5)

    def xtestAlardLupton1(self):
        self.policy.set('kernelBasisSet', 'alard-lupton')
        self.policy.set('useRegularization', False)

        nGauss = 1
        wGauss = [2.5,]
        self.policy.set('alardNGauss', nGauss)
        self.policy.set('alardSigGauss', wGauss[0])
        self.policy.set('alardDegGauss', 0)
        for i in range(1, len(wGauss)):
            self.policy.add('alardSigGauss', wGauss[i])
            self.policy.add('alardDegGauss', 0)

        kset = self.runGaussianField(0, title = 'AL1')
        diffimTools.displayKernelMosaic(kset, 6)


    def xtestAlardLupton2(self):
        self.policy.set('kernelBasisSet', 'alard-lupton')
        self.policy.set('useRegularization', False)

        wGauss = [1.0, 2.5, 5.0]
        dGauss = [4, 3, 2]
        self.policy.set('alardSigGauss', wGauss[0])
        self.policy.set('alardDegGauss', dGauss[0])
        for i in range(1, len(wGauss)):
            self.policy.add('alardSigGauss', wGauss[i])
            self.policy.add('alardDegGauss', dGauss[i])

        kset = self.runGaussianField(0, title = 'AL2')
        diffimTools.displayKernelMosaic(kset, 7)

    def xtestAlardLupton3(self):
        self.policy.set('kernelBasisSet', 'alard-lupton')
        self.policy.set('useRegularization', False)

        wGauss = [2.5, 1.0, 5.0]
        dGauss = [3, 4, 2]
        self.policy.set('alardSigGauss', wGauss[0])
        self.policy.set('alardDegGauss', dGauss[0])
        for i in range(1, len(wGauss)):
            self.policy.add('alardSigGauss', wGauss[i])
            self.policy.add('alardDegGauss', dGauss[i])

        kset = self.runGaussianField(0, title = 'AL3')
        diffimTools.displayKernelMosaic(kset, 8)

    def xtestRegularization1(self):
        self.policy.set('kernelBasisSet', 'delta-function')
        self.policy.set('useRegularization', True)
        self.policy.set('lambdaValue', 100000.)
        kset = self.runGaussianField(0, title = 'DFr1')
        #diffimTools.displayKernelMosaic(kset, 9)

    def testRegularization2(self):
        self.policy.set('kernelBasisSet', 'delta-function')
        self.policy.set('useRegularization', True)
        self.policy.set('lambdaValue', 1000000000.)
        kset = self.runGaussianField(0, title = 'DFr2')
        #diffimTools.displayKernelMosaic(kset, 10)

    def xtestShow(self):
        pylab.show()  
    
    def runGaussianField(self, order, title, writeIm = True, writeKern = True):
        self.policy.set('spatialKernelOrder', order)
        
        # set up basis list
        nGauss = 1
        wGauss = [2.5,]
        #nGauss = self.policy.get('alardNGauss')
        #wGauss = self.policy.getDoubleArray('alardSigGauss')
        basisList = afwMath.KernelList()
        for i in range(nGauss):
            gaussFunction = afwMath.GaussianFunction2D(wGauss[i], wGauss[i])
            gaussKernel   = afwMath.AnalyticKernel(self.kSize, self.kSize, gaussFunction)
            basisList.append(gaussKernel)

        # first kernel has no spatial variation; the others have no kernel sum
        basisList = afwMath.KernelList(ipDiffim.renormalizeKernelList(basisList))

        tMi, sMi, sKernel, kernelCellSet = diffimTools.makeFakeKernelSet(self.policy,
                                                                         basisList,
                                                                         nCell = 50,
                                                                         addNoise=True,
                                                                         bgValue = 1.e2)

        # single kernel visitor
        basisListToFit = ipDiffim.makeKernelBasisList(self.policy)
        
        if self.policy.get("useRegularization"):
            hMat = ipDiffim.makeRegularizationMatrix(self.policy)
            bskv = ipDiffim.BuildSingleKernelVisitorF(basisListToFit, self.policy, hMat)
        else:
            bskv = ipDiffim.BuildSingleKernelVisitorF(basisListToFit, self.policy)

        # visit candidates by hand
        frame = 7
        for cell in kernelCellSet.getCellList():
            for cand in cell.begin(False): # False = include bad candidates
                cand = ipDiffim.cast_KernelCandidateF(cand)
                bskv.processCandidate(cand)
                if writeIm:
                    (cand.getMiToConvolvePtr()).writeFits('%s/tmi/tmi_%s_%d.fits' %
                                                          (title, title, cand.getId()))
                    (cand.getMiToNotConvolvePtr()).writeFits('%s/smi/smi_%s_%d.fits' %
                                                             (title, title, cand.getId()))
                if writeKern:
                    (cand.getImage()).writeFits('%s/k/k_%s_%d.fits' % (title, title, cand.getId()))

                # this thing is a memory hog!  see if this helps.
                del cand

        #s1o, s1i = self.accumulateDiffimStats(kernelCellSet, title)
        #s2o, s2i = self.compareNeighbors(kernelCellSet, title)
        #print '# CAW', title, s1i.mean(), s1i.std(), s1o.mean(), s1o.std(), \
        #      s2i.mean(), s2i.std(), s2o.mean(), s2o.std()

        return kernelCellSet

    def compareNeighbors(self, kernelCellSet, title):
        # lotsa permutations here...
        allStats1 = []
        allStats2 = []
        for cell1 in kernelCellSet.getCellList():
            for cand1 in cell1.begin(False):
                cand1   = ipDiffim.cast_KernelCandidateF(cand1)
                kernel1 = cand1.getKernel(ipDiffim.KernelCandidateF.ORIG)

                for cell2 in kernelCellSet.getCellList():
                    for cand2 in cell2.begin(False):
                        cand2  = ipDiffim.cast_KernelCandidateF(cand2)
                        if cand1.getId() == cand2.getId():
                            continue

                        diffimFull = cand2.getDifferenceImage(kernel1, 0.0)
                        itcvFull   = cand2.getMiToConvolvePtr().getVariance()
                        itncvFull  = cand2.getMiToNotConvolvePtr().getVariance()

                        bbox = kernel1.shrinkBBox(diffim.getBBox(afwImage.LOCAL))
                        diffim = afwImage.MaskedImageF(diffimFull, bbox, afwImage.LOCAL)
                        itcv   = afwImage.ImageF(itcvFull, bbox)
                        itncv  = afwImage.ImageF(itncvFull, bbox)
                
                        pval   = diffimTools.vectorFromImage(diffim.getImage())
                        vval   = diffimTools.vectorFromImage(diffim.getVariance())
                        itcv   = diffimTools.vectorFromImage(itcv)
                        itncv  = diffimTools.vectorFromImage(itncv)
                        
                        allStats1.append( pval / numpy.sqrt(vval) )
                        allStats2.append( pval / numpy.sqrt(itcv + itncv) )
                        


        allStats1 = numpy.ravel(allStats1) 
        pylab.figure()
        n, b, p = pylab.hist(allStats1, bins=50, normed=True)
        pylab.title('%s : For neighbors (output variance)' % (title))
        pylab.xlabel('N Sigma')
        pylab.ylabel('N Pix')
        
        allStats2 = numpy.ravel(allStats2) 
        pylab.figure()
        n, b, p = pylab.hist(allStats2, bins=50, normed=True)
        pylab.title('%s : For neighbors (input variance)' % (title))
        pylab.xlabel('N Sigma')
        pylab.ylabel('N Pix')
        return allStats1, allStats2

    def accumulateDiffimStats(self, kernelCellSet, title):
        allStats1 = []
        allStats2 = []
        for cell in kernelCellSet.getCellList():
            for cand in cell.begin(False): # False = include bad candidates
                cand = ipDiffim.cast_KernelCandidateF(cand)
                kernel = cand.getKernel(ipDiffim.KernelCandidateF.ORIG)
                diffim = cand.getDifferenceImage(ipDiffim.KernelCandidateF.ORIG)
                itcv   = cand.getMiToConvolvePtr().getVariance()
                itncv  = cand.getMiToNotConvolvePtr().getVariance()
                
                
                bbox = kernel.shrinkBBox(diffim.getBBox(afwImage.LOCAL))
                diffim = afwImage.MaskedImageF(diffim, bbox, afwImage.LOCAL)
                itcv   = afwImage.ImageF(itcv, bbox)
                itncv  = afwImage.ImageF(itncv, bbox)
                
                pval   = diffimTools.vectorFromImage(diffim.getImage())
                vval   = diffimTools.vectorFromImage(diffim.getVariance())
                itcv   = diffimTools.vectorFromImage(itcv)
                itncv  = diffimTools.vectorFromImage(itncv)

                allStats1.append( pval / numpy.sqrt(vval) )
                allStats2.append( pval / numpy.sqrt(itcv + itncv) )

        allStats1 = numpy.ravel(allStats1) 
        pylab.figure()
        n, b, p = pylab.hist(allStats1, bins=50, normed=True)
        pylab.title('%s : For itself (output variance)' % (title))
        pylab.xlabel('N Sigma')
        pylab.ylabel('N Pix')
        
        allStats2 = numpy.ravel(allStats2) 
        pylab.figure()
        n, b, p = pylab.hist(allStats2, bins=50, normed=True)
        pylab.title('%s : For itself (input variance)' % (title))
        pylab.xlabel('N Sigma')
        pylab.ylabel('N Pix')
        return allStats1, allStats2
        

 
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
