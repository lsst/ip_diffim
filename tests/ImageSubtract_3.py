import unittest
import numpy as num
import pdb
import os

import eups
import lsst.pex.policy
import lsst.utils.tests   as tests
import lsst.afw.math      as afwMath
import lsst.afw.image     as afwImage
import lsst.afw.detection as afwDetection
import lsst.ip.diffim     as ipDiffim

import lsst.ip.diffim.diffimTools as ipDiffimTools
import lsst.ip.diffim.runPca      as runPca

try:
    type(verbosity)
except NameError:
    verbosity = 5
logging.Trace.setVerbosity('lsst.ip.diffim', verbosity)

ipDiffimDir = eups.productDir("ip_diffim")
if not ipDiffimDir:
    raise RuntimeError("Could not get path to ip_diffim")

policyPath = os.path.join(ipDiffimDir, "pipeline", "ImageSubtractStageDictionary.paf")
policy = lsst.pex.policy.Policy.createPolicy(policyPath)

class fitFunctionUT(unittest.TestCase):
    def testOrder2(self, width1=1.0, width2=2.7, npts=100):
        # (order + 1) * (order + 2) / 2 coefficients
        inPars = num.array( (1, 2, 3, 4, 5.5, 7.9) )
        f      = afwMath.PolynomialFunction2D(2)
        f.setParameters(inPars)

        values = num.zeros(npts)
        errors = num.zeros(npts)
        cols   = num.zeros(npts)
        rows   = num.zeros(npts)

        idx = 0
        for i in num.arange(-5, 5, 1):
            for j in num.arange(-5, 5, 1):
                values[idx] = f(i, j)
                errors[idx] = num.sqrt(values[idx])
                cols[idx]   = i
                rows[idx]   = j

                idx        += 1

        pars     = num.zeros(2)
        stepsize = 0.1 * num.ones(2)
        fit      = ipDiffimTools.fitFunction(afwMath.PolynomialFunction2D(2),
                                             values, errors,
                                             cols, rows,
                                             policy)
        for i in range(len(inPars)):
            self.assertAlmostEqual(inPars[i], fit.parameterList[i])
    
class rejectKernelSumOutliersUT(unittest.TestCase):
    def testDefault(self):
        pass
    
class createSpatialModelKernelCellsUT(unittest.TestCase):
    def testDefault(self):
        pass
    
class runPcaUT(unittest.TestCase):
    def testDefault(self):
        pass
    
class spatialModelKernelPcaUT(unittest.TestCase):
    def testDefault(self):
        pass
    
class spatialModelByPixelUT(unittest.TestCase):
    def testDefault(self):
        pass
    
class evaluateModelByPixelUT(unittest.TestCase):
    def testDefault(self):
        pass
    
class spatialModelByPcaUT(unittest.TestCase):
    def testDefault(self):
        pass
    
class evaluateModelByPcaUT(unittest.TestCase):
    def testDefault(self):
        pass
    

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []

    # diffimTools.py
    suites += unittest.makeSuite(fitFunctionUT)
    suites += unittest.makeSuite(rejectKernelSumOutliersUT)
    suites += unittest.makeSuite(createSpatialModelKernelCellsUT)

    # runPca.py
    suites += unittest.makeSuite(runPcaUT)

    # spatialModelKernelFit.py
    suites += unittest.makeSuite(spatialModelKernelPcaUT)
    suites += unittest.makeSuite(spatialModelByPixelUT)
    suites += unittest.makeSuite(evaluateModelByPixelUT)
    suites += unittest.makeSuite(spatialModelByPcaUT)
    suites += unittest.makeSuite(evaluateModelByPcaUT)
    
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    """Tests each unit of python code"""
    run(True)
