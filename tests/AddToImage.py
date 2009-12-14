#!/usr/bin/env python
import unittest
import lsst.utils.tests as tests
import lsst.ip.diffim as ipDiffim
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage

class DiffimTestCases(unittest.TestCase):
    def setUp(self):
        self.value    = 100.
        self.function = afwMath.PolynomialFunction2D(2)
        nParams       = self.function.getNParameters()
        params        = [i for i in range(nParams)]
        self.function.setParameters(params)
    
    def tearDown(self):
        del self.value
        del self.function

    def testDouble(self):
        img = afwImage.ImageF(10,10,0)
        img += self.value
        for j in range(img.getHeight()):
            for i in range(img.getWidth()):
                self.assertEqual(img.get(i,j), self.value)

    def testFunction(self):
        img = afwImage.ImageF(10,10,0)
        img += self.function
        for j in range(img.getHeight()):
            for i in range(img.getWidth()):
                self.assertEqual(img.get(i,j), self.function(i,j))

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
    run(True)
    
        


     
