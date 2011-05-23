#!/usr/bin/env python
import unittest
import lsst.utils.tests as tests

class DiffimTestCases(unittest.TestCase):
    # To be written
    def setUp(self):
        pass
    def tearDown(self):
        pass

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
