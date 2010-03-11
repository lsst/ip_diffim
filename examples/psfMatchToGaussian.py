#!/usr/bin/env python
import unittest
import eups
import os
import sys
import lsst.utils.tests as tests
import lsst.ip.diffim as ipDiffim
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)

class DiffimTestCases(unittest.TestCase):
    def setUp(self):
        self.diffimDir    = eups.productDir('ip_diffim')
        self.diffimPolicy = os.path.join(self.diffimDir, 'pipeline', 'ImageSubtractStageDictionary.paf')
        self.policy       = ipDiffim.generateDefaultPolicy(self.diffimPolicy, fwhm=5)
        
        self.defDataDir = eups.productDir('afwdata')
        if self.defDataDir:
            defSciencePath = os.path.join(self.defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                           "v5-e0-c011-a00.sci")
            self.scienceImage   = afwImage.MaskedImageF(defSciencePath)

            # BY HAND STARS IN THE IMAGE
            fp1  = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(414, 1755), 1, 1))
            fp2  = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(134, 1956), 1, 1))
            fp3  = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(234, 1594), 1, 1))
            fp4  = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(470, 1654), 1, 1))
            fp5  = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(42, 1520), 1, 1))
            fp6  = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(283, 1340), 1, 1))
            fp8  = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(448, 1254), 1, 1))
            fp9  = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(114, 933), 1, 1))
            fp10 = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(104, 706), 1, 1))
            fp11 = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(50, 637), 1, 1))
            fp12 = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(341, 762), 1, 1))
            fp13 = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(407, 579), 1, 1))
            fp14 = afwDetection.Footprint(afwImage.BBox(afwImage.PointI(236, 324), 1, 1))

            self.footprints = [afwDetection.growFootprint(fp1, 30, False),
                               afwDetection.growFootprint(fp2, 30, False),
                               afwDetection.growFootprint(fp3, 30, False),
                               afwDetection.growFootprint(fp4, 30, False),
                               afwDetection.growFootprint(fp5, 30, False),
                               afwDetection.growFootprint(fp6, 30, False),
                               afwDetection.growFootprint(fp8, 30, False),
                               afwDetection.growFootprint(fp9, 30, False),
                               afwDetection.growFootprint(fp10, 30, False),
                               afwDetection.growFootprint(fp11, 30, False),
                               afwDetection.growFootprint(fp12, 30, False),
                               afwDetection.growFootprint(fp13, 30, False),
                               afwDetection.growFootprint(fp14, 30, False)]

    def testGaussian(self, sigGauss=5.):
        if not self.defDataDir:
            print >> sys.stderr, "Warning: afwdata is not set up; not running PsfMatchToGaussian.py"
            return

        # NOTE - you need to subtract off background from the image
        # you run detection on.  Here it is the template.
        diffimTools.backgroundSubtract(self.policy, [self.scienceImage,])

        ipDiffim.makePsfMatchingKernelToGaussian(self.scienceImage,
                                                 sigGauss,
                                                 self.policy,
                                                 self.footprints)
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
    
        
