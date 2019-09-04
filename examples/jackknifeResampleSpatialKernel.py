#!/usr/bin/env python

# This file is part of ip_diffim.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import sys
import unittest
import lsst.utils.tests as tests

import lsst.utils
import lsst.afw.display as afwDisplay
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
from lsst.log import Log
import lsst.log.utils as logUtils
import lsst.pex.config as pexConfig

afwDisplay.setDefaultMaskTransparency(75)

verbosity = 2
logUtils.traceSetAt("ip.diffim", verbosity)

logger = Log.getLogger("ip.diffim.JackknifeResampleKernel")

display = False
writefits = False

try:
    defDataDir = lsst.utils.getPackageDir('afwdata')
except Exception:
    defDataDir = None

if defDataDir:
    defTemplatePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                   "v5-e0-c011-a00.sci.fits")
    defSciencePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                  "v26-e0-c011-a00.sci.fits")

# THIS IS "LEAVE-ONE-OUT" CROSS VALIDATION OF THE SPATIAL KERNEL


class DiffimTestCases(unittest.TestCase):

    def setUp(self):
        if not defDataDir:
            return

        self.configAL = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.configAL.kernel.name = "AL"
        self.subconfigAL = self.configAL.kernel.active

        self.configDF = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.configDF.kernel.name = "DF"
        self.subconfigDF = self.configDF.kernel.active

        self.configDFr = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.configDFr.kernel.name = "DF"
        self.subconfigDFr = self.configDFr.kernel.active

        self.subconfigDF.useRegularization = False
        self.subconfigDFr.useRegularization = True

        self.scienceExposure = afwImage.ExposureF(defSciencePath)
        self.templateExposure = afwImage.ExposureF(defTemplatePath)

        warper = afwMath.Warper.fromConfig(self.subconfigAL.warpingConfig)
        self.templateExposure = warper.warpExposure(self.scienceExposure.getWcs(), self.templateExposure,
                                                    destBBox=self.scienceExposure.getBBox())

        self.scienceMaskedImage = self.scienceExposure.getMaskedImage()
        self.templateMaskedImage = self.templateExposure.getMaskedImage()

        bgConfig = self.subconfigAL.afwBackgroundConfig
        diffimTools.backgroundSubtract(bgConfig, [self.templateMaskedImage,
                                                  self.scienceMaskedImage])

    def stats(self, cid, diffim, core=5):
        ps = pexConfig.makePropertySet(self.config)
        dStats = ipDiffim.ImageStatisticsF(ps)
        dStats.apply(diffim)
        logger.debug("Candidate %d : Residuals all (%d px): %.3f +/- %.3f",
                     cid, dStats.getNpix(), dStats.getMean(), dStats.getRms())

        dStats.apply(diffim, core)
        logger.debug("Candidate %d : Residuals core (%d px): %.3f +/- %.3f",
                     cid, dStats.getNpix(), dStats.getMean(), dStats.getRms())

    def assess(self, cand, kFn1, bgFn1, kFn2, bgFn2, frame0):
        tmi = cand.getTemplateMaskedImage()
        smi = cand.getScienceMaskedImage()

        im1 = afwImage.ImageD(kFn1.getDimensions())
        kFn1.computeImage(im1, False,
                          afwImage.indexToPosition(int(cand.getXCenter())),
                          afwImage.indexToPosition(int(cand.getYCenter())))
        fk1 = afwMath.FixedKernel(im1)
        bg1 = bgFn1(afwImage.indexToPosition(int(cand.getXCenter())),
                    afwImage.indexToPosition(int(cand.getYCenter())))
        d1 = ipDiffim.convolveAndSubtract(tmi, smi, fk1, bg1)

        ####

        im2 = afwImage.ImageD(kFn2.getDimensions())
        kFn2.computeImage(im2, False,
                          afwImage.indexToPosition(int(cand.getXCenter())),
                          afwImage.indexToPosition(int(cand.getYCenter())))
        fk2 = afwMath.FixedKernel(im2)
        bg2 = bgFn2(afwImage.indexToPosition(int(cand.getXCenter())),
                    afwImage.indexToPosition(int(cand.getYCenter())))
        d2 = ipDiffim.convolveAndSubtract(tmi, smi, fk2, bg2)

        if display:
            disp = afwDisplay.Display(frame=frame0)
            disp.mtv(tmi, title="Template Masked Image")
            disp.dot("Cand %d" % (cand.getId()), 0, 0)

            afwDisplay.Display(frame=frame0 + 1).mtv(smi, title="Science Masked Image")
            afwDisplay.Display(frame=frame0 + 2).mtv(im1, title="Masked Image: 1")
            afwDisplay.Display(frame=frame0 + 3).mtv(d1, title="Difference Image: 1")
            afwDisplay.Display(frame=frame0 + 4).mtv(im2, title="Masked Image: 2")
            afwDisplay.Display(frame=frame0 + 5).mtv(d2, title="Difference Image: 2")

        logger.debug("Full Spatial Model")
        self.stats(cand.getId(), d1)

        logger.debug("N-1 Spatial Model")
        self.stats(cand.getId(), d2)

    def setStatus(self, cellSet, cid, value):
        # ideally
        # cellSet.getCandidateById(id).setStatus(value)
        for cell in cellSet.getCellList():
            for cand in cell.begin(False):
                if (cand.getId() == cid):
                    cand.setStatus(value)
                    return cand

    def jackknifeResample(self, psfmatch, results):

        kernel = results.psfMatchingKernel
        bg = results.backgroundModel
        cellSet = results.kernelCellSet

        goodList = []
        for cell in cellSet.getCellList():
            print()
            for cand in cell.begin(False):
                if cand.getStatus() == afwMath.SpatialCellCandidate.GOOD:
                    goodList.append(cand.getId())
                else:
                    # This is so that UNKNOWNs are not processed
                    cand.setStatus(afwMath.SpatialCellCandidate.BAD)

        nStarPerCell = self.config.nStarPerCell
        ps = pexConfig.makePropertySet(self.config)
        for idx in range(len(goodList)):
            cid = goodList[idx]

            print()  # clear the screen
            logger.debug("Removing candidate %d", cid)

            cand = self.setStatus(cellSet, cid, afwMath.SpatialCellCandidate.BAD)

            # From _solve
            regionBBox = cellSet.getBBox()
            spatialkv = ipDiffim.BuildSpatialKernelVisitorF(kernel.getKernelList(), regionBBox, ps)
            cellSet.visitCandidates(spatialkv, nStarPerCell)
            spatialkv.solveLinearEquation()
            jkKernel, jkBg = spatialkv.getSolutionPair()

            # jkResults = psfmatch._solve(cellSet, kernel.getKernelList())
            # jkKernel  = jkResults[1]
            # jkBg      = jkResults[2]

            # lots of windows
            # self.assess(cand, kernel, bg, jkKernel, jkBg, 6*idx+1)

            # only 6 windows
            self.assess(cand, kernel, bg, jkKernel, jkBg, 1)

            self.setStatus(cellSet, cid, afwMath.SpatialCellCandidate.GOOD)

    def runTest(self, mode):
        logger.debug("Mode %s", mode)
        if mode == "DF":
            self.config = self.configDF
        elif mode == "DFr":
            self.config = self.configDFr
        elif mode == "AL":
            self.config = self.configAL
        else:
            raise
        subconfig = self.config.kernel.active
        psfmatch = ipDiffim.ImagePsfMatchTask(self.config)
        candidateList = []
        if self.scienceExposure.getPsf():
            candidateList = psfmatch.makeCandidateList(self.templateExposure, self.scienceExposure,
                                                       subconfig.kernelSize)
            results = psfmatch.subtractMaskedImages(self.templateMaskedImage, self.scienceMaskedImage,
                                                    candidateList)
            self.jackknifeResample(psfmatch, results)

    def test(self):
        if not defDataDir:
            print("Warning: afwdata not set up; not running JackknifeResampleSpatialKernel.py",
                  file=sys.stderr)
            return

        self.runTest(mode="AL")

    def tearDown(self):
        del self.subconfigAL
        del self.subconfigDF
        del self.subconfigDFr
        del self.scienceExposure
        del self.templateExposure
        del self.scienceMaskedImage
        del self.templateMaskedImage
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
    if len(sys.argv) > 3:
        defTemplatePath = sys.argv[1]
        defSciencePath = sys.argv[2]

    if '-d' in sys.argv:
        display = True

    run(True)

# python tests/JackknifeResampleSpatialKernel.py -d
# python tests/JackknifeResampleSpatialKernel.py $AFWDATA_DIR/CFHT/D4/cal-53535-i-797722_1_tmpl
# ... $AFWDATA_DIR/CFHT/D4/cal-53535-i-797722_1 -d
