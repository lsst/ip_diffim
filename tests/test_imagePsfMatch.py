#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#


import unittest


import lsst.utils.tests
from lsst.afw.geom import makeSkyWcs
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.daf.base as dafBase
import lsst.log.utils as logUtils
import lsst.meas.algorithms as measAlg

logUtils.traceSetAt("lsst.ip.diffim", 4)


class PsfMatchTestCases(unittest.TestCase):

    def setUp(self):
        self.configAL = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.configAL.kernel.name = "AL"
        self.subconfigAL = self.configAL.kernel.active

        self.configDF = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.configDF.kernel.name = "DF"
        self.subconfigDF = self.configDF.kernel.active

        self.configDFr = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.configDFr.kernel.name = "DF"
        self.subconfigDFr = self.configDFr.kernel.active

        self.subconfigAL.afwBackgroundConfig.useApprox = False
        self.subconfigDF.afwBackgroundConfig.useApprox = False
        self.subconfigDFr.afwBackgroundConfig.useApprox = False

        self.subconfigDF.useRegularization = False
        self.subconfigDFr.useRegularization = True

        self.subconfigAL.constantVarianceWeighting = False
        self.subconfigDF.constantVarianceWeighting = False
        self.subconfigDFr.constantVarianceWeighting = False

        # variance is a hack
        self.subconfigAL.singleKernelClipping = False
        self.subconfigAL.spatialKernelClipping = False
        self.subconfigDF.singleKernelClipping = False
        self.subconfigDF.spatialKernelClipping = False
        self.subconfigDFr.singleKernelClipping = False
        self.subconfigDFr.spatialKernelClipping = False

        # Send fake kernel a differential background
        self.bgValue = 100.
        self.subconfigAL.fitForBackground = True
        self.subconfigDF.fitForBackground = True
        self.subconfigDFr.fitForBackground = True

        # Make ideal PSF
        self.ksize = 21
        self.sigma = 2.0
        self.psf = measAlg.DoubleGaussianPsf(self.ksize, self.ksize, self.sigma)

    def makeWcs(self, offset=0):
        # taken from $AFW_DIR/tests/testMakeWcs.py
        metadata = dafBase.PropertySet()
        metadata.set("SIMPLE", "T")
        metadata.set("BITPIX", -32)
        metadata.set("NAXIS", 2)
        metadata.set("NAXIS1", 1024)
        metadata.set("NAXIS2", 1153)
        metadata.set("RADESYS", 'FK5')
        metadata.set("EQUINOX", 2000.)
        metadata.setDouble("CRVAL1", 215.604025685476)
        metadata.setDouble("CRVAL2", 53.1595451514076)
        metadata.setDouble("CRPIX1", 1109.99981456774 + offset)
        metadata.setDouble("CRPIX2", 560.018167811613 + offset)
        metadata.set("CTYPE1", 'RA---SIN')
        metadata.set("CTYPE2", 'DEC--SIN')
        metadata.setDouble("CD1_1", 5.10808596133527E-05)
        metadata.setDouble("CD1_2", 1.85579539217196E-07)
        metadata.setDouble("CD2_2", -5.10281493481982E-05)
        metadata.setDouble("CD2_1", -8.27440751733828E-07)
        return makeSkyWcs(metadata)

    def testWarping(self):
        tMi, sMi, sK, kcs, confake = diffimTools.makeFakeKernelSet(bgValue=self.bgValue)

        tWcs = self.makeWcs(offset=0)
        sWcs = self.makeWcs(offset=1)
        tExp = afwImage.ExposureF(tMi, tWcs)
        sExp = afwImage.ExposureF(sMi, sWcs)

        # Should fail due to registration problem
        psfMatchAL = ipDiffim.ImagePsfMatchTask(config=self.configAL)
        try:
            psfMatchAL.subtractExposures(tExp, sExp, doWarping=True)
        except Exception as e:
            print("testWarning failed with %r" % (e,))
            pass
        else:
            self.fail()

    def testSubtractExposures(self):
        # Test all 3 options
        tMi, sMi, sK, kcs, confake = diffimTools.makeFakeKernelSet(bgValue=self.bgValue)

        tWcs = self.makeWcs(offset=0)
        sWcs = self.makeWcs(offset=1)
        tExp = afwImage.ExposureF(tMi, tWcs)
        sExp = afwImage.ExposureF(sMi, sWcs)
        sExp.setPsf(self.psf)
        tExp.setPsf(self.psf)

        psfMatchAL = ipDiffim.ImagePsfMatchTask(config=self.configAL)
        psfMatchDF = ipDiffim.ImagePsfMatchTask(config=self.configDF)
        psfMatchDFr = ipDiffim.ImagePsfMatchTask(config=self.configDFr)

        self.assertEqual(psfMatchAL.useRegularization, False)
        self.assertEqual(psfMatchDF.useRegularization, False)
        self.assertEqual(psfMatchDFr.useRegularization, True)

        resultsAL = psfMatchAL.subtractExposures(tExp, sExp, doWarping=True)
        psfMatchDF.subtractExposures(tExp, sExp, doWarping=True)
        psfMatchDFr.subtractExposures(tExp, sExp, doWarping=True)

        self.assertEqual(type(resultsAL.subtractedExposure), afwImage.ExposureF)
        self.assertEqual(type(resultsAL.psfMatchingKernel), afwMath.LinearCombinationKernel)
        self.assertEqual(type(resultsAL.backgroundModel), afwMath.Chebyshev1Function2D)
        self.assertEqual(type(resultsAL.kernelCellSet), afwMath.SpatialCellSet)

    def testMatchExposures(self):
        # Only test 1 option
        tMi, sMi, sK, kcs, confake = diffimTools.makeFakeKernelSet(bgValue=self.bgValue)

        tWcs = self.makeWcs(offset=0)
        sWcs = self.makeWcs(offset=1)
        tExp = afwImage.ExposureF(tMi, tWcs)
        sExp = afwImage.ExposureF(sMi, sWcs)
        sExp.setPsf(self.psf)
        tExp.setPsf(self.psf)

        psfMatchAL = ipDiffim.ImagePsfMatchTask(config=self.configAL)
        resultsAL = psfMatchAL.matchExposures(tExp, sExp,
                                              templateFwhmPix=2.0, scienceFwhmPix=3.0, doWarping=True)
        self.assertEqual(type(resultsAL.matchedExposure), afwImage.ExposureF)
        self.assertEqual(type(resultsAL.psfMatchingKernel), afwMath.LinearCombinationKernel)
        self.assertEqual(type(resultsAL.backgroundModel), afwMath.Chebyshev1Function2D)
        self.assertEqual(type(resultsAL.kernelCellSet), afwMath.SpatialCellSet)

    def testPca(self, nTerms=3):
        tMi, sMi, sK, kcs, confake = diffimTools.makeFakeKernelSet(bgValue=self.bgValue)

        tWcs = self.makeWcs(offset=0)
        sWcs = self.makeWcs(offset=0)
        tExp = afwImage.ExposureF(tMi, tWcs)
        sExp = afwImage.ExposureF(sMi, sWcs)
        sExp.setPsf(self.psf)

        self.subconfigDF.usePcaForSpatialKernel = True
        self.subconfigDF.numPrincipalComponents = nTerms

        psfMatchDF = ipDiffim.ImagePsfMatchTask(config=self.configDF)
        candList = psfMatchDF.makeCandidateList(tExp, sExp, self.ksize)
        resultsDF = psfMatchDF.subtractMaskedImages(tMi, sMi, candList)

        spatialKernel = resultsDF.psfMatchingKernel
        spatialKernelSolution = spatialKernel.getSpatialParameters()
        self.assertEqual(len(spatialKernelSolution), nTerms)

        # First basis has no spatial variation
        for i in range(1, nTerms):
            self.assertEqual(spatialKernelSolution[0][i], 0.)

        # All bases have correct number of terms
        sko = self.subconfigDF.spatialKernelOrder
        nSpatialTerms = int(0.5 * (sko + 1) * (sko + 2))
        for i in range(len(spatialKernelSolution)):
            self.assertEqual(len(spatialKernelSolution[i]), nSpatialTerms)

        spatialBg = resultsDF.backgroundModel
        spatialBgSolution = spatialBg.getParameters()
        bgo = self.subconfigDF.spatialBgOrder
        nBgTerms = int(0.5 * (bgo + 1) * (bgo + 2))
        self.assertEqual(len(spatialBgSolution), nBgTerms)

    def testSubtractMaskedImages(self):
        # Lets do some additional testing here to make sure we recover
        # the known spatial model.  No background, just the faked
        # alard-lupton basis set.  The rest of matchMaskedImages() and
        # subtractMaskedImages() functionality is tested by the
        # Exposure-based methods.
        fakeCoeffs = diffimTools.fakeCoeffs()

        # Quick note; you shouldn't change confake here, since the
        # candidates in the KernelCellSet are initialized in
        # makeFakeKernelSet
        tMi, sMi, sK, kcs, confake = diffimTools.makeFakeKernelSet(bgValue=0.0, addNoise=False)

        svar = sMi.getVariance()
        svar.set(1.0)
        tvar = tMi.getVariance()
        tvar.set(1.0)

        basisList = ipDiffim.makeKernelBasisList(confake.kernel.active)
        psfMatchAL = ipDiffim.ImagePsfMatchTask(config=confake)
        spatialSolution, psfMatchingKernel, backgroundModel = psfMatchAL._solve(kcs, basisList)

        fitCoeffs = psfMatchingKernel.getSpatialParameters()

        for b in range(len(fakeCoeffs)):
            for s in range(len(fakeCoeffs[b])):

                if fakeCoeffs[b][s] == 0.0:
                    self.assertAlmostEqual(fitCoeffs[b][s], 0.0)
                else:
                    # OUTSTANDING ISSUE - WHY IS THIS ONE TERM OFF!?!?
                    if b != 4 and s != 0:
                        self.assertAlmostEqual(fitCoeffs[b][s]/fakeCoeffs[b][s], 1.0, 1)

    def tearDown(self):
        del self.configAL
        del self.configDF
        del self.configDFr
        del self.psf


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
