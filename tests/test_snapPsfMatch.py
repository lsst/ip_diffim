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
import lsst.afw.image as afwImage
from lsst.afw.geom import makeSkyWcs
import lsst.meas.algorithms as measAlg
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.daf.base as dafBase
import lsst.utils.logging as logUtils

logUtils.trace_set_at("lsst.ip.diffim", 4)


class PsfMatchTestCases(lsst.utils.tests.TestCase):

    def setUp(self):
        self.configAL = ipDiffim.SnapPsfMatchTask.ConfigClass()
        self.configAL.kernel.name = "AL"
        self.subconfigAL = self.configAL.kernel.active

        self.configDF = ipDiffim.SnapPsfMatchTask.ConfigClass()
        self.configDF.kernel.name = "DF"
        self.subconfigDF = self.configDF.kernel.active

        self.configDFr = ipDiffim.SnapPsfMatchTask.ConfigClass()
        self.configDFr.kernel.name = "DF"
        self.subconfigDFr = self.configDFr.kernel.active

        self.subconfigDF.useRegularization = False
        self.subconfigDFr.useRegularization = True

        self.subconfigAL.afwBackgroundConfig.useApprox = False
        self.subconfigDF.afwBackgroundConfig.useApprox = False
        self.subconfigDFr.afwBackgroundConfig.useApprox = False

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

    def testSnap(self):
        tMi, sMi, sK, kcs, confake = diffimTools.makeFakeKernelSet(bgValue=self.bgValue)

        tWcs = self.makeWcs(offset=0)
        sWcs = self.makeWcs(offset=0)
        tExp = afwImage.ExposureF(tMi, tWcs)
        sExp = afwImage.ExposureF(sMi, sWcs)
        sExp.setPsf(self.psf)
        psfMatchAL = ipDiffim.SnapPsfMatchTask(config=self.configAL)
        psfMatchDF = ipDiffim.SnapPsfMatchTask(config=self.configDF)
        psfMatchDFr = ipDiffim.SnapPsfMatchTask(config=self.configDFr)
        psfMatchAL.subtractMaskedImages(tMi, sMi, psfMatchAL.makeCandidateList(tExp, sExp, self.ksize))
        psfMatchDF.subtractMaskedImages(tMi, sMi, psfMatchDF.makeCandidateList(tExp, sExp, self.ksize))
        psfMatchDFr.subtractMaskedImages(tMi, sMi, psfMatchDFr.makeCandidateList(tExp, sExp, self.ksize))

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
