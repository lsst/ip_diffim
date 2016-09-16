from builtins import range
#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008-2016 LSST Corporation.
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
import unittest
import numpy as np
import lsst.utils.tests
import lsst.afw.table as afwTable
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
from lsst.meas.algorithms import LoadReferenceObjectsTask, getRefFluxField
import lsst.ip.diffim as ipDiffim


class DiaCatalogSourceSelectorTest(lsst.utils.tests.TestCase):

    def setUp(self):
        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField("test_flux", type=float)
        schema.addField("test_fluxSigma", type=float)
        self.sourceSelector = ipDiffim.DiaCatalogSourceSelectorTask(schema=schema)
        for flag in self.sourceSelector.config.badFlags:
            schema.addField(flag, type="Flag")
        table = afwTable.SourceTable.make(schema)
        table.definePsfFlux("test")
        self.srcCat = afwTable.SourceCatalog(table)
        self.exposure = afwImage.ExposureF()

    def tearDown(self):
        del self.sourceSelector
        del self.exposure
        del self.srcCat

    def makeRefCatalog(self):
        schema = LoadReferenceObjectsTask.makeMinimalSchema(filterNameList=["g", "r"],
                                                            addFluxSigma=False, addIsPhotometric=True, addIsResolved=True)
        catalog = afwTable.SimpleCatalog(schema)
        return catalog

    def makeMatches(self, refCat, srcCat, nSrc):
        for i in range(nSrc):

            refSrc = refCat.addNew()
            srcSrc = srcCat.addNew()

            coord = afwCoord.Coord(afwGeom.Point2D(*np.random.randn(2)), afwGeom.degrees)

            refSrc.set("g_flux", 10**(-0.4*18))
            refSrc.set("r_flux", 10**(-0.4*18))
            refSrc.set("resolved", False)
            refSrc.set("photometric", True)
            refSrc.setCoord(coord)

            srcSrc.setCoord(coord)
            srcSrc.set(srcSrc.getTable().getPsfFluxKey(), 10.)
            srcSrc.set(srcSrc.getTable().getPsfFluxErrKey(), 1.)
            for flag in self.sourceSelector.config.badFlags:
                srcSrc.set(flag, False)

        mat = afwTable.matchRaDec(refCat, srcCat, 1.0 * afwGeom.arcseconds, False)
        self.assertEqual(len(mat), nSrc)
        return mat

    def testCuts(self):
        nSrc = 5

        refCat = self.makeRefCatalog()

        matches = self.makeMatches(refCat, self.srcCat, nSrc)
        sources = self.sourceSelector.selectStars(self.exposure, self.srcCat, matches).starCat
        self.assertEqual(len(sources), nSrc)

        # Set one of the source flags to be bad
        matches[0].second.set(self.sourceSelector.config.badFlags[0], True)
        sources = self.sourceSelector.selectStars(self.exposure, self.srcCat, matches).starCat
        self.assertEqual(len(sources), nSrc-1)

        # Set one of the ref flags to be bad
        matches[1].first.set("photometric", False)
        sources = self.sourceSelector.selectStars(self.exposure, self.srcCat, matches).starCat
        self.assertEqual(len(sources), nSrc-2)

        # Set one of the colors to be bad
        grMin = self.sourceSelector.config.grMin
        rFluxField = getRefFluxField(refCat.schema, "r")
        gFluxField = getRefFluxField(refCat.schema, "g")
        gFlux = 10**(-0.4 * (grMin - 0.1)) * matches[2].first.get(rFluxField)
        matches[2].first.set(gFluxField, gFlux)
        sources = self.sourceSelector.selectStars(self.exposure, self.srcCat, matches).starCat
        self.assertEqual(len(sources), nSrc-3)

        # Set one of the types to be bad
        if self.sourceSelector.config.selectStar and not self.sourceSelector.config.selectGalaxy:
            matches[3].first.set("resolved", True)
            sources = self.sourceSelector.selectStars(self.exposure, self.srcCat, matches).starCat
            self.assertEqual(len(sources), nSrc-4)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
