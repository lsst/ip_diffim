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

import unittest

import numpy as np

import lsst.afw.geom
import lsst.afw.image
import lsst.afw.table
import lsst.geom
import lsst.ip.diffim
import lsst.meas.algorithms
import lsst.utils.tests


def _makeWcs(center=(50, 50)):
    return lsst.afw.geom.makeSkyWcs(
        lsst.geom.Point2D(*center),
        lsst.geom.SpherePoint(0, 0, lsst.geom.degrees),
        lsst.afw.geom.makeCdMatrix(0.2*lsst.geom.arcseconds, 0*lsst.geom.degrees))


def _makeCcds(wcs, specs):
    """Build a per-CCD coadd input catalog from (visit, ccd, Box2D) specs."""
    schema = lsst.afw.table.ExposureTable.makeMinimalSchema()
    schema.addField("visit", type=np.int64, doc="visit id")
    schema.addField("ccd", type=np.int32, doc="ccd id")
    ccds = lsst.afw.table.ExposureCatalog(schema)
    for visit, ccd, polyBox in specs:
        record = ccds.addNew()
        record.setWcs(wcs)
        record.setBBox(lsst.geom.Box2I(polyBox))
        record.setValidPolygon(lsst.afw.geom.Polygon(polyBox.getCorners()))
        record["visit"] = visit
        record["ccd"] = ccd
    return ccds


class CheckTemplateQualityTestCase(lsst.utils.tests.TestCase):
    """Test the template-quality mask-setting checks."""

    def testComputeDepthMap(self):
        """_computeDepthMap counts the distinct input images covering each pixel
        and deduplicates inputs on (visit, ccd).
        """
        wcs = _makeWcs()
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(100, 100))
        fullBox = lsst.geom.Box2D(box)
        leftBox = lsst.geom.Box2D(
            lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(50, 100)))

        # The left half is covered by visits {1, 2, 3} and the right half by {1};
        # the duplicate (visit, ccd) must not inflate the depth.
        ccds = _makeCcds(wcs, [(1, 1, fullBox), (2, 1, leftBox), (3, 1, leftBox),
                               (1, 1, fullBox)])
        task = lsst.ip.diffim.CheckTemplateQualityTask()
        depth = task._computeDepthMap([ccds], wcs, box)
        self.assertTrue((depth[:, :50] == 3).all())
        self.assertTrue((depth[:, 50:] == 1).all())

    def testMaskShallowCoverage(self):
        """A pixel is flagged HIGH_VARIANCE only where every contributing tract
        is shallow; a deep tract filling another's gap keeps it unflagged.
        """
        wcs = _makeWcs()
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(100, 100))
        fullBox = lsst.geom.Box2D(box)
        rightBox = lsst.geom.Box2D(
            lsst.geom.Box2I(lsst.geom.Point2I(50, 0), lsst.geom.Extent2I(50, 100)))

        # Tract 0 covers everything but is shallow (depth 1); tract 1 covers the
        # right half deeply (depth 2). So the left half is shallow with no deep
        # tract -> flagged; the right half is filled by tract 1 -> not flagged.
        rightHasData = np.zeros((100, 100), dtype=bool)
        rightHasData[:, 50:] = True
        hasDataByTract = {0: np.ones((100, 100), dtype=bool), 1: rightHasData}
        ccdInputCatalogsByTract = {
            0: [_makeCcds(wcs, [(1, 1, fullBox)])],
            1: [_makeCcds(wcs, [(10, 1, rightBox), (11, 1, rightBox)])],
        }

        def flaggedMask(growRadius):
            template = lsst.afw.image.ExposureF(box)
            template.setWcs(wcs)
            config = lsst.ip.diffim.CheckTemplateQualityConfig()
            config.minNumberOfInputImages = 2
            config.coverageGrowRadius = growRadius
            task = lsst.ip.diffim.CheckTemplateQualityTask(config=config)
            task.maskShallowCoverage(template, hasDataByTract, ccdInputCatalogsByTract)
            return (template.mask.array
                    & template.mask.getPlaneBitMask("HIGH_VARIANCE")) != 0

        flagged = flaggedMask(0)
        self.assertTrue(flagged[:, :50].all())
        self.assertEqual(flagged[:, 50:].sum(), 0)

        # Growing by 5 expands the masked region by 5 columns into the deep half.
        self.assertEqual(flaggedMask(5).any(axis=0).sum(), 55)

    def testLabelConstantDepthRegions(self):
        """Connected regions of constant depth get distinct labels, including
        two disconnected regions that share the same depth value.
        """
        depth = np.zeros((10, 30), dtype=np.int32)
        depth[:, 0:10] = 3
        depth[:, 10:20] = 5
        depth[:, 20:30] = 3
        hasData = np.ones(depth.shape, dtype=bool)
        labels = lsst.ip.diffim.CheckTemplateQualityTask._labelConstantDepthRegions(depth, hasData)
        # Two depth-3 blocks (disconnected) and one depth-5 block -> 3 regions.
        self.assertEqual(len(np.unique(labels)), 3)
        self.assertNotEqual(labels[5, 5], labels[5, 25])   # same depth, disconnected
        self.assertNotEqual(labels[5, 5], labels[5, 15])   # different depth

    def testResolveNarrowRegions(self):
        """Narrow shallow regions are masked; narrow deep regions are merged
        into their largest neighbor.
        """
        task = lsst.ip.diffim.CheckTemplateQualityTask()  # minRegionWidth=10, minNumberOfInputImages=2
        depth = np.full((100, 100), 5, dtype=np.int32)
        depth[:, 0:4] = 8       # narrow deep strip at the left edge (one neighbor)
        depth[:, 50:54] = 1     # narrow shallow strip, splits the depth-5 background
        hasData = np.ones(depth.shape, dtype=bool)
        labels = task._labelConstantDepthRegions(depth, hasData)
        maskNarrow = task._resolveNarrowRegions(labels, depth)

        # The narrow shallow strip is masked and removed from the labels.
        self.assertEqual(maskNarrow.sum(), 4*100)
        self.assertTrue(maskNarrow[:, 50:54].all())
        self.assertTrue((labels[:, 50:54] == 0).all())
        # The narrow deep strip is merged into its neighbor and not masked.
        self.assertFalse(maskNarrow[:, 0:4].any())
        self.assertTrue((labels[:, 0:4] == labels[5, 10]).all())

    def testFlagDiscontinuousRegions(self):
        """A region whose PSF size deviates from the dominant model is flagged,
        while a smooth gradient consistent with the polynomial order is not.
        """
        task = lsst.ip.diffim.CheckTemplateQualityTask()  # order 1, threshold 0.05

        # Two large regions at radius 5 and a small deviant region at radius 8.
        x1 = np.array([-0.8, -0.6, -0.4, -0.8, -0.6, -0.4])
        y1 = np.array([-0.5, -0.5, -0.5, 0.5, 0.5, 0.5])
        x2 = np.array([0.4, 0.6, 0.8, 0.4, 0.6, 0.8])
        y2 = y1
        x3 = np.array([0.9, 0.85])
        y3 = np.array([0.9, 0.95])
        samples = {
            "x": np.concatenate([x1, x2, x3]),
            "y": np.concatenate([y1, y2, y3]),
            "radius": np.concatenate([np.full(6, 5.0), np.full(6, 5.0), np.full(2, 8.0)]),
            "weight": np.concatenate([np.full(6, 4e4), np.full(6, 4e4), np.full(2, 1e3)]),
            "label": np.concatenate([np.full(6, 1), np.full(6, 2), np.full(2, 3)]),
        }
        self.assertEqual(task._flagDiscontinuousRegions(samples), {3})

        # A linear gradient in PSF size is fully absorbed by the order-1 fit.
        xs = np.linspace(-1, 1, 6)
        gradient = {
            "x": xs,
            "y": np.zeros(6),
            "radius": 5 + 2*xs,
            "weight": np.full(6, 4e4),
            "label": np.arange(1, 7),
        }
        self.assertEqual(task._flagDiscontinuousRegions(gradient), set())

    def testMaskPsfDiscontinuity(self):
        """End-to-end: a spatially constant PSF over two depth regions is not
        flagged, while masking that leaves too few usable pixels raises.
        """
        wcs = _makeWcs(center=(100, 100))
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(200, 200))

        def makeTemplate():
            template = lsst.afw.image.ExposureF(bbox)
            template.setWcs(wcs)
            template.setPsf(lsst.meas.algorithms.SingleGaussianPsf(25, 25, 2.0))
            return template

        fullBox = lsst.geom.Box2D(bbox)
        leftHalf = lsst.geom.Box2D(
            lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(100, 200)))
        # Depth 2 on the left, 1 on the right, but a spatially constant PSF: no
        # discontinuity, so nothing should be flagged.
        template = makeTemplate()
        lsst.ip.diffim.CheckTemplateQualityTask().maskPsfDiscontinuity(
            template, [_makeCcds(wcs, [(1, 1, fullBox), (2, 1, leftHalf)])])
        if "PSF_DISCONTINUITY" in template.mask.getMaskPlaneDict():
            bit = template.mask.getPlaneBitMask("PSF_DISCONTINUITY")
            self.assertEqual(np.count_nonzero(template.mask.array & bit), 0)

        # A narrow shallow strip (cols 196-199, depth 1) is masked; requiring a
        # high usable fraction then rejects the template.
        almostFull = lsst.geom.Box2D(
            lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(196, 200)))
        config = lsst.ip.diffim.CheckTemplateQualityConfig()
        config.psfDiscontinuityMinUsableFraction = 0.99
        with self.assertRaises(lsst.ip.diffim.PsfDiscontinuityError):
            lsst.ip.diffim.CheckTemplateQualityTask(config=config).maskPsfDiscontinuity(
                makeTemplate(), [_makeCcds(wcs, [(1, 1, almostFull), (2, 1, fullBox)])])


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
