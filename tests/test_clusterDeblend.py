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

import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.geom as geom
import lsst.utils.tests
from lsst.ip.diffim.clusterDeblend import deblend_clustered_peaks
from lsst.ip.diffim.peakClustering import find_peak_clusters

SIZE = 40


def makeSchemas():
    schema = afwTable.SourceTable.makeMinimalSchema()
    mergeList = afwDetection.FootprintMergeList(schema, ["positive", "negative"])
    return schema, mergeList.getPeakSchema()


def makeFootprintAndImage(peakSchema, peaks, ampValue=100.0, badPeaks=()):
    """Build a full-frame footprint with the given peaks, and a matching image.

    Parameters
    ----------
    peaks : `list` of `tuple`
        Each ``(ix, iy, isNegative)``.
    ampValue : `float`
        Magnitude placed in the 3x3 block around each peak (negated for
        negative peaks) in the difference image.
    badPeaks : `set` of `int`
        Indices of peaks whose 3x3 block is flagged bad in the mask.
    """
    box = geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(SIZE, SIZE))
    footprint = afwDetection.Footprint(afwGeom.SpanSet(box), peakSchema)
    negativeKey = peakSchema.find("merge_peak_negative").key
    positiveKey = peakSchema.find("merge_peak_positive").key

    image = np.zeros((SIZE, SIZE), dtype=np.float32)
    mask = np.zeros((SIZE, SIZE), dtype=np.int32)
    badBit = 1
    for i, (ix, iy, isNegative) in enumerate(peaks):
        peak = footprint.addPeak(ix, iy, -ampValue if isNegative else ampValue)
        peak.set(negativeKey, isNegative)
        peak.set(positiveKey, not isNegative)
        block = (slice(iy - 1, iy + 2), slice(ix - 1, ix + 2))
        image[block] = -ampValue if isNegative else ampValue
        if i in badPeaks:
            mask[block] = badBit
    return footprint, image, mask, badBit


def makeClusters(schema, footprint, linkingRadius):
    catalog = afwTable.SourceCatalog(schema)
    record = catalog.addNew()
    record.setId(1)
    record.setFootprint(footprint)
    results = find_peak_clusters(catalog, linkingRadius)
    return results[0] if results else None


class ClusterDeblendTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.schema, self.peakSchema = makeSchemas()

    def test_extract_clusters_from_footprint(self):
        """A footprint that splits into two clusters yields two sources."""
        # Lone peak at (8, 20); dipole at (30, 20)/(32, 20).
        footprint, image, mask, _ = makeFootprintAndImage(
            self.peakSchema, [(8, 20, False), (30, 20, False), (32, 20, True)])
        fc = makeClusters(self.schema, footprint, linkingRadius=5.0)
        self.assertTrue(fc.is_deblendable)

        result = deblend_clustered_peaks(fc, image, mask, (0, 0), sigma=4.0, bad_bitmask=0)
        self.assertIsNotNone(result)
        self.assertEqual(len(result.children), 1)

        # The first cluster (the lone peak at (8, 20)) stays on the parent.
        self.assertEqual(len(result.parent.getPeaks()), 1)
        self.assertEqual((result.parent.getPeaks()[0].getIx(),
                          result.parent.getPeaks()[0].getIy()), (8, 20))

        # The dipole is the other cluster, keeping both its peaks.
        child = result.children[0]
        self.assertEqual(len(child.getPeaks()), 2)
        self.assertTrue(child.getSpans().contains(geom.Point2I(30, 20)))
        self.assertTrue(child.getSpans().contains(geom.Point2I(32, 20)))

        # The extracted pixels exactly tile the original footprint.
        self.assertEqual(result.parent.getSpans().getArea() + child.getSpans().getArea(),
                         footprint.getSpans().getArea())
        # Peak schema (with merge flags) is preserved on the outputs.
        self.assertIn("merge_peak_negative", child.getPeaks().schema.getNames())

    def test_splits_beyond_radius_positives(self):
        """Two positive peaks beyond the linking radius split into two sources
        even with no remainder cluster.
        """
        footprint, image, mask, _ = makeFootprintAndImage(
            self.peakSchema, [(8, 20, False), (32, 20, False)])
        fc = makeClusters(self.schema, footprint, linkingRadius=5.0)
        self.assertTrue(fc.is_deblendable)

        result = deblend_clustered_peaks(fc, image, mask, (0, 0), sigma=4.0, bad_bitmask=0)
        self.assertIsNotNone(result)
        # One peak stays on the parent, the other becomes a new source.
        self.assertEqual(len(result.children), 1)
        self.assertEqual(len(result.parent.getPeaks()), 1)
        self.assertEqual(len(result.children[0].getPeaks()), 1)
        self.assertEqual(result.parent.getSpans().getArea()
                         + result.children[0].getSpans().getArea(),
                         footprint.getSpans().getArea())

    def test_within_radius_positive_blend_stays(self):
        """Two positive peaks within the linking radius are one cluster and
        are left as a single source.
        """
        footprint, image, mask, _ = makeFootprintAndImage(
            self.peakSchema, [(20, 20, False), (22, 20, False)])
        fc = makeClusters(self.schema, footprint, linkingRadius=5.0)
        self.assertFalse(fc.is_deblendable)
        result = deblend_clustered_peaks(fc, image, mask, (0, 0), sigma=4.0, bad_bitmask=0)
        self.assertIsNone(result)

    def test_single_cluster_dipole_returns_none(self):
        """A footprint that is a single cluster (here a dipole) is untouched."""
        footprint, image, mask, _ = makeFootprintAndImage(
            self.peakSchema, [(20, 20, False), (22, 20, True)])
        fc = makeClusters(self.schema, footprint, linkingRadius=5.0)
        self.assertFalse(fc.is_deblendable)
        result = deblend_clustered_peaks(fc, image, mask, (0, 0), sigma=4.0, bad_bitmask=0)
        self.assertIsNone(result)

    def test_extract_negative_cluster(self):
        """Negative and positive clusters are both extracted (polarity is not
        a constraint in this version).
        """
        # Lone negative at (8, 20); lone positive at (32, 20); beyond radius.
        footprint, image, mask, _ = makeFootprintAndImage(
            self.peakSchema, [(8, 20, True), (32, 20, False)])
        fc = makeClusters(self.schema, footprint, linkingRadius=5.0)
        self.assertTrue(fc.is_deblendable)

        result = deblend_clustered_peaks(fc, image, mask, (0, 0), sigma=4.0, bad_bitmask=0)
        self.assertIsNotNone(result)
        self.assertEqual(len(result.children), 1)
        # Both peaks are extracted, one per output footprint.
        outputs = [result.parent] + result.children
        extracted = {(fp.getPeaks()[0].getIx(), fp.getPeaks()[0].getIy()) for fp in outputs}
        self.assertEqual(extracted, {(8, 20), (32, 20)})

    def test_dominated_isolated_peak_is_merged(self):
        """A faint 'isolated' peak dominated by a bright neighbour is merged
        back in, leaving nothing to extract.
        """
        # Faint positive at (12, 20); bright dipole at (16, 20)/(18, 20).
        peaks = [(12, 20, False), (16, 20, False), (18, 20, True)]
        footprint, image, mask, _ = makeFootprintAndImage(self.peakSchema, peaks, ampValue=100.0)
        # Make the first peak faint: overwrite its 3x3 block with a small value.
        image[19:22, 11:14] = 5.0
        fc = makeClusters(self.schema, footprint, linkingRadius=3.0)
        self.assertTrue(fc.is_deblendable)
        result = deblend_clustered_peaks(fc, image, mask, (0, 0), sigma=2.0, bad_bitmask=0)
        self.assertIsNone(result)

    def test_bad_masked_amplitude_excluded(self):
        """A peak whose 3x3 block is entirely bad-masked gets zero amplitude
        and is merged away (proving bad pixels are excluded from the sum).
        """
        footprint, image, mask, badBit = makeFootprintAndImage(
            self.peakSchema, [(8, 20, False), (30, 20, False), (32, 20, True)],
            ampValue=100.0, badPeaks={0})
        fc = makeClusters(self.schema, footprint, linkingRadius=5.0)
        self.assertTrue(fc.is_deblendable)
        # With the bright isolated peak's flux excluded, it has no influence
        # and cannot be extracted.
        result = deblend_clustered_peaks(fc, image, mask, (0, 0), sigma=4.0, bad_bitmask=badBit)
        self.assertIsNone(result)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
