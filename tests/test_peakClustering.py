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

import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.geom as geom
import lsst.utils.tests
from lsst.ip.diffim.peakClustering import cluster_peaks, find_peak_clusters


def makePeakSchema():
    """Build a source schema and the matching merged-peak schema.

    The peak schema carries the ``merge_peak_positive``/``merge_peak_negative``
    flags, mimicking the footprints produced by ``FootprintMergeList`` in
    ``DetectAndMeasureTask``.
    """
    schema = afwTable.SourceTable.makeMinimalSchema()
    mergeList = afwDetection.FootprintMergeList(schema, ["positive", "negative"])
    return schema, mergeList.getPeakSchema()


def makeFootprint(peakSchema, peaks, size=100):
    """Make a footprint spanning a ``size`` x ``size`` box with the given peaks.

    Parameters
    ----------
    peakSchema : `lsst.afw.table.Schema`
        Schema carrying the merge-peak flags.
    peaks : `list` [`tuple`]
        Each entry is ``(fx, fy, value, isNegative)``.
    """
    spans = afwGeom.SpanSet(geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(size, size)))
    footprint = afwDetection.Footprint(spans, peakSchema)
    negativeKey = peakSchema.find("merge_peak_negative").key
    positiveKey = peakSchema.find("merge_peak_positive").key
    for fx, fy, value, isNegative in peaks:
        peak = footprint.addPeak(fx, fy, value)
        peak.set(negativeKey, isNegative)
        peak.set(positiveKey, not isNegative)
    return footprint


def clusterContaining(clusters, peakIndex):
    """Return the single cluster whose members include ``peakIndex``."""
    matches = [c for c in clusters if peakIndex in c.peak_indices]
    assert len(matches) == 1, f"peak {peakIndex} found in {len(matches)} clusters"
    return matches[0]


class PeakClusteringTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.schema, self.peakSchema = makePeakSchema()
        self.radius = 5.0

    def test_two_peaks_far_apart(self):
        """Peaks beyond the linking radius form separate singleton clusters."""
        footprint = makeFootprint(self.peakSchema,
                                  [(10, 10, 100, False), (30, 30, 120, False)])
        clusters = cluster_peaks(footprint.getPeaks(), self.radius)
        self.assertEqual(len(clusters), 2)
        for cluster in clusters:
            self.assertTrue(cluster.is_singleton)
            self.assertTrue(cluster.is_isolated_positive)

    def test_two_peaks_close(self):
        """Peaks within the linking radius merge into one cluster."""
        footprint = makeFootprint(self.peakSchema,
                                  [(10, 10, 100, False), (12, 10, 90, False)])
        clusters = cluster_peaks(footprint.getPeaks(), self.radius)
        self.assertEqual(len(clusters), 1)
        self.assertFalse(clusters[0].is_singleton)
        self.assertFalse(clusters[0].is_isolated_positive)
        self.assertEqual(sorted(clusters[0].peak_indices), [0, 1])

    def test_transitive_chain(self):
        """Single-linkage chains A-B-C even when A and C exceed the radius."""
        # 10->14->18: adjacent gaps are 4 (< 5), but 10->18 is 8 (> 5).
        footprint = makeFootprint(self.peakSchema,
                                  [(10, 10, 100, False), (14, 10, 100, False),
                                   (18, 10, 100, False)])
        clusters = cluster_peaks(footprint.getPeaks(), self.radius)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(sorted(clusters[0].peak_indices), [0, 1, 2])

    def test_dipole_single_cluster(self):
        """A close positive/negative pair clusters together, not isolated."""
        footprint = makeFootprint(self.peakSchema,
                                  [(10, 10, 100, False), (12, 10, -90, True)])
        clusters = cluster_peaks(footprint.getPeaks(), self.radius)
        self.assertEqual(len(clusters), 1)
        cluster = clusters[0]
        self.assertFalse(cluster.is_isolated_positive)
        self.assertEqual(cluster.is_negative.tolist(), [False, True])
        self.assertEqual(cluster.peak_values.tolist(), [100, -90])

    def test_isolated_positive_beside_dipole(self):
        """An isolated positive peak is separated from a nearby dipole."""
        footprint = makeFootprint(self.peakSchema,
                                  [(50, 50, 200, False),           # isolated positive
                                   (10, 10, 100, False), (12, 10, -90, True)])  # dipole
        clusters = cluster_peaks(footprint.getPeaks(), self.radius)
        self.assertEqual(len(clusters), 2)
        isolated = clusterContaining(clusters, 0)
        dipole = clusterContaining(clusters, 1)
        self.assertTrue(isolated.is_isolated_positive)
        self.assertFalse(dipole.is_isolated_positive)
        self.assertEqual(sorted(dipole.peak_indices), [1, 2])

    def test_isolated_negative_not_extracted(self):
        """A lone negative peak is a singleton but not an isolated positive."""
        footprint = makeFootprint(self.peakSchema,
                                  [(50, 50, 200, False), (10, 10, -80, True)])
        clusters = cluster_peaks(footprint.getPeaks(), self.radius)
        negative = clusterContaining(clusters, 1)
        self.assertTrue(negative.is_singleton)
        self.assertFalse(negative.is_isolated_positive)

    def test_positions_and_values(self):
        """Cluster arrays carry the peaks' positions and signed values."""
        footprint = makeFootprint(self.peakSchema, [(10, 20, 100, False)])
        (cluster,) = cluster_peaks(footprint.getPeaks(), self.radius)
        self.assertEqual(cluster.positions.tolist(), [[10, 20]])
        self.assertEqual(cluster.peak_values.tolist(), [100])
        self.assertEqual(cluster.is_negative.tolist(), [False])


class FindPeakClustersTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.schema, self.peakSchema = makePeakSchema()
        self.radius = 5.0

    def makeCatalog(self, footprints):
        catalog = afwTable.SourceCatalog(self.schema)
        for i, footprint in enumerate(footprints):
            record = catalog.addNew()
            record.setId(i + 1)
            record.setFootprint(footprint)
        return catalog

    def test_skips_single_peak_footprints(self):
        """Footprints with fewer than min_peaks peaks are not returned."""
        single = makeFootprint(self.peakSchema, [(10, 10, 100, False)])
        multi = makeFootprint(self.peakSchema,
                              [(10, 10, 100, False), (30, 30, 120, False)])
        catalog = self.makeCatalog([single, multi])
        results = find_peak_clusters(catalog, self.radius)
        self.assertEqual(len(results), 1)
        # Only the second (id=2) footprint qualifies.
        self.assertEqual(results[0].source_id, 2)
        self.assertEqual(results[0].n_peaks, 2)

    def test_handles_missing_footprint(self):
        """Sources with no footprint are skipped rather than raising."""
        multi = makeFootprint(self.peakSchema,
                              [(10, 10, 100, False), (30, 30, 120, False)])
        catalog = afwTable.SourceCatalog(self.schema)
        catalog.addNew().setId(1)  # no footprint attached
        record = catalog.addNew()
        record.setId(2)
        record.setFootprint(multi)
        results = find_peak_clusters(catalog, self.radius)
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].source_id, 2)

    def test_deblendable_classification(self):
        """is_deblendable requires both an extraction candidate and a
        remainder.
        """
        # Isolated positive + dipole: deblendable.
        mixed = makeFootprint(self.peakSchema,
                              [(50, 50, 200, False),
                               (10, 10, 100, False), (12, 10, -90, True)])
        # Two isolated positives beyond the radius: deblendable into two
        # sources even though there is no remainder cluster.
        allIsolated = makeFootprint(self.peakSchema,
                                    [(10, 10, 100, False), (50, 50, 120, False)])
        # One tight blend: a remainder but nothing isolated to extract.
        blend = makeFootprint(self.peakSchema,
                              [(10, 10, 100, False), (12, 10, 90, False)])
        catalog = self.makeCatalog([mixed, allIsolated, blend])
        results = find_peak_clusters(catalog, self.radius)
        byId = {r.source_id: r for r in results}

        self.assertTrue(byId[1].is_deblendable)
        self.assertEqual(len(byId[1].isolated_positive), 1)
        self.assertEqual(len(byId[1].retained), 1)

        self.assertTrue(byId[2].is_deblendable)
        self.assertEqual(len(byId[2].isolated_positive), 2)
        self.assertEqual(len(byId[2].retained), 0)

        self.assertFalse(byId[3].is_deblendable)
        self.assertEqual(len(byId[3].isolated_positive), 0)
        self.assertEqual(len(byId[3].retained), 1)

    def test_source_handle_is_live(self):
        """The stored source is the live catalog record, so deblender edits
        made through it land on the catalog.
        """
        multi = makeFootprint(self.peakSchema,
                              [(10, 10, 100, False), (30, 30, 120, False)])
        catalog = self.makeCatalog([multi])
        results = find_peak_clusters(catalog, self.radius)
        # Mutating through the stored handle must be visible in the catalog.
        results[0].source.setId(4242)
        self.assertEqual(catalog[0].getId(), 4242)
        self.assertEqual(results[0].footprint.getBBox(), catalog[0].getFootprint().getBBox())


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
