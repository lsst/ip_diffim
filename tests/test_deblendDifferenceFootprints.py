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
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.geom as geom
import lsst.utils.tests
from lsst.ip.diffim.deblendDifferenceFootprints import (
    DeblendDifferenceFootprintsTask,
    _cluster_peak_indices,
    _deblend_footprint,
)


def makeSchemas():
    """Build a source schema with the deblending flags and the matching merged
    peak schema.
    """
    schema = afwTable.SourceTable.makeMinimalSchema()
    schema.addField("is_negative", type="Flag", doc="Source is negative.")
    schema.addField("is_deblended", type="Flag", doc="Source was deblended.")
    mergeList = afwDetection.FootprintMergeList(schema, ["positive", "negative"])
    peakSchema = mergeList.getPeakSchema()
    peakSchema.addField("significance", type=np.float64,
                        doc="Ratio of peak value to the standard deviation.")
    return schema, peakSchema


def makeFootprint(peakSchema, peaks):
    """Build a full-frame footprint carrying the given peaks.

    Parameters
    ----------
    peaks : `list` of `tuple`
        Each ``(ix, iy, isNegative, significance)``.
    """
    box = geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(64, 64))
    footprint = afwDetection.Footprint(afwGeom.SpanSet(box), peakSchema)
    negativeKey = peakSchema.find("merge_peak_negative").key
    positiveKey = peakSchema.find("merge_peak_positive").key
    significanceKey = peakSchema.find("significance").key
    for ix, iy, isNegative, significance in peaks:
        peak = footprint.addPeak(ix, iy, -significance if isNegative else significance)
        peak.set(negativeKey, isNegative)
        peak.set(positiveKey, not isNegative)
        peak.set(significanceKey, significance)
    return footprint


class ClusterPeaksTestCase(lsst.utils.tests.TestCase):
    """Tests of the initial single-linkage peak clustering."""

    def setUp(self):
        self.schema, self.peakSchema = makeSchemas()
        self.radius = 5.0

    def _cluster(self, peaks):
        # peaks: (ix, iy, isNegative); significance is irrelevant to clustering.
        footprint = makeFootprint(self.peakSchema, [(ix, iy, neg, 100.0) for ix, iy, neg in peaks])
        return _cluster_peak_indices(footprint.getPeaks(), self.radius)

    def test_two_peaks_far_apart(self):
        """Peaks beyond the linking radius form separate singleton clusters."""
        clusters = self._cluster([(10, 10, False), (30, 30, False)])
        self.assertEqual(sorted(clusters), [[0], [1]])

    def test_two_peaks_close(self):
        """Peaks within the linking radius merge into one cluster."""
        clusters = self._cluster([(10, 10, False), (12, 10, False)])
        self.assertEqual(len(clusters), 1)
        self.assertEqual(sorted(clusters[0]), [0, 1])

    def test_transitive_chain(self):
        """Single-linkage chains A-B-C even when A and C exceed the radius."""
        # 10->14->18: adjacent gaps are 4 (< 5), but 10->18 is 8 (> 5).
        A = (10, 10, False)
        B = (14, 10, False)
        C = (18, 10, False)

        # A and C alone should not merge
        clustersAC = self._cluster([A, C])
        self.assertEqual(len(clustersAC), 2)
        clustersABC = self._cluster([A, B, C])
        self.assertEqual(len(clustersABC), 1)
        self.assertEqual(sorted(clustersABC[0]), [0, 1, 2])

    def test_dipole_single_cluster(self):
        """A close positive/negative pair clusters together."""
        clusters = self._cluster([(10, 10, False), (12, 10, True)])
        self.assertEqual(len(clusters), 1)
        self.assertEqual(sorted(clusters[0]), [0, 1])

    def test_lone_peak_beside_dipole(self):
        """A lone peak is separated from a nearby dipole into two clusters."""
        clusters = self._cluster([(50, 50, False),           # lone peak
                                  (10, 10, False), (15, 10, True)])  # dipole
        lone = [c for c in clusters if 0 in c]
        dipole = [c for c in clusters if 1 in c]
        self.assertEqual(len(clusters), 2)
        self.assertEqual(lone[0], [0])
        self.assertEqual(sorted(dipole[0]), [1, 2])


class DeblendFootprintTestCase(lsst.utils.tests.TestCase):
    """Tests of the per-footprint deblending routine."""

    def setUp(self):
        self.schema, self.peakSchema = makeSchemas()

    def _deblend(self, peaks, linkingRadius=5.0, sigma=4.0, max_contamination=0.05):
        footprint = makeFootprint(self.peakSchema, peaks)

        clusters = _cluster_peak_indices(footprint.getPeaks(), linkingRadius)
        result = _deblend_footprint(footprint, clusters, sigma=sigma,
                                    max_contamination=max_contamination)
        return footprint, clusters, result

    def test_extract_clusters_from_footprint(self):
        """A footprint that splits into two clusters yields two footprints."""
        # Lone peak at (8, 20); dipole at (30, 20)/(35, 20).
        footprint, clusters, result = self._deblend(
            [(8, 20, False, 100.0), (30, 20, False, 100.0), (35, 20, True, 100.0)])
        self.assertEqual(len(clusters), 2)
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

    def test_splits_beyond_radius(self):
        """Two peaks beyond the linking radius split into two footprints."""
        _, _, result = self._deblend([(8, 20, False, 100.0), (32, 20, False, 100.0)])
        self.assertIsNotNone(result)
        self.assertEqual(len(result.children), 1)
        self.assertEqual(len(result.parent.getPeaks()), 1)
        self.assertEqual(len(result.children[0].getPeaks()), 1)

    def test_within_radius_blend_stays(self):
        """Two peaks within the linking radius are one cluster: no split."""
        footprint, clusters, result = self._deblend([(20, 20, False, 100.0), (22, 20, False, 100.0)])
        self.assertEqual(len(clusters), 1)
        self.assertIsNone(result)

    def test_single_cluster_dipole_returns_none(self):
        """A footprint that is a single cluster (a dipole) is untouched."""
        footprint, clusters, result = self._deblend([(20, 20, False, 100.0), (22, 20, True, 100.0)])
        self.assertEqual(len(clusters), 1)
        self.assertIsNone(result)

    def test_extract_negative_cluster(self):
        """Negative and positive clusters are both extracted."""
        _, _, result = self._deblend([(8, 20, True, 100.0), (32, 20, False, 100.0)])
        self.assertIsNotNone(result)
        self.assertEqual(len(result.children), 1)
        outputs = [result.parent] + result.children
        extracted = {(fp.getPeaks()[0].getIx(), fp.getPeaks()[0].getIy()) for fp in outputs}
        self.assertEqual(extracted, {(8, 20), (32, 20)})

    def test_low_significance_peak_is_merged(self):
        """A faint (low-significance) peak heavily contaminated by a bright
        neighbour is merged back in, leaving nothing to extract.
        """
        # Faint positive at (12, 20); bright dipole at (16, 20)/(18, 20).
        _, _, result = self._deblend(
            [(12, 20, False, 5.0), (16, 20, False, 100.0), (18, 20, True, 100.0)],
            linkingRadius=3.0, sigma=2.0)
        self.assertIsNone(result)

    def test_contamination_merges_beyond_cluster_radius(self):
        """Peaks the KD step leaves separate are still merged when their PSF
        flux contamination exceeds ``max_contamination``.
        """
        # d=5 px, sigma=2 -> overlap exp(-25/16)=0.21 > 0.05, but the peaks are
        # beyond linkingRadius=3 px, so the contamination test drives the merge.
        _, clusters, result = self._deblend(
            [(20, 20, False, 100.0), (25, 20, False, 100.0)],
            linkingRadius=3.0, sigma=2.0, max_contamination=0.05)
        self.assertEqual(len(clusters), 2)   # KD kept them separate
        self.assertIsNone(result)            # contamination merged them

    def test_no_merge_when_contamination_below_limit(self):
        """Peaks whose mutual contamination is below the limit stay separate."""
        # d=10 px, sigma=2 -> overlap exp(-100/16)=0.0019 < 0.05.
        _, clusters, result = self._deblend(
            [(20, 20, False, 100.0), (30, 20, False, 100.0)],
            linkingRadius=3.0, sigma=2.0, max_contamination=0.05)
        self.assertEqual(len(clusters), 2)
        self.assertIsNotNone(result)
        self.assertEqual(len(result.children), 1)


class RunTestCase(lsst.utils.tests.TestCase):
    """Tests of the task-level ``run`` over a diaSource catalog."""

    def setUp(self):
        self.schema, self.peakSchema = makeSchemas()

    def makeExposure(self, psfSigma=2.0):
        # The PSF sets the clustering scale; the mask is updated in place with
        # the NOT_DEBLENDED plane for any guardrail-skipped footprints.
        bbox = geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(64, 64))
        exposure = afwImage.ExposureF(bbox)
        exposure.setPsf(afwDetection.GaussianPsf(11, 11, psfSigma))
        return exposure

    def makeCatalog(self, footprints):
        catalog = afwTable.SourceCatalog(self.schema)
        for footprint in footprints:
            catalog.addNew().setFootprint(footprint)
        return catalog

    def makeTask(self, **configKwargs):
        config = DeblendDifferenceFootprintsTask.ConfigClass()
        # sigma ~= 2 px, FWHM ~= 4.7 px, so clusterRadius=2 links within ~9 px.
        config.clusterRadius = 2.0
        for key, value in configKwargs.items():
            setattr(config, key, value)
        return DeblendDifferenceFootprintsTask(config=config)

    def test_run_splits_multicluster_footprint(self):
        """A footprint with two well-separated peaks becomes two diaSources."""
        footprint = makeFootprint(self.peakSchema, [(10, 20, False, 100.0), (50, 20, False, 100.0)])
        catalog = self.makeCatalog([footprint])
        task = self.makeTask()
        difference = self.makeExposure()

        result = task.run(catalog, difference)

        self.assertEqual(len(result.diaSources), 2)
        self.assertTrue(all(record["is_deblended"] for record in result.diaSources))
        self.assertEqual(task.metadata["nDeblendedFootprints"], 1)
        self.assertEqual(task.metadata["nDeblendedDiaSources"], 1)
        self.assertEqual(task.metadata["nSkipped"], 0)
        # Nothing was skipped, so no pixels are flagged NOT_DEBLENDED.
        notDeblended = difference.mask.getPlaneBitMask("NOT_DEBLENDED")
        self.assertTrue(np.all((difference.mask.array & notDeblended) == 0))

    def test_run_leaves_single_cluster_untouched(self):
        """A footprint that is one cluster is not modified or flagged."""
        footprint = makeFootprint(self.peakSchema, [(20, 20, False, 100.0), (22, 20, False, 100.0)])
        catalog = self.makeCatalog([footprint])
        task = self.makeTask()

        result = task.run(catalog, self.makeExposure())

        self.assertEqual(len(result.diaSources), 1)
        self.assertFalse(result.diaSources[0]["is_deblended"])
        self.assertEqual(task.metadata["nDeblendedFootprints"], 0)

    def test_run_sets_is_negative(self):
        """An extracted purely-negative cluster is flagged is_negative."""
        footprint = makeFootprint(self.peakSchema, [(10, 20, False, 100.0), (50, 20, True, 100.0)])
        catalog = self.makeCatalog([footprint])
        task = self.makeTask()

        result = task.run(catalog, self.makeExposure())

        self.assertEqual(len(result.diaSources), 2)
        byNegative = sorted(record["is_negative"] for record in result.diaSources)
        self.assertEqual(byNegative, [False, True])

    def test_run_skips_too_many_peaks(self):
        """A deblendable footprint with more peaks than maxPeaks is skipped and
        its pixels are flagged NOT_DEBLENDED.
        """
        footprint = makeFootprint(
            self.peakSchema, [(10, 20, False, 100.0), (30, 20, False, 100.0), (50, 20, False, 100.0)])
        catalog = self.makeCatalog([footprint])
        task = self.makeTask(maxPeaks=2)
        difference = self.makeExposure()

        result = task.run(catalog, difference)

        self.assertEqual(len(result.diaSources), 1)
        self.assertFalse(result.diaSources[0]["is_deblended"])
        self.assertEqual(task.metadata["nSkipped"], 1)
        # The whole (full-frame) footprint is flagged NOT_DEBLENDED.
        notDeblended = difference.mask.getPlaneBitMask("NOT_DEBLENDED")
        self.assertTrue(np.all((difference.mask.array & notDeblended) != 0))

    def test_run_skips_too_large_area(self):
        """A deblendable footprint bigger than maxFootprintArea is skipped and
        its pixels are flagged NOT_DEBLENDED.
        """
        footprint = makeFootprint(self.peakSchema, [(10, 20, False, 100.0), (50, 20, False, 100.0)])
        catalog = self.makeCatalog([footprint])
        task = self.makeTask(maxFootprintArea=100)  # footprint area is 64**2
        difference = self.makeExposure()

        result = task.run(catalog, difference)

        self.assertEqual(len(result.diaSources), 1)
        self.assertFalse(result.diaSources[0]["is_deblended"])
        self.assertEqual(task.metadata["nSkipped"], 1)
        notDeblended = difference.mask.getPlaneBitMask("NOT_DEBLENDED")
        self.assertTrue(np.all((difference.mask.array & notDeblended) != 0))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
