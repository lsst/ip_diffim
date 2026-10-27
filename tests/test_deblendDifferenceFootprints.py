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
    _merge_trail_clusters,
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


def makeDiagonalFootprint(peakSchema, peaks, width=3):
    """Build a thin diagonal footprint carrying the given peaks.
    """
    mask = afwImage.Mask(geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(64, 64)))
    bit = mask.getPlaneBitMask("DETECTED")
    rows, columns = np.mgrid[0:64, 0:64]
    mask.array[np.abs(rows - columns) <= width] = bit
    footprint = afwDetection.Footprint(afwGeom.SpanSet.fromMask(mask, bit), peakSchema)
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

    def test_trail_survives_contamination_merge(self):
        """A cluster that absorbs a trailed source through the contamination
        merge is still a trailed source, and is left as detected.
        """
        footprint = makeFootprint(self.peakSchema, [(10, 20, False, 100.0), (22, 20, False, 100.0),
                                                    (34, 20, False, 100.0), (38, 20, False, 100.0)])
        # 4 px from the trail's last peak, with sigma=2: overlap exp(-1) > 0.05.
        result = _deblend_footprint(footprint, [[3], [0, 1, 2]], sigma=2.0, trails=[False, True])
        self.assertIsNotNone(result)
        self.assertEqual(result.children, [])
        self.assertEqual(result.trails, [True])
        self.assertIs(result.parent, footprint)


class RunTestCase(lsst.utils.tests.TestCase):
    """Tests of the task-level ``run`` over a diaSource catalog."""

    def setUp(self):
        self.schema, self.peakSchema = makeSchemas()

    def makeExposure(self, psfSigma=2.0):
        # The PSF sets the clustering scale; the mask is updated in place with
        # the NOT_DEBLENDED plane for any guardrail-skipped footprints.
        bbox = geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(64, 64))
        exposure = afwImage.ExposureF(bbox)
        exposure.setPsf(afwDetection.GaussianPsf(21, 21, psfSigma))
        exposure.variance.array[:] = 1.0
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
        # The size limits apply before the collinear peaks could form a trail.
        self.assertEqual(task.metadata["nTrailedDiaSources"], 0)
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


class FootprintAreaTestCase(lsst.utils.tests.TestCase):
    """Tests for the limits on footprint size."""

    def setUp(self):
        self.schema, self.peakSchema = makeSchemas()
        # Well beyond the KD tree linking radius, and both on the diagonal.
        self.peaks = [(10, 10, False, 100.0), (50, 50, False, 100.0)]

    def makeExposure(self):
        bbox = geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(64, 64))
        exposure = afwImage.ExposureF(bbox)
        exposure.setPsf(afwDetection.GaussianPsf(21, 21, 2.0))
        exposure.variance.array[:] = 1.0
        return exposure

    def runTask(self, footprint, **configKwargs):
        catalog = afwTable.SourceCatalog(self.schema)
        catalog.addNew().setFootprint(footprint)
        config = DeblendDifferenceFootprintsTask.ConfigClass()
        config.clusterRadius = 2.0
        for key, value in configKwargs.items():
            setattr(config, key, value)
        task = DeblendDifferenceFootprintsTask(config=config)
        return task, task.run(catalog, self.makeExposure()).diaSources

    def test_diagonal_footprint_is_deblended(self):
        """A narrow but long diagonal footprint should be deblended."""
        # The footprint should have ~400 pixels
        footprint = makeDiagonalFootprint(self.peakSchema, self.peaks)
        self.assertLess(footprint.getArea(), 1000)
        self.assertGreater(footprint.getBBox().getArea(), 1000)

        task, diaSources = self.runTask(footprint, maxFootprintArea=1000)

        self.assertEqual(len(diaSources), 2)
        self.assertEqual(task.metadata["nSkipped"], 0)

    def test_compact_footprint_of_same_bbox_is_skipped(self):
        """An equivalent filled footprint will be skipped."""
        footprint = makeFootprint(self.peakSchema, self.peaks)
        self.assertGreater(footprint.getArea(), 1000)
        self.assertGreater(footprint.getBBox().getArea(), 1000)

        task, diaSources = self.runTask(footprint, maxFootprintArea=1000)

        self.assertEqual(len(diaSources), 1)
        self.assertEqual(task.metadata["nSkipped"], 1)

    def test_bbox_backstop_skips_widely_spread_footprint(self):
        """A narrow diagonal footprint will be skipped if the bbox is large
        enough.
        """
        footprint = makeDiagonalFootprint(self.peakSchema, self.peaks)
        self.assertLess(footprint.getArea(), 1000)
        self.assertGreater(footprint.getBBox().getArea(), 1000)

        task, diaSources = self.runTask(footprint, maxFootprintArea=1000,
                                        maxFootprintBBoxArea=1000)

        self.assertEqual(len(diaSources), 1)
        self.assertEqual(task.metadata["nSkipped"], 1)


class TrailTestCase(lsst.utils.tests.TestCase):
    """Tests that clusters of positive peaks along one line are kept together
    as a trailed source.

    With a 2 pixel Gaussian PSF the FWHM is ~4.7 pixels, so peaks 12 pixels
    apart are beyond the ~9.4 pixel linking radius: each forms its own cluster,
    and only the trail stage can keep them together.
    """

    def setUp(self):
        self.schema, self.peakSchema = makeSchemas()

    def runTask(self, peaks, footprint=None, **configKwargs):
        """Deblend one full-frame footprint carrying the given peaks.

        Parameters
        ----------
        peaks : `list` of `tuple`
            Each ``(ix, iy)`` for a positive peak of significance 20, or
            ``(ix, iy, isNegative[, significance])``.
        footprint : `lsst.afw.detection.Footprint`, optional
            Footprint to deblend instead, carrying its own peaks.
        """
        bbox = geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(64, 64))
        self.difference = afwImage.ExposureF(bbox)
        self.difference.setPsf(afwDetection.GaussianPsf(21, 21, 2.0))
        config = DeblendDifferenceFootprintsTask.ConfigClass()
        config.clusterRadius = 2.0
        for key, value in configKwargs.items():
            setattr(config, key, value)
        catalog = afwTable.SourceCatalog(self.schema)
        if footprint is None:
            footprintPeaks = [(peak[0], peak[1], len(peak) > 2 and peak[2],
                               peak[3] if len(peak) > 3 else 20.0) for peak in peaks]
            footprint = makeFootprint(self.peakSchema, footprintPeaks)
        catalog.addNew().setFootprint(footprint)
        self.task = DeblendDifferenceFootprintsTask(config=config)
        return self.task.run(catalog, self.difference).diaSources

    def assertPeakCounts(self, diaSources, counts):
        self.assertEqual(sorted(len(record.getFootprint().getPeaks()) for record in diaSources), counts)

    def notDeblendedPixels(self):
        bit = self.difference.mask.getPlaneBitMask("NOT_DEBLENDED")
        return (self.difference.mask.array & bit) != 0

    def test_isolated_peaks_on_a_line(self):
        """A line of isolated peaks stays one diaSource, left as detected."""
        diaSources = self.runTask([(14, 32), (26, 32), (38, 32), (50, 32)])

        self.assertPeakCounts(diaSources, [4])
        self.assertFalse(diaSources[0]["is_deblended"])
        self.assertTrue(np.all(self.notDeblendedPixels()))
        self.assertEqual(self.task.metadata["nTrailedDiaSources"], 1)

    def test_three_isolated_peaks_on_a_line(self):
        """A line of exactly the minimum number of isolated peaks is found."""
        diaSources = self.runTask([(20, 32), (32, 32), (44, 32)])

        self.assertPeakCounts(diaSources, [3])

    def test_vertical_line(self):
        """A vertical line of isolated peaks stays one diaSource."""
        diaSources = self.runTask([(32, 14), (32, 26), (32, 38), (32, 50)])

        self.assertPeakCounts(diaSources, [4])

    def test_two_isolated_peaks_split(self):
        """Two isolated peaks are too few to make a line, so they deblend."""
        diaSources = self.runTask([(20, 32), (32, 32)])

        self.assertPeakCounts(diaSources, [1, 1])
        self.assertTrue(all(record["is_deblended"] for record in diaSources))
        self.assertFalse(np.any(self.notDeblendedPixels()))
        self.assertEqual(self.task.metadata["nTrailedDiaSources"], 0)

    def test_cluster_extends_along_its_line(self):
        """A cluster of collinear positive peaks picks up isolated peaks on its
        line beyond the linking radius.

        Two isolated peaks are too few for a line of isolated peaks, so only
        the cluster's own line can pick them up.
        """
        # The cluster's peaks are 8 pixels (1.7 FWHM) apart, enough to fix a
        # direction.
        diaSources = self.runTask([(10, 32), (18, 32), (30, 32), (42, 32)])

        self.assertPeakCounts(diaSources, [4])

    def test_cluster_and_one_peak_split(self):
        """A cluster and one isolated peak on its line are too few clusters to
        make a trail, so they deblend.
        """
        diaSources = self.runTask([(10, 32), (18, 32), (30, 32)])

        self.assertPeakCounts(diaSources, [1, 2])
        self.assertEqual(self.task.metadata["nTrailedDiaSources"], 0)

    def test_cluster_absorbs_collinear_clusters(self):
        """A cluster absorbs whole clusters of positive peaks on its line."""
        diaSources = self.runTask([(4, 32), (12, 32), (24, 32), (32, 32), (44, 32), (52, 32)])

        self.assertPeakCounts(diaSources, [6])

    def test_negative_peak_blocks_absorption(self):
        """A cluster on the line that holds a negative peak stays separate,
        and the trail deblended from it is not flagged NOT_DEBLENDED.
        """
        # The negative peak is on the line, within the linking radius of the
        # last crest, so only its sign keeps the cluster out of the trail.
        diaSources = self.runTask([(8, 32), (20, 32), (32, 32), (44, 32), (56, 32), (59, 32, True)])

        self.assertPeakCounts(diaSources, [2, 4])
        self.assertTrue(all(record["is_deblended"] for record in diaSources))
        self.assertFalse(np.any(self.notDeblendedPixels()))
        self.assertEqual(self.task.metadata["nTrailedDiaSources"], 1)

    def test_isolated_negative_peak_not_in_line(self):
        """An isolated negative peak on the line is not part of the trail."""
        diaSources = self.runTask([(8, 32), (20, 32), (32, 32), (44, 32), (56, 32, True)])

        self.assertPeakCounts(diaSources, [1, 4])

    def test_two_peak_cluster_refit(self):
        """A two-peak cluster picks up peaks that lie within tolerance of the
        line refit to the trail, though not of the pair's own axis.
        """
        # (40, 34) is 0.42 FWHM off the pair's axis, and fits best of the two.
        diaSources = self.runTask([(16, 32), (24, 32), (40, 34), (56, 36)])

        self.assertPeakCounts(diaSources, [4])

    def test_two_peak_cluster_off_axis(self):
        """A two-peak cluster does not pick up a peak far off its axis."""
        diaSources = self.runTask([(16, 32), (24, 32), (28, 46)])

        self.assertPeakCounts(diaSources, [1, 2])

    def test_compact_cluster_joins_line_as_a_point(self):
        """A cluster spanning less than 1.5 FWHM joins a line of isolated peaks
        at its mean position.
        """
        diaSources = self.runTask([(8, 32), (20, 32), (31, 31), (33, 33), (44, 32)])

        self.assertPeakCounts(diaSources, [5])

    def test_compact_cluster_span(self):
        """A cluster spread across a line joins it as a point only when its
        span is below ``compactClusterSpan``.
        """
        # The cluster's peaks are 6 pixels (1.3 FWHM) apart, each 0.64 FWHM
        # off the line, but their mean is on it.
        peaks = [(20, 32), (32, 29), (32, 35), (44, 32), (56, 32)]
        self.assertPeakCounts(self.runTask(peaks), [5])
        self.assertPeakCounts(self.runTask(peaks, compactClusterSpan=1.0), [2, 3])

    def test_compact_cluster_absorbed_at_its_mean(self):
        """A trail grown from a wider cluster absorbs a compact cluster whose
        mean lies on its line, though its peaks do not.
        """
        # The middle pair's peaks are 6 pixels (1.3 FWHM) apart and 0.64 FWHM
        # off the line; the first pair's are 8 pixels (1.7 FWHM) apart, so only
        # the first fixes a direction.
        peaks = [(8, 32), (16, 32), (28, 29), (28, 35), (40, 32), (52, 32)]
        self.assertPeakCounts(self.runTask(peaks), [6])
        self.assertPeakCounts(self.runTask(peaks, compactClusterSpan=1.0), [2, 4])

    def test_tolerance_in_fwhm(self):
        """A peak 1 pixel (0.21 FWHM) off the line joins the trail, and one 2
        pixels (0.42 FWHM) off does not.
        """
        line = [(12, 32), (32, 32), (42, 32), (52, 32)]
        self.assertPeakCounts(self.runTask(line + [(22, 33)]), [5])
        self.assertPeakCounts(self.runTask(line + [(22, 34)]), [1, 4])

    def test_off_line_peak_blocks_absorption(self):
        """A cluster with a peak off the line stays separate."""
        # 8 pixels from the line: within the linking radius of the last crest,
        # but well beyond the 0.25 FWHM tolerance.
        diaSources = self.runTask([(8, 32), (20, 32), (32, 32), (44, 32), (56, 32), (56, 40)])

        self.assertPeakCounts(diaSources, [2, 4])

    def test_wide_gap_joins_line(self):
        """A collinear peak far along the line still joins the trail when the
        footprint covers the space between them.
        """
        # 32 pixels, ~6.8 FWHM, from the nearest other peak.
        diaSources = self.runTask([(4, 32), (16, 32), (28, 32), (60, 32)])

        self.assertPeakCounts(diaSources, [4])

    def notchedFootprint(self, columns, peaks):
        """Build a footprint whose band along row 32 is cut by a notch, with
        the two sides joined below it.

        Parameters
        ----------
        columns : `tuple` [`int`, `int`]
            First and last column of the notch.
        peaks : `list` of `tuple`
            Each ``(ix, iy)`` for a positive peak of significance 20.
        """
        mask = afwImage.Mask(geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(64, 64)))
        bit = mask.getPlaneBitMask("DETECTED")
        rows, cols = np.mgrid[0:64, 0:64]
        notch = (rows <= 39) & (cols >= columns[0]) & (cols <= columns[1])
        mask.array[(rows >= 30) & (rows <= 44) & ~notch] = bit
        footprint = afwDetection.Footprint(afwGeom.SpanSet.fromMask(mask, bit), self.peakSchema)
        for ix, iy in peaks:
            footprint.addPeak(ix, iy, 20.0)["significance"] = 20.0
        return footprint

    def test_footprint_gap_splits_line(self):
        """A gap in the footprint splits a line of isolated peaks into a trail
        on each side.
        """
        footprint = self.notchedFootprint((30, 36), [(3, 32), (15, 32), (27, 32),
                                                     (39, 32), (51, 32), (63, 32)])

        diaSources = self.runTask([], footprint=footprint)

        self.assertPeakCounts(diaSources, [3, 3])
        self.assertEqual(self.task.metadata["nTrailedDiaSources"], 2)

    def test_footprint_gap_stops_cluster_trail(self):
        """A cluster does not extend its line across a gap in the footprint."""
        footprint = self.notchedFootprint((22, 26), [(10, 32), (18, 32), (30, 32), (42, 32)])

        diaSources = self.runTask([], footprint=footprint)

        self.assertPeakCounts(diaSources, [1, 1, 2])
        self.assertEqual(self.task.metadata["nTrailedDiaSources"], 0)

    def test_trail_footprint_is_contiguous(self):
        """A neighbor beside a gap between a trail's crests does not cut the
        trail's footprint in two.
        """
        # A narrow band holding the trail, with a bulge for the neighbor.
        mask = afwImage.Mask(geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(64, 64)))
        bit = mask.getPlaneBitMask("DETECTED")
        rows, columns = np.mgrid[0:64, 0:64]
        mask.array[(np.abs(rows - 32) <= 2) & (columns >= 2) & (columns <= 61)] = bit
        mask.array[(rows - 37)**2 + (columns - 32)**2 <= 9] = bit
        footprint = afwDetection.Footprint(afwGeom.SpanSet.fromMask(mask, bit), self.peakSchema)
        for ix, iy in [(8, 32), (20, 32), (44, 32), (56, 32), (32, 37)]:
            footprint.addPeak(ix, iy, 20.0)["significance"] = 20.0

        diaSources = self.runTask([], footprint=footprint)

        self.assertPeakCounts(diaSources, [1, 4])
        trail = [record for record in diaSources if len(record.getFootprint().getPeaks()) == 4][0]
        self.assertEqual(len(trail.getFootprint().getSpans().split()), 1)

    def test_split_line_keeps_every_peak(self):
        """Every peak ends in exactly one cluster when the trail grown from a
        line also absorbs a peak just off that line.
        """
        # Positions in FWHM: the peak at 9.5 is off the line the finder
        # returns, but within tolerance of the line refit to the trail.
        positions = np.array([(0, 0), (3, 0), (6, 0), (9.5, 0.27), (13, 0), (16, 0), (19, 0), (22, -0.2)])

        def findLines(points):
            return [[0, 1, 2, 4, 5, 6, 7]]

        def connected(start, end):
            return True

        clusters, trails = _merge_trail_clusters([[index] for index in range(8)], positions,
                                                 np.zeros(8, dtype=bool), findLines, connected,
                                                 min_points=3, tolerance=0.25, compact_span=1.5)

        self.assertEqual(sorted(peak for cluster in clusters for peak in cluster), list(range(8)))
        self.assertEqual(sum(trails), 1)

    def test_short_trail_does_not_pivot(self):
        """A short trail absorbs the far peaks of its own line before a nearer
        peak off the line that a refit could also reach.
        """
        # Positions in FWHM along and across a trail; the last peak is a
        # bright neighbor about 1.2 FWHM off the line.
        positions = np.array([(2.85, 0.15), (4.06, -0.15), (5.56, -0.15),
                              (-2.25, -0.15), (-4.5, 0.0), (0.0, 1.21)])

        def findLines(points):
            return []

        def connected(start, end):
            return True

        clusters, trails = _merge_trail_clusters([[0, 1, 2], [3], [4], [5]], positions,
                                                 np.zeros(6, dtype=bool), findLines, connected,
                                                 min_points=3, tolerance=0.25, compact_span=1.5)

        flags = {tuple(sorted(cluster)): trail for cluster, trail in zip(clusters, trails)}
        self.assertEqual(flags, {(0, 1, 2, 3, 4): True, (5,): False})

    def test_contamination_merge_into_trail(self):
        """A faint peak merged into a trail by the contamination merge leaves
        the footprint whole, as one trailed source.
        """
        # 10 pixels off the line, beyond the linking radius, but 40 times
        # fainter, so its PSF flux is contaminated by more than 5%.
        diaSources = self.runTask([(14, 32), (26, 32), (38, 32), (50, 32), (26, 42, False, 0.5)])

        self.assertPeakCounts(diaSources, [5])
        self.assertFalse(diaSources[0]["is_deblended"])
        self.assertTrue(np.all(self.notDeblendedPixels()))
        self.assertEqual(self.task.metadata["nTrailedDiaSources"], 1)

    def test_two_trails_in_one_footprint(self):
        """Two separate lines of peaks become two trailed sources."""
        diaSources = self.runTask([(8, 12), (20, 12), (32, 12), (44, 12),
                                   (8, 52), (20, 52), (32, 52), (44, 52)])

        self.assertPeakCounts(diaSources, [4, 4])
        self.assertFalse(np.any(self.notDeblendedPixels()))
        self.assertEqual(self.task.metadata["nTrailedDiaSources"], 2)
        self.assertEqual(self.task.metadata["nDeblendedFootprints"], 1)
        self.assertEqual(self.task.metadata["nDeblendedDiaSources"], 1)

    def test_trail_and_separate_transient(self):
        """A trail is kept together while an unrelated transient is deblended
        from it, and neither is flagged NOT_DEBLENDED.
        """
        diaSources = self.runTask([(14, 32), (26, 32), (38, 32), (50, 32), (32, 50)])

        self.assertPeakCounts(diaSources, [1, 4])
        self.assertTrue(all(record["is_deblended"] for record in diaSources))
        self.assertFalse(np.any(self.notDeblendedPixels()))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
