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


def addGaussian(image, x, y, amplitude, psfSigma):
    """Add one Gaussian of the given peak amplitude to an image array."""
    yy, xx = np.mgrid[0:image.shape[0], 0:image.shape[1]]
    image += amplitude*np.exp(-((xx - x)**2 + (yy - y)**2)/(2.0*psfSigma**2))


def addPointSource(image, x, y, significance, psfSigma):
    """Add a point source of the given detection significance.

    Detection smooths with a unit-sum Gaussian of the PSF width, which halves a
    point source's peak and reduces the noise by ``sqrt(4*pi*psfSigma**2)``, so
    a raw peak of ``significance/(psfSigma*sqrt(pi))`` is detected at
    ``significance`` in an image with unit variance.
    """
    addGaussian(image, x, y, significance/(psfSigma*np.sqrt(np.pi)), psfSigma)


def addTrail(image, start, end, significance, psfSigma, step=0.25):
    """Add a PSF-convolved trail of the given ridge significance.

    A line of unit flux per pixel of length has a smoothed ridge of
    ``1/(2*psfSigma*sqrt(pi))``, the same as the smoothed image noise, so the
    linear flux density equals the ridge significance.
    """
    length = np.hypot(end[0] - start[0], end[1] - start[1])
    amplitude = significance*step/(2.0*np.pi*psfSigma**2)
    for index in range(int(length/step) + 1):
        fraction = index*step/length
        addGaussian(image, start[0] + fraction*(end[0] - start[0]),
                    start[1] + fraction*(end[1] - start[1]), amplitude, psfSigma)


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
    """Tests that trailed sources are not split into several diaSources.

    Every case uses the same two crests, far enough apart that the clustering
    and contamination steps both leave them separate, so only the trail test can
    keep them together.  The pixels then decide: a trail merges, an equivalent
    pair of point sources does not.
    """

    def setUp(self):
        self.schema, self.peakSchema = makeSchemas()
        self.psfSigma = 2.0
        # 12 pixels apart, beyond the ~9.4 pixel (2 FWHM) linking radius.
        self.crests = [(20, 32), (32, 32)]
        self.significance = 20.0

    def makeExposure(self):
        bbox = geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(64, 64))
        exposure = afwImage.ExposureF(bbox)
        exposure.setPsf(afwDetection.GaussianPsf(21, 21, self.psfSigma))
        exposure.variance.array[:] = 1.0
        return exposure

    def makePairExposure(self):
        """An exposure with flux only at the two crests."""
        exposure = self.makeExposure()
        for x, y in self.crests:
            addPointSource(exposure.image.array, x, y, self.significance, self.psfSigma)
        return exposure

    def makeTrailExposure(self):
        """An exposure with a continuous trail running through both crests."""
        exposure = self.makeExposure()
        # Extend past both crests so their smoothed flux is not cut by the ends.
        addTrail(exposure.image.array, (14, 32), (38, 32), self.significance, self.psfSigma)
        return exposure

    def runTask(self, exposure, peaks, **configKwargs):
        config = DeblendDifferenceFootprintsTask.ConfigClass()
        config.clusterRadius = 2.0
        for key, value in configKwargs.items():
            setattr(config, key, value)
        catalog = afwTable.SourceCatalog(self.schema)
        catalog.addNew().setFootprint(makeFootprint(self.peakSchema, peaks))
        task = DeblendDifferenceFootprintsTask(config=config)
        return task.run(catalog, exposure).diaSources

    def crestPeaks(self, extra=()):
        peaks = [(x, y, False, self.significance) for x, y in self.crests]
        return peaks + list(extra)

    def test_trail_stays_one_source(self):
        """Two crests joined by a continuous trail stay a single diaSource.

        The control run pins down what did the merging: the same footprint,
        peaks and config over an image with flux only at the crests deblend into
        two, so it is the flux between them that keeps them together.  A size
        guardrail would have skipped both runs alike, and leaves a footprint
        looking just like a merged trail.
        """
        self.assertEqual(len(self.runTask(self.makePairExposure(), self.crestPeaks())), 2)

        exposure = self.makeTrailExposure()
        diaSources = self.runTask(exposure, self.crestPeaks())

        self.assertEqual(len(diaSources), 1)
        self.assertFalse(diaSources[0]["is_deblended"])
        self.assertEqual(len(diaSources[0].getFootprint().getPeaks()), 2)
        # The trail is flagged NOT_DEBLENDED for the trailed source measurement.
        notDeblended = exposure.mask.getPlaneBitMask("NOT_DEBLENDED")
        self.assertTrue(np.all((exposure.mask.array & notDeblended) != 0))

    def test_point_source_pair_still_splits(self):
        """The same peaks, with no flux between them, still deblend."""
        exposure = self.makePairExposure()

        diaSources = self.runTask(exposure, self.crestPeaks())

        self.assertEqual(len(diaSources), 2)
        self.assertTrue(all(record["is_deblended"] for record in diaSources))
        notDeblended = exposure.mask.getPlaneBitMask("NOT_DEBLENDED")
        self.assertTrue(np.all((exposure.mask.array & notDeblended) == 0))

    def test_trail_and_separate_transient(self):
        """A trail merges while an unrelated transient stays its own diaSource.

        The source count cannot tell this from an ordinary two-way split, so
        check instead that NOT_DEBLENDED covers the trail's pixels and not the
        transient's.
        """
        exposure = self.makeTrailExposure()
        # Well off the trail line, so the collinearity test keeps it separate.
        transient = (32, 50)
        addPointSource(exposure.image.array, *transient, self.significance, self.psfSigma)
        peaks = self.crestPeaks([(*transient, False, self.significance)])

        diaSources = self.runTask(exposure, peaks)

        # Both crests on one diaSource, the transient on its own.
        self.assertEqual(len(diaSources), 2)
        self.assertEqual(sorted(len(record.getFootprint().getPeaks()) for record in diaSources),
                         [1, 2])
        notDeblended = exposure.mask.getPlaneBitMask("NOT_DEBLENDED")
        for x, y in self.crests:
            self.assertTrue(exposure.mask.array[y, x] & notDeblended)
        self.assertFalse(exposure.mask.array[transient[1], transient[0]] & notDeblended)

    def test_negative_peak_blocks_trail(self):
        """A cluster holding a negative peak is never merged as a trail.

        Trailed sources appear only in the science image, so a negative peak
        means the cluster is not part of one, however continuous the flux is.
        """
        exposure = self.makeTrailExposure()
        # 3 pixels from the first crest, so it joins that cluster, and on the
        # far side so it does not chain the two crests together.
        negative = [(17, 32, True, self.significance)]

        diaSources = self.runTask(exposure, self.crestPeaks(negative))

        self.assertEqual(len(diaSources), 2)

    def test_neighbor_wings_do_not_make_a_trail(self):
        """Flux from a source that is not one of the footprint's peaks does not
        merge two point sources.

        The wings of a bright neighbor leave a residual between the two peaks
        that is significant against the noise but only a small fraction of the
        peaks themselves, so only the fraction test can reject it.
        """
        exposure = self.makePairExposure()
        # Below the midpoint, and never added as a peak, so it cannot be
        # modelled: its wings leave ~3 sigma between the crests, well above
        # minTrailSignificance but far below minTrailFraction of a crest.
        addPointSource(exposure.image.array, 26, 26, 2.0*self.significance, self.psfSigma)

        self.assertEqual(len(self.runTask(exposure, self.crestPeaks())), 2)

        # Without the fraction test the same peaks merge, so it is what rejects
        # them, not the noise threshold.
        self.assertEqual(len(self.runTask(exposure, self.crestPeaks(), minTrailFraction=0.0)), 1)

    def test_off_line_peak_blocks_trail(self):
        """A cluster holding a peak off the trail line is not merged."""
        exposure = self.makeTrailExposure()
        # 8 pixels perpendicular: within the linking radius of the first crest,
        # so it joins that cluster, but well beyond trailWidth of the line.
        offLine = [(20, 40, False, self.significance)]
        addPointSource(exposure.image.array, 20, 40, self.significance, self.psfSigma)

        diaSources = self.runTask(exposure, self.crestPeaks(offLine))

        self.assertEqual(len(diaSources), 2)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
