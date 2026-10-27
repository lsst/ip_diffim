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

"""Deblend merged difference-image footprints by clustering their peaks.
"""

__all__ = ["DeblendDifferenceFootprintsConfig", "DeblendDifferenceFootprintsTask"]

import dataclasses
import functools

from astropy.stats import gaussian_sigma_to_fwhm
import numpy as np
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

import lsst.afw.detection
import lsst.afw.geom
import lsst.afw.image
import lsst.afw.table
import lsst.geom
from lsst.meas.algorithms import FindGlintTrailsTask, fit_line
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

# Per-peak flag set by `FootprintMergeList` marking an individual peak as
# originating from the negative-polarity detection footprint set.
_NEGATIVE_PEAK_FLAG = "merge_peak_negative"


@dataclasses.dataclass
class _DeblendedFootprint:
    """The footprints resulting from deblending one multi-peak footprint.

    Parameters
    ----------
    parent : `lsst.afw.detection.Footprint`
        Footprint of the first cluster, set back on the original source record.
    children : `list` [`lsst.afw.detection.Footprint`]
        One new footprint per remaining cluster, each carrying that cluster's
        peaks.
    trails : `list` [`bool`]
        Whether each footprint of ``[parent] + children`` holds clusters that
        were merged as one trailed source.
    """

    parent: lsst.afw.detection.Footprint
    children: list
    trails: list


def _cluster_peak_indices(peaks, linking_radius):
    """Group a footprint's peaks into spatial clusters.

    Peaks are clustered with a KD tree with spacing equal to the linking radius.

    Parameters
    ----------
    peaks : `lsst.afw.detection.PeakCatalog`
        Peaks of a single (merged) footprint.
    linking_radius : `float`
        Maximum distance, in pixels, between two peaks for them to be
        clustered together.

    Returns
    -------
    clusters : `list` [`list` [`int`]]
        One entry per cluster, each a list of rows into ``peaks``.
    """
    positions = np.array([(peak.getIx(), peak.getIy()) for peak in peaks], dtype=float)
    n = len(positions)
    if n <= 1:
        return [list(range(n))]
    pairs = cKDTree(positions).query_pairs(r=linking_radius, output_type="ndarray")
    # Fill the sparse matrix in "coordinate" (COO) format
    if len(pairs) == 0:
        graph = coo_matrix((n, n))
    else:
        graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    return [np.nonzero(labels == label)[0].tolist() for label in np.unique(labels)]


def _segment_on_mask(mask, origin, scale, start, end):
    """Test whether the straight segment between two points lies entirely on
    set pixels of a mask.

    Parameters
    ----------
    mask : `numpy.ndarray` [`bool`]
        Pixels to test against.
    origin : `numpy.ndarray`
        ``(x, y)`` pixel position of ``mask[0, 0]``.
    scale : `float`
        Pixels per unit of ``start`` and ``end``.
    start, end : `numpy.ndarray`
        ``(x, y)`` positions of the ends of the segment, as pixel coordinates
        divided by ``scale``.

    Returns
    -------
    on_mask : `bool`
        Whether every pixel on an unbroken path of pixels along the segment
        is set in ``mask``.
    """
    start = start*scale - origin
    end = end*scale - origin
    # Samples at most half a pixel apart round to the same or adjacent
    # (possibly diagonal) pixels, so the pixels tested form an unbroken path.
    # A pixel whose corner the segment barely clips may not be tested.
    n_samples = int(np.ceil(2.0*np.hypot(*(end - start)))) + 1
    x, y = np.linspace(start, end, n_samples).T
    return bool(np.all(mask[np.rint(y).astype(int), np.rint(x).astype(int)]))


def _is_continuous(line, trail, candidate, connected):
    """Test whether a candidate cluster joins a trail without crossing a gap
    in its footprint.

    Parameters
    ----------
    line : `lsst.meas.algorithms.FittedLine`
        Line fit to the points of ``trail`` and ``candidate`` together.
    trail : `numpy.ndarray`
        ``(n, 2)`` array of the trail's points, as pixel coordinates divided
        by the PSF FWHM.
    candidate : `numpy.ndarray`
        ``(m, 2)`` array of the candidate cluster's points, in the same units
        as ``trail``.
    connected : callable
        Takes two positions and returns whether the segment between them lies
        on the footprint.

    Returns
    -------
    continuous : `bool`
        Whether each candidate point is joined to its neighbors along the line.
    """
    points = np.concatenate((trail, candidate))
    n_old = len(trail)
    # Sort the points along the line. Each point is projected onto the line's
    # direction, giving its distance along the trail from the line's center.
    order = np.argsort((points - line.center) @ line.direction)
    return all(connected(points[a], points[b]) for a, b in zip(order[:-1], order[1:])
               if a >= n_old or b >= n_old)


def _split_at_gaps(members, points, connected):
    """Split a line of single-point clusters wherever the footprint has a gap
    between neighbors along the line.

    Parameters
    ----------
    members : `list` [`int`]
        Clusters on the line, each standing as one point.
    points : `list` [`numpy.ndarray`]
        ``(1, 2)`` array of each cluster's point, as pixel coordinates divided
        by the PSF FWHM.
    connected : callable
        Takes two positions and returns whether the segment between them lies
        on the footprint.

    Returns
    -------
    runs : `list` [`list` [`int`]]
        The members of each continuous piece of the line, in order along it.
    """
    xy = np.concatenate([points[index] for index in members])
    line = fit_line(*xy.T)
    order = np.argsort((xy - line.center) @ line.direction)
    runs = [[members[order[0]]]]
    for a, b in zip(order[:-1], order[1:]):
        if not connected(xy[a], xy[b]):
            runs.append([])
        runs[-1].append(members[b])
    return runs


def _grow_trail(members, points, absorbed, positive, tolerance, connected):
    """Find the positive clusters that extend a trail while all its points
    stay collinear and joined by the footprint.

    A cluster is added when every point of the trail and the cluster lies
    within ``tolerance`` of the line fit to them all, and the footprint has no
    gaps between the cluster's points and its neighbors. Clusters are added one
    at a time, best fit first, so that clusters on the line fix its direction
    before a nearby peak off the line is tested; such a peak then no longer
    fits and is left out of the trail.

    Parameters
    ----------
    members : `list` [`int`]
        Clusters already in the trail.
    points : `list` [`numpy.ndarray`]
        ``(n, 2)`` array of the points that stand for each cluster on a line,
        as pixel coordinates divided by the PSF FWHM.
    absorbed : `list` [`bool`]
        Whether each cluster already belongs to another.
    positive : `list` [`bool`]
        Whether all of each cluster's peaks are positive.
    tolerance : `float`
        Maximum distance, in PSF FWHM, of any point from the trail's line.
    connected : callable
        Takes two positions and returns whether the segment between them lies
        on the footprint.

    Returns
    -------
    members : `list` [`int`]
        Clusters of the grown trail, starting with the original ``members``.
    """
    members = list(members)
    trail = np.concatenate([points[index] for index in members])
    while True:
        best = None
        best_residual = tolerance
        for index in range(len(points)):
            if index in members or absorbed[index] or not positive[index]:
                continue
            # Fit a new line to the trail and this cluster. A line fit to only
            # a few closely spaced points on a trail can miss a distant cluster
            # on the same trail, even if they have a good simultaneous fit.
            xy = np.concatenate((trail, points[index]))
            line = fit_line(*xy.T)
            residual = np.max(line.distance(*xy.T))
            if residual <= best_residual and _is_continuous(line, trail, points[index], connected):
                best = index
                best_residual = residual
        if best is None:
            return members
        members.append(best)
        trail = np.concatenate((trail, points[best]))


def _merge_trail_clusters(clusters, positions, negative, find_lines, connected, min_points, tolerance,
                          compact_span):
    """Merge clusters of positive peaks that lie along one line.

    Find lines among the positive clusters of peaks, and grow each into a
    trail that absorbs the positive clusters on its line (see `_grow_trail`).
    A cluster whose peaks span less than ``compact_span`` is treated as a
    single point at its mean position, while a wider cluster starts a line
    along its own peaks. A trail is kept only if it holds at least
    ``min_points`` clusters. Clusters with a negative peak are never part of a
    trail, and a trail never spans a gap in the footprint.

    Parameters
    ----------
    clusters : `list` [`list` [`int`]]
        Clusters of peaks within one footprint, as rows into its peaks.
    positions : `numpy.ndarray`
        ``(n_peaks, 2)`` array of the peaks' positions, as pixel coordinates
        divided by the PSF FWHM.
    negative : `numpy.ndarray` [`bool`]
        Whether each peak is negative.
    find_lines : callable
        Takes an ``(n, 2)`` array of positions and returns a list of lines,
        each a list of rows into that array.
    connected : callable
        Takes two ``(x, y)`` positions, in the same units as ``positions``,
        and returns whether the straight segment between them lies on the
        footprint.
    min_points : `int`
        Minimum number of clusters in a trail.
    tolerance : `float`
        Maximum distance, in PSF FWHM, of any point from a trail's line.
    compact_span : `float`
        Clusters whose peaks span less than this distance (in PSF FWHM) will
        be considered as a single point for line fitting.

    Returns
    -------
    clusters : `list` [`list` [`int`]]
        The merged clusters.
    trails : `list` [`bool`]
        Whether each returned cluster was merged as a trailed source.
    """
    clusters = [list(cluster) for cluster in clusters]
    positive = [not np.any(negative[cluster]) for cluster in clusters]
    absorbed = [False]*len(clusters)
    trails = [False]*len(clusters)
    compact = [np.max(np.hypot(*(positions[cluster][:, np.newaxis] - positions[cluster]).T)) < compact_span
               for cluster in clusters]
    points = [positions[cluster].mean(axis=0, keepdims=True) if small else positions[cluster]
              for cluster, small in zip(clusters, compact)]

    def join(members):
        """Merge the clusters of a trail into its first."""
        seed = members[0]
        for index in members[1:]:
            clusters[seed] = clusters[seed] + clusters[index]
            points[seed] = np.concatenate((points[seed], points[index]))
            absorbed[index] = True
        trails[seed] = True

    # Larger clusters fix their line better, so they seed first.
    for index in np.argsort([-len(cluster) for cluster in clusters], kind="stable"):
        if not compact[index] and positive[index] and not absorbed[index]:
            members = _grow_trail([index], points, absorbed, positive, tolerance, connected)
            if len(members) >= min_points:
                join(members)

    candidates = [index for index in range(len(clusters))
                  if compact[index] and positive[index] and not absorbed[index]]
    if len(candidates) >= min_points:
        for line in find_lines(np.concatenate([points[index] for index in candidates])):
            # Lines may share points, and a trail grown from an earlier line
            # may already hold some of this one.
            members = [candidates[row] for row in line
                       if not absorbed[candidates[row]] and not trails[candidates[row]]]
            if len(members) < min_points:
                continue
            for run in _split_at_gaps(members, points, connected):
                # A trail grown from an earlier run may hold some of this one.
                run = [index for index in run if not absorbed[index]]
                if len(run) >= min_points:
                    join(_grow_trail(run, points, absorbed, positive, tolerance, connected))

    keep = [index for index in range(len(clusters)) if not absorbed[index]]
    return [clusters[index] for index in keep], [trails[index] for index in keep]


def _merge_clusters(clusters, contamination, n_peaks, max_contamination):
    """Merge clusters whose mutual PSF-flux contamination exceeds the limit.

    Two clusters are merged when, at one cluster's peak, the summed
    contamination from the other cluster exceeds a ``max_contamination``
    fraction of that peak's own PSF flux.  Merging repeats until every
    surviving cluster's peaks are contaminated by less than the limit.

    Parameters
    ----------
    clusters : `list` [`list` [`int`]]
        Initial clusters, each a list of peak indices.
    contamination : `numpy.ndarray`
        ``(n_peaks, n_peaks)`` array where ``contamination[p, j]`` is the
        fractional PSF-flux contamination peak ``p`` imposes on peak ``j``.
    n_peaks : `int`
        Number of peaks.
    max_contamination : `float`
        Maximum fractional contamination a peak may receive from another
        cluster before the two clusters are merged.

    Returns
    -------
    clusters : `list` [`list` [`int`]]
        The merged clusters.
    """
    clusters = [list(cluster) for cluster in clusters]
    while True:
        owner = [0]*n_peaks
        for index, cluster in enumerate(clusters):
            for peak in cluster:
                owner[peak] = index
        merged = False
        for j in range(n_peaks):
            # Total contamination of peak j from each cluster; a peak is not
            # contaminated by its own cluster.
            received = np.array([contamination[cluster, j].sum() for cluster in clusters])
            received[owner[j]] = 0.0
            worst = int(np.argmax(received))
            if received[worst] > max_contamination:
                low, high = sorted((owner[j], worst))
                clusters[low] = clusters[low] + clusters[high]
                del clusters[high]
                merged = True
                break
        if not merged:
            return clusters


def _spans_to_bool(span_set, bbox):
    """Rasterize a SpanSet to a boolean array over ``bbox``."""
    array = np.zeros((bbox.getHeight(), bbox.getWidth()), dtype=bool)
    min_x, min_y = bbox.getMinX(), bbox.getMinY()
    for span in span_set:
        array[span.getY() - min_y, span.getX0() - min_x:span.getX1() - min_x + 1] = True
    return array


def _footprint_from_bool(pixels, bbox, peak_schema, region):
    """Build a peak-schema-preserving footprint from a boolean pixel array."""
    mask = lsst.afw.image.Mask(bbox)
    bit = mask.getPlaneBitMask("DETECTED")
    mask.array[:] = 0
    mask.array[pixels] = bit
    footprint = lsst.afw.detection.Footprint(lsst.afw.geom.SpanSet.fromMask(mask, bit), peak_schema)
    footprint.setRegion(region)
    return footprint


def _deblend_footprint(footprint, clusters, sigma, max_contamination=0.05, trails=None):
    """Extract every cluster of peaks from a footprint into new footprints.

    Nearby clusters are first merged when the PSF-flux contamination one would
    impose on the other exceeds ``max_contamination``; each remaining cluster
    is then extracted into its own footprint.  Peak amplitudes are the
    detection significance (``peak['significance']``), used as a proxy for
    relative flux.  A cluster that absorbs a trailed source is itself a
    trailed source.

    Parameters
    ----------
    footprint : `lsst.afw.detection.Footprint`
        The multi-peak footprint to deblend.
    clusters : `list` [`list` [`int`]]
        Clusters of peaks within the footprint, as rows into
        ``footprint.getPeaks()``.
    sigma : `float`
        Width (in pixels) of the 2D Gaussian used to weight each peak's
        influence on the footprint pixels.
    max_contamination : `float`, optional
        Maximum fractional PSF-flux contamination a source may receive from a
        neighboring cluster before the two are merged.
    trails : `list` [`bool`], optional
        Whether each of ``clusters`` is a trailed source; none are if `None`.

    Returns
    -------
    result : `_DeblendedFootprint` or `None`
        The parent footprint (kept on the original record) and the new
        footprints, or `None` if the footprint did not split into more than
        one cluster and holds no trailed source.  A footprint whose clusters
        all merged into one trailed source is returned unchanged as the
        parent, with no children.
    """
    if trails is None:
        trails = [False]*len(clusters)
    trail_peaks = {peak for cluster, trail in zip(clusters, trails) if trail for peak in cluster}
    peaks = footprint.getPeaks()
    n_peaks = len(peaks)
    peak_schema = peaks.getSchema()

    positions = np.array([(peak.getIx(), peak.getIy()) for peak in peaks], dtype=int)
    amplitudes = np.array([peak["significance"] for peak in peaks], dtype=float)

    # inv_two_sigma_sq weights the PSF *value* at a pixel (for pixel
    # assignment); inv_four_sigma_sq is the PSF-flux *overlap* of two peaks
    # (for the contamination merge).
    inv_two_sigma_sq = 1.0/(2.0*sigma*sigma)
    inv_four_sigma_sq = 1.0/(4.0*sigma*sigma)
    delta = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]
    dist_sq = np.sum(delta*delta, axis=2)
    # contamination[p, j]: fractional PSF-flux bias peak p imposes on peak j --
    # the Gaussian-PSF flux overlap scaled by the significance (flux) ratio.
    contamination = (amplitudes[:, np.newaxis]/amplitudes[np.newaxis, :])*np.exp(-dist_sq*inv_four_sigma_sq)

    clusters = _merge_clusters(clusters, contamination, n_peaks, max_contamination)
    trails = [any(peak in trail_peaks for peak in cluster) for cluster in clusters]

    # Extract every cluster into its own source; skip footprints that did not
    # split into more than one cluster.
    if len(clusters) < 2:
        if trails[0]:
            return _DeblendedFootprint(parent=footprint, children=[], trails=[True])
        return None

    # Assign every footprint pixel to the cluster of highest influence.
    bbox = footprint.getBBox()
    min_x, min_y = bbox.getMinX(), bbox.getMinY()
    footprint_mask = _spans_to_bool(footprint.getSpans(), bbox)
    xx, yy = np.meshgrid(np.arange(min_x, bbox.getMaxX() + 1),
                         np.arange(min_y, bbox.getMaxY() + 1))
    influence = np.empty((len(clusters), bbox.getHeight(), bbox.getWidth()), dtype=float)
    for index, cluster in enumerate(clusters):
        accumulator = np.zeros(footprint_mask.shape, dtype=float)
        for peak in cluster:
            r_sq = (xx - positions[peak, 0])**2 + (yy - positions[peak, 1])**2
            accumulator += amplitudes[peak]*np.exp(-r_sq*inv_two_sigma_sq)
        if trails[index]:
            # Spread a trail's influence along the segment between its end
            # peaks, so that a neighbor beside a gap between its crests does
            # not cut its footprint in two.
            rows = [peak for peak in cluster if peak in trail_peaks]
            line = fit_line(positions[rows, 0], positions[rows, 1])
            along = (positions[rows] - line.center) @ line.direction
            dx = xx - line.center[0]
            dy = yy - line.center[1]
            t = np.clip(dx*line.direction[0] + dy*line.direction[1], along.min(), along.max())
            d_sq = (dx - t*line.direction[0])**2 + (dy - t*line.direction[1])**2
            accumulator = np.maximum(accumulator,
                                     np.median(amplitudes[rows])*np.exp(-d_sq*inv_two_sigma_sq))
        influence[index] = accumulator
    labels = np.argmax(influence, axis=0)
    # Guarantee each cluster owns its own peak
    for index, cluster in enumerate(clusters):
        for peak in cluster:
            labels[positions[peak, 1] - min_y, positions[peak, 0] - min_x] = index
    # Far from every peak the influences are tiny and can cross over, leaving
    # a trail stray pixels cut off from its line.  Give any piece of a trail
    # holding none of its peaks to the next most influential cluster.
    for index, cluster in enumerate(clusters):
        if not trails[index]:
            continue
        pieces, _ = ndimage.label(footprint_mask & (labels == index))
        held = [pieces[positions[peak, 1] - min_y, positions[peak, 0] - min_x] for peak in cluster]
        stray = (pieces > 0) & ~np.isin(pieces, held)
        if np.any(stray):
            influence[index][stray] = -np.inf
            labels[stray] = np.argmax(influence[:, stray], axis=0)

    # One footprint per cluster, carrying that cluster's peaks.
    region = footprint.getRegion()
    footprints = []
    for index, cluster in enumerate(clusters):
        fp = _footprint_from_bool(footprint_mask & (labels == index), bbox, peak_schema, region)
        for peak in cluster:
            fp.getPeaks().append(peaks[peak])
        footprints.append(fp)

    # The first cluster stays on the original record; the rest are new records.
    return _DeblendedFootprint(parent=footprints[0], children=footprints[1:], trails=trails)


class DeblendDifferenceFootprintsConfig(pexConfig.Config):
    """Configuration for `DeblendDifferenceFootprintsTask`."""

    clusterRadius = pexConfig.Field(
        dtype=float,
        default=2.,
        doc="Initial clustering radius (in PSF FWHM) to group tightly-packed "
            "peaks. Clusters may ultimately merge with significantly "
            "contaminating neighbors. Each final cluster is extracted into its "
            "own diaSource."
    )
    maxFootprintArea = pexConfig.Field(
        dtype=int,
        default=100000,
        doc="Skip deblending footprints larger than this, in pixels."
    )
    maxFootprintBBoxArea = pexConfig.Field(
        dtype=int,
        default=600000,
        doc="Additional check to skip deblending footprints with an enclosing "
            "bounding box larger than this, in pixels."
    )
    maxPeaks = pexConfig.Field(
        dtype=int,
        default=200,
        doc="Footprints with more than this many peaks are left un-deblended."
    )
    maxContamination = pexConfig.Field(
        dtype=float,
        default=0.05,
        doc="Maximum fractional PSF-flux contamination a diaSource may receive "
            "from a neighboring cluster before the two are merged into one "
            "(un-deblended) source. For example, 0.05 keeps each deblended "
            "source's PSF flux biased by less than 5% by its neighbors. "
            "Contamination is the Gaussian-PSF flux overlap scaled by the "
            "neighbor-to-source significance ratio."
    )
    compactClusterSpan = pexConfig.Field(
        dtype=float,
        default=1.5,
        doc="Clusters whose peaks span less than this distance (in PSF FWHM) "
            "will be treated as one point at their mean position."
    )
    findTrails = pexConfig.ConfigurableField(
        target=FindGlintTrailsTask,
        doc="Subtask to find lines of isolated positive peaks within a "
            "footprint, which are kept together as one trailed source. Peak "
            "positions are passed in units of the PSF FWHM, so its ``radius`` "
            "and ``threshold`` are in FWHM rather than pixels; ``threshold`` "
            "also limits how far any peak of a trail may lie from its line."
    )

    def setDefaults(self):
        super().setDefaults()
        self.findTrails.bad_flags = []
        # Larger than any footprint the size limits allow.
        self.findTrails.radius = 1000.0
        self.findTrails.min_points = 3
        # Peaks of detected trails lie within ~0.25 FWHM of the trail's line.
        self.findTrails.threshold = 0.25


class DeblendDifferenceFootprintsTask(pipeBase.Task):
    """Deblend merged difference-image footprints by clustering their peaks.

    Each footprint whose peaks split into more than one spatial cluster is
    replaced by one footprint per cluster: the first cluster stays on the
    original record and the rest become new diaSources.  Footprint pixels are
    assigned to the cluster of highest Gaussian-weighted influence.  Footprints
    that are too large or have too many peaks are left un-deblended.

    Clusters of positive peaks lying along one line are merged first, so that
    a trailed source is not split into several diaSources: the trailed source
    measurement handles them as long as their peaks share one footprint.
    """

    ConfigClass = DeblendDifferenceFootprintsConfig
    _DefaultName = "deblendDifferenceFootprints"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask("findTrails")
        self._lineSchema = lsst.afw.table.SourceTable.makeMinimalSchema()
        self._lineCentroidKey = lsst.afw.table.Point2DKey.addFields(
            self._lineSchema, "centroid", "Peak position, as pixel coordinates divided by the PSF FWHM", "")
        self._lineSchema.getAliasMap().set("slot_Centroid", "centroid")

    def run(self, diaSources, difference):
        """Deblend the multi-cluster footprints of a diaSource catalog.

        Parameters
        ----------
        diaSources : `lsst.afw.table.SourceCatalog`
            The merged diaSources; modified in place and possibly extended.
            The schema must carry the ``is_deblended`` and ``is_negative``
            flags, and footprint peak catalogs the ``merge_peak_negative`` and
            ``significance`` fields.
        difference : `lsst.afw.image.ExposureF`
            Difference image the footprints were detected on.  Its PSF sets the
            clustering and influence scale, and its mask is updated in place:
            the ``NOT_DEBLENDED`` plane is set for footprints that can't be
            deblended, and for those left as detected because all their peaks
            belong to one trailed source.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results as a struct with attributes:

            ``diaSources``
                The updated catalog (`lsst.afw.table.SourceCatalog`).
        """
        sigma = difference.psf.computeShape(difference.psf.getAveragePosition()).getDeterminantRadius()
        fwhm = sigma*gaussian_sigma_to_fwhm
        linkingRadius = self.config.clusterRadius*fwhm
        isDeblendedKey = diaSources.schema.find("is_deblended").key
        isNegativeKey = diaSources.schema.find("is_negative").key
        difference.mask.addMaskPlane("NOT_DEBLENDED")
        notDeblendedBitmask = difference.mask.getPlaneBitMask("NOT_DEBLENDED")

        nFootprints = 0
        nChildren = 0
        nSkipped = 0
        nTrails = 0
        # Snapshot the records up front so appending children below does not
        # disturb this iteration.
        for source in list(diaSources):
            footprint = source.getFootprint()
            if footprint is None or len(footprint.getPeaks()) < 2:
                continue
            clusters = _cluster_peak_indices(footprint.getPeaks(), linkingRadius)
            if len(clusters) < 2:
                continue
            # Leave very large footprints un-deblended
            if (footprint.getArea() > self.config.maxFootprintArea
                    or footprint.getBBox().getArea() > self.config.maxFootprintBBoxArea
                    or len(footprint.getPeaks()) > self.config.maxPeaks):
                nSkipped += 1
                self.log.debug("Skipped footprint with area %d, bounding box area %d, and %d peaks.",
                               footprint.getArea(), footprint.getBBox().getArea(),
                               len(footprint.getPeaks()))
                footprint.getSpans().setMask(difference.mask, notDeblendedBitmask)
                continue
            peaks = footprint.getPeaks()
            negativeKey = peaks.getSchema().find(_NEGATIVE_PEAK_FLAG).key
            positions = np.array([(peak.getIx(), peak.getIy()) for peak in peaks], dtype=float)/fwhm
            negative = np.array([peak.get(negativeKey) for peak in peaks], dtype=bool)
            bbox = footprint.getBBox()
            connected = functools.partial(_segment_on_mask, _spans_to_bool(footprint.getSpans(), bbox),
                                          np.array([bbox.getMinX(), bbox.getMinY()]), fwhm)
            clusters, trails = _merge_trail_clusters(clusters, positions, negative, self._findLines,
                                                     connected, self.findTrails.config.min_points,
                                                     self.findTrails.config.threshold,
                                                     self.config.compactClusterSpan)
            result = _deblend_footprint(footprint, clusters, sigma, self.config.maxContamination, trails)
            if result is None:
                continue
            nTrails += sum(result.trails)
            if result.children:
                self._extractRecords(diaSources, source, result, isDeblendedKey, isNegativeKey)
                nFootprints += 1
                nChildren += len(result.children)
            else:
                # All of the peaks belong to one trailed source, left as
                # detected for the trailed source measurement.
                footprint.getSpans().setMask(difference.mask, notDeblendedBitmask)

        self.metadata["nDeblendedFootprints"] = nFootprints
        self.metadata["nDeblendedDiaSources"] = nChildren
        self.metadata["nSkipped"] = nSkipped
        self.metadata["nTrailedDiaSources"] = nTrails
        self.log.info("Deblended %d footprints into %d additional diaSources "
                      "(%d skipped, %d trailed sources kept together)",
                      nFootprints, nChildren, nSkipped, nTrails)
        return pipeBase.Struct(diaSources=diaSources)

    def _findLines(self, positions):
        """Find lines of points with the ``findTrails`` subtask.

        Parameters
        ----------
        positions : `numpy.ndarray`
            ``(n, 2)`` array of the points' positions, as pixel coordinates
            divided by the PSF FWHM.

        Returns
        -------
        lines : `list` [`list` [`int`]]
            Rows of ``positions`` on each line found.  A row may appear on
            more than one line.
        """
        catalog = lsst.afw.table.SourceCatalog(self._lineSchema)
        for x, y in positions:
            catalog.addNew().set(self._lineCentroidKey, lsst.geom.Point2D(x, y))
        rows = {record.getId(): row for row, record in enumerate(catalog)}
        trails = self.findTrails.run(catalog).trails
        return [[rows[sourceId] for sourceId in trail["id"]] for trail in trails]

    def _extractRecords(self, diaSources, source, result, isDeblendedKey, isNegativeKey):
        """Write one deblended footprint per cluster back to the catalog.

        The first cluster stays on ``source``; the rest become new records
        appended to ``diaSources``.  Every resulting record is flagged
        ``is_deblended``, and ``is_negative`` is set when all of that
        footprint's peaks are negative.
        """
        records = [source] + [diaSources.addNew() for _ in result.children]
        negativePeakKey = result.parent.getPeaks().getSchema().find(_NEGATIVE_PEAK_FLAG).key
        for record, newFootprint in zip(records, [result.parent] + result.children):
            record.setFootprint(newFootprint)
            record.set(isDeblendedKey, True)
            record.set(isNegativeKey,
                       all(peak.get(negativePeakKey) for peak in newFootprint.getPeaks()))
