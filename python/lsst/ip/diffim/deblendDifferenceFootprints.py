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
import itertools

import numpy as np
from scipy.signal import fftconvolve
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

import lsst.afw.detection
import lsst.afw.geom
import lsst.afw.image
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

# Per-peak flag set by `FootprintMergeList` marking an individual peak as
# originating from the negative-polarity detection footprint set.
_NEGATIVE_PEAK_FLAG = "merge_peak_negative"

# Conversion from a Gaussian sigma to its full width at half maximum.
_SIGMA_TO_FWHM = 2.0*np.sqrt(2.0*np.log(2.0))


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
        Whether each footprint of ``[parent] + children`` holds peaks that were
        kept together as one trailed source.
    """

    parent: lsst.afw.detection.Footprint
    children: list
    trails: list


@dataclasses.dataclass
class _TrailModel:
    """Image and PSF quantities used to test a candidate trail for continuity.

    Peak amplitudes are sampled from the smoothed image rather than recovered
    from ``peak['significance']``, so the model and the data it is subtracted
    from are on the same scale by construction, with no need to invert the local
    noise that detection divided each significance by.

    Parameters
    ----------
    image : `numpy.ndarray`
        Pixels of the difference image.
    xy0 : `tuple` [`int`]
        Pixel coordinates of ``image[0, 0]``.
    smoothing : `numpy.ndarray`
        Unit-sum Gaussian kernel matching the one detection smoothed the image
        with, so that sampled values are values of the same smoothed image the
        peaks were detected on.
    profile : `numpy.ndarray`
        Profile of a point source in the smoothed image (the PSF convolved with
        ``smoothing``), scaled to unit peak, so that scaling it by a peak's
        sampled amplitude models that peak.
    noise : `float`
        One overall standard deviation of the smoothed image, which scales the
        significance threshold.
    spacing : `float`
        Spacing, in pixels, between samples along a candidate trail.
    width : `float`
        Maximum distance, in pixels, a peak may lie from the line joining the
        ends of a candidate trail.
    threshold : `float`
        Flux, in units of ``noise``, that must be left unexplained by the point
        source model at every sample along a candidate trail.
    fraction : `float`
        Flux, as a fraction of the fainter of the two peaks being bridged, that
        must be left unexplained at every sample.
    """

    image: np.ndarray
    xy0: tuple
    smoothing: np.ndarray
    profile: np.ndarray
    noise: float
    spacing: float
    width: float
    threshold: float
    fraction: float


def _make_trail_model(difference, sigma, spacing, width, threshold, fraction):
    """Precompute the image and PSF quantities the trail test needs.

    Parameters
    ----------
    difference : `lsst.afw.image.ExposureF`
        Difference image the footprints were detected on.
    sigma : `float`
        Width (in pixels) of the PSF, which is also the width of the Gaussian
        detection smoothed the image with.
    spacing, width, threshold, fraction : `float`
        Values for the matching `_TrailModel` fields.

    Returns
    -------
    trail : `_TrailModel`
        Quantities used by `_merge_trail_clusters`.
    """
    # Detection smooths with a separable unit-sum Gaussian of the PSF width, so
    # peak significances are values of that smoothed image divided by its noise.
    offsets = np.arange(-int(3.0*sigma + 0.5), int(3.0*sigma + 0.5) + 1)
    gaussian = np.exp(-0.5*(offsets/sigma)**2)
    smoothing = np.outer(gaussian, gaussian)
    smoothing /= smoothing.sum()

    # A point source in the smoothed image has the profile of the PSF convolved
    # with that Gaussian; unit peak lets it scale by the peak significance.
    psf = difference.psf
    profile = fftconvolve(psf.computeKernelImage(psf.getAveragePosition()).array, smoothing, mode="full")
    profile /= profile.max()

    # Scale the significance threshold by one overall noise estimate rather than
    # the local variance, which is inflated around bright sources -- exactly
    # where a trail may cross one, and where a local threshold would lose it.
    variance = difference.variance.array
    good = np.isfinite(variance) & (variance > 0)
    noise = np.sqrt(np.median(variance[good])*np.sum(smoothing*smoothing)) if np.any(good) else 0.0

    bbox = difference.getBBox()
    return _TrailModel(image=difference.image.array, xy0=(bbox.getMinX(), bbox.getMinY()),
                       smoothing=smoothing, profile=profile, noise=float(noise),
                       spacing=spacing, width=width, threshold=threshold, fraction=fraction)


def _sampled_value(trail, x, y):
    """Value of the smoothed difference image at one pixel.

    Sampling a peak's own position gives its amplitude, on exactly the scale the
    model is subtracted from.  Off-image and non-finite pixels are dropped and
    the remaining kernel weights renormalized, so the estimate stays usable at
    the image edge.
    """
    half_y = (trail.smoothing.shape[0] - 1)//2
    half_x = (trail.smoothing.shape[1] - 1)//2
    row = y - trail.xy0[1]
    column = x - trail.xy0[0]
    top = max(row - half_y, 0)
    bottom = min(row + half_y + 1, trail.image.shape[0])
    left = max(column - half_x, 0)
    right = min(column + half_x + 1, trail.image.shape[1])
    if top >= bottom or left >= right:
        return 0.0
    kernel = trail.smoothing[top - row + half_y:bottom - row + half_y,
                             left - column + half_x:right - column + half_x]
    pixels = trail.image[top:bottom, left:right]
    good = np.isfinite(pixels)
    weight = kernel[good].sum()
    if weight <= 0:
        return 0.0
    return float((kernel[good]*pixels[good]).sum()/weight)


def _model_value(trail, positions, amplitudes, x, y):
    """Value at one pixel of the point source model of every peak.

    Every peak of the footprint is included: a neighbor left out of the model
    would look like flux the peaks cannot explain.
    """
    half_y = (trail.profile.shape[0] - 1)//2
    half_x = (trail.profile.shape[1] - 1)//2
    total = 0.0
    for row in range(len(positions)):
        offset_y = y - positions[row, 1] + half_y
        offset_x = x - positions[row, 0] + half_x
        if 0 <= offset_y < trail.profile.shape[0] and 0 <= offset_x < trail.profile.shape[1]:
            total += amplitudes[row]*trail.profile[offset_y, offset_x]
    return total


def _peak_pair(positions, rows_a, rows_b, farthest=False):
    """Rows of the closest (or farthest) pair of peaks between two sets."""
    delta = positions[rows_a][:, np.newaxis, :] - positions[rows_b][np.newaxis, :, :]
    dist_sq = np.sum(delta*delta, axis=2)
    flat = np.argmax(dist_sq) if farthest else np.argmin(dist_sq)
    first, second = np.unravel_index(flat, dist_sq.shape)
    return rows_a[first], rows_b[second]


def _collinear(positions, rows, first, last, tolerance):
    """Do all ``rows`` lie within ``tolerance`` of the line through two peaks?"""
    direction = (positions[last] - positions[first]).astype(float)
    length = np.hypot(*direction)
    if length == 0:
        return False
    normal = np.array([-direction[1], direction[0]])/length
    return bool(np.all(np.abs((positions[rows] - positions[first]) @ normal) <= tolerance))


def _continuous_trail(trail, positions, amplitudes, first, last):
    """Is there flux the peaks cannot explain all along the gap between them?

    The smoothed difference image is sampled at roughly ``trail.spacing``
    intervals between the two peaks and the point source model is subtracted.  A
    trailed source leaves a residual at every sample, since its ridge is flat
    while the model of its crests falls off between them; two separate point
    sources leave none.  Every sample must clear the threshold: a single filled
    sample beside a gap is not a trail.

    The residual must be both significant against the noise and a large enough
    fraction of the peaks it bridges.  The two guard different regimes: without
    the noise term, faint pairs merge on noise alone; without the fraction term,
    bright pairs merge on the wings of a neighbor too far away to be one of
    ``positions``, or on the difference between the real PSF and this model.  A
    trail leaves at least half of its crest amplitude unexplained, whereas such
    contamination leaves only a few percent.
    """
    start = positions[first].astype(float)
    end = positions[last].astype(float)
    gap = np.hypot(*(end - start))
    minimum = max(trail.threshold*trail.noise,
                  trail.fraction*min(amplitudes[first], amplitudes[last]))
    # Single-linkage clustering leaves the clusters more than the linking
    # radius apart, so there is always at least one sample between them.
    n_samples = max(1, int(round(gap/trail.spacing)) - 1)
    for step in range(1, n_samples + 1):
        point = start + (end - start)*step/(n_samples + 1)
        x = int(round(point[0]))
        y = int(round(point[1]))
        residual = (_sampled_value(trail, x, y)
                    - _model_value(trail, positions, amplitudes, x, y))
        if residual < minimum:
            return False
    return True


def _is_trail_pair(trail, positions, amplitudes, negative, cluster_a, cluster_b):
    """Do the peaks of two clusters belong to one trailed source?"""
    rows = list(cluster_a) + list(cluster_b)
    # Trailed sources appear only in the science image, so a cluster holding a
    # negative peak cannot be part of one.  A dipole's negative lobe sits on
    # top of its positive lobe, so the collinearity test below cannot catch it.
    if negative[rows].any():
        return False
    # The two most widely separated peaks give the best conditioned line;
    # taking the closest pair instead lets a cluster's own axis disagree with
    # it and trigger a spurious rejection here.
    ends = _peak_pair(positions, rows, rows, farthest=True)
    if not _collinear(positions, rows, ends[0], ends[1], trail.width):
        return False
    # Sample across the gap between the clusters, the only place a break in the
    # trail can hide; the peaks within a cluster are already close-packed.
    near = _peak_pair(positions, list(cluster_a), list(cluster_b))
    return _continuous_trail(trail, positions, amplitudes, near[0], near[1])


def _merge_trail_clusters(clusters, positions, amplitudes, negative, trail):
    """Merge clusters whose peaks lie along one continuous trailed source.

    Trailed sources are of high scientific value and are handled downstream by
    a dedicated measurement as long as their peaks share one footprint, so this
    errs towards merging.  Merging repeats until no further pair of clusters
    looks like one trail.

    Parameters
    ----------
    clusters : `list` [`list` [`int`]]
        Clusters surviving the contamination merge, each a list of peak
        indices.
    positions : `numpy.ndarray`
        ``(n_peaks, 2)`` array of peak pixel coordinates.
    amplitudes : `numpy.ndarray`
        Smoothed image value at each peak, from `_sampled_value`.
    negative : `numpy.ndarray`
        Boolean array flagging the negative-polarity peaks.
    trail : `_TrailModel`
        Image and PSF quantities for the continuity test.

    Returns
    -------
    clusters : `list` [`list` [`int`]]
        The merged clusters.
    trails : `list` [`bool`]
        Whether each returned cluster was merged as a trailed source.
    """
    clusters = [list(cluster) for cluster in clusters]
    trails = [False]*len(clusters)
    while True:
        pair = next(((first, second)
                     for first, second in itertools.combinations(range(len(clusters)), 2)
                     if _is_trail_pair(trail, positions, amplitudes, negative,
                                       clusters[first], clusters[second])), None)
        if pair is None:
            return clusters, trails
        first, second = pair
        clusters[first] = clusters[first] + clusters[second]
        trails[first] = True
        del clusters[second]
        del trails[second]


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
    if len(pairs) == 0:
        graph = coo_matrix((n, n))
    else:
        graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    return [np.nonzero(labels == label)[0].tolist() for label in np.unique(labels)]


def _merge_clusters(clusters, contamination, n_peaks, max_contamination):
    """Merge clusters whose mutual PSF-flux contamination exceeds the limit.

    Two clusters are merged when, at one cluster's peak, the summed
    contamination from the other cluster exceeds ``max_contamination`` (a
    fraction of that peak's own PSF flux).  Merging repeats until every
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


def _deblend_footprint(footprint, clusters, sigma, max_contamination=0.05, trail=None):
    """Extract every cluster of peaks from a footprint into new footprints.

    Nearby clusters are first merged when the PSF-flux contamination one would
    impose on the other exceeds ``max_contamination``, then clusters lying
    along one continuous trail are merged; each remaining cluster is finally
    extracted into its own footprint.  Peak amplitudes are the detection
    significance (``peak['significance']``), used as a proxy for relative flux.

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
    trail : `_TrailModel`, optional
        Image and PSF quantities for the trail continuity test.  Clusters lying
        along one continuous trail are kept together; if `None`, no trail test
        is applied.

    Returns
    -------
    result : `_DeblendedFootprint` or `None`
        The parent footprint (kept on the original record) and the new
        footprints, or `None` if the footprint did not split into more than one
        cluster and holds no trail.
    """
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

    if trail is None or len(clusters) < 2:
        trails = [False]*len(clusters)
    else:
        negative = np.array([peak.get(peak_schema.find(_NEGATIVE_PEAK_FLAG).key) for peak in peaks],
                            dtype=bool)
        # The trail test compares the image against a model of its peaks, so it
        # needs their amplitudes on the scale of the image rather than the
        # significances used above.
        trail_amplitudes = np.array([_sampled_value(trail, x, y) for x, y in positions])
        clusters, trails = _merge_trail_clusters(clusters, positions, trail_amplitudes, negative, trail)

    # Extract every cluster into its own source; skip footprints that did not
    # split into more than one cluster.  A footprint whose clusters all merged
    # into one trail is left exactly as detected, but still reported so that it
    # can be flagged.
    if len(clusters) < 2:
        if any(trails):
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
        influence[index] = accumulator
    labels = np.argmax(influence, axis=0)
    # Guarantee each cluster owns its own peak
    for index, cluster in enumerate(clusters):
        for peak in cluster:
            labels[positions[peak, 1] - min_y, positions[peak, 0] - min_x] = index

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
        default=400000,
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
    trailWidth = pexConfig.Field(
        dtype=float,
        default=0.5,
        doc="Maximum distance (in PSF FWHM) a peak may lie from the line "
            "joining the ends of a candidate trail. Clusters holding a peak "
            "further from that line than this are not merged as a trail."
    )
    minTrailSignificance = pexConfig.Field(
        dtype=float,
        default=2.0,
        doc="Flux, in units of the smoothed image noise, that must be left "
            "unexplained by the point source model of a candidate trail's "
            "peaks at every sample point along it. A trailed source leaves at "
            "least half its peak significance unexplained, since its ridge is "
            "flat while the model of its crests falls off between them, so any "
            "detectable trail leaves at least 2.5; two separate point sources "
            "leave nothing. Lower values keep more trails together."
    )
    minTrailFraction = pexConfig.Field(
        dtype=float,
        default=0.3,
        doc="Flux, as a fraction of the fainter of the two peaks it bridges, "
            "that must be left unexplained by the point source model at every "
            "sample point along a candidate trail. A trailed source leaves at "
            "least half of its crest amplitude unexplained, while the wings of "
            "a bright neighbor outside the footprint leave only a few percent, "
            "so this rejects bright pairs that the noise threshold alone would "
            "merge. Lower values keep more trails together."
    )

    def validate(self):
        super().validate()
        # Peaks of separate clusters are more than clusterRadius FWHM apart, so
        # sampling a candidate trail every FWHM only guarantees a sample
        # between two clusters if the clusters are at least that far apart.
        if self.clusterRadius < 1.0:
            raise pexConfig.FieldValidationError(
                DeblendDifferenceFootprintsConfig.clusterRadius, self,
                "clusterRadius must be at least 1.0 PSF FWHM, so that a candidate trail is "
                "sampled at least once between two clusters.")


class DeblendDifferenceFootprintsTask(pipeBase.Task):
    """Deblend merged difference-image footprints by clustering their peaks.

    Each footprint whose peaks split into more than one spatial cluster is
    replaced by one footprint per cluster: the first cluster stays on the
    original record and the rest become new diaSources.  Footprint pixels are
    assigned to the cluster of highest Gaussian-weighted influence.  Footprints
    that are too large or have too many peaks are left un-deblended.

    Clusters lying along one continuous trail are kept together, so that a
    trailed source is not split into several diaSources: they are of high
    scientific value, and the trailed source measurement handles them as long
    as their peaks share one footprint.
    """

    ConfigClass = DeblendDifferenceFootprintsConfig
    _DefaultName = "deblendDifferenceFootprints"

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
            clustering and influence scale, its pixels are sampled to test
            candidate trails for continuity, and its mask is updated in place:
            the ``NOT_DEBLENDED`` plane is set for footprints that can't be
            deblended and for those holding a trailed source.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results as a struct with attributes:

            ``diaSources``
                The updated catalog (`lsst.afw.table.SourceCatalog`).
        """
        sigma = difference.psf.computeShape(difference.psf.getAveragePosition()).getDeterminantRadius()
        fwhm = sigma*_SIGMA_TO_FWHM
        linkingRadius = self.config.clusterRadius*fwhm
        # Sample candidate trails every FWHM: samples closer together than that
        # are strongly correlated by the smoothing, so they add cost but no
        # independent information.
        trailModel = _make_trail_model(difference, sigma, spacing=fwhm,
                                       width=self.config.trailWidth*fwhm,
                                       threshold=self.config.minTrailSignificance,
                                       fraction=self.config.minTrailFraction)
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
            result = _deblend_footprint(footprint, clusters, sigma, self.config.maxContamination,
                                        trailModel)
            if result is None:
                continue
            if result.children:
                self._extractRecords(diaSources, source, result, isDeblendedKey, isNegativeKey)
                nFootprints += 1
                nChildren += len(result.children)
            # Flag the trailed sources, which the trailed source measurement
            # handles now that their peaks were kept in one footprint.
            for newFootprint, isTrail in zip([result.parent] + result.children, result.trails):
                if isTrail:
                    newFootprint.getSpans().setMask(difference.mask, notDeblendedBitmask)
                    nTrails += 1

        self.metadata["nDeblendedFootprints"] = nFootprints
        self.metadata["nDeblendedDiaSources"] = nChildren
        self.metadata["nSkipped"] = nSkipped
        self.metadata["nTrailedDiaSources"] = nTrails
        self.log.info("Deblended %d footprints into %d additional diaSources "
                      "(%d skipped, %d trailed sources kept together)",
                      nFootprints, nChildren, nSkipped, nTrails)
        return pipeBase.Struct(diaSources=diaSources)

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
