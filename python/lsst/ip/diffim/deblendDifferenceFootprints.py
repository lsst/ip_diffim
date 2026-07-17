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

import numpy as np
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
_SIGMA_TO_FWHM = np.sqrt(8.0*np.log(2.0))


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
    """

    parent: lsst.afw.detection.Footprint
    children: list


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


def _merge_dominated_clusters(clusters, weighted, n_peaks):
    """Merge clusters until every peak owns (dominates) its own peak pixel.

    A cluster is absorbed into another when the other has a strictly greater
    summed influence at one of the first cluster's own peak pixels.

    Parameters
    ----------
    clusters : `list` [`list` [`int`]]
        Initial clusters, each a list of peak indices.
    weighted : `numpy.ndarray`
        ``(n_peaks, n_peaks)`` array where ``weighted[p, j]`` is peak ``p``'s
        influence evaluated at peak ``j``'s pixel.
    n_peaks : `int`
        Number of peaks.

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
            # Check whether the owning cluster has the most significant
            # influence at peak j's pixel; if not, merge it into the one that
            # does.
            influence = np.array([weighted[cluster, j].sum() for cluster in clusters])
            best = int(np.argmax(influence))
            if influence[best] > influence[owner[j]]:
                low, high = sorted((owner[j], best))
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


def _deblend_footprint(footprint, clusters, sigma):
    """Extract every cluster of peaks from a footprint into new footprints.

    Each peak's influence is calculated as a gaussian-approximation of the PSF
    with amplitude equal to the absolute value of the signal to noise of the
    peak on the detection image (```peak['significance']```).

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

    Returns
    -------
    result : `_DeblendedFootprint` or `None`
        The parent footprint (kept on the original record) and the new
        footprints, or `None` if the footprint did not split into more than
        one cluster.
    """
    peaks = footprint.getPeaks()
    n_peaks = len(peaks)
    peak_schema = peaks.getSchema()

    positions = np.array([(peak.getIx(), peak.getIy()) for peak in peaks], dtype=int)
    amplitudes = np.array([peak["significance"] for peak in peaks], dtype=float)

    # weighted[p, j]: influence of peak p evaluated at peak j's pixel.
    inv_two_sigma_sq = 1.0/(2.0*sigma*sigma)
    delta = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]
    dist_sq = np.sum(delta*delta, axis=2)
    weighted = amplitudes[:, np.newaxis]*np.exp(-dist_sq*inv_two_sigma_sq)

    clusters = _merge_dominated_clusters(clusters, weighted, n_peaks)

    # Extract every cluster into its own source; skip footprints that did not
    # split into more than one cluster.
    if len(clusters) < 2:
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
    return _DeblendedFootprint(parent=footprints[0], children=footprints[1:])


class DeblendDifferenceFootprintsConfig(pexConfig.Config):
    """Configuration for `DeblendDifferenceFootprintsTask`."""

    clusterRadius = pexConfig.Field(
        dtype=float,
        default=3.0,
        doc="Radius (in PSF FWHM) within which detected peaks are grouped into "
            "one cluster. Each cluster is extracted into its own diaSource."
    )
    maxFootprintArea = pexConfig.Field(
        dtype=int,
        default=50000,
        doc="Footprints whose bounding-box area (pixels) exceeds this are left "
            "un-deblended, bounding the per-footprint cost so a single large "
            "footprint cannot dominate the processing time."
    )
    maxPeaks = pexConfig.Field(
        dtype=int,
        default=40,
        doc="Footprints with more than this many peaks are left un-deblended, "
            "bounding the per-footprint deblending cost."
    )


class DeblendDifferenceFootprintsTask(pipeBase.Task):
    """Deblend merged difference-image footprints by clustering their peaks.

    Each footprint whose peaks split into more than one spatial cluster is
    replaced by one footprint per cluster: the first cluster stays on the
    original record and the rest become new diaSources.  Footprint pixels are
    assigned to the cluster of highest Gaussian-weighted influence.  Footprints
    that are too large or have too many peaks are left un-deblended.
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
            Difference image the footprints were detected on.  Only its PSF is
            used, to set the clustering and influence scale.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results as a struct with attributes:

            ``diaSources``
                The updated catalog (`lsst.afw.table.SourceCatalog`).
        """
        sigma = difference.psf.computeShape(difference.psf.getAveragePosition()).getDeterminantRadius()
        linkingRadius = self.config.clusterRadius*sigma*_SIGMA_TO_FWHM
        isDeblendedKey = diaSources.schema.find("is_deblended").key
        isNegativeKey = diaSources.schema.find("is_negative").key

        nFootprints = 0
        nChildren = 0
        nSkipped = 0
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
            if (footprint.getBBox().getArea() > self.config.maxFootprintArea
                    or len(footprint.getPeaks()) > self.config.maxPeaks):
                nSkipped += 1
                continue
            result = _deblend_footprint(footprint, clusters, sigma)
            if result is None:
                continue
            self._extractRecords(diaSources, source, result, isDeblendedKey, isNegativeKey)
            nFootprints += 1
            nChildren += len(result.children)

        self.metadata["nDeblendedFootprints"] = nFootprints
        self.metadata["nDeblendedDiaSources"] = nChildren
        self.metadata["nSkipped"] = nSkipped
        self.log.info("Deblended %d footprints into %d additional diaSources (%d skipped)",
                      nFootprints, nChildren, nSkipped)
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
