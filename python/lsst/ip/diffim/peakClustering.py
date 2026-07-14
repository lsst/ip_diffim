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

"""Identify spatial clusters of peaks within difference image footprints.
"""

__all__ = ["PeakCluster", "FootprintPeakClusters", "cluster_peaks",
           "find_peak_clusters"]

import dataclasses

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

import lsst.afw.detection
import lsst.afw.table


# Per-peak flag set by `FootprintMergeList` marking an individual peak as
# originating from the negative-polarity detection footprint set.
_NEGATIVE_PEAK_FLAG = "merge_peak_negative"


@dataclasses.dataclass
class PeakCluster:
    """A group of peaks within a single footprint that lie close together.

    Peaks are grouped by single-linkage clustering on their positions,
    ignoring polarity, so a dipole's positive and negative lobes should fall
    into the same cluster.

    Parameters
    ----------
    peak_indices : `list` [`int`]
        Rows into the parent footprint's peak catalog
        (``footprint.getPeaks()``) that belong to this cluster.
    positions : `numpy.ndarray`
        ``(n, 2)`` array of integer peak-pixel positions ``(i_x, i_y)`` in
        the parent (image) coordinate frame.
    peak_values : `numpy.ndarray`
        ``(n,)`` array of the signed peak values.
    is_negative : `numpy.ndarray`
        ``(n,)`` boolean array, `True` for peaks from the negative detection
        set (from the ``merge_peak_negative`` flag).
    """

    peak_indices: list
    positions: np.ndarray
    peak_values: np.ndarray
    is_negative: np.ndarray


@dataclasses.dataclass
class FootprintPeakClusters:
    """The peak clusters found within one multi-peak footprint.

    Parameters
    ----------
    source : `lsst.afw.table.SourceRecord`
        The source record whose footprint was clustered.  This is the live
        record from the catalog, so a downstream deblender can modify it (and
        its footprint) in place without a separate lookup.
    clusters : `list` [`PeakCluster`]
        The clusters of peaks within the footprint.
    """

    source: lsst.afw.table.SourceRecord
    clusters: list

    @property
    def footprint(self):
        """The footprint whose peaks were clustered
        (`lsst.afw.detection.Footprint`).
        """
        return self.source.getFootprint()

    @property
    def source_id(self):
        """Id of the source record (`int`).

        Only meaningful once ids have been assigned to the catalog; while
        clustering runs inside the merge step, ids are not yet set.
        """
        return self.source.getId()

    @property
    def n_peaks(self):
        """Total number of peaks in the footprint (`int`)."""
        return sum(len(cluster.peak_indices) for cluster in self.clusters)

    @property
    def is_deblendable(self):
        """Whether the footprint splits into more than one cluster, so that
        deblending would extract more than one source (`bool`).
        """
        return len(self.clusters) >= 2


def _single_linkage_labels(positions, linking_radius):
    """Assign a cluster label to each position by single-linkage clustering.

    Two positions are linked when they lie within ``linking_radius`` of each
    other; clusters are the connected components of the resulting graph.

    Parameters
    ----------
    positions : `numpy.ndarray`
        ``(n, 2)`` array of positions.
    linking_radius : `float`
        Maximum distance, in pixels, between two peaks for them to be linked
        into the same cluster.

    Returns
    -------
    labels : `numpy.ndarray`
        ``(n,)`` integer array of cluster labels.
    """
    n = len(positions)
    if n <= 1:
        return np.zeros(n, dtype=int)
    pairs = cKDTree(positions).query_pairs(r=linking_radius, output_type="ndarray")
    if len(pairs) == 0:
        graph = coo_matrix((n, n))
    else:
        graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    return labels


def cluster_peaks(peaks, linking_radius):
    """Group the peaks of one footprint into spatial clusters.

    Peaks are clustered by single-linkage on their positions, ignoring
    polarity: two peaks within ``linking_radius`` of each other join the same
    cluster, and clustering is transitive (chaining).  A peak with no
    neighbour within the radius forms its own singleton cluster.

    Parameters
    ----------
    peaks : `lsst.afw.detection.PeakCatalog`
        Peaks of a single (merged) footprint.  The schema must carry the
        ``merge_peak_negative`` flag.
    linking_radius : `float`
        Maximum distance, in pixels, between two peaks for them to be
        clustered together.

    Returns
    -------
    clusters : `list` [`PeakCluster`]
        One entry per cluster.  Empty if ``peaks`` is empty.
    """
    n = len(peaks)
    positions = np.empty((n, 2), dtype=float)
    peak_values = np.empty(n, dtype=float)
    is_negative = np.empty(n, dtype=bool)
    negative_key = peaks.getSchema().find(_NEGATIVE_PEAK_FLAG).key
    for i, peak in enumerate(peaks):
        # Detection peaks store the integer local-maximum pixel; f_x/f_y equal
        # i_x/i_y here, so use the integer accessors to be explicit.
        positions[i] = (peak.getIx(), peak.getIy())
        peak_values[i] = peak.getPeakValue()
        is_negative[i] = peak.get(negative_key)

    labels = _single_linkage_labels(positions, linking_radius)
    clusters = []
    for label in np.unique(labels):
        (members,) = np.nonzero(labels == label)
        clusters.append(PeakCluster(
            peak_indices=members.tolist(),
            positions=positions[members],
            peak_values=peak_values[members],
            is_negative=is_negative[members],
        ))
    return clusters


def find_peak_clusters(catalog, linking_radius, *, min_peaks=2):
    """Find peak clusters in every multi-peak footprint of a catalog.

    Parameters
    ----------
    catalog : `lsst.afw.table.SourceCatalog`
        Sources whose footprints are to be searched, e.g. the merged
        ``initialDiaSources``.  Footprint peak catalogs must carry the
        ``merge_peak_negative`` flag.
    linking_radius : `float`
        Maximum distance, in pixels, between two peaks for them to be
        clustered together.
    min_peaks : `int`, optional
        Only footprints with at least this many peaks are considered; those
        with fewer are not deblending candidates and are skipped.

    Returns
    -------
    results : `list` [`FootprintPeakClusters`]
        One entry per footprint with at least ``min_peaks`` peaks.
    """
    results = []
    for source in catalog:
        footprint = source.getFootprint()
        if footprint is None:
            continue
        peaks = footprint.getPeaks()
        if len(peaks) < min_peaks:
            continue
        results.append(FootprintPeakClusters(
            source=source,
            clusters=cluster_peaks(peaks, linking_radius),
        ))
    return results
