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

"""Simple footprint deblender for difference images.
"""

__all__ = ["DeblendedPeaks", "deblend_clustered_peaks"]

import dataclasses

import numpy as np

import lsst.afw.detection
import lsst.afw.geom
import lsst.afw.image


@dataclasses.dataclass
class DeblendedPeaks:
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


def _peak_amplitude(image, mask, col, row, bad_bitmask):
    """Summed absolute image flux in the 3x3 box around a peak.

    Parameters
    ----------
    image, mask : `numpy.ndarray`
        Difference-image and mask arrays.
    col, row : `int`
        Peak location in array coordinates (i.e. image coordinates minus the
        image origin).
    bad_bitmask : `int`
        Bitmask of mask planes to exclude from the sum.

    Returns
    -------
    amplitude : `float`
    """
    height, width = image.shape
    c0, c1 = max(col - 1, 0), min(col + 2, width)
    r0, r1 = max(row - 1, 0), min(row + 2, height)
    sub_image = image[r0:r1, c0:c1]
    good = (mask[r0:r1, c0:c1] & bad_bitmask) == 0
    return float(np.sum(np.abs(sub_image[good])))


def _merge_peak_clusters(peak_clusters, weighted, n_peaks):
    """Merge clusters of peaks until every peak owns (dominates) its own peak
    pixel.

    Parameters
    ----------
    peaks : `list` [`list` [`int`]]
        Initial peaks, each a list of peak indices.
    weighted : `numpy.ndarray`
        ``(n_peaks, n_peaks)`` array where ``weighted[p, j]`` is peak ``p``'s
        influence evaluated at peak ``j``'s pixel.
    n_peaks : `int`
        Number of peaks.

    Returns
    -------
    peaks : `list` [`list` [`int`]]
        The merged peaks.
    """
    peak_clusters = [list(peak_cluster) for peak_cluster in peak_clusters]
    while True:
        owner = [0]*n_peaks
        for index, peak_cluster in enumerate(peak_clusters):
            for peak in peak_cluster:
                owner[peak] = index
        merged = False
        for j in range(n_peaks):
            # Check if the current peak is the most significant influence at
            # its location. If not, add it to the more significant cluster.
            influence = np.array([weighted[peak_cluster, j].sum() for peak_cluster in peak_clusters])
            best = int(np.argmax(influence))
            # Merge only if another cluster of peaks *strictly* dominates the
            # owner at this peak pixel.
            if influence[best] > influence[owner[j]]:
                low, high = sorted((owner[j], best))
                peak_clusters[low] = peak_clusters[low] + peak_clusters[high]
                del peak_clusters[high]
                merged = True
                break
        if not merged:
            return peak_clusters


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


def deblend_clustered_peaks(footprint_clusters, image, mask, image_xy0, sigma, bad_bitmask):
    """Extract every cluster of peaks from a footprint into new footprints.

    Parameters
    ----------
    footprint_clusters : `~lsst.ip.diffim.peakClustering.FootprintPeakClusters`
        Clustering result for one multi-peak footprint.
    image, mask : `numpy.ndarray`
        Difference-image and mask arrays used to weight the peaks.
    image_xy0 : `tuple` [`int`, `int`]
        Image origin ``(x0, y0)``, so peak image coordinates can index the
        arrays.
    sigma : `float`
        Width (in pixels) of the 2D Gaussian to use to determine the relative
        significance of each peak to the pixels in the footprint.
    bad_bitmask : `int`
        Bitmask of mask planes to exclude from the peak amplitude sums.

    Returns
    -------
    result : `DeblendedPeaks` or `None`
        The parent footprint (kept on the original record) and the new
        footprints, or `None` if the footprint did not split into more than
        one cluster.
    """
    footprint = footprint_clusters.footprint
    peaks = footprint.getPeaks()
    n_peaks = len(peaks)
    peak_schema = peaks.getSchema()
    x0, y0 = image_xy0

    positions = np.zeros((n_peaks, 2), dtype=int)
    amplitudes = np.zeros(n_peaks, dtype=float)
    for i, peak in enumerate(peaks):
        ix, iy = peak.getIx(), peak.getIy()
        positions[i] = (ix, iy)
        amplitudes[i] = _peak_amplitude(image, mask, ix - x0, iy - y0, bad_bitmask)

    # weighted[p, j]: influence of peak p evaluated at peak j's pixel.
    inv_two_sigma_sq = 1.0/(2.0*sigma*sigma)
    delta = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]
    dist_sq = np.sum(delta*delta, axis=2)
    weighted = amplitudes[:, np.newaxis]*np.exp(-dist_sq*inv_two_sigma_sq)

    peak_clusters = _merge_peak_clusters(
        [cluster.peak_indices for cluster in footprint_clusters.clusters], weighted, n_peaks)

    # Extract every cluster into its own source; skip footprints that did not
    # split into more than one cluster.
    if len(peak_clusters) < 2:
        return None

    # Assign every footprint pixel to the cluster of highest influence.
    bbox = footprint.getBBox()
    min_x, min_y = bbox.getMinX(), bbox.getMinY()
    footprint_mask = _spans_to_bool(footprint.getSpans(), bbox)
    xx, yy = np.meshgrid(np.arange(min_x, bbox.getMaxX() + 1),
                         np.arange(min_y, bbox.getMaxY() + 1))
    influence = np.empty((len(peak_clusters), bbox.getHeight(), bbox.getWidth()), dtype=float)
    for index, peak_cluster in enumerate(peak_clusters):
        accumulator = np.zeros(footprint_mask.shape, dtype=float)
        for peak in peak_cluster:
            r_sq = (xx - positions[peak, 0])**2 + (yy - positions[peak, 1])**2
            accumulator += amplitudes[peak]*np.exp(-r_sq*inv_two_sigma_sq)
        influence[index] = accumulator
    labels = np.argmax(influence, axis=0)
    # Guarantee each cluster owns its own peak pixels (breaks any influence
    # ties in the owner's favour, so no cluster ends up empty).
    for index, peak_cluster in enumerate(peak_clusters):
        for peak in peak_cluster:
            labels[positions[peak, 1] - min_y, positions[peak, 0] - min_x] = index

    # One footprint per cluster, carrying that cluster's peaks.
    region = footprint.getRegion()
    footprints = []
    for index, peak_cluster in enumerate(peak_clusters):
        fp = _footprint_from_bool(footprint_mask & (labels == index), bbox, peak_schema, region)
        for peak in peak_cluster:
            fp.getPeaks().append(peaks[peak])
        footprints.append(fp)

    # The first cluster stays on the original record; the rest are new records.
    return DeblendedPeaks(parent=footprints[0], children=footprints[1:])
