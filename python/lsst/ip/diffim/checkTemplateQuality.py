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
import numpy as np
import scipy.ndimage

import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.afw.geom as afwGeom
import lsst.pex.config as pexConfig
import lsst.pex.exceptions as pexExceptions
import lsst.pipe.base as pipeBase

from lsst.meas.algorithms import SubtractBackgroundTask
from lsst.utils.timer import timeMethod

__all__ = [
    "CheckTemplateQualityTask",
    "CheckTemplateQualityConfig",
    "PsfDiscontinuityError",
]


class PsfDiscontinuityError(pipeBase.AlgorithmError):
    """Raised when too little of the template has a usable, smoothly-varying
    PSF for image differencing.

    The template is partitioned into regions of constant input-image depth,
    within which the PSF varies smoothly. Regions whose PSF differs from the
    dominant smooth model, and narrow shallow regions, are masked. If the
    remaining usable fraction is too small, PSF matching the template to a
    science image would leave structured residuals and produce large numbers
    of junk diaSources, so the template is rejected instead.

    Parameters
    ----------
    usableFraction : `float`
        Fraction of the template's data-bearing pixels that remain usable
        (not flagged as ``PSF_DISCONTINUITY``).
    threshold : `float`
        Configured minimum usable fraction below which the template is
        rejected.
    """
    def __init__(self, *, usableFraction, threshold):
        super().__init__(
            f"Template PSF is too discontinuous: only {usableFraction:.1%} of the data pixels"
            f" have a usable PSF (require at least {threshold:.1%}); image differencing on this"
            " template would produce a large number of junk diaSources."
        )
        self.usableFraction = usableFraction
        self.threshold = threshold

    @property
    def metadata(self) -> dict:
        return {
            "usableFraction": float(self.usableFraction),
            "threshold": float(self.threshold),
        }


class CheckTemplateQualityConfig(pexConfig.Config):
    varianceBackground = pexConfig.ConfigurableField(
        target=SubtractBackgroundTask,
        doc="Task to estimate the background variance.",
    )
    highVarianceThreshold = pexConfig.RangeField(
        dtype=float,
        default=4,
        min=1,
        doc="Set the HIGH_VARIANCE mask plane for regions with variance"
        " greater than the median by this factor.",
    )
    highVarianceMaskFraction = pexConfig.Field(
        dtype=float,
        default=0.1,
        doc="Minimum fraction of unmasked pixels needed to set the"
        " HIGH_VARIANCE mask plane.",
    )
    doMaskShallowCoverage = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Set the HIGH_VARIANCE mask plane on template regions covered by"
        " fewer input images than ``minNumberOfInputImages`` in every"
        " contributing tract. This flags shallow regions (e.g. chip gaps) where"
        " the template PSF varies discontinuously, which breaks the assumption"
        " of a smoothly-varying PSF when PSF-matching the template to a science"
        " image.",
    )
    minNumberOfInputImages = pexConfig.RangeField(
        dtype=int,
        default=2,
        min=1,
        doc="A tract is considered to cover a pixel deeply if at least this many"
        " of its distinct input images (deduplicated across overlapping patches)"
        " cover the pixel."
        " Only used if ``doMaskShallowCoverage`` is set.",
    )
    coverageGrowRadius = pexConfig.RangeField(
        dtype=int,
        default=0,
        min=0,
        doc="Grow the shallow-coverage mask region by this many pixels, to"
        " account for the spread of the PSF discontinuity by the PSF-matching"
        " convolution kernel."
        " Only used if ``doMaskShallowCoverage`` is set.",
    )
    maxInputVisitsForCoverageCheck = pexConfig.RangeField(
        dtype=int,
        default=6,
        min=1,
        doc="Skip the (expensive) shallow-coverage depth computation for a tract"
        " whose shallowest input coadd has at least this many input visits."
        " Only used if ``doMaskShallowCoverage`` is set.",
    )
    doMaskPsfDiscontinuity = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Set the PSF_DISCONTINUITY mask plane on template regions whose PSF"
        " differs substantially from the dominant smooth PSF model (e.g. across"
        " coadd patch or chip-gap boundaries), which the AL PSF-matching kernel"
        " cannot absorb. Such regions are excluded from PSF-matching kernel"
        " candidates and from difference-image source detection and measurement.",
    )
    psfDiscontinuityPolyOrder = pexConfig.RangeField(
        dtype=int,
        default=1,
        min=0,
        doc="Order of the 2D polynomial fit to the per-region PSF size, modeling"
        " the smoothly-varying component of the template PSF that AL PSF matching"
        " can absorb. Should match the ``spatialKernelOrder`` of the downstream"
        " AlardLuptonSubtractTask (default 1)."
        " Only used if ``doMaskPsfDiscontinuity`` is set.",
    )
    psfDiscontinuityThreshold = pexConfig.RangeField(
        dtype=float,
        default=0.05,
        min=0.0,
        doc="Flag a region as discontinuous if the residual of its PSF determinant"
        " radius from the smooth polynomial fit, as a fraction of the median"
        " radius, exceeds this value."
        " Only used if ``doMaskPsfDiscontinuity`` is set.",
    )
    psfDiscontinuitySigmaClip = pexConfig.RangeField(
        dtype=float,
        default=3.0,
        min=0.0,
        doc="Sigma-clipping threshold for the robust polynomial fit to the"
        " per-region PSF sizes, so the fit tracks the dominant smooth PSF."
        " Only used if ``doMaskPsfDiscontinuity`` is set.",
    )
    psfDiscontinuityNumSigmaClipIter = pexConfig.RangeField(
        dtype=int,
        default=3,
        min=0,
        doc="Number of sigma-clipping iterations for the robust polynomial fit."
        " Only used if ``doMaskPsfDiscontinuity`` is set.",
    )
    psfDiscontinuityMinUsableFraction = pexConfig.RangeField(
        dtype=float,
        default=0.1,
        min=0.0,
        max=1.0,
        doc="Raise `PsfDiscontinuityError` if the fraction of data pixels that"
        " remain usable (not flagged PSF_DISCONTINUITY) falls below this value."
        " Only used if ``doMaskPsfDiscontinuity`` is set.",
    )
    psfDiscontinuityMinRegionWidth = pexConfig.RangeField(
        dtype=int,
        default=10,
        min=1,
        doc="Depth regions narrower than this (maximum inscribed diameter, in"
        " pixels) are merged into their largest neighbor if their depth exceeds"
        " ``minNumberOfInputImages``, and masked otherwise."
        " Only used if ``doMaskPsfDiscontinuity`` is set.",
    )
    psfDiscontinuitySampleArea = pexConfig.RangeField(
        dtype=int,
        default=40000,
        min=1,
        doc="Approximate number of pixels represented by each PSF sample point"
        " in the polynomial fit. A region contributes one sample per this many"
        " pixels of its area (all sharing the region's single PSF size), and each"
        " sample is weighted by the pixels it represents (capped at this value)."
        " Only used if ``doMaskPsfDiscontinuity`` is set.",
    )

    def setDefaults(self):
        # Background subtraction of the variance plane
        self.varianceBackground.algorithm = "LINEAR"
        self.varianceBackground.binSize = 32
        self.varianceBackground.useApprox = False
        self.varianceBackground.statisticsProperty = "MEDIAN"
        self.varianceBackground.doFilterSuperPixels = True
        self.varianceBackground.ignoredPixelMask = [
            "BAD",
            "EDGE",
            "DETECTED",
            "DETECTED_NEGATIVE",
            "NO_DATA",
        ]


class CheckTemplateQualityTask(pipeBase.Task):
    """Set mask planes on an assembled template flagging regions that would
    degrade image differencing.

    Three checks are applied, each setting a mask plane that the difference
    imaging task already honors:

    - high-variance regions (``HIGH_VARIANCE``);
    - shallow-coverage regions, e.g. chip gaps, built from too few input images
      in every contributing tract (``HIGH_VARIANCE``);
    - regions whose PSF varies discontinuously and cannot be absorbed by the AL
      PSF-matching kernel (``PSF_DISCONTINUITY``).
    """
    ConfigClass = CheckTemplateQualityConfig
    _DefaultName = "checkTemplateQuality"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.makeSubtask("varianceBackground")

    @timeMethod
    def run(self, template, hasDataByTract, ccdInputCatalogsByTract):
        """Set template-quality mask planes on the template in place.

        Parameters
        ----------
        template : `lsst.afw.image.Exposure`
            The assembled template, with PSF set; modified in place.
        hasDataByTract : `dict` [`int`, `numpy.ndarray` [`bool`]]
            Per-tract boolean array (template frame) marking where each tract has
            data, captured before the tracts were merged. Used by the
            shallow-coverage check.
        ccdInputCatalogsByTract : `dict` [`int`, `list` \
                [`lsst.afw.table.ExposureCatalog`]]
            Per-tract per-CCD coadd input catalogs (``CoaddInputs.ccds``).

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            A struct with attribute:

            ``template``
                The same exposure, with quality mask planes set
                (`lsst.afw.image.Exposure`).

        Raises
        ------
        PsfDiscontinuityError
            If the PSF discontinuity check leaves too little usable template.
        """
        self.checkHighVariance(template)
        if self.config.doMaskShallowCoverage:
            self.maskShallowCoverage(template, hasDataByTract, ccdInputCatalogsByTract)
        if self.config.doMaskPsfDiscontinuity:
            allCcdInputs = [catalog for catalogs in ccdInputCatalogsByTract.values()
                            for catalog in catalogs]
            self.maskPsfDiscontinuity(template, allCcdInputs)
        return pipeBase.Struct(template=template)

    def checkHighVariance(self, template):
        """Set a mask plane for regions with unusually high variance.

        Parameters
        ----------
        template : `lsst.afw.image.Exposure`
            The warped template exposure, which will be modified in place.
        """
        highVarianceMaskPlaneBit = template.mask.addMaskPlane("HIGH_VARIANCE")
        ignoredPixelBits = template.mask.getPlaneBitMask(self.varianceBackground.config.ignoredPixelMask)
        goodMask = (template.mask.array & ignoredPixelBits) == 0
        goodFraction = np.count_nonzero(goodMask)/template.mask.array.size
        if goodFraction < self.config.highVarianceMaskFraction:
            self.log.info("Not setting HIGH_VARIANCE mask plane, only %2.1f%% of"
                          " pixels were unmasked for background estimation, but"
                          " %2.1f%% are required", 100*goodFraction, 100*self.config.highVarianceMaskFraction)
        else:
            varianceExposure = template.clone()
            varianceExposure.image.array = varianceExposure.variance.array
            varianceBackground = self.varianceBackground.run(varianceExposure).background.getImage().array
            threshold = self.config.highVarianceThreshold*np.nanmedian(varianceBackground)
            highVariancePix = varianceBackground > threshold
            template.mask.array[highVariancePix] |= 2**highVarianceMaskPlaneBit

    def maskShallowCoverage(self, template, hasDataByTract, ccdInputCatalogsByTract):
        """Flag template regions built from too few input images.

        Shallow coverage is assessed per-tract and combined: a pixel is flagged
        only if it has template data and no contributing tract covers it deeply.
        A single tract covering a pixel with sufficient depth is enough to keep
        it unflagged, so depth contributed by one tract filling another's gap in
        an overlap region is respected.

        Parameters
        ----------
        template : `lsst.afw.image.Exposure`
            The warped template exposure, modified in place.
        hasDataByTract : `dict` [`int`, `numpy.ndarray` [`bool`]]
            Per-tract boolean array marking where each tract has data.
        ccdInputCatalogsByTract : `dict` [`int`, `list` \
                [`lsst.afw.table.ExposureCatalog`]]
            Per-tract per-CCD coadd input catalogs.
        """
        bbox = template.getBBox()
        wcs = template.getWcs()
        anyData = np.zeros((bbox.getHeight(), bbox.getWidth()), dtype=bool)
        anyDeep = np.zeros((bbox.getHeight(), bbox.getWidth()), dtype=bool)
        for tract, hasData in hasDataByTract.items():
            self._accumulateShallowCoverage(tract, hasData, ccdInputCatalogsByTract[tract],
                                            wcs, bbox, anyData, anyDeep)

        noData = (template.mask.array & template.mask.getPlaneBitMask("NO_DATA")) != 0
        shallow = anyData & ~anyDeep & ~noData
        nShallow = int(np.count_nonzero(shallow))
        if nShallow == 0:
            return
        self.log.info("Flagging %d pixels (%.1f%%) covered by fewer than %d input"
                      " images in every contributing tract as HIGH_VARIANCE.",
                      nShallow, 100*nShallow/shallow.size, self.config.minNumberOfInputImages)

        highVarianceBitMask = 2**template.mask.addMaskPlane("HIGH_VARIANCE")
        shallowMask = afwImage.Mask(template.getBBox())
        shallowMask.array[shallow] = 1
        spanSet = afwGeom.SpanSet.fromMask(shallowMask, 1)
        if self.config.coverageGrowRadius > 0:
            spanSet = spanSet.dilated(self.config.coverageGrowRadius).clippedTo(template.getBBox())
        spanSet.setMask(template.mask, highVarianceBitMask)

    def _accumulateShallowCoverage(self, tract, hasData, ccdInputCatalogs, wcs, bbox,
                                   anyData, anyDeep):
        """Accumulate one tract's contribution to the shallow-coverage check.

        Parameters
        ----------
        tract : `int`
            Identifier of the tract being accumulated, for logging.
        hasData : `numpy.ndarray` [`bool`]
            Template-frame mask of where this tract has data.
        ccdInputCatalogs : `list` [`lsst.afw.table.ExposureCatalog`]
            The per-CCD coadd input catalogs (``CoaddInputs.ccds``) of this
            tract's patches.
        wcs : `lsst.afw.geom.SkyWcs`
            WCS of the template frame.
        bbox : `lsst.geom.Box2I`
            Bounding box of the template frame.
        anyData : `numpy.ndarray` [`bool`]
            Accumulator of pixels covered by any tract; updated in place.
        anyDeep : `numpy.ndarray` [`bool`]
            Accumulator of pixels covered deeply by any tract; updated in place.
        """
        anyData |= hasData

        if not ccdInputCatalogs:
            # No coverage information for this tract, so we cannot assess its
            # depth; do not flag it (treat its whole data region as deep).
            self.log.warning("No coadd input coverage available for tract %s; not"
                             " checking it for shallow coverage.", tract)
            anyDeep |= hasData
            return

        # Cheap test: the per-pixel depth cannot exceed the number of input
        # visits, so a tract whose shallowest patch has many visits is very
        # unlikely to have shallow regions. Skip the depth computation for it.
        visitCounts = [len(set(ccds["visit"])) for ccds in ccdInputCatalogs
                       if "visit" in set(ccds.schema.getNames())]
        if visitCounts and min(visitCounts) >= self.config.maxInputVisitsForCoverageCheck:
            self.log.info("Skipping shallow-coverage depth computation for tract %s:"
                          " shallowest input coadd has %d visits (>= %d).",
                          tract, min(visitCounts), self.config.maxInputVisitsForCoverageCheck)
            anyDeep |= hasData
            return

        depth = self._computeDepthMap(ccdInputCatalogs, wcs, bbox)
        anyDeep |= hasData & (depth >= self.config.minNumberOfInputImages)

    def _computeDepthMap(self, ccdInputCatalogs, wcs, bbox):
        """Compute the per-pixel number of input images covering the template.

        Each input CCD footprint is warped to the template frame and
        accumulated. Inputs are deduplicated on ``(visit, ccd)`` so that the
        same exposure contributing to overlapping patches is only counted once.

        Note that this is a purely geometric measure of coverage: it does not
        account for pixels individually rejected during coadd assembly (cosmic
        rays, saturation, etc.), so it captures missing-footprint shallowness
        (chip gaps) but not shallowness caused by per-pixel masking.

        Parameters
        ----------
        ccdInputCatalogs : `list` [`lsst.afw.table.ExposureCatalog`]
            Per-CCD coadd input catalogs (``CoaddInputs.ccds``) to accumulate.
        wcs : `lsst.afw.geom.SkyWcs`
            WCS of the template frame.
        bbox : `lsst.geom.Box2I`
            Bounding box of the template frame.

        Returns
        -------
        depth : `numpy.ndarray` [`numpy.int32`]
            Per-pixel count of distinct input images covering each pixel.
        """
        x0, y0 = bbox.getMinX(), bbox.getMinY()
        depth = np.zeros((bbox.getHeight(), bbox.getWidth()), dtype=np.int32)
        seen = set()
        for ccds in ccdInputCatalogs:
            # Deduplicate on the physical image; fall back to the record id if
            # the visit/ccd fields are somehow absent.
            useVisitCcd = {"visit", "ccd"} <= set(ccds.schema.getNames())
            for record in ccds:
                key = (record["visit"], record["ccd"]) if useVisitCcd else record.getId()
                if key in seen:
                    continue
                seen.add(key)
                recordWcs = record.getWcs()
                if recordWcs is None:
                    continue
                polygon = record.getValidPolygon()
                if polygon is None:
                    polygon = afwGeom.Polygon(geom.Box2D(record.getBBox()))
                # Inputs are stored in their own (CCD) pixel frame, so warp the
                # footprint to the template frame before rasterizing.
                transform = afwGeom.makeWcsPairTransform(recordWcs, wcs)
                polygon = polygon.transform(transform)
                polyBox = geom.Box2I(polygon.getBBox())
                polyBox.clip(bbox)
                if polyBox.isEmpty():
                    continue
                covered = polygon.createImage(polyBox).array > 0
                ymin, xmin = polyBox.getMinY() - y0, polyBox.getMinX() - x0
                sub = depth[ymin:ymin + polyBox.getHeight(), xmin:xmin + polyBox.getWidth()]
                sub[covered] += 1
        return depth

    def maskPsfDiscontinuity(self, template, ccdInputCatalogs):
        """Flag template regions whose PSF cannot be matched by a smooth kernel.

        The template is partitioned into connected regions of constant input
        image depth; within each the input set, and so the PSF, is continuous.
        Narrow regions are merged into their largest neighbor or masked. The PSF
        determinant radius is sampled once per remaining region and a smooth
        polynomial of the AL spatial-kernel order is fit (area-weighted, robust);
        regions whose PSF deviates from that model are flagged with the
        ``PSF_DISCONTINUITY`` mask plane, which the difference imaging task
        excludes from PSF-matching kernel candidates and from difference-image
        source detection and measurement.

        Parameters
        ----------
        template : `lsst.afw.image.Exposure`
            The assembled template, with PSF set; modified in place.
        ccdInputCatalogs : `list` [`lsst.afw.table.ExposureCatalog`]
            Per-CCD coadd input catalogs (``CoaddInputs.ccds``) of every coadd
            that contributed to the template.

        Raises
        ------
        PsfDiscontinuityError
            If less than ``config.psfDiscontinuityMinUsableFraction`` of the data
            pixels remain usable after masking.
        """
        if not ccdInputCatalogs:
            self.log.warning("No coadd input coverage available; skipping the PSF"
                             " discontinuity check.")
            return
        bbox = template.getBBox()
        depth = self._computeDepthMap(ccdInputCatalogs, template.getWcs(), bbox)
        noData = (template.mask.array & template.mask.getPlaneBitMask("NO_DATA")) != 0
        hasData = (depth > 0) & ~noData
        nData = int(np.count_nonzero(hasData))
        if nData == 0:
            return

        labels = self._labelConstantDepthRegions(depth, hasData)
        maskNarrow = self._resolveNarrowRegions(labels, depth)
        samples = self._samplePsfRegions(template.getPsf(), labels, bbox)
        flagged = self._flagDiscontinuousRegions(samples) if samples is not None else set()

        discontinuous = maskNarrow.copy()
        if flagged:
            discontinuous |= np.isin(labels, list(flagged))
        discontinuous &= hasData
        nMasked = int(np.count_nonzero(discontinuous))
        usableFraction = (nData - nMasked)/nData
        self.log.info("PSF discontinuity check: masking %d of %d data pixels"
                      " (%.1f%% remain usable).", nMasked, nData, 100*usableFraction)
        if usableFraction < self.config.psfDiscontinuityMinUsableFraction:
            raise PsfDiscontinuityError(usableFraction=usableFraction,
                                        threshold=self.config.psfDiscontinuityMinUsableFraction)
        if nMasked == 0:
            return
        discontinuityBit = 2**template.mask.addMaskPlane("PSF_DISCONTINUITY")
        template.mask.array[discontinuous] |= discontinuityBit

    @staticmethod
    def _labelConstantDepthRegions(depth, hasData):
        """Label connected regions of constant input-image depth.

        Within a connected region of constant depth the set of contributing
        input images -- and therefore the PSF -- is continuous, so the regions
        are the units across which the PSF may vary discontinuously.

        Parameters
        ----------
        depth : `numpy.ndarray` [`numpy.int32`]
            Per-pixel input image count.
        hasData : `numpy.ndarray` [`bool`]
            Pixels with valid template data.

        Returns
        -------
        labels : `numpy.ndarray` [`numpy.int32`]
            Region label image (0 where there is no data).
        """
        labels = np.zeros(depth.shape, dtype=np.int32)
        nextLabel = 0
        for value in np.unique(depth[hasData]):
            valueMask = hasData & (depth == value)
            components, nComponents = scipy.ndimage.label(valueMask)
            labels[valueMask] = components[valueMask] + nextLabel
            nextLabel += nComponents
        return labels

    def _resolveNarrowRegions(self, labels, depth):
        """Merge or mask depth regions too narrow to characterize.

        Regions narrower than ``config.psfDiscontinuityMinRegionWidth`` (the
        maximum inscribed diameter) are handled specially: those deeper than
        ``config.minNumberOfInputImages`` are merged into their largest
        neighboring region (``labels`` is updated in place), and shallower ones
        are returned for masking. Narrow regions are processed smallest-first so
        that slivers cascade into the dominant region.

        Parameters
        ----------
        labels : `numpy.ndarray` [`numpy.int32`]
            Region label image; modified in place by merges.
        depth : `numpy.ndarray` [`numpy.int32`]
            Per-pixel input image count.

        Returns
        -------
        maskNarrow : `numpy.ndarray` [`bool`]
            Pixels in narrow, shallow regions that should be masked directly.
        """
        minWidth = self.config.psfDiscontinuityMinRegionWidth
        minDepth = self.config.minNumberOfInputImages
        maskNarrow = np.zeros(labels.shape, dtype=bool)
        areas = np.bincount(labels.ravel())
        slices = scipy.ndimage.find_objects(labels)

        narrowDepths = {}
        for label, slc in enumerate(slices, start=1):
            if slc is None:
                continue
            regionMask = labels[slc] == label
            # Pad with a background border so the distance transform sees the
            # region edges even when the region fills its bounding box.
            distance = scipy.ndimage.distance_transform_edt(np.pad(regionMask, 1))
            if 2*distance.max() < minWidth:
                narrowDepths[label] = int(depth[slc][regionMask][0])

        for label in sorted(narrowDepths, key=lambda key: areas[key]):
            if narrowDepths[label] <= minDepth:
                maskNarrow[labels == label] = True
                labels[labels == label] = 0
                continue
            neighbor = self._largestNeighbor(labels, slices[label - 1], label, areas)
            if neighbor is not None:
                areas[neighbor] += areas[label]
                labels[labels == label] = neighbor
        return maskNarrow

    @staticmethod
    def _largestNeighbor(labels, slc, label, areas):
        """Return the label of the largest region adjacent to ``label``, or
        `None` if it has no labeled neighbor."""
        rows, cols = slc
        row0 = max(0, rows.start - 1)
        row1 = min(labels.shape[0], rows.stop + 1)
        col0 = max(0, cols.start - 1)
        col1 = min(labels.shape[1], cols.stop + 1)
        sub = labels[row0:row1, col0:col1]
        isRegion = sub == label
        ring = scipy.ndimage.binary_dilation(isRegion) & ~isRegion
        neighbors = np.unique(sub[ring])
        neighbors = neighbors[neighbors > 0]
        if len(neighbors) == 0:
            return None
        return int(neighbors[np.argmax(areas[neighbors])])

    def _samplePsfRegions(self, psf, labels, bbox):
        """Sample the template PSF size once per region.

        The PSF is assumed constant within each region, so it is evaluated once,
        at the in-region pixel nearest the centroid (guaranteed inside the
        region, which may be concave). That value is then replicated onto
        spatially-spread sample points (a lattice spaced ``~sqrt(sampleArea)``)
        so a large region constrains the polynomial fit over its extent, and
        each sample is weighted by the number of pixels it represents.

        Parameters
        ----------
        psf : `lsst.afw.detection.Psf`
            The template PSF.
        labels : `numpy.ndarray` [`numpy.int32`]
            Region label image.
        bbox : `lsst.geom.Box2I`
            Template bounding box.

        Returns
        -------
        samples : `dict` [`str`, `numpy.ndarray`] or `None`
            Arrays ``x``, ``y`` (normalized to [-1, 1] over ``bbox``),
            ``radius``, ``weight`` and ``label`` for the fit, or `None` if no
            region could be sampled.
        """
        sampleArea = self.config.psfDiscontinuitySampleArea
        step = max(1, int(round(np.sqrt(sampleArea))))
        x0, y0 = bbox.getMinX(), bbox.getMinY()
        height, width = labels.shape
        xNorm = max(width - 1, 1)
        yNorm = max(height - 1, 1)
        slices = scipy.ndimage.find_objects(labels)

        xList, yList, radiusList, weightList, labelList = [], [], [], [], []
        for label, slc in enumerate(slices, start=1):
            if slc is None:
                continue
            regionMask = labels[slc] == label
            area = int(np.count_nonzero(regionMask))
            if area == 0:
                continue
            rows, cols = np.nonzero(regionMask)
            rows = rows + slc[0].start
            cols = cols + slc[1].start
            # Evaluate the PSF once, at the in-region pixel nearest the centroid.
            nearest = np.argmin((rows - rows.mean())**2 + (cols - cols.mean())**2)
            try:
                shape = psf.computeShape(geom.Point2D(cols[nearest] + x0, rows[nearest] + y0))
            except pexExceptions.Exception:
                self.log.debug("Could not evaluate the template PSF in region %d;"
                               " excluding it from the discontinuity fit.", label)
                continue
            radius = shape.getDeterminantRadius()
            if not np.isfinite(radius):
                continue
            # Replicate the (constant) PSF value onto spread sample locations.
            lattice = np.zeros(regionMask.shape, dtype=bool)
            lattice[::step, ::step] = True
            spread = regionMask & lattice
            if spread.any():
                sRows, sCols = np.nonzero(spread)
                sRows = sRows + slc[0].start
                sCols = sCols + slc[1].start
            else:
                sRows = rows[nearest:nearest + 1]
                sCols = cols[nearest:nearest + 1]
            weight = min(sampleArea, area/len(sRows))
            xList.append(2*sCols/xNorm - 1)
            yList.append(2*sRows/yNorm - 1)
            radiusList.append(np.full(len(sRows), radius))
            weightList.append(np.full(len(sRows), weight))
            labelList.append(np.full(len(sRows), label))
        if not radiusList:
            return None
        return {
            "x": np.concatenate(xList),
            "y": np.concatenate(yList),
            "radius": np.concatenate(radiusList),
            "weight": np.concatenate(weightList),
            "label": np.concatenate(labelList),
        }

    def _flagDiscontinuousRegions(self, samples):
        """Fit a smooth polynomial to the per-region PSF sizes and flag regions
        whose PSF deviates beyond the configured threshold.

        Parameters
        ----------
        samples : `dict` [`str`, `numpy.ndarray`]
            Sample arrays from `_samplePsfRegions`.

        Returns
        -------
        flagged : `set` [`int`]
            Labels of regions to mask as PSF discontinuities.
        """
        order = self.config.psfDiscontinuityPolyOrder
        x = samples["x"]
        y = samples["y"]
        radius = samples["radius"]
        weight = samples["weight"]
        label = samples["label"]

        terms = [(p, q) for p in range(order + 1) for q in range(order + 1 - p)]
        basis = np.stack([x**p * y**q for p, q in terms], axis=-1)
        nCoeffs = len(terms)
        if len(x) <= nCoeffs:
            self.log.info("Too few region samples (%d) to fit an order-%d polynomial;"
                          " skipping PSF discontinuity residual flagging.", len(x), order)
            return set()

        good = np.ones(len(x), dtype=bool)
        coeffs = self._weightedPolyFit(basis[good], radius[good], weight[good])
        for _ in range(self.config.psfDiscontinuityNumSigmaClipIter):
            residual = radius - basis @ coeffs
            sigma = np.sqrt(np.average(residual[good]**2, weights=weight[good]))
            if sigma <= 0:
                break
            updated = good & (np.abs(residual) <= self.config.psfDiscontinuitySigmaClip*sigma)
            if np.array_equal(updated, good) or np.count_nonzero(updated) <= nCoeffs:
                break
            good = updated
            coeffs = self._weightedPolyFit(basis[good], radius[good], weight[good])

        residual = radius - basis @ coeffs
        regionLabels = np.unique(label)
        medianRadius = np.median([radius[label == value][0] for value in regionLabels])
        threshold = self.config.psfDiscontinuityThreshold
        flagged = set()
        for value in regionLabels:
            regionResidual = np.median(np.abs(residual[label == value]))
            if regionResidual/medianRadius > threshold:
                flagged.add(int(value))
        return flagged

    @staticmethod
    def _weightedPolyFit(basis, values, weights):
        """Weighted linear least-squares fit, returning the coefficients."""
        sqrtWeight = np.sqrt(weights)
        coeffs, *_ = np.linalg.lstsq(basis*sqrtWeight[:, np.newaxis], values*sqrtWeight, rcond=None)
        return coeffs
