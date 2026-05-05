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
import scipy.signal

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.geom as geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.afw.math._warper import computeWarpedBBox
from lsst.meas.algorithms import KernelPsf
from lsst.utils.timer import timeMethod

from .getTemplate import (
    GetTemplateConfig,
    GetTemplateConnections,
    GetTemplateTask,
)

__all__ = [
    "GetPsfMatchedTemplateTask",
    "GetPsfMatchedTemplateConfig",
]


class GetPsfMatchedTemplateConfig(GetTemplateConfig, pipelineConnections=GetTemplateConnections):
    targetPsfSigmaPad = pexConfig.Field(
        dtype=float,
        default=0.05,
        doc="Fractional tolerance for treating a patch as already at the "
        "target PSF width: any patch whose Gaussian-equivalent sigma is "
        "within this fraction of the maximum accepted sigma is passed "
        "through without convolution.",
    )
    outlierRejectionSigma = pexConfig.Field(
        dtype=float,
        default=3.0,
        doc="Reject input patches whose Gaussian PSF sigma is more than "
        "this many robust standard deviations *above* the median. "
        "Below-median sigmas are always kept, since a narrow-PSF outlier "
        "is harmless under the matching convolution.",
    )
    outlierRejectionMaxIter = pexConfig.Field(
        dtype=int,
        default=3,
        doc="Maximum number of upper-outlier rejection iterations.",
    )
    minAcceptedPatches = pexConfig.RangeField(
        dtype=int,
        default=1,
        min=1,
        doc="Minimum number of input patches that must survive PSF "
        "outlier rejection for a template to be produced.",
    )
    matchingKernelSigmas = pexConfig.RangeField(
        dtype=float,
        default=4.0,
        min=1.0,
        doc="Half-width of the Gaussian matching kernel, in units of the "
        "kernel's own sigma. The kernel will span this many sigma on each "
        "side of the center.",
    )
    minMatchingKernelSigma = pexConfig.Field(
        dtype=float,
        default=5e-2,
        doc="Floor on the Gaussian matching kernel sigma in pixels. Patches "
        "whose intrinsic sigma is already within this distance of the "
        "target are passed through without convolution.",
    )


class GetPsfMatchedTemplateTask(GetTemplateTask):
    """Build a template by convolving each input *patch* with a Gaussian
    matching kernel so that all accepted patches share a common
    Gaussian-approximated PSF width before being merged and warped.

    The flow extends `GetTemplateTask` by inserting a per-patch
    PSF-matching stage *before* the per-tract merge:

    1. *Metadata pass* — for every input patch (across all tracts),
       fetch only the ``psf``, ``wcs``, ``photoCalib``, and ``bbox``
       components from the butler (no maskedImage I/O). Estimate each
       patch's Gaussian-equivalent PSF sigma from this PSF.
    2. Apply one-sided sigma clipping to the global distribution of
       patch sigmas: drop non-finite/non-positive entries, and reject
       any patches whose sigma is more than
       ``outlierRejectionSigma`` robust standard deviations *above*
       the median. Below-median sigmas are always kept, since a
       narrow-PSF outlier is harmless under the matching convolution.
    3. Choose a single target Gaussian sigma equal to the max accepted
       sigma. Patches whose sigma is within ``targetPsfSigmaPad``
       (fractionally) of the target are passed through unconvolved.
    4. *Per-tract data pass* — for each tract, load the full Exposure
       only for accepted patches, convolve each `MaskedImage` with a 2D
       Gaussian kernel of sigma ``sqrt(target**2 - intrinsic**2)``, and
       replace the patch's catalog record PSF with a `KernelPsf`
       representing the original PSF convolved with the same Gaussian.
       The matched patches are immediately merged and warped, then
       freed before moving to the next tract.

    Splitting the metadata and pixel-data loads keeps peak memory at
    one tract's worth of accepted patches (matching the parent task's
    profile) and avoids reading rejected patches' pixels at all.

    The output template's PSF is a `CoaddPsf` built from the matched
    per-patch PSFs, so it represents the actual convolved PSF rather
    than just its Gaussian-equivalent width. The matching kernel
    correlates the per-patch noise; the merged template's variance
    plane reports per-pixel marginal variance honestly, and the
    representative matching-kernel parameters are stamped onto the
    template's metadata so a downstream
    `lsst.ip.diffim.AlardLuptonSubtractTask` can compose them with its
    own A&L kernel and decorrelate the diffim with the right effective
    kernel + pre-convolution variance — see
    ``MATCHING_KERNEL_SIGMA`` / ``MATCHING_KERNEL_NORMSQ`` below.
    """
    ConfigClass = GetPsfMatchedTemplateConfig
    _DefaultName = "getPsfMatchedTemplate"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @timeMethod
    def run(self, *, coaddExposureHandles, bbox, wcs, dataIds, physical_filter):
        """Build a PSF-matched template; parameters mirror
        `GetTemplateTask.run`.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
           A struct with attributes:

           ``template`` : `lsst.afw.image.ExposureF`
               The PSF-matched template coadd

        Raises
        ------
        NoWorkFound
            If no patches overlap the science exposure, or if fewer than
            ``minAcceptedPatches`` survive upper-outlier rejection.
        """
        tractData = self._loadMetadata(coaddExposureHandles, dataIds)
        if not tractData:
            raise pipeBase.NoWorkFound("No input patches found.")
        band, photoCalib = self._validateBandAndCalib(tractData, dataIds)

        bbox.grow(self.config.templateBorderSize)

        targetSigma, numRejected = self._selectAndMatchPatches(tractData)

        warped, catalogs, includedKernelSigmas = self._processTracts(
            tractData, bbox, wcs,
        )
        if not warped:
            raise pipeBase.NoWorkFound(
                "No tracts overlap science exposure after PSF matching."
            )

        template, count, _ = self._merge(warped, bbox, wcs)
        if count == 0:
            raise pipeBase.NoWorkFound("No valid pixels in PSF-matched template.")

        finalCatalog = afwTable.ExposureCatalog(self.schema)
        finalCatalog.reserve(sum(len(c) for c in catalogs))
        for c in catalogs:
            finalCatalog.extend(c)

        self._finalizeTemplate(
            template, finalCatalog, wcs,
            band=band, photoCalib=photoCalib, physical_filter=physical_filter,
            targetSigma=targetSigma, numRejected=numRejected,
            includedKernelSigmas=includedKernelSigmas,
        )
        return pipeBase.Struct(template=template)

    def _loadMetadata(self, coaddExposureHandles, dataIds):
        """Phase 1: pull only the per-patch ``psf``/``wcs``/``photoCalib``/
        ``bbox`` components from the butler (no maskedImage I/O) and
        record the Gaussian-equivalent PSF sigma for each patch.

        Returns
        -------
        tractData : `dict` [`int`, `dict`]
            Per-tract bundle with keys ``patches`` (list of per-patch
            dicts carrying ``dataId``, ``ref``, ``record``, ``sigma``),
            ``wcs`` (the tract's shared WCS, taken from the first
            patch), and ``totalBox`` (union of patch bboxes). Tracts
            with no input patches are omitted.
        """
        tractData = {}
        for tract in coaddExposureHandles:
            catalog, totalBox, refByDataId = self._makeMetadataCatalog(
                coaddExposureHandles[tract], dataIds[tract]
            )
            if len(catalog) == 0:
                continue
            tractData[tract] = {
                "patches": [
                    {
                        "dataId": dataId,
                        "ref": refByDataId[dataId],
                        "record": catalog[i],
                        "sigma": self._patchPsfSigma(catalog[i]),
                    }
                    for i, dataId in enumerate(dataIds[tract])
                ],
                "wcs": catalog[0].wcs,
                "pixelScale": catalog[0].wcs.getPixelScale().asArcseconds(),
                "totalBox": totalBox,
            }
        return tractData

    def _validateBandAndCalib(self, tractData, dataIds):
        """Verify that all input patches share a single band and a single
        photometric calibration.

        photoCalibs are read out of the per-patch records loaded in
        Phase 1 rather than re-fetched from the butler.
        """
        bands = {dataId["band"] for tract in dataIds for dataId in dataIds[tract]}
        if len(bands) > 1:
            raise RuntimeError(
                f"GetPsfMatchedTemplateTask called with multiple bands: {bands}"
            )
        photoCalibs = [
            p["record"].getPhotoCalib()
            for info in tractData.values() for p in info["patches"]
        ]
        if not all(photoCalibs[0] == c for c in photoCalibs):
            raise RuntimeError(
                "GetPsfMatchedTemplateTask called with exposures with different"
                f" photoCalibs: {photoCalibs}"
            )
        return bands.pop(), photoCalibs[0]

    def _selectAndMatchPatches(self, tractData):
        """Phase 2 + 3: one-sided sigma-clip the global PSF-sigma
        distribution, choose a single target sigma, and assign each
        accepted patch its matching-kernel sigma.

        Mutates each per-patch dict in ``tractData`` in place, adding
        ``accepted`` (`bool`) and, for accepted patches, ``kernelSigma``
        (`float` or `None` for passthroughs).

        Returns
        -------
        targetSigma : `float`
            Common Gaussian sigma all accepted patches will be matched to.
        numRejected : `int`
            Number of patches rejected by the upper-outlier clip.
        """
        flatPatches = []
        flatScales = []
        for info in tractData.values():
            flatPatches.extend(info["patches"])
            flatScales.extend([info["pixelScale"]]*len(info["patches"]))
        sigmas = np.array([p["sigma"] for p in flatPatches], dtype=float)
        keep = self._sigmaClip(sigmas)

        if keep.sum() < self.config.minAcceptedPatches:
            raise pipeBase.NoWorkFound(
                f"Only {int(keep.sum())} patch(es) survived PSF outlier rejection; "
                f"require at least {self.config.minAcceptedPatches}."
            )

        for p, k, sigma, scale in zip(flatPatches, keep, sigmas, flatScales):
            p["accepted"] = bool(k)
            if not k:
                self.log.info(
                    "Rejecting patch tract=%s patch=%s from PSF-matched template:"
                    " PSF FWHM=%.3f arcsec.",
                    p["dataId"].get("tract"), p["dataId"].get("patch"),
                    self._psfFwhmArcsec(sigma, scale),
                )

        targetIdx = int(np.argmax(np.where(keep, sigmas, -np.inf)))
        targetSigma = float(sigmas[targetIdx])
        self.log.info(
            "Selected target PSF FWHM %.3f arcsec (passthrough tolerance %.3f).",
            self._psfFwhmArcsec(targetSigma, flatScales[targetIdx]),
            self.config.targetPsfSigmaPad,
        )
        for p, k, scale in zip(flatPatches, keep, flatScales):
            if k:
                p["kernelSigma"] = self._matchingKernelSigma(
                    p["sigma"], targetSigma, pixelScale=scale,
                )
        return targetSigma, int((~keep).sum())

    def _processTracts(self, tractData, bbox, wcs):
        """Phase 4: process each tract independently, freeing its data
        before moving on so peak memory stays at one tract's worth of
        patches. Empty tracts are silently skipped.

        Returns
        -------
        warped : `dict` [`int`, `lsst.afw.image.MaskedImage`]
            Warped, PSF-matched MaskedImage per surviving tract.
        catalogs : `list` [`lsst.afw.table.ExposureCatalog`]
            Per-tract catalogs of records that contributed to ``warped``.
        includedKernelSigmas : `list`
            Per-included-patch matching-kernel sigmas, parallel to the
            concatenation of ``catalogs``. Used by
            ``_finalizeTemplate`` to derive the representative kernel
            metadata.
        """
        warped = {}
        catalogs = []
        includedKernelSigmas = []
        for tract in list(tractData):
            info = tractData.pop(tract)
            result = self._processTract(tract, info, bbox, wcs)
            if result is None:
                continue
            warpedImage, tempCatalog, kernelSigmas = result
            warped[tract] = warpedImage
            catalogs.append(tempCatalog)
            includedKernelSigmas.extend(kernelSigmas)
        return warped, catalogs, includedKernelSigmas

    def _processTract(self, tract, info, bbox, wcs):
        """Run Phase 4 for a single tract: load accepted patches,
        convolve to the common Gaussian width, merge per-tract, and warp
        onto the science WCS.

        Returns
        -------
        result : `tuple` or `None`
            ``(warpedMaskedImage, tempCatalog, kernelSigmas)`` if the
            tract contributes data, or `None` if it should be skipped
            (no accepted patches, no good pixels after merge, or no
            overlap with the science bbox).
        """
        accepted = [p for p in info["patches"] if p["accepted"]]
        if not accepted:
            self.log.info(
                "All patches in tract %s rejected; not including in output.", tract,
            )
            return None

        keptImages = self._loadAndMatchPatches(accepted)

        warpedBox = computeWarpedBBox(info["wcs"], bbox, wcs)
        warpedBox.grow(5)  # to ensure we catch all relevant input pixels

        unwarped, count, included = self._merge(keptImages, warpedBox, info["wcs"])
        del keptImages
        if count == 0:
            self.log.info(
                "No valid pixels from PSF-matched patches in tract %s;"
                " not including in output.", tract,
            )
            return None

        warpedBox.clip(info["totalBox"])
        potentialInput = self.warper.warpExposure(
            wcs, unwarped.subset(warpedBox), destBBox=bbox
        )
        del unwarped
        if np.all(
            potentialInput.mask.array
            & potentialInput.mask.getPlaneBitMask("NO_DATA")
        ):
            self.log.info(
                "No overlap from PSF-matched patches in tract %s;"
                " not including in output.", tract,
            )
            return None

        # `_merge` enumerates `keptImages` in insertion order, which
        # mirrors `accepted`, so `accepted[i]` recovers each included
        # patch's record and matching-kernel sigma.
        tempCatalog = afwTable.ExposureCatalog(self.schema)
        tempCatalog.reserve(len(included))
        kernelSigmas = []
        for i in included:
            tempCatalog.append(accepted[i]["record"])
            kernelSigmas.append(accepted[i]["kernelSigma"])
        return potentialInput.maskedImage, tempCatalog, kernelSigmas

    def _loadAndMatchPatches(self, accepted):
        """Load full Exposures for one tract's accepted patches and
        convolve image + PSF to the common matching width. Patches with
        ``kernelSigma is None`` are passed through untouched.
        """
        keptImages = {}
        for p in accepted:
            coadd = p["ref"].get()
            if p["kernelSigma"] is None:
                keptImages[p["dataId"]] = coadd.maskedImage
            else:
                matchingKernel = self._buildMatchingKernel(p["kernelSigma"])
                keptImages[p["dataId"]] = self._convolveImage(
                    coadd.maskedImage, matchingKernel,
                )
                p["record"].setPsf(
                    self._convolvePsf(p["record"].getPsf(), p["kernelSigma"])
                )
            del coadd
        return keptImages

    def _finalizeTemplate(self, template, catalog, wcs, *, band, photoCalib,
                          physical_filter, targetSigma, numRejected,
                          includedKernelSigmas):
        """Apply finishing touches to the merged template: HIGH_VARIANCE
        masking, NO_DATA variance sentinel, CoaddPsf, filter, photoCalib,
        and PSF-match metadata.
        """
        self.checkHighVariance(template)
        # Replace NaN / non-positive variance pixels with a large finite
        # sentinel so downstream detection doesn't divide by zero or
        # NaN-poison SNRs at NO_DATA regions. After `checkHighVariance`
        # so its background statistics aren't skewed by the sentinel.
        varArr = template.variance.array
        varArr[~(np.isfinite(varArr) & (varArr > 0))] = self._noDataVariance

        # Each accepted record carries the convolved KernelPsf, so the
        # resulting CoaddPsf reflects the matching, not just the
        # Gaussian-equivalent width.
        template.setPsf(self._makePsf(template, catalog, wcs))
        template.setFilter(afwImage.FilterLabel(band, physical_filter))
        template.setPhotoCalib(photoCalib)

        template.metadata["PSF_MATCHED_TARGET_SIGMA"] = float(targetSigma)
        template.metadata["PSF_MATCHED_NUM_PATCHES"] = int(len(catalog))
        template.metadata["PSF_MATCHED_NUM_REJECTED"] = int(numRejected)

        sigmaRepr, normSq = self._representativeKernel(includedKernelSigmas)
        template.metadata["MATCHING_KERNEL_SIGMA"] = sigmaRepr
        template.metadata["MATCHING_KERNEL_NORMSQ"] = normSq

    def _representativeKernel(self, includedKernelSigmas):
        """Compute the representative matching-kernel scalars stamped
        into template metadata so a downstream
        `AlardLuptonSubtractTask` can compose them with its own A&L
        kernel and decorrelate the diffim with the right effective
        kernel and pre-convolution variance.

        Two scalars are returned: the average matching-kernel sigma
        over included patches (px), and Σ K² for that kernel (the
        noise-variance attenuation factor, used to recover the original
        ``tvar`` from the post-convolution variance plane). Patches that
        were passed through without convolution
        (``kernelSigma is None``) contribute σ=0 to the average — they
        didn't correlate the noise — and pull ``normSq`` toward 1.
        """
        sigmasForMean = [
            float(s) if s is not None else 0.0
            for s in includedKernelSigmas
        ]
        sigmaRepr = float(np.mean(sigmasForMean)) if sigmasForMean else 0.0
        if sigmaRepr <= 0:
            return sigmaRepr, 1.0
        reprKernel = self._buildMatchingKernel(sigmaRepr)
        reprImage = afwImage.ImageD(reprKernel.getDimensions())
        reprKernel.computeImage(reprImage, True)
        return sigmaRepr, float((reprImage.array**2).sum())

    def _patchPsfSigma(self, record):
        """Estimate the Gaussian-equivalent PSF sigma (in pixels) for one
        input patch from its catalog record.

        Returns `numpy.nan` if the PSF cannot be evaluated; such patches
        are dropped before sigma clipping.
        """
        psf = record.getPsf()
        if psf is None:
            return float("nan")
        try:
            position = psf.getAveragePosition()
        except Exception:
            position = geom.Box2D(record.getBBox()).getCenter()
        try:
            shape = psf.computeShape(position)
        except Exception as err:
            self.log.debug(
                "Failed to evaluate patch PSF shape (tract=%s patch=%s): %s",
                record.get("tract"), record.get("patch"), err,
            )
            return float("nan")
        return float(shape.getDeterminantRadius())

    def _makeMetadataCatalog(self, exposureRefs, dataIds):
        """Build an `ExposureCatalog` from per-patch *component* loads,
        without reading any patch ``maskedImage`` from the butler.

        Mirrors `GetTemplateTask._makeExposureCatalog` but fetches only
        the small ``psf``, ``wcs``, ``photoCalib``, and ``bbox``
        components so that the expensive pixel-data load can be deferred
        to Phase 4 and skipped entirely for rejected patches.

        Parameters
        ----------
        exposureRefs : iterable of `lsst.daf.butler.DeferredDatasetHandle`
            Handles for the patches to enumerate.
        dataIds : `list` [`lsst.daf.butler.DataCoordinate`]
            Data IDs aligned with ``exposureRefs``.

        Returns
        -------
        catalog : `lsst.afw.table.ExposureCatalog`
            One record per input patch, populated as
            ``_makeExposureCatalog`` would have but without touching the
            pixel data.
        totalBox : `lsst.geom.Box2I`
            Union of the per-patch bounding boxes.
        refByDataId : `dict` [`lsst.daf.butler.DataCoordinate`,
                              `lsst.daf.butler.DeferredDatasetHandle`]
            Mapping from each input ``dataId`` back to its butler handle,
            used by the per-tract pass to load the full ``Exposure``.
        """
        catalog = afwTable.ExposureCatalog(self.schema)
        catalog.reserve(len(exposureRefs))
        totalBox = geom.Box2I()
        refByDataId = {}

        for ref, dataId in zip(exposureRefs, dataIds):
            psf = ref.get(component="psf")
            wcs = ref.get(component="wcs")
            photoCalib = ref.get(component="photoCalib")
            bbox = ref.get(component="bbox")
            totalBox = totalBox.expandedTo(bbox)
            record = catalog.addNew()
            record.setPsf(psf)
            record.setWcs(wcs)
            record.setPhotoCalib(photoCalib)
            record.setBBox(bbox)
            record.setValidPolygon(
                afwGeom.Polygon(geom.Box2D(bbox).getCorners())
            )
            record.set("tract", dataId["tract"])
            record.set("patch", dataId["patch"])
            record.set("weight", 1)
            refByDataId[dataId] = ref

        return catalog, totalBox, refByDataId

    def _sigmaClip(self, sigmas):
        """Return a mask of accepted PSF sigmas using one-sided
        (upper) outlier rejection.

        Every finite, positive sigma at or below the running median is
        always kept — a narrow-PSF outlier is harmless because it just
        gets convolved up to the target. A finite, positive sigma
        *above* the running median is rejected only if its excess
        ``(sigma − median)`` is more than
        ``config.outlierRejectionSigma`` robust standard deviations
        (1.4826·MAD, falling back to stddev when MAD is zero); a single
        wide-PSF patch would otherwise pull the global target high
        enough to broaden every other patch dramatically. The clipping
        iterates up to ``config.outlierRejectionMaxIter`` times.
        """
        keep = np.isfinite(sigmas) & (sigmas > 0)
        if not keep.any():
            return keep

        for _ in range(self.config.outlierRejectionMaxIter):
            if keep.sum() <= 2:
                break
            median = np.median(sigmas[keep])
            mad = np.median(np.abs(sigmas[keep] - median))
            # 1.4826 makes MAD a consistent estimator of the stddev for
            # Gaussian-distributed inputs; fall back to stddev if MAD is 0.
            scale = 1.4826*mad if mad > 0 else float(np.std(sigmas[keep]))
            if scale <= 0:
                break
            threshold = self.config.outlierRejectionSigma*scale
            # Asymmetric: only sigmas more than `threshold` *above* the
            # median are rejected; below-median sigmas always stay.
            new_keep = keep & (sigmas - median <= threshold)
            if np.array_equal(new_keep, keep):
                break
            keep = new_keep

        self.log.info(
            "Accepted %d of %d input patches after upper-outlier rejection.",
            int(keep.sum()), int(np.isfinite(sigmas).sum()),
        )
        return keep

    # FWHM = 2 * sqrt(2 * ln(2)) * sigma for a Gaussian.
    _SIGMA_TO_FWHM = 2.0*np.sqrt(2.0*np.log(2.0))

    @classmethod
    def _psfFwhmArcsec(cls, sigmaPixels, pixelScaleArcsec):
        """Convert a Gaussian-equivalent PSF sigma in pixels to FWHM in
        arcseconds, given the local pixel scale.
        """
        return cls._SIGMA_TO_FWHM*sigmaPixels*pixelScaleArcsec

    def _matchingKernelSigma(self, sigma, targetSigma, pixelScale=None):
        """Return the matching-Gaussian sigma in pixels for a patch with
        intrinsic Gaussian-equivalent ``sigma``, or `None` if no
        convolution should be applied.

        Returns `None` when the patch's sigma is within
        ``config.targetPsfSigmaPad`` (fractionally) of the target — the
        max accepted sigma always satisfies this — or when the implied
        kernel sigma is below ``config.minMatchingKernelSigma``.

        ``pixelScale`` (arcsec/pixel) is used only to format the debug
        log in arcseconds; pass `None` to skip the FWHM annotation.
        """
        passthroughSigma = (1.0 - self.config.targetPsfSigmaPad)*targetSigma
        if sigma >= passthroughSigma:
            if pixelScale is not None:
                self.log.debug(
                    "Skipping PSF match: input FWHM %.3f arcsec within pad"
                    " of target FWHM %.3f arcsec.",
                    self._psfFwhmArcsec(sigma, pixelScale),
                    self._psfFwhmArcsec(targetSigma, pixelScale),
                )
            else:
                self.log.debug(
                    "Skipping PSF match: input sigma within passthrough pad of target."
                )
            return None

        kernelSigma = float(np.sqrt(targetSigma**2 - sigma**2))
        if kernelSigma < self.config.minMatchingKernelSigma:
            self.log.debug(
                "Skipping PSF match: kernel sigma %.3g below floor %.3g.",
                kernelSigma, self.config.minMatchingKernelSigma,
            )
            return None
        return kernelSigma

    def _buildMatchingKernel(self, kernelSigma):
        """Return a normalized 2D Gaussian `AnalyticKernel` of width
        ``kernelSigma`` sized to ``config.matchingKernelSigmas`` half-widths.
        """
        halfWidth = int(np.ceil(self.config.matchingKernelSigmas*kernelSigma))
        ksize = 2*halfWidth + 1
        gaussian = afwMath.GaussianFunction2D(kernelSigma, kernelSigma)
        return afwMath.AnalyticKernel(ksize, ksize, gaussian)

    def _convolveImage(self, maskedImage, kernel):
        """Convolve a patch `MaskedImage` with the supplied matching
        kernel.

        Zero out non-finite image pixels so convolution doesn't smear a
        NaN halo over a kernel footprint. The variance plane still carries
        NaN/zero at those pixels, so `_merge` rejects them anyway.
        """
        sanitized = maskedImage.clone()
        bad = (~np.isfinite(sanitized.image.array)
               | ~np.isfinite(sanitized.variance.array)
               | (sanitized.variance.array <= 0))
        if bad.any():
            sanitized.image.array[bad] = 0.0
            sanitized.variance.array[bad] = self._noDataVariance
            sanitized.mask.array[bad] |= sanitized.mask.getPlaneBitMask("NO_DATA")

        out = afwImage.MaskedImageF(sanitized.getBBox())
        control = afwMath.ConvolutionControl(doNormalize=True)
        afwMath.convolve(out, sanitized, kernel, control)
        return out

    def _convolvePsf(self, originalPsf, kernelSigma):
        """Return a `KernelPsf` representing
        ``originalPsf`` ⊛ ``Gaussian(kernelSigma)``.

        The original PSF is sampled at its average position; if that
        evaluation fails, the original PSF is returned unchanged so the
        record at least keeps a valid PSF for the final ``CoaddPsf``.
        """
        try:
            position = originalPsf.getAveragePosition()
        except Exception:
            position = geom.Point2D(0.0, 0.0)

        try:
            psfImage = originalPsf.computeKernelImage(position)
        except Exception as err:
            self.log.warning(
                "Could not compute PSF kernel image; leaving record PSF"
                " unconvolved: %s", err,
            )
            return originalPsf

        psfArr = np.asarray(psfImage.array, dtype=np.float64)

        halfWidth = int(np.ceil(self.config.matchingKernelSigmas*kernelSigma))
        coords = np.arange(-halfWidth, halfWidth + 1, dtype=np.float64)
        yy, xx = np.meshgrid(coords, coords, indexing="ij")
        gauss = np.exp(-0.5*(xx**2 + yy**2)/(kernelSigma**2))
        gauss /= gauss.sum()

        # 'full' so we don't truncate the wings of either kernel.
        matchedArr = scipy.signal.fftconvolve(psfArr, gauss, mode="full")
        total = matchedArr.sum()
        if total > 0:
            matchedArr /= total

        ny, nx = matchedArr.shape
        matchedImage = afwImage.ImageD(nx, ny)
        matchedImage.array[:, :] = matchedArr
        # Match the `psf.computeKernelImage` convention: xy0 puts (0, 0)
        # at the central pixel of the kernel.
        matchedImage.setXY0(geom.Point2I(-(nx//2), -(ny//2)))
        return KernelPsf(afwMath.FixedKernel(matchedImage))
