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
import collections

import numpy as np
from scipy import ndimage

import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.afw.geom as afwGeom
from lsst.afw.image import VisitInfo
import lsst.afw.table as afwTable
from lsst.afw.math._warper import computeWarpedBBox
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.skymap import BaseSkyMap
from lsst.ip.diffim.dcrModel import calculateDcr, fitThroughput
from lsst.meas.algorithms import CoaddPsf, CoaddPsfConfig, SubtractBackgroundTask
from lsst.utils.timer import timeMethod

__all__ = [
    "GetTemplateTask",
    "GetTemplateConfig",
]


class GetTemplateConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("instrument", "visit", "detector"),
    defaultTemplates={"coaddName": "goodSeeing", "warpTypeSuffix": "", "fakesType": ""},
):
    bbox = pipeBase.connectionTypes.Input(
        doc="Bounding box of exposure to determine the geometry of the output template.",
        name="{fakesType}calexp.bbox",
        storageClass="Box2I",
        dimensions=("instrument", "visit", "detector"),
    )
    wcs = pipeBase.connectionTypes.Input(
        doc="WCS of the exposure that we will construct the template for.",
        name="{fakesType}calexp.wcs",
        storageClass="Wcs",
        dimensions=("instrument", "visit", "detector"),
    )
    skyMap = pipeBase.connectionTypes.Input(
        doc="Geometry of the tracts and patches that the coadds are defined on.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        dimensions=("skymap",),
        storageClass="SkyMap",
    )
    coaddExposures = pipeBase.connectionTypes.Input(
        doc="Coadds that may overlap the desired region, as possible inputs to the template."
        " Will be restricted to those that directly overlap the projected bounding box.",
        dimensions=("tract", "patch", "skymap", "band"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Coadd{warpTypeSuffix}",
        multiple=True,
        deferLoad=True,
        deferGraphConstraint=True,
    )
    visitInfo = pipeBase.connectionTypes.Input(
        doc="VisitInfo of the exposure that we will construct the template for."
        " Used to determine the observing conditions for DCR.",
        name="{fakesType}calexp.visitInfo",
        storageClass="VisitInfo",
        dimensions=("instrument", "visit", "detector"),
    )
    throughput = pipeBase.connectionTypes.Input(
        doc="Bandpass of the filter used for the observation.",
        name="standard_passband",
        storageClass="ArrowAstropy",
        dimensions=("band", "instrument"),
    )
    dcrCorrectionCatalogs = pipeBase.connectionTypes.Input(
        doc="Catalog of sub-band fluxes and footprints for moderately bright sources.",
        storageClass="SourceCatalog",
        dimensions=("tract", "patch", "skymap", "band"),
        name="dcr_correction_catalog",
        multiple=True,
        deferLoad=True,
        deferGraphConstraint=True,
    )

    template = pipeBase.connectionTypes.Output(
        doc="Warped template, pixel matched to the bounding box and WCS.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Diff_templateExp{warpTypeSuffix}",
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not config.useDcrCorrection:
            self.inputs.remove("visitInfo")
            self.inputs.remove("dcrCorrectionCatalogs")
            self.inputs.remove("throughput")


class GetTemplateConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=GetTemplateConnections
):
    templateBorderSize = pexConfig.Field(
        dtype=int,
        default=20,
        doc="Number of pixels to grow the requested template image to account for warping",
    )
    warp = pexConfig.ConfigField(
        dtype=afwMath.Warper.ConfigClass,
        doc="warper configuration",
    )
    coaddPsf = pexConfig.ConfigField(
        doc="Configuration for CoaddPsf",
        dtype=CoaddPsfConfig,
    )
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
    useDcrCorrection = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Use the DCR catalog to correct the shape of included sources",
    )
    dcrModelScale = pexConfig.Field(
        dtype=float,
        default=1,
        doc="Fudge scaling factor of the model",
    )
    dcrScale = pexConfig.Field(
        dtype=float,
        default=1,
        doc="Fudge scaling factor of the shift",
    )
    dcrWavelengthShift = pexConfig.Field(
        dtype=float,
        default=0,
        doc="Fudge shift of the effective wavelength",
    )

    def setDefaults(self):
        # Use a smaller cache: per SeparableKernel.computeCache, this should
        # give a warping error of a fraction of a count (these must match).
        self.warp.cacheSize = 100000
        self.coaddPsf.cacheSize = self.warp.cacheSize
        # The WCS for LSST should be smoothly varying, so we can use a longer
        # interpolation length for WCS evaluations.
        self.warp.interpLength = 100
        self.warp.warpingKernelName = "lanczos3"
        self.coaddPsf.warpingKernelName = self.warp.warpingKernelName

        # Background subtraction of the variance plane
        self.varianceBackground.algorithm = "LINEAR"
        self.varianceBackground.binSize = 32
        self.varianceBackground.useApprox = False
        self.varianceBackground.statisticsProperty = "MEDIAN"
        self.varianceBackground.doFilterSuperPixels = True
        self.varianceBackground.ignoredPixelMask = ["BAD",
                                                    "EDGE",
                                                    "DETECTED",
                                                    "DETECTED_NEGATIVE",
                                                    "NO_DATA",
                                                    ]


class GetTemplateTask(pipeBase.PipelineTask):
    ConfigClass = GetTemplateConfig
    _DefaultName = "getTemplate"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.warper = afwMath.Warper.fromConfig(self.config.warp)
        self.schema = afwTable.ExposureTable.makeMinimalSchema()
        self.schema.addField(
            "tract", type=np.int32, doc="Which tract this exposure came from."
        )
        self.schema.addField(
            "patch",
            type=np.int32,
            doc="Which patch in the tract this exposure came from.",
        )
        self.schema.addField(
            "weight",
            type=float,
            doc="Weight for each exposure, used to make the CoaddPsf; should always be 1.",
        )
        self.makeSubtask("varianceBackground")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        bbox = inputs.pop("bbox")
        wcs = inputs.pop("wcs")
        coaddExposures = inputs.pop("coaddExposures")
        skymap = inputs.pop("skyMap")
        if self.config.useDcrCorrection:
            visitInfo = inputs.pop("visitInfo")
            dcrCorrectionCatalogs = inputs.pop("dcrCorrectionCatalogs")
            throughput = fitThroughput(inputs.pop("throughput"))
            self.effectiveWavelength = throughput.effectiveWavelength
            self.bandwidth = throughput.bandwidth
        else:
            visitInfo = None
            dcrCorrectionCatalogs = None
            self.effectiveWavelength = None
            self.bandwidth = None

        # This should not happen with a properly configured execution context.
        assert not inputs, "runQuantum got more inputs than expected"

        results = self.getExposures(coaddExposures, bbox, skymap, wcs, dcrCorrectionCatalogs)
        physical_filter = butlerQC.quantum.dataId["physical_filter"]
        outputs = self.run(
            coaddExposureHandles=results.coaddExposures,
            bbox=bbox,
            wcs=wcs,
            dataIds=results.dataIds,
            physical_filter=physical_filter,
            visit=outputRefs.template.dataId["visit"],
            visitInfo=visitInfo,
            dcrCorrectionHandles=results.dcrCorrectionCatalogs,
        )
        butlerQC.put(outputs, outputRefs)

    def getExposures(self, coaddExposureHandles, bbox, skymap, wcs, dcrCorrectionCatalogHandles):

        """Return a data structure containing the coadds that overlap the
        specified bbox projected onto the sky, and a corresponding data
        structure of their dataIds.
        These are the appropriate inputs to this task's `run` method.

        The spatial index in the butler registry has generous padding and often
        supplies patches near, but not directly overlapping the desired region.
        This method filters the inputs so that `run` does not have to read in
        all possibly-matching coadd exposures.

        Parameters
        ----------
        coaddExposureHandles : `iterable` \
                          [`lsst.daf.butler.DeferredDatasetHandle` of \
                           `lsst.afw.image.Exposure`]
            Dataset handles to exposures that might overlap the desired
            region.
        bbox : `lsst.geom.Box2I`
            Template bounding box of the pixel geometry onto which the
            coaddExposures will be resampled.
        skymap : `lsst.skymap.SkyMap`
            Geometry of the tracts and patches the coadds are defined on.
        wcs : `lsst.afw.geom.SkyWcs`
            Template WCS onto which the coadds will be resampled.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
           A struct with attributes:

           ``coaddExposures``
               Dict of coadd exposures that overlap the projected bbox,
               indexed on tract id
               (`dict` [`int`, `list` [`lsst.daf.butler.DeferredDatasetHandle` of
                                       `lsst.afw.image.Exposure`] ]).
           ``dataIds``
               Dict of data IDs of the coadd exposures that overlap the
               projected bbox, indexed on tract id
               (`dict` [`int`, `list [`lsst.daf.butler.DataCoordinate`] ]).

        Raises
        ------
        NoWorkFound
            Raised if no patches overlap the input detector bbox, or the input
            WCS is None.
        """
        if wcs is None:
            raise pipeBase.NoWorkFound(
                "WCS is None; cannot find overlapping exposures."
            )

        # Exposure's validPolygon would be more accurate
        detectorPolygon = geom.Box2D(bbox)
        detectorCorners = wcs.pixelToSky(detectorPolygon.getCorners())
        overlappingArea = 0
        coaddExposures = collections.defaultdict(list)
        dataIds = collections.defaultdict(list)
        if self.config.useDcrCorrection:
            # Create a temporary lookup table of the DCR catalogs, so that the list
            # per tract will have the same order as the coadds and dataIds
            dcrCorrectionCatalogs = collections.defaultdict(list)
            dcrCorrectionTempRefs = collections.defaultdict(dict)
            for dcrRef in dcrCorrectionCatalogHandles:
                dataId = dcrRef.dataId
                dcrCorrectionTempRefs[dataId["tract"]][dataId["patch"]] = dcrRef
        else:
            dcrCorrectionCatalogs = None
        for coaddRef in coaddExposureHandles:
            dataId = coaddRef.dataId
            patchWcs = skymap[dataId["tract"]].getWcs()
            patchBBox = skymap[dataId["tract"]][dataId["patch"]].getOuterBBox()
            patchPolygon = afwGeom.Polygon(geom.Box2D(patchBBox))
            # Calculate detector/patch overlap in patch coordinates rather than
            # detector coordinates because the skymap's inverse mapping
            # (patchWcs.skyToPixel()) is more stable than the detector's for
            # arbitrary sky coordinates.
            detectorInPatchCoordinates = afwGeom.Polygon(patchWcs.skyToPixel(detectorCorners))
            if patchPolygon.intersection(detectorInPatchCoordinates):
                overlappingArea += patchPolygon.intersectionSingle(
                    detectorInPatchCoordinates
                ).calculateArea()
                self.log.info(
                    "Using template input tract=%s, patch=%s",
                    dataId["tract"],
                    dataId["patch"],
                )
                coaddExposures[dataId["tract"]].append(coaddRef)
                dataIds[dataId["tract"]].append(dataId)

                if self.config.useDcrCorrection:
                    # A patch with no DCR catalog is still usable, it just
                    # cannot be corrected. Append `None` to keep the list
                    # aligned with the coadds and dataIds.
                    dcrRef = dcrCorrectionTempRefs[dataId["tract"]].get(dataId["patch"])
                    if dcrRef is None:
                        self.log.warning(
                            "No DCR correction catalog for tract=%s, patch=%s;"
                            " including this patch without a DCR correction.",
                            dataId["tract"],
                            dataId["patch"],
                        )
                    dcrCorrectionCatalogs[dataId["tract"]].append(dcrRef)

        if not overlappingArea:
            raise pipeBase.NoWorkFound("No patches overlap detector")

        return pipeBase.Struct(coaddExposures=coaddExposures, dataIds=dataIds,
                               dcrCorrectionCatalogs=dcrCorrectionCatalogs)

    @timeMethod
    def run(self, *, coaddExposureHandles, bbox, wcs, dataIds, physical_filter, visit=None,
            visitInfo=None, dcrCorrectionHandles=None):
        """Warp coadds from multiple tracts and patches to form a template to
        subtract from a science image.

        Tract and patch overlap regions are combined by a variance-weighted
        average, and the variance planes are combined with the same weights,
        not added in quadrature; the overlap regions are not statistically
        independent, because they're derived from the same original data.
        The PSF on the template is created by combining the CoaddPsf on each
        template image into a meta-CoaddPsf.

        Parameters
        ----------
        coaddExposureHandles : `dict` [`int`,  `list` of \
                          [`lsst.daf.butler.DeferredDatasetHandle` of \
                           `lsst.afw.image.Exposure`]]
            Coadds to be mosaicked, indexed on tract id.
        bbox : `lsst.geom.Box2I`
            Template Bounding box of the detector geometry onto which to
            resample the ``coaddExposureHandles``. Modified in-place to include the
            template border.
        wcs : `lsst.afw.geom.SkyWcs`
            Template WCS onto which to resample the ``coaddExposureHandles``.
        dataIds : `dict` [`int`, `list` [`lsst.daf.butler.DataCoordinate`]]
            Record of the tract and patch of each coaddExposure, indexed on
            tract id.
        physical_filter : `str`
            Physical filter of the science image.
        visit : `int`, optional
            If supplied, over-write the visit ID in the template's visitInfo
            so that downstream source injection tasks can link the template and
            science image for the visit.
        visitInfo : `lsst.afw.image.VisitInfo`, optional
            VisitInfo of exposure used for applying the DCR correction catalog

        Returns
        -------
        result : `lsst.pipe.base.Struct`
           A struct with attributes:

           ``template``
               A template coadd exposure assembled out of patches
               (`lsst.afw.image.ExposureF`).

        Raises
        ------
        NoWorkFound
            If no coadds are found with sufficient un-masked pixels.
        """
        band, photoCalib = self._checkInputs(dataIds, coaddExposureHandles)

        bbox.grow(self.config.templateBorderSize)

        warped = {}
        catalogs = []
        for tract in coaddExposureHandles:
            dcrCorrectionHandlesTract = dcrCorrectionHandles[tract] if self.config.useDcrCorrection else None
            maskedImages, catalog, totalBox = self._makeExposureCatalog(
                coaddExposureHandles[tract], dataIds[tract], visitInfo, dcrCorrectionHandlesTract,
            )
            warpedBox = computeWarpedBBox(catalog[0].wcs, bbox, wcs)
            warpedBox.grow(5)  # to ensure we catch all relevant input pixels
            # Combine images from individual patches together.
            unwarped, count, included = self._merge(
                maskedImages, warpedBox, catalog[0].wcs
            )
            # Delete `maskedImages` after combining into one large image to reduce peak memory use
            del maskedImages
            if count == 0:
                self.log.info(
                    "No valid pixels from coadd patches in tract %s; not including in output.",
                    tract,
                )
                continue
            warpedBox.clip(totalBox)
            potentialInput = self.warper.warpExposure(
                wcs, unwarped.subset(warpedBox), destBBox=bbox
            )

            # Delete the single large `unwarped` image after warping to reduce peak memory use
            del unwarped
            if np.all(
                potentialInput.mask.array
                & potentialInput.mask.getPlaneBitMask("NO_DATA")
            ):
                self.log.info(
                    "No overlap from coadd patches in tract %s; not including in output.",
                    tract,
                )
                continue

            # Trim the exposure catalog to just the patches that were used.
            tempCatalog = afwTable.ExposureCatalog(self.schema)
            tempCatalog.reserve(len(included))
            for i in included:
                tempCatalog.append(catalog[i])
            catalogs.append(tempCatalog)
            warped[tract] = potentialInput.maskedImage

        if len(warped) == 0:
            raise pipeBase.NoWorkFound("No patches found to overlap science exposure.")
        # At this point, all entries will be valid, so we can ignore included.
        template, count, _ = self._merge(warped, bbox, wcs)
        if count == 0:
            raise pipeBase.NoWorkFound("No valid pixels in warped template.")

        # Make a single catalog containing all the inputs that were accepted.
        catalog = afwTable.ExposureCatalog(self.schema)
        catalog.reserve(sum([len(c) for c in catalogs]))
        for c in catalogs:
            catalog.extend(c)

        # Set a mask plane for any regions with exceptionally high variance.
        self.checkHighVariance(template)
        if visit is not None:
            template.getInfo().setVisitInfo(VisitInfo(id=visit))
        template.setFilter(afwImage.FilterLabel(band, physical_filter))
        template.setPhotoCalib(photoCalib)
        template.setPsf(self._makePsf(template, catalog, wcs))
        # Record the input coadd patches as the template's coadd inputs.
        coaddInputs = afwImage.CoaddInputs(afwTable.ExposureTable.makeMinimalSchema(), self.schema)
        coaddInputs.ccds.extend(catalog, deep=True)
        template.getInfo().setCoaddInputs(coaddInputs)
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

    @staticmethod
    def _checkInputs(dataIds, coaddExposures):
        """Check that the all the dataIds are from the same band and that
        the exposures all have the same photometric calibration.

        Parameters
        ----------
        dataIds : `dict` [`int`, `list` [`lsst.daf.butler.DataCoordinate`]]
            Record of the tract and patch of each coaddExposure.
        coaddExposures : `dict` [`int`,  `list` of \
                          [`lsst.daf.butler.DeferredDatasetHandle` of \
                           `lsst.afw.image.Exposure` or
                           `lsst.afw.image.Exposure`]]
            Coadds to be mosaicked.

        Returns
        -------
        band : `str`
            Filter band of all the input exposures.
        photoCalib : `lsst.afw.image.PhotoCalib`
            Photometric calibration of all of the input exposures.

        Raises
        ------
        RuntimeError
            Raised if the bands or calibrations of the input exposures are not
            all the same.
        """
        bands = set(dataId["band"] for tract in dataIds for dataId in dataIds[tract])
        if len(bands) > 1:
            raise RuntimeError(f"GetTemplateTask called with multiple bands: {bands}")
        band = bands.pop()
        photoCalibs = [
            exposure.get(component="photoCalib")
            for exposures in coaddExposures.values()
            for exposure in exposures
        ]
        if not all([photoCalibs[0] == x for x in photoCalibs]):
            msg = f"GetTemplateTask called with exposures with different photoCalibs: {photoCalibs}"
            raise RuntimeError(msg)
        photoCalib = photoCalibs[0]
        return band, photoCalib

    def _makeExposureCatalog(self, exposureRefs, dataIds, visitInfo=None, dcrCorrectionRefs=None):
        """Make an exposure catalog for one tract.

        Parameters
        ----------
        exposureRefs : `list` of [`lsst.daf.butler.DeferredDatasetHandle` of \
                        `lsst.afw.image.Exposure`]
            Exposures to include in the catalog.
        dataIds : `list` [`lsst.daf.butler.DataCoordinate`]
            Data ids of each of the included exposures; must have "tract" and
            "patch" entries.
        dcrCorrectionRefs : `list` of [`lsst.daf.butler.DeferredDatasetHandle` of \
                        `lsst.afw.table.SourceCatalog`], optional
            Catalogs with heavy footprints of DCR models for sources in each
            patch, in the same order as ``exposureRefs``. Entries may be `None`
            for patches that have no DCR correction catalog, which are then
            used without a DCR correction.

        Returns
        -------
        images : `dict` [`lsst.afw.image.MaskedImage`]
            MaskedImages of each of the input exposures, for warping.
        catalog : `lsst.afw.table.ExposureCatalog`
            Catalog of metadata for each exposure
        totalBox : `lsst.geom.Box2I`
            The union of the bounding boxes of all the input exposures.
        """
        catalog = afwTable.ExposureCatalog(self.schema)
        catalog.reserve(len(exposureRefs))
        exposures = (exposureRef.get() for exposureRef in exposureRefs)
        if dcrCorrectionRefs is None:
            dcrCatalogs = [None, ]*len(exposureRefs)
        else:
            # `getExposures` builds this list alongside the coadds, so a length
            # mismatch means they are no longer aligned and the corrections
            # would be applied to the wrong patches.
            if len(dcrCorrectionRefs) != len(exposureRefs):
                raise RuntimeError(f"Got {len(dcrCorrectionRefs)} DCR correction catalogs for "
                                   f"{len(exposureRefs)} coadds; these must correspond one to one.")
            dcrCatalogs = (None if ref is None else ref.get() for ref in dcrCorrectionRefs)
        images = {}
        totalBox = geom.Box2I()

        for coadd, dataId, dcrCatalog in zip(exposures, dataIds, dcrCatalogs):
            if dcrCatalog is None:
                images[dataId] = coadd.maskedImage
            else:
                images[dataId] = self.applyDcr(coadd, visitInfo, dcrCatalog).maskedImage
            bbox = coadd.getBBox()
            totalBox = totalBox.expandedTo(bbox)
            record = catalog.addNew()
            record.setPsf(coadd.psf)
            record.setWcs(coadd.wcs)
            record.setPhotoCalib(coadd.photoCalib)
            record.setBBox(bbox)
            record.setValidPolygon(afwGeom.Polygon(geom.Box2D(bbox).getCorners()))
            record.set("tract", dataId["tract"])
            record.set("patch", dataId["patch"])
            # Weight is used by CoaddPsf, but the PSFs from overlapping patches
            # should be very similar, so this value mostly shouldn't matter.
            record.set("weight", 1)

        return images, catalog, totalBox

    def _merge(self, maskedImages, bbox, wcs):
        """Merge the images that came from one tract into one larger image,
        ignoring NaN pixels and non-finite variance pixels from individual
        exposures.

        Parameters
        ----------
        maskedImages : `dict` [`lsst.afw.image.MaskedImage` or
                               `lsst.afw.image.Exposure`]
            Images to be merged into one larger bounding box.
        bbox : `lsst.geom.Box2I`
            Bounding box defining the image to merge into.
        wcs : `lsst.afw.geom.SkyWcs`
            WCS of all of the input images to set on the output image.

        Returns
        -------
        merged : `lsst.afw.image.MaskedImage`
            Merged image with all of the inputs at their respective bbox
            positions.
        count : `int`
            Count of the number of good pixels (those with positive weights)
            in the merged image.
        included : `list` [`int`]
            List of indexes of patches that were included in the merged
            result, to be used to trim the exposure catalog.
        """
        merged = afwImage.ExposureF(bbox, wcs)
        weights = afwImage.ImageF(bbox)
        included = []  # which patches were included in the result
        for i, (dataId, maskedImage) in enumerate(maskedImages.items()):
            # Only merge into the trimmed box, to save memory
            clippedBox = geom.Box2I(maskedImage.getBBox())
            clippedBox.clip(bbox)
            if clippedBox.area == 0:
                self.log.debug("%s does not overlap template region.", dataId)
                continue  # nothing in this image overlaps the output
            maskedImage = maskedImage.subset(clippedBox)
            # Catch both zero-value and NaN variance plane pixels
            good = (maskedImage.variance.array > 0) & (
                np.isfinite(maskedImage.variance.array)
            )
            weight = maskedImage.variance.array[good] ** (-0.5)
            bad = np.isnan(maskedImage.image.array) | ~good
            # Note that modifying the patch MaskedImage in place is fine;
            # we're throwing it away at the end anyway.
            maskedImage.image.array[bad] = 0.0
            maskedImage.variance.array[bad] = 0.0
            # Reset mask, too, since these pixels don't contribute to sum.
            maskedImage.mask.array[bad] = 0
            # Cannot use `merged.maskedImage *= weight` because that operator
            # multiplies the variance by the weight twice; in this case
            # `weight` are the exact values we want to scale by.
            maskedImage.image.array[good] *= weight
            maskedImage.variance.array[good] *= weight
            weights[clippedBox].array[good] += weight
            # Free memory before creating new large arrays
            del weight
            merged.maskedImage[clippedBox] += maskedImage
            included.append(i)

        good = weights.array > 0

        # Cannot use `merged.maskedImage /= weights` because that
        # operator divides the variance by the weight twice; in this case
        # `weights` are the exact values we want to scale by.
        weights = weights.array[good]
        merged.image.array[good] /= weights
        merged.variance.array[good] /= weights

        merged.mask.array[~good] |= merged.mask.getPlaneBitMask("NO_DATA")

        return merged, good.sum(), included

    def _makePsf(self, template, catalog, wcs):
        """Return a PSF containing the PSF at each of the input regions.

        Note that although this includes all the exposures from the catalog,
        the PSF knows which part of the template the inputs came from, so when
        evaluated at a given position it will not include inputs that never
        went in to those pixels.

        Parameters
        ----------
        template : `lsst.afw.image.Exposure`
            Generated template the PSF is for.
        catalog : `lsst.afw.table.ExposureCatalog`
            Catalog of exposures that went into the template that contains all
            of the input PSFs.
        wcs : `lsst.afw.geom.SkyWcs`
            WCS of the template, to warp the PSFs to.

        Returns
        -------
        coaddPsf : `lsst.meas.algorithms.CoaddPsf`
            The meta-psf constructed from all of the input catalogs.
        """
        # CoaddPsf centroid not only must overlap image, but must overlap the
        # part of image with data. Use centroid of region with data.
        boolmask = template.mask.array & template.mask.getPlaneBitMask("NO_DATA") == 0
        maskx = afwImage.makeMaskFromArray(boolmask.astype(afwImage.MaskPixel))
        centerCoord = afwGeom.SpanSet.fromMask(maskx, 1).computeCentroid()

        ctrl = self.config.coaddPsf.makeControl()
        coaddPsf = CoaddPsf(
            catalog, wcs, centerCoord, ctrl.warpingKernelName, ctrl.cacheSize
        )
        return coaddPsf

    def applyDcr(self, coadd, visitInfo, dcrCatalog):
        """Summary

        Parameters
        ----------
        coadd : TYPE
            Description
        dcrCatalog : TYPE
            Description
        """
        nSubfilters = None
        coaddBBox = coadd.getBBox()
        for dcrRecord in dcrCatalog:
            if nSubfilters is None:
                nSubfilters = int(dcrRecord['numSubfilters'])
            wl_use = self.effectiveWavelength - self.config.dcrWavelengthShift
            dcrShift = calculateDcr(visitInfo, coadd.wcs, wl_use, self.bandwidth,
                                    nSubfilters, bbox=coaddBBox)
            footprint = dcrRecord.getFootprint()
            fpBBox = footprint.getBBox()
            # The footprints were not necessarily measured on this coadd, so
            # they may extend past its edge. Only correct the pixels we have;
            # the rest of the source is handled by the neighboring patch.
            bbox = fpBBox.clippedTo(coaddBBox)
            if bbox.isEmpty():
                self.log.debug("DCR correction footprint %s does not overlap the coadd; skipping.",
                               dcrRecord.getId())
                continue
            # flux = dcrRecord['modelFlux']
            # ``fill`` must be set: it defaults to NaN, which would poison
            # every pixel of the bbox that is outside the footprint.
            modelImage = footprint.extractImage(fill=0.)
            modelImage.array *= self.config.dcrModelScale
            # Shift the full model and clip afterwards, so that flux moving in
            # from beyond the coadd edge is not lost.
            shiftedImage = afwImage.ImageF(fpBBox)

            coadd[bbox].image.array -= modelImage[bbox].array
            for subfilter, shift in enumerate(dcrShift):
                shift2 = tuple(sh*self.config.dcrScale for sh in shift)
                subFlux = dcrRecord[f'subfilterWeight_{subfilter}']
                shiftedImage.array[:] = ndimage.shift(modelImage.array, shift2)
                coadd[bbox].image.array += subFlux*shiftedImage[bbox].array
        return coadd
