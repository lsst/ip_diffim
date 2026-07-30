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
        doc="Fraction of the DCR correction to apply to each source. Zero leaves"
        " the coadd unchanged, and one applies the full correction.",
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
    dcrWeightTolerance = pexConfig.Field(
        dtype=float,
        default=1e-3,
        doc="Maximum amount by which the sum of a source's subfilter weights may"
        " differ from one. The DCR correction conserves the flux of a source"
        " only if its weights sum to one, so sources that fail this check are"
        " left uncorrected.",
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
            effectiveWavelength = throughput.effectiveWavelength
            bandwidth = throughput.bandwidth
        else:
            visitInfo = None
            dcrCorrectionCatalogs = None
            effectiveWavelength = None
            bandwidth = None

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
            effectiveWavelength=effectiveWavelength,
            bandwidth=bandwidth,
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
            visitInfo=None, dcrCorrectionHandles=None, effectiveWavelength=None, bandwidth=None):
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
            VisitInfo of the science exposure, used to calculate the DCR of each
            subfilter. Required if ``dcrCorrectionHandles`` is supplied.
        dcrCorrectionHandles : `dict` [`int`, `list` \
                [`lsst.daf.butler.DeferredDatasetHandle`]], optional
            DCR correction catalogs of each patch, indexed on tract id and in the
            same order as ``coaddExposureHandles``.
        effectiveWavelength : `float`, optional
            Effective wavelength of the filter, in nm. Required if
            ``dcrCorrectionHandles`` is supplied.
        bandwidth : `float`, optional
            Bandwidth of the filter, in nm. Required if
            ``dcrCorrectionHandles`` is supplied.

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
        ValueError
            If ``dcrCorrectionHandles`` is supplied without the filter
            parameters needed to calculate DCR.
        """
        if dcrCorrectionHandles is not None and (visitInfo is None or effectiveWavelength is None
                                                 or bandwidth is None):
            raise ValueError("`visitInfo`, `effectiveWavelength` and `bandwidth` must all be supplied "
                             "in order to apply a DCR correction.")
        band, photoCalib = self._checkInputs(dataIds, coaddExposureHandles)

        bbox.grow(self.config.templateBorderSize)

        warped = {}
        catalogs = []
        for tract in coaddExposureHandles:
            dcrCorrectionHandlesTract = (dcrCorrectionHandles[tract]
                                         if dcrCorrectionHandles is not None else None)
            maskedImages, catalog, totalBox = self._makeExposureCatalog(
                coaddExposureHandles[tract], dataIds[tract], visitInfo, dcrCorrectionHandlesTract,
                effectiveWavelength, bandwidth,
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

    def _makeExposureCatalog(self, exposureRefs, dataIds, visitInfo=None, dcrCorrectionRefs=None,
                             effectiveWavelength=None, bandwidth=None):
        """Make an exposure catalog for one tract.

        Parameters
        ----------
        exposureRefs : `list` of [`lsst.daf.butler.DeferredDatasetHandle` of \
                        `lsst.afw.image.Exposure`]
            Exposures to include in the catalog.
        dataIds : `list` [`lsst.daf.butler.DataCoordinate`]
            Data ids of each of the included exposures; must have "tract" and
            "patch" entries.
        visitInfo : `lsst.afw.image.VisitInfo`, optional
            VisitInfo of the science exposure, used to calculate the DCR of each
            subfilter.
        dcrCorrectionRefs : `list` of [`lsst.daf.butler.DeferredDatasetHandle` of \
                        `lsst.afw.table.SourceCatalog`], optional
            Catalogs with heavy footprints of DCR models for sources in each
            patch, in the same order as ``exposureRefs``. Entries may be `None`
            for patches that have no DCR correction catalog, which are then
            used without a DCR correction.
        effectiveWavelength : `float`, optional
            Effective wavelength of the filter, in nm.
        bandwidth : `float`, optional
            Bandwidth of the filter, in nm.

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
                images[dataId] = self.applyDcr(coadd, visitInfo, dcrCatalog,
                                               effectiveWavelength, bandwidth).maskedImage
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

    def applyDcr(self, coadd, visitInfo, dcrCatalog, effectiveWavelength, bandwidth):
        """Replace the DCR of the coadd with the DCR of the science visit, for
        each source in a DCR correction catalog.

        The sources in a coadd are smeared along the parallactic direction by
        the average DCR of all of the visits that went into it, while in the
        science image they are shifted by the DCR of that observation alone. For
        each modeled source this subtracts the model of the source as it appears
        in the coadd, and adds back the un-shifted model of the source shifted to
        the position of each subfilter in this observation.

        Parameters
        ----------
        coadd : `lsst.afw.image.Exposure`
            One patch of the coadd, modified in place. Only the image plane is
            changed; the mask and variance planes are left alone.
        visitInfo : `lsst.afw.image.VisitInfo`
            Metadata of the science exposure the template is being built for,
            used to calculate the DCR of each subfilter.
        dcrCatalog : `lsst.afw.table.SourceCatalog`
            Catalog of sub-band fluxes and model footprints for this patch. Each
            source is represented by a pair of records sharing the same id: one
            with ``isCoaddModel`` set, whose footprint is the DCR-smeared model
            of the source in the coadd, and one without, whose footprint is the
            un-shifted model.
        effectiveWavelength : `float`
            Effective wavelength of the filter, in nm.
        bandwidth : `float`
            Bandwidth of the filter, in nm.

        Returns
        -------
        coadd : `lsst.afw.image.Exposure`
            The same coadd that was supplied, with the corrections applied.
        """
        if len(dcrCatalog) == 0:
            return coadd
        coaddBBox = coadd.getBBox()
        nSubfilters = self._getNumSubfilters(dcrCatalog)
        # The shift depends only on the observation and the coadd, so calculate
        # it once for the whole catalog.
        shiftedWavelength = effectiveWavelength - self.config.dcrWavelengthShift
        dcrShift = [tuple(sh*self.config.dcrScale for sh in shift)
                    for shift in calculateDcr(visitInfo, coadd.wcs, shiftedWavelength,
                                              bandwidth, nSubfilters, bbox=coaddBBox)]
        if not np.all(np.isfinite(dcrShift)):
            self.log.warning("DCR shifts are not all finite (%s); no DCR correction can be applied to"
                             " this patch. Check the filter throughput and the exposure metadata.",
                             dcrShift)
            return coadd

        # Each source has a second record holding the model of the source as it
        # appears in the coadd, which is what has to be subtracted.
        coaddModels = {record.getId(): record.getFootprint()
                       for record in dcrCatalog if record['isCoaddModel']}

        nBadWeights = 0
        worstWeightSum = 1.
        worstDeviation = 0.
        nUnpaired = 0
        nNotFinite = 0
        for dcrRecord in dcrCatalog:
            if dcrRecord['isCoaddModel']:
                continue
            subfilterWeights = [dcrRecord[f'subfilterWeight_{subfilter}']
                                for subfilter in range(nSubfilters)]
            # The correction only redistributes a source's flux between the
            # subfilters if the weights sum to one. If they do not, applying it
            # would change the flux of the source in the template, so leave the
            # source uncorrected instead.
            # Note the inverted comparison: a NaN weight must fail this check,
            # and every direct comparison with NaN is False.
            deviation = abs(sum(subfilterWeights) - 1)
            if not (deviation <= self.config.dcrWeightTolerance):
                nBadWeights += 1
                # Rank NaN weights above any finite deviation, so that they are
                # the ones reported as the worst offender.
                rank = np.inf if np.isnan(deviation) else deviation
                if rank > worstDeviation:
                    worstDeviation = rank
                    worstWeightSum = sum(subfilterWeights)
                continue
            footprint = dcrRecord.getFootprint()
            fpBBox = footprint.getBBox()
            coaddFootprint = coaddModels.get(dcrRecord.getId())
            if coaddFootprint is None or coaddFootprint.getBBox() != fpBBox:
                nUnpaired += 1
                continue
            # The footprints were not necessarily measured on this coadd, so
            # they may extend past its edge. Only correct the pixels we have;
            # the rest of the source is handled by the neighboring patch.
            bbox = fpBBox.clippedTo(coaddBBox)
            if bbox.isEmpty():
                self.log.debug("DCR correction footprint %s does not overlap the coadd; skipping.",
                               dcrRecord.getId())
                continue
            # ``fill`` must be set: it defaults to NaN, which would poison
            # every pixel of the bbox that is outside the footprint.
            # Build the whole correction on the footprint's own bbox and clip
            # once at the end, so that flux moving in from beyond the coadd edge
            # is not lost.
            correction = afwImage.ImageF(fpBBox)
            correction.array -= coaddFootprint.extractImage(fill=0.).array
            unshiftedModel = footprint.extractImage(fill=0.).array
            for shift, subFlux in zip(dcrShift, subfilterWeights):
                correction.array += subFlux*ndimage.shift(unshiftedModel, shift)
            # Last check before the correction reaches the template: a single
            # non-finite pixel in either stored model would spread over the
            # whole footprint when it is shifted, and would then be averaged
            # into the template with full weight by `_merge`.
            if not np.all(np.isfinite(correction.array)):
                nNotFinite += 1
                continue
            coadd[bbox].image.array += self.config.dcrModelScale*correction[bbox].array

        if nNotFinite:
            self.log.warning(
                "Left %d sources uncorrected because their DCR correction was not finite; the model"
                " footprints in the DCR correction catalog contain NaN or infinite pixels.",
                nNotFinite,
            )
        if nUnpaired:
            self.log.warning(
                "Left %d sources uncorrected because they have no matching coadd model record;"
                " the DCR correction catalog should contain one record with `isCoaddModel` set"
                " and one without, on the same bbox, for every source.",
                nUnpaired,
            )
        if nBadWeights:
            self.log.warning(
                "Left %d of %d sources uncorrected because their subfilter weights do not sum to"
                " one (worst sum %g); the DCR correction would have changed their flux.",
                nBadWeights,
                len(dcrCatalog) - len(coaddModels),
                worstWeightSum,
            )
        return coadd

    @staticmethod
    def _getNumSubfilters(dcrCatalog):
        """Return the number of subfilters that a DCR correction catalog was
        modeled with.

        Parameters
        ----------
        dcrCatalog : `lsst.afw.table.SourceCatalog`
            Catalog of sub-band fluxes and footprints for one patch.

        Returns
        -------
        nSubfilters : `int`
            Number of subfilters used to model the sources in the catalog.

        Raises
        ------
        RuntimeError
            If the records do not all use the same number of subfilters, or if
            the catalog is missing any of the corresponding weight fields.
        """
        nSubfilters = {int(record['numSubfilters']) for record in dcrCatalog}
        if len(nSubfilters) > 1:
            raise RuntimeError("DCR correction catalog mixes records with different numbers of "
                               f"subfilters: {sorted(nSubfilters)}.")
        nSubfilters = nSubfilters.pop()

        names = dcrCatalog.schema.getNames()
        missing = [f"subfilterWeight_{subfilter}" for subfilter in range(nSubfilters)
                   if f"subfilterWeight_{subfilter}" not in names]
        if missing:
            raise RuntimeError(f"DCR correction catalog has numSubfilters={nSubfilters}, but is "
                               f"missing the weight fields {missing}.")
        return nSubfilters
