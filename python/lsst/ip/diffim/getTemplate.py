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
from deprecated.sphinx import deprecated

import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
from lsst.afw.math._warper import computeWarpedBBox
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from lsst.skymap import BaseSkyMap
from lsst.ip.diffim.dcrModel import DcrModel
from lsst.meas.algorithms import CoaddPsf, CoaddPsfConfig
from lsst.utils.timer import timeMethod

__all__ = [
    "GetTemplateTask",
    "GetTemplateConfig",
    "GetDcrTemplateTask",
    "GetDcrTemplateConfig",
]


class GetTemplateConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("instrument", "visit", "detector", "skymap"),
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

    template = pipeBase.connectionTypes.Output(
        doc="Warped template, pixel matched to the bounding box and WCS.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Diff_templateExp{warpTypeSuffix}",
    )


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

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        bbox = inputs.pop("bbox")
        wcs = inputs.pop("wcs")
        coaddExposures = inputs.pop("coaddExposures")
        skymap = inputs.pop("skyMap")

        # This should not happen with a properly configured execution context.
        assert not inputs, "runQuantum got more inputs than expected"

        results = self.getExposures(coaddExposures, bbox, skymap, wcs)
        physical_filter = butlerQC.quantum.dataId["physical_filter"]
        outputs = self.run(
            coaddExposureHandles=results.coaddExposures,
            bbox=bbox,
            wcs=wcs,
            dataIds=results.dataIds,
            physical_filter=physical_filter,
        )
        butlerQC.put(outputs, outputRefs)

    @deprecated(
        reason="Replaced by getExposures, which uses explicit arguments instead of a kwargs dict. "
        "This method will be removed after v29.",
        version="v29.0",
        category=FutureWarning,
    )
    def getOverlappingExposures(self, inputs):
        return self.getExposures(
            inputs["coaddExposures"], inputs["bbox"], inputs["skyMap"], inputs["wcs"]
        )

    def getExposures(self, coaddExposureHandles, bbox, skymap, wcs):
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
        overlappingArea = 0
        coaddExposures = collections.defaultdict(list)
        dataIds = collections.defaultdict(list)

        for coaddRef in coaddExposureHandles:
            dataId = coaddRef.dataId
            patchWcs = skymap[dataId["tract"]].getWcs()
            patchBBox = skymap[dataId["tract"]][dataId["patch"]].getOuterBBox()
            patchCorners = patchWcs.pixelToSky(geom.Box2D(patchBBox).getCorners())
            patchPolygon = afwGeom.Polygon(wcs.skyToPixel(patchCorners))
            if patchPolygon.intersection(detectorPolygon):
                overlappingArea += patchPolygon.intersectionSingle(
                    detectorPolygon
                ).calculateArea()
                self.log.info(
                    "Using template input tract=%s, patch=%s",
                    dataId["tract"],
                    dataId["patch"],
                )
                coaddExposures[dataId["tract"]].append(coaddRef)
                dataIds[dataId["tract"]].append(dataId)

        if not overlappingArea:
            raise pipeBase.NoWorkFound("No patches overlap detector")

        return pipeBase.Struct(coaddExposures=coaddExposures, dataIds=dataIds)

    @timeMethod
    def run(self, *, coaddExposureHandles, bbox, wcs, dataIds, physical_filter):
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
            maskedImages, catalog, totalBox = self._makeExposureCatalog(
                coaddExposureHandles[tract], dataIds[tract]
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

        template.setPsf(self._makePsf(template, catalog, wcs))
        template.setFilter(afwImage.FilterLabel(band, physical_filter))
        template.setPhotoCalib(photoCalib)
        return pipeBase.Struct(template=template)

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

    def _makeExposureCatalog(self, exposureRefs, dataIds):
        """Make an exposure catalog for one tract.

        Parameters
        ----------
        exposureRefs : `list` of [`lsst.daf.butler.DeferredDatasetHandle` of \
                        `lsst.afw.image.Exposure`]
            Exposures to include in the catalog.
        dataIds : `list` [`lsst.daf.butler.DataCoordinate`]
            Data ids of each of the included exposures; must have "tract" and
            "patch" entries.

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
        images = {}
        totalBox = geom.Box2I()

        for coadd, dataId in zip(exposures, dataIds):
            images[dataId] = coadd.maskedImage
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


class GetDcrTemplateConnections(
    GetTemplateConnections,
    dimensions=("instrument", "visit", "detector", "skymap"),
    defaultTemplates={"coaddName": "dcr", "warpTypeSuffix": "", "fakesType": ""},
):
    visitInfo = pipeBase.connectionTypes.Input(
        doc="VisitInfo of calexp used to determine observing conditions.",
        name="{fakesType}calexp.visitInfo",
        storageClass="VisitInfo",
        dimensions=("instrument", "visit", "detector"),
    )
    dcrCoadds = pipeBase.connectionTypes.Input(
        doc="Input DCR template to match and subtract from the exposure",
        name="{fakesType}dcrCoadd{warpTypeSuffix}",
        storageClass="ExposureF",
        dimensions=("tract", "patch", "skymap", "band", "subfilter"),
        multiple=True,
        deferLoad=True,
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        self.inputs.remove("coaddExposures")


class GetDcrTemplateConfig(
    GetTemplateConfig, pipelineConnections=GetDcrTemplateConnections
):
    numSubfilters = pexConfig.Field(
        doc="Number of subfilters in the DcrCoadd.",
        dtype=int,
        default=3,
    )
    effectiveWavelength = pexConfig.Field(
        doc="Effective wavelength of the filter in nm.",
        optional=False,
        dtype=float,
    )
    bandwidth = pexConfig.Field(
        doc="Bandwidth of the physical filter.",
        optional=False,
        dtype=float,
    )

    def validate(self):
        if self.effectiveWavelength is None or self.bandwidth is None:
            raise ValueError(
                "The effective wavelength and bandwidth of the physical filter "
                "must be set in the getTemplate config for DCR coadds. "
                "Required until transmission curves are used in DM-13668."
            )


class GetDcrTemplateTask(GetTemplateTask):
    ConfigClass = GetDcrTemplateConfig
    _DefaultName = "getDcrTemplate"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        bbox = inputs.pop("bbox")
        wcs = inputs.pop("wcs")
        dcrCoaddExposureHandles = inputs.pop("dcrCoadds")
        skymap = inputs.pop("skyMap")
        visitInfo = inputs.pop("visitInfo")

        # This should not happen with a properly configured execution context.
        assert not inputs, "runQuantum got more inputs than expected"

        results = self.getExposures(
            dcrCoaddExposureHandles, bbox, skymap, wcs, visitInfo
        )
        physical_filter = butlerQC.quantum.dataId["physical_filter"]
        outputs = self.run(
            coaddExposureHandles=results.coaddExposures,
            bbox=bbox,
            wcs=wcs,
            dataIds=results.dataIds,
            physical_filter=physical_filter,
        )
        butlerQC.put(outputs, outputRefs)

    @deprecated(
        reason="Replaced by getExposures, which uses explicit arguments instead of a kwargs dict. "
        "This method will be removed after v29.",
        version="v29.0",
        category=FutureWarning,
    )
    def getOverlappingExposures(self, inputs):
        return self.getExposures(
            inputs["dcrCoadds"],
            inputs["bbox"],
            inputs["skyMap"],
            inputs["wcs"],
            inputs["visitInfo"],
        )

    def getExposures(self, dcrCoaddExposureHandles, bbox, skymap, wcs, visitInfo):
        """Return lists of coadds and their corresponding dataIds that overlap
        the detector.

        The spatial index in the registry has generous padding and often
        supplies patches near, but not directly overlapping the detector.
        Filters inputs so that we don't have to read in all input coadds.

        Parameters
        ----------
        dcrCoaddExposureHandles :  `list` \
                                  [`lsst.daf.butler.DeferredDatasetHandle` of \
                                  `lsst.afw.image.Exposure`]
            Data references to exposures that might overlap the detector.
        bbox : `lsst.geom.Box2I`
            Template Bounding box of the detector geometry onto which to
            resample the coaddExposures.
        skymap : `lsst.skymap.SkyMap`
            Input definition of geometry/bbox and projection/wcs for
            template exposures.
        wcs : `lsst.afw.geom.SkyWcs`
            Template WCS onto which to resample the coaddExposures.
        visitInfo : `lsst.afw.image.VisitInfo`
            Metadata for the science image.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
           A struct with attibutes:

           ``coaddExposures``
               Dict of coadd exposures that overlap the projected bbox,
               indexed on tract id
               (`dict` [`int`, `list` [`lsst.afw.image.Exposure`] ]).
           ``dataIds``
               Dict of data IDs of the coadd exposures that overlap the
               projected bbox, indexed on tract id
               (`dict` [`int`, `list [`lsst.daf.butler.DataCoordinate`] ]).

        Raises
        ------
        pipeBase.NoWorkFound
            Raised if no patches overlatp the input detector bbox.
        """
        # Check that the patches actually overlap the detector
        # Exposure's validPolygon would be more accurate
        if wcs is None:
            raise pipeBase.NoWorkFound("Exposure has no WCS; cannot create a template.")

        detectorPolygon = geom.Box2D(bbox)
        overlappingArea = 0
        dataIds = collections.defaultdict(list)
        patchList = dict()
        for coaddRef in dcrCoaddExposureHandles:
            dataId = coaddRef.dataId
            subfilter = dataId["subfilter"]
            patchWcs = skymap[dataId["tract"]].getWcs()
            patchBBox = skymap[dataId["tract"]][dataId["patch"]].getOuterBBox()
            patchCorners = patchWcs.pixelToSky(geom.Box2D(patchBBox).getCorners())
            patchPolygon = afwGeom.Polygon(wcs.skyToPixel(patchCorners))
            if patchPolygon.intersection(detectorPolygon):
                overlappingArea += patchPolygon.intersectionSingle(
                    detectorPolygon
                ).calculateArea()
                self.log.info(
                    "Using template input tract=%s, patch=%s, subfilter=%s"
                    % (dataId["tract"], dataId["patch"], dataId["subfilter"])
                )
                if dataId["tract"] in patchList:
                    patchList[dataId["tract"]].append(dataId["patch"])
                else:
                    patchList[dataId["tract"]] = [
                        dataId["patch"],
                    ]
                if subfilter == 0:
                    dataIds[dataId["tract"]].append(dataId)

        if not overlappingArea:
            raise pipeBase.NoWorkFound("No patches overlap detector")

        self.checkPatchList(patchList)

        coaddExposures = self.getDcrModel(patchList, dcrCoaddExposureHandles, visitInfo)
        return pipeBase.Struct(coaddExposures=coaddExposures, dataIds=dataIds)

    def checkPatchList(self, patchList):
        """Check that all of the DcrModel subfilters are present for each
        patch.

        Parameters
        ----------
        patchList : `dict`
            Dict of the patches containing valid data for each tract.

        Raises
        ------
        RuntimeError
            If the number of exposures found for a patch does not match the
            number of subfilters.
        """
        for tract in patchList:
            for patch in set(patchList[tract]):
                if patchList[tract].count(patch) != self.config.numSubfilters:
                    raise RuntimeError(
                        "Invalid number of DcrModel subfilters found: %d vs %d expected",
                        patchList[tract].count(patch),
                        self.config.numSubfilters,
                    )

    def getDcrModel(self, patchList, coaddRefs, visitInfo):
        """Build DCR-matched coadds from a list of exposure references.

        Parameters
        ----------
        patchList : `dict`
            Dict of the patches containing valid data for each tract.
        coaddRefs : `list` [`lsst.daf.butler.DeferredDatasetHandle`]
            Data references to `~lsst.afw.image.Exposure` representing
            DcrModels that overlap the detector.
        visitInfo : `lsst.afw.image.VisitInfo`
            Metadata for the science image.

        Returns
        -------
        coaddExposures : `list` [`lsst.afw.image.Exposure`]
            Coadd exposures that overlap the detector.
        """
        coaddExposures = collections.defaultdict(list)
        for tract in patchList:
            for patch in set(patchList[tract]):
                coaddRefList = [
                    coaddRef
                    for coaddRef in coaddRefs
                    if _selectDataRef(coaddRef, tract, patch)
                ]

                dcrModel = DcrModel.fromQuantum(
                    coaddRefList,
                    self.config.effectiveWavelength,
                    self.config.bandwidth,
                    self.config.numSubfilters,
                )
                coaddExposures[tract].append(dcrModel.buildMatchedExposureHandle(visitInfo=visitInfo))
        return coaddExposures


def _selectDataRef(coaddRef, tract, patch):
    condition = (coaddRef.dataId["tract"] == tract) & (
        coaddRef.dataId["patch"] == patch
    )
    return condition
