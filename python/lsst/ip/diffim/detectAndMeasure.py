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

import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
import lsst.geom
from lsst.ip.diffim.utils import evaluateMaskFraction
from lsst.meas.algorithms import SkyObjectsTask, SourceDetectionTask, SetPrimaryFlagsTask, MaskStreaksTask
from lsst.meas.base import ForcedMeasurementTask, ApplyApCorrTask, DetectorVisitIdGeneratorConfig
import lsst.meas.deblender
import lsst.meas.extensions.trailedSources  # noqa: F401
import lsst.meas.extensions.shapeHSM
import lsst.pex.config as pexConfig
from lsst.pex.exceptions import InvalidParameterError
import lsst.pipe.base as pipeBase
import lsst.utils
from lsst.utils.timer import timeMethod

from . import DipoleFitTask

__all__ = ["DetectAndMeasureConfig", "DetectAndMeasureTask",
           "DetectAndMeasureScoreConfig", "DetectAndMeasureScoreTask"]


class DetectAndMeasureConnections(pipeBase.PipelineTaskConnections,
                                  dimensions=("instrument", "visit", "detector"),
                                  defaultTemplates={"coaddName": "deep",
                                                    "warpTypeSuffix": "",
                                                    "fakesType": ""}):
    science = pipeBase.connectionTypes.Input(
        doc="Input science exposure.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}calexp"
    )
    matchedTemplate = pipeBase.connectionTypes.Input(
        doc="Warped and PSF-matched template used to create the difference image.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Diff_matchedExp",
    )
    difference = pipeBase.connectionTypes.Input(
        doc="Result of subtracting template from science.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Diff_differenceTempExp",
    )
    outputSchema = pipeBase.connectionTypes.InitOutput(
        doc="Schema (as an example catalog) for output DIASource catalog.",
        storageClass="SourceCatalog",
        name="{fakesType}{coaddName}Diff_diaSrc_schema",
    )
    diaSources = pipeBase.connectionTypes.Output(
        doc="Detected diaSources on the difference image.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="SourceCatalog",
        name="{fakesType}{coaddName}Diff_diaSrc",
    )
    subtractedMeasuredExposure = pipeBase.connectionTypes.Output(
        doc="Difference image with detection mask plane filled in.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Diff_differenceExp",
    )
    maskedStreaks = pipeBase.connectionTypes.Output(
        doc='Streak profile information.',
        storageClass="ArrowNumpyDict",
        dimensions=("instrument", "visit", "detector"),
        name="{fakesType}{coaddName}Diff_streaks",
    )

    def __init__(self, *, config):
        super().__init__(config=config)
        if not (self.config.writeStreakInfo and self.config.doMaskStreaks):
            self.outputs.remove("maskedStreaks")


class DetectAndMeasureConfig(pipeBase.PipelineTaskConfig,
                             pipelineConnections=DetectAndMeasureConnections):
    """Config for DetectAndMeasureTask
    """
    doMerge = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Merge positive and negative diaSources with grow radius "
            "set by growFootprint",
        deprecated="This field is no longer used and will be removed after v28."
    )
    doForcedMeasurement = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Force photometer diaSource locations on PVI?")
    doAddMetrics = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Add columns to the source table to hold analysis metrics?"
    )
    detection = pexConfig.ConfigurableField(
        target=SourceDetectionTask,
        doc="Final source detection for diaSource measurement",
    )
    deblend = pexConfig.ConfigurableField(
        target=lsst.meas.deblender.SourceDeblendTask,
        doc="Task to split blended sources into their components."
    )
    measurement = pexConfig.ConfigurableField(
        target=DipoleFitTask,
        doc="Task to measure sources on the difference image.",
    )
    doApCorr = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        doc="Run subtask to apply aperture corrections"
    )
    applyApCorr = lsst.pex.config.ConfigurableField(
        target=ApplyApCorrTask,
        doc="Task to apply aperture corrections"
    )
    forcedMeasurement = pexConfig.ConfigurableField(
        target=ForcedMeasurementTask,
        doc="Task to force photometer science image at diaSource locations.",
    )
    growFootprint = pexConfig.Field(
        dtype=int,
        default=2,
        doc="Grow positive and negative footprints by this many pixels before merging",
        deprecated="This field is no longer used and will be removed after v28."
    )
    diaSourceMatchRadius = pexConfig.Field(
        dtype=float,
        default=0.5,
        doc="Match radius (in arcseconds) for DiaSource to Source association"
    )
    doSkySources = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Generate sky sources?",
    )
    skySources = pexConfig.ConfigurableField(
        target=SkyObjectsTask,
        doc="Generate sky sources",
    )
    doMaskStreaks = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Turn on streak masking",
    )
    maskStreaks = pexConfig.ConfigurableField(
        target=MaskStreaksTask,
        doc="Subtask for masking streaks. Only used if doMaskStreaks is True. "
            "Adds a mask plane to an exposure, with the mask plane name set by streakMaskName.",
    )
    writeStreakInfo = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Record the parameters of any detected streaks. For LSST, this should be turned off except for "
            "development work."
    )
    setPrimaryFlags = pexConfig.ConfigurableField(
        target=SetPrimaryFlagsTask,
        doc="Task to add isPrimary and deblending-related flags to the catalog."
    )
    badSourceFlags = lsst.pex.config.ListField(
        dtype=str,
        doc="Sources with any of these flags set are removed before writing the output catalog.",
        default=("base_PixelFlags_flag_offimage",
                 "base_PixelFlags_flag_interpolatedCenterAll",
                 "base_PixelFlags_flag_badCenterAll",
                 "base_PixelFlags_flag_edgeCenterAll",
                 "base_PixelFlags_flag_saturatedCenterAll",
                 ),
    )
    idGenerator = DetectorVisitIdGeneratorConfig.make_field()

    def setDefaults(self):
        # DiaSource Detection
        self.detection.thresholdPolarity = "both"
        self.detection.thresholdValue = 5.0
        self.detection.reEstimateBackground = False
        self.detection.thresholdType = "pixel_stdev"
        self.detection.excludeMaskPlanes = ["EDGE"]

        # Add filtered flux measurement, the correct measurement for pre-convolved images.
        self.measurement.algorithms.names.add("base_PeakLikelihoodFlux")
        self.measurement.plugins.names |= ["ext_trailedSources_Naive",
                                           "base_LocalPhotoCalib",
                                           "base_LocalWcs",
                                           "ext_shapeHSM_HsmSourceMoments",
                                           "ext_shapeHSM_HsmPsfMoments",
                                           ]
        self.measurement.slots.psfShape = "ext_shapeHSM_HsmPsfMoments"
        self.measurement.slots.shape = "ext_shapeHSM_HsmSourceMoments"
        self.measurement.plugins["base_SdssCentroid"].maxDistToPeak = 5.0
        self.forcedMeasurement.plugins = ["base_TransformedCentroid", "base_PsfFlux"]
        self.forcedMeasurement.copyColumns = {
            "id": "objectId", "parent": "parentObjectId", "coord_ra": "coord_ra", "coord_dec": "coord_dec"}
        self.forcedMeasurement.slots.centroid = "base_TransformedCentroid"
        self.forcedMeasurement.slots.shape = None

        # Keep track of which footprints contain streaks
        self.measurement.plugins["base_PixelFlags"].masksFpAnywhere = [
            "STREAK", "INJECTED", "INJECTED_TEMPLATE"]
        self.measurement.plugins["base_PixelFlags"].masksFpCenter = [
            "STREAK", "INJECTED", "INJECTED_TEMPLATE"]
        self.skySources.avoidMask = ["DETECTED", "DETECTED_NEGATIVE", "BAD", "NO_DATA", "EDGE"]


class DetectAndMeasureTask(lsst.pipe.base.PipelineTask):
    """Detect and measure sources on a difference image.
    """
    ConfigClass = DetectAndMeasureConfig
    _DefaultName = "detectAndMeasure"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        # Add coordinate error fields:
        afwTable.CoordKey.addErrorFields(self.schema)

        self.algMetadata = dafBase.PropertyList()
        self.makeSubtask("detection", schema=self.schema)
        self.makeSubtask("deblend", schema=self.schema)
        self.makeSubtask("setPrimaryFlags", schema=self.schema, isSingleFrame=True)
        self.makeSubtask("measurement", schema=self.schema,
                         algMetadata=self.algMetadata)
        if self.config.doApCorr:
            self.makeSubtask("applyApCorr", schema=self.measurement.schema)
        if self.config.doForcedMeasurement:
            self.schema.addField(
                "ip_diffim_forced_PsfFlux_instFlux", "D",
                "Forced PSF flux measured on the direct image.",
                units="count")
            self.schema.addField(
                "ip_diffim_forced_PsfFlux_instFluxErr", "D",
                "Forced PSF flux error measured on the direct image.",
                units="count")
            self.schema.addField(
                "ip_diffim_forced_PsfFlux_area", "F",
                "Forced PSF flux effective area of PSF.",
                units="pixel")
            self.schema.addField(
                "ip_diffim_forced_PsfFlux_flag", "Flag",
                "Forced PSF flux general failure flag.")
            self.schema.addField(
                "ip_diffim_forced_PsfFlux_flag_noGoodPixels", "Flag",
                "Forced PSF flux not enough non-rejected pixels in data to attempt the fit.")
            self.schema.addField(
                "ip_diffim_forced_PsfFlux_flag_edge", "Flag",
                "Forced PSF flux object was too close to the edge of the image to use the full PSF model.")
            self.makeSubtask("forcedMeasurement", refSchema=self.schema)

        self.schema.addField("refMatchId", "L", "unique id of reference catalog match")
        self.schema.addField("srcMatchId", "L", "unique id of source match")
        if self.config.doSkySources:
            self.makeSubtask("skySources", schema=self.schema)
        if self.config.doMaskStreaks:
            self.makeSubtask("maskStreaks")

        # Check that the schema and config are consistent
        for flag in self.config.badSourceFlags:
            if flag not in self.schema:
                raise pipeBase.InvalidQuantumError("Field %s not in schema" % flag)

        # initialize InitOutputs
        self.outputSchema = afwTable.SourceCatalog(self.schema)
        self.outputSchema.getTable().setMetadata(self.algMetadata)

    def runQuantum(self, butlerQC: pipeBase.QuantumContext,
                   inputRefs: pipeBase.InputQuantizedConnection,
                   outputRefs: pipeBase.OutputQuantizedConnection):
        inputs = butlerQC.get(inputRefs)
        idGenerator = self.config.idGenerator.apply(butlerQC.quantum.dataId)
        idFactory = idGenerator.make_table_id_factory()
        outputs = self.run(**inputs, idFactory=idFactory)
        butlerQC.put(outputs, outputRefs)

    @timeMethod
    def run(self, science, matchedTemplate, difference,
            idFactory=None):
        """Detect and measure sources on a difference image.

        The difference image will be convolved with a gaussian approximation of
        the PSF to form a maximum likelihood image for detection.
        Close positive and negative detections will optionally be merged into
        dipole diaSources.
        Sky sources, or forced detections in background regions, will optionally
        be added, and the configured measurement algorithm will be run on all
        detections.

        Parameters
        ----------
        science : `lsst.afw.image.ExposureF`
            Science exposure that the template was subtracted from.
        matchedTemplate : `lsst.afw.image.ExposureF`
            Warped and PSF-matched template that was used produce the
            difference image.
        difference : `lsst.afw.image.ExposureF`
            Result of subtracting template from the science image.
        idFactory : `lsst.afw.table.IdFactory`, optional
            Generator object used to assign ids to detected sources in the
            difference image. Ids from this generator are not set until after
            deblending and merging positive/negative peaks.

        Returns
        -------
        measurementResults : `lsst.pipe.base.Struct`

            ``subtractedMeasuredExposure`` : `lsst.afw.image.ExposureF`
                Subtracted exposure with detection mask applied.
            ``diaSources``  : `lsst.afw.table.SourceCatalog`
                The catalog of detected sources.
        """
        if idFactory is None:
            idFactory = lsst.meas.base.IdGenerator().make_table_id_factory()

        # Ensure that we start with an empty detection and deblended mask.
        mask = difference.mask
        clearMaskPlanes = ["DETECTED", "DETECTED_NEGATIVE", "NOT_DEBLENDED"]
        for mp in clearMaskPlanes:
            if mp not in mask.getMaskPlaneDict():
                mask.addMaskPlane(mp)
        mask &= ~mask.getPlaneBitMask(clearMaskPlanes)

        # Don't use the idFactory until after deblend+merge, so that we aren't
        # generating ids that just get thrown away (footprint merge doesn't
        # know about past ids).
        table = afwTable.SourceTable.make(self.schema)
        results = self.detection.run(
            table=table,
            exposure=difference,
            doSmooth=True,
        )

        sources = self._buildCatalogAndDeblend(difference, results.positive, results.negative, idFactory)

        return self._measureSources(science, matchedTemplate, difference, sources)

    def _measureSources(self, science, matchedTemplate, difference, initialDiaSources):
        """Measure and process the results of source detection.

        Parameters
        ----------
        science : `lsst.afw.image.ExposureF`
            Science exposure that the template was subtracted from.
        matchedTemplate : `lsst.afw.image.ExposureF`
            Warped and PSF-matched template that was used produce the
            difference image.
        difference : `lsst.afw.image.ExposureF`
            Result of subtracting template from the science image.
        initialDiaSources : `lsst.afw.table.SourceCatalog`
            Detected sources on the difference exposure.

        Returns
        -------
        measurementResults : `lsst.pipe.base.Struct`

            ``subtractedMeasuredExposure`` : `lsst.afw.image.ExposureF`
                Subtracted exposure with detection mask applied.
            ``diaSources``  : `lsst.afw.table.SourceCatalog`
                The catalog of detected sources.
        """
        self.metadata.add("nDiaSources", len(initialDiaSources))
        initialDiaSources.setMetadata(self.algMetadata)

        if self.config.doMaskStreaks:
            streakInfo = self._runStreakMasking(difference.maskedImage)

        self.measureDiaSources(initialDiaSources, science, difference, matchedTemplate)
        diaSources = self._removeBadSources(initialDiaSources)

        if self.config.doForcedMeasurement:
            self.measureForcedSources(diaSources, science, difference.getWcs())

        self.calculateMetrics(difference)

        measurementResults = pipeBase.Struct(
            subtractedMeasuredExposure=difference,
            diaSources=diaSources,
        )
        if self.config.doMaskStreaks and self.config.writeStreakInfo:
            measurementResults.mergeItems(streakInfo, 'maskedStreaks')

        return measurementResults

    def _buildCatalogAndDeblend(self, difference, positiveFootprints, negativeFootprints, idFactory):
        """Create a source catalog and deblend when possible.

        This method creates a source catalog from the positive and negative
        footprints, and then deblends the footprints that have only positive
        peaks. This is because `lsst.meas.deblender.SourceDeblendTask` is not
        designed to handle footprints that have negative flux, resulting in
        garbage output.

        Parameters
        ----------
        difference : `lsst.afw.image.Exposure`
            Result of subtracting template from the science image.
        positiveFootprints, negativeFootprints : `lsst.afw.detection.FootprintSet`
            Positive and negative polarity footprints measured on
            ``difference`` to be deblended separately.
        idFactory : `lsst.afw.table.IdFactory`
            Generator object used to assign ids to detected sources in the
            difference image.

        Returns
        -------
        sources : `lsst.afw.table.SourceCatalog`
            Positive and negative deblended children.
        positives, negatives : `lsst.afw.detection.FootprintSet`
            Deblended positive and negative polarity footprints measured on
            ``difference``.
        """
        # Merge the positive and negative footprints.
        # The original detection FootprintSets already grew each detection,
        # so there is no need to grow them again.
        merged_footprints = positiveFootprints
        merged_footprints.merge(negativeFootprints, 0, 0, False)

        # Create a source catalog from the footprints.
        table = afwTable.SourceTable.make(self.schema, idFactory)
        sources = afwTable.SourceCatalog(table)
        merged_footprints.makeSources(sources)

        # Sky sources must be added before deblending, otherwise the
        # sources with parent == 0 will be out of order and
        # downstream measurement tasks cannot run.
        if self.config.doSkySources:
            self.addSkySources(sources, difference.mask, difference.info.id)

        # Find the footprints with only positive peaks and no negative peaks.
        footprints = [src.getFootprint() for src in sources]
        nPeaks = np.array([len(fp.peaks) for fp in footprints])
        blend_ids = []
        skipped_ids = []
        for src in sources[nPeaks > 1]:
            peaks = src.getFootprint().peaks
            peak_flux = np.array([peak.getPeakValue() for peak in peaks])
            if np.all(peak_flux > 0) or np.all(peak_flux < 0):
                blend_ids.append(src.getId())
            else:
                skipped_ids.append(src.getId())

        # Deblend the footprints that can be deblended with SourceDeblendTask.
        temp_cat = afwTable.SourceCatalog(sources.table)
        for sid in blend_ids:
            temp_cat.append(sources.find(sid))
        first_child_index = len(temp_cat)
        self.deblend.run(exposure=difference, sources=temp_cat)

        # Update the deblended parents and add their children
        # to the sources catalog.
        sources.extend(temp_cat[first_child_index:], deep=False)

        # Set deblending flags.
        # Since SourceDeblendTask already has a method for setting
        # flags for skipped parents, we just call that method here
        # instead of setting the flags manually.
        for sid in skipped_ids:
            src = sources.find(sid)
            self.deblend.skipParent(src, difference.mask)

        # Set detection and primary flags
        self.setPrimaryFlags.run(sources)

        table = afwTable.SourceTable.make(self.schema, idFactory)
        _sources = afwTable.SourceCatalog(table)
        _sources.extend(sources, deep=True)

        return _sources

    def _removeBadSources(self, diaSources):
        """Remove unphysical diaSources from the catalog.

        Parameters
        ----------
        diaSources : `lsst.afw.table.SourceCatalog`
            The catalog of detected sources.

        Returns
        -------
        diaSources : `lsst.afw.table.SourceCatalog`
            The updated catalog of detected sources, with any source that has a
            flag in ``config.badSourceFlags`` set removed.
        """
        selector = np.ones(len(diaSources), dtype=bool)
        for flag in self.config.badSourceFlags:
            flags = diaSources[flag]
            nBad = np.count_nonzero(flags)
            if nBad > 0:
                self.log.debug("Found %d unphysical sources with flag %s.", nBad, flag)
                selector &= ~flags
        nBadTotal = np.count_nonzero(~selector)
        self.metadata.add("nRemovedBadFlaggedSources", nBadTotal)
        self.log.info("Removed %d unphysical sources.", nBadTotal)
        return diaSources[selector].copy(deep=True)

    def addSkySources(self, diaSources, mask, seed,
                      subtask=None):
        """Add sources in empty regions of the difference image
        for measuring the background.

        Parameters
        ----------
        diaSources : `lsst.afw.table.SourceCatalog`
            The catalog of detected sources.
        mask : `lsst.afw.image.Mask`
            Mask plane for determining regions where Sky sources can be added.
        seed : `int`
            Seed value to initialize the random number generator.
        """
        if subtask is None:
            subtask = self.skySources
        skySourceFootprints = subtask.run(mask=mask, seed=seed, catalog=diaSources)
        self.metadata.add(f"n_{subtask.getName()}", len(skySourceFootprints))

    def measureDiaSources(self, diaSources, science, difference, matchedTemplate):
        """Use (matched) template and science image to constrain dipole fitting.

        Parameters
        ----------
        diaSources : `lsst.afw.table.SourceCatalog`
            The catalog of detected sources.
        science : `lsst.afw.image.ExposureF`
            Science exposure that the template was subtracted from.
        difference : `lsst.afw.image.ExposureF`
            Result of subtracting template from the science image.
        matchedTemplate : `lsst.afw.image.ExposureF`
            Warped and PSF-matched template that was used produce the
            difference image.
        """
        # Ensure that the required mask planes are present
        for mp in self.config.measurement.plugins["base_PixelFlags"].masksFpAnywhere:
            difference.mask.addMaskPlane(mp)
        # Note that this may not be correct if we convolved the science image.
        # In the future we may wish to persist the matchedScience image.
        self.measurement.run(diaSources, difference, science, matchedTemplate)
        if self.config.doApCorr:
            apCorrMap = difference.getInfo().getApCorrMap()
            if apCorrMap is None:
                self.log.warning("Difference image does not have valid aperture correction; skipping.")
            else:
                self.applyApCorr.run(
                    catalog=diaSources,
                    apCorrMap=apCorrMap,
                )

    def measureForcedSources(self, diaSources, science, wcs):
        """Perform forced measurement of the diaSources on the science image.

        Parameters
        ----------
        diaSources : `lsst.afw.table.SourceCatalog`
            The catalog of detected sources.
        science : `lsst.afw.image.ExposureF`
            Science exposure that the template was subtracted from.
        wcs : `lsst.afw.geom.SkyWcs`
            Coordinate system definition (wcs) for the exposure.
        """
        # Run forced psf photometry on the PVI at the diaSource locations.
        # Copy the measured flux and error into the diaSource.
        forcedSources = self.forcedMeasurement.generateMeasCat(science, diaSources, wcs)
        self.forcedMeasurement.run(forcedSources, science, diaSources, wcs)
        mapper = afwTable.SchemaMapper(forcedSources.schema, diaSources.schema)
        mapper.addMapping(forcedSources.schema.find("base_PsfFlux_instFlux")[0],
                          "ip_diffim_forced_PsfFlux_instFlux", True)
        mapper.addMapping(forcedSources.schema.find("base_PsfFlux_instFluxErr")[0],
                          "ip_diffim_forced_PsfFlux_instFluxErr", True)
        mapper.addMapping(forcedSources.schema.find("base_PsfFlux_area")[0],
                          "ip_diffim_forced_PsfFlux_area", True)
        mapper.addMapping(forcedSources.schema.find("base_PsfFlux_flag")[0],
                          "ip_diffim_forced_PsfFlux_flag", True)
        mapper.addMapping(forcedSources.schema.find("base_PsfFlux_flag_noGoodPixels")[0],
                          "ip_diffim_forced_PsfFlux_flag_noGoodPixels", True)
        mapper.addMapping(forcedSources.schema.find("base_PsfFlux_flag_edge")[0],
                          "ip_diffim_forced_PsfFlux_flag_edge", True)
        for diaSource, forcedSource in zip(diaSources, forcedSources):
            diaSource.assign(forcedSource, mapper)

    def calculateMetrics(self, difference):
        """Add image QA metrics to the Task metadata.

        Parameters
        ----------
        difference : `lsst.afw.image.Exposure`
            The target image to calculate metrics for.

        """
        mask = difference.mask
        badPix = (mask.array & mask.getPlaneBitMask(self.config.detection.excludeMaskPlanes)) > 0
        self.metadata.add("nGoodPixels", np.sum(~badPix))
        self.metadata.add("nBadPixels", np.sum(badPix))
        detPosPix = (mask.array & mask.getPlaneBitMask("DETECTED")) > 0
        detNegPix = (mask.array & mask.getPlaneBitMask("DETECTED_NEGATIVE")) > 0
        self.metadata.add("nPixelsDetectedPositive", np.sum(detPosPix))
        self.metadata.add("nPixelsDetectedNegative", np.sum(detNegPix))
        detPosPix &= badPix
        detNegPix &= badPix
        self.metadata.add("nBadPixelsDetectedPositive", np.sum(detPosPix))
        self.metadata.add("nBadPixelsDetectedNegative", np.sum(detNegPix))

        metricsMaskPlanes = list(mask.getMaskPlaneDict().keys())
        for maskPlane in metricsMaskPlanes:
            try:
                self.metadata.add("%s_mask_fraction"%maskPlane.lower(), evaluateMaskFraction(mask, maskPlane))
            except InvalidParameterError:
                self.metadata.add("%s_mask_fraction"%maskPlane.lower(), -1)
                self.log.info("Unable to calculate metrics for mask plane %s: not in image"%maskPlane)

    def _runStreakMasking(self, maskedImage):
        """Do streak masking at put results into catalog.

        Parameters
        ----------
        maskedImage: `lsst.afw.image.maskedImage`
            The image in which to search for streaks. Must have a detection
            mask.

        Returns
        -------
        streakInfo: `lsst.pipe.base.Struct`
            ``rho`` : `np.ndarray`
                Angle of detected streak.
            ``theta`` : `np.ndarray`
                Distance from center of detected streak.
            ``sigma`` : `np.ndarray`
                Width of streak profile.
        """
        streaks = self.maskStreaks.run(maskedImage)
        if self.config.writeStreakInfo:
            rhos = np.array([line.rho for line in streaks.lines])
            thetas = np.array([line.theta for line in streaks.lines])
            sigmas = np.array([line.sigma for line in streaks.lines])
            streakInfo = {'rho': rhos, 'theta': thetas, 'sigma': sigmas}
        else:
            streakInfo = {'rho': np.array([]), 'theta': np.array([]), 'sigma': np.array([])}
        return pipeBase.Struct(maskedStreaks=streakInfo)


class DetectAndMeasureScoreConnections(DetectAndMeasureConnections):
    scoreExposure = pipeBase.connectionTypes.Input(
        doc="Maximum likelihood image for detection.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Diff_scoreExp",
    )


class DetectAndMeasureScoreConfig(DetectAndMeasureConfig,
                                  pipelineConnections=DetectAndMeasureScoreConnections):
    pass


class DetectAndMeasureScoreTask(DetectAndMeasureTask):
    """Detect DIA sources using a score image,
    and measure the detections on the difference image.

    Source detection is run on the supplied score, or maximum likelihood,
    image. Note that no additional convolution will be done in this case.
    Close positive and negative detections will optionally be merged into
    dipole diaSources.
    Sky sources, or forced detections in background regions, will optionally
    be added, and the configured measurement algorithm will be run on all
    detections.
    """
    ConfigClass = DetectAndMeasureScoreConfig
    _DefaultName = "detectAndMeasureScore"

    @timeMethod
    def run(self, science, matchedTemplate, difference, scoreExposure,
            idFactory=None):
        """Detect and measure sources on a score image.

        Parameters
        ----------
        science : `lsst.afw.image.ExposureF`
            Science exposure that the template was subtracted from.
        matchedTemplate : `lsst.afw.image.ExposureF`
            Warped and PSF-matched template that was used produce the
            difference image.
        difference : `lsst.afw.image.ExposureF`
            Result of subtracting template from the science image.
        scoreExposure : `lsst.afw.image.ExposureF`
            Score or maximum likelihood difference image
        idFactory : `lsst.afw.table.IdFactory`, optional
            Generator object used to assign ids to detected sources in the
            difference image. Ids from this generator are not set until after
            deblending and merging positive/negative peaks.

        Returns
        -------
        measurementResults : `lsst.pipe.base.Struct`

            ``subtractedMeasuredExposure`` : `lsst.afw.image.ExposureF`
                Subtracted exposure with detection mask applied.
            ``diaSources``  : `lsst.afw.table.SourceCatalog`
                The catalog of detected sources.
        """
        if idFactory is None:
            idFactory = lsst.meas.base.IdGenerator().make_table_id_factory()

        # Ensure that we start with an empty detection mask.
        mask = scoreExposure.mask
        mask &= ~(mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE"))

        # Don't use the idFactory until after deblend+merge, so that we aren't
        # generating ids that just get thrown away (footprint merge doesn't
        # know about past ids).
        table = afwTable.SourceTable.make(self.schema)
        results = self.detection.run(
            table=table,
            exposure=scoreExposure,
            doSmooth=False,
        )
        # Copy the detection mask from the Score image to the difference image
        difference.mask.assign(scoreExposure.mask, scoreExposure.getBBox())

        sources = self._buildCatalogAndDeblend(difference, results.positive, results.negative, idFactory)

        return self._measureSources(science, matchedTemplate, difference, sources)
