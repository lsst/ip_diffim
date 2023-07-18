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

from deprecated.sphinx import deprecated

import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
from lsst.meas.algorithms import SkyObjectsTask, SourceDetectionTask
from lsst.meas.base import ForcedMeasurementTask, ApplyApCorrTask, DetectorVisitIdGeneratorConfig
import lsst.meas.extensions.trailedSources  # noqa: F401
import lsst.meas.extensions.shapeHSM
from lsst.obs.base import ExposureIdInfo
import lsst.pex.config as pexConfig
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


class DetectAndMeasureConfig(pipeBase.PipelineTaskConfig,
                             pipelineConnections=DetectAndMeasureConnections):
    """Config for DetectAndMeasureTask
    """
    doMerge = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Merge positive and negative diaSources with grow radius "
            "set by growFootprint"
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
        doc="Grow positive and negative footprints by this many pixels before merging"
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
    idGenerator = DetectorVisitIdGeneratorConfig.make_field()

    def setDefaults(self):
        # DiaSource Detection
        self.detection.thresholdPolarity = "both"
        self.detection.thresholdValue = 5.0
        self.detection.reEstimateBackground = False
        self.detection.thresholdType = "pixel_stdev"
        self.detection.excludeMaskPlanes = ["EDGE"]

        # Add filtered flux measurement, the correct measurement for pre-convolved images.
        self.measurement.algorithms.names.add('base_PeakLikelihoodFlux')
        self.measurement.plugins.names |= ['ext_trailedSources_Naive',
                                           'base_LocalPhotoCalib',
                                           'base_LocalWcs',
                                           'ext_shapeHSM_HsmSourceMoments',
                                           'ext_shapeHSM_HsmPsfMoments',
                                           ]
        self.measurement.slots.psfShape = "ext_shapeHSM_HsmPsfMoments"
        self.measurement.slots.shape = "ext_shapeHSM_HsmSourceMoments"
        self.measurement.plugins["base_NaiveCentroid"].maxDistToPeak = 5.0
        self.measurement.plugins["base_SdssCentroid"].maxDistToPeak = 5.0
        self.forcedMeasurement.plugins = ["base_TransformedCentroid", "base_PsfFlux"]
        self.forcedMeasurement.copyColumns = {
            "id": "objectId", "parent": "parentObjectId", "coord_ra": "coord_ra", "coord_dec": "coord_dec"}
        self.forcedMeasurement.slots.centroid = "base_TransformedCentroid"
        self.forcedMeasurement.slots.shape = None


class DetectAndMeasureTask(lsst.pipe.base.PipelineTask):
    """Detect and measure sources on a difference image.
    """
    ConfigClass = DetectAndMeasureConfig
    _DefaultName = "detectAndMeasure"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.schema = afwTable.SourceTable.makeMinimalSchema()

        self.algMetadata = dafBase.PropertyList()
        self.makeSubtask("detection", schema=self.schema)
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
            self.makeSubtask("skySources")
            self.skySourceKey = self.schema.addField("sky_source", type="Flag", doc="Sky objects.")

        # initialize InitOutputs
        self.outputSchema = afwTable.SourceCatalog(self.schema)
        self.outputSchema.getTable().setMetadata(self.algMetadata)

    @staticmethod
    @deprecated(
        reason=(
            "ID factory construction now depends on configuration; use the "
            "idGenerator config field. Will be removed after v27."
        ),
        version="v26.0",
        category=FutureWarning,
    )
    def makeIdFactory(expId, expBits):
        """Create IdFactory instance for unique 64 bit diaSource id-s.

        Parameters
        ----------
        expId : `int`
            Exposure id.

        expBits: `int`
            Number of used bits in ``expId``.

        Notes
        -----
        The diasource id-s consists of the ``expId`` stored fixed in the highest value
        ``expBits`` of the 64-bit integer plus (bitwise or) a generated sequence number in the
        low value end of the integer.

        Returns
        -------
        idFactory: `lsst.afw.table.IdFactory`
        """
        return ExposureIdInfo(expId, expBits).makeSourceIdFactory()

    def runQuantum(self, butlerQC: pipeBase.ButlerQuantumContext,
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
            Generator object to assign ids to detected sources in the difference image.

        Returns
        -------
        measurementResults : `lsst.pipe.base.Struct`

            ``subtractedMeasuredExposure`` : `lsst.afw.image.ExposureF`
                Subtracted exposure with detection mask applied.
            ``diaSources``  : `lsst.afw.table.SourceCatalog`
                The catalog of detected sources.
        """
        # Ensure that we start with an empty detection mask.
        mask = difference.mask
        mask &= ~(mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE"))

        table = afwTable.SourceTable.make(self.schema, idFactory)
        table.setMetadata(self.algMetadata)
        results = self.detection.run(
            table=table,
            exposure=difference,
            doSmooth=True,
        )

        return self.processResults(science, matchedTemplate, difference, results.sources, table,
                                   positiveFootprints=results.positive, negativeFootprints=results.negative)

    def processResults(self, science, matchedTemplate, difference, sources, table,
                       positiveFootprints=None, negativeFootprints=None,):
        """Measure and process the results of source detection.

        Parameters
        ----------
        sources : `lsst.afw.table.SourceCatalog`
            Detected sources on the difference exposure.
        positiveFootprints : `lsst.afw.detection.FootprintSet`, optional
            Positive polarity footprints.
        negativeFootprints : `lsst.afw.detection.FootprintSet`, optional
            Negative polarity footprints.
        table : `lsst.afw.table.SourceTable`
            Table object that will be used to create the SourceCatalog.
        science : `lsst.afw.image.ExposureF`
            Science exposure that the template was subtracted from.
        matchedTemplate : `lsst.afw.image.ExposureF`
            Warped and PSF-matched template that was used produce the
            difference image.
        difference : `lsst.afw.image.ExposureF`
            Result of subtracting template from the science image.

        Returns
        -------
        measurementResults : `lsst.pipe.base.Struct`

            ``subtractedMeasuredExposure`` : `lsst.afw.image.ExposureF`
                Subtracted exposure with detection mask applied.
            ``diaSources``  : `lsst.afw.table.SourceCatalog`
                The catalog of detected sources.
        """
        if self.config.doMerge:
            fpSet = positiveFootprints
            fpSet.merge(negativeFootprints, self.config.growFootprint,
                        self.config.growFootprint, False)
            diaSources = afwTable.SourceCatalog(table)
            fpSet.makeSources(diaSources)
            self.log.info("Merging detections into %d sources", len(diaSources))
        else:
            diaSources = sources

        if self.config.doSkySources:
            self.addSkySources(diaSources, difference.mask, difference.info.id)

        self.measureDiaSources(diaSources, science, difference, matchedTemplate)

        if self.config.doForcedMeasurement:
            self.measureForcedSources(diaSources, science, difference.getWcs())

        measurementResults = pipeBase.Struct(
            subtractedMeasuredExposure=difference,
            diaSources=diaSources,
        )

        return measurementResults

    def addSkySources(self, diaSources, mask, seed):
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
        skySourceFootprints = self.skySources.run(mask=mask, seed=seed)
        if skySourceFootprints:
            for foot in skySourceFootprints:
                s = diaSources.addNew()
                s.setFootprint(foot)
                s.set(self.skySourceKey, True)

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
        forcedSources = self.forcedMeasurement.generateMeasCat(
            science, diaSources, wcs)
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
            Generator object to assign ids to detected sources in the difference image.

        Returns
        -------
        measurementResults : `lsst.pipe.base.Struct`

            ``subtractedMeasuredExposure`` : `lsst.afw.image.ExposureF`
                Subtracted exposure with detection mask applied.
            ``diaSources``  : `lsst.afw.table.SourceCatalog`
                The catalog of detected sources.
        """
        # Ensure that we start with an empty detection mask.
        mask = scoreExposure.mask
        mask &= ~(mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE"))

        table = afwTable.SourceTable.make(self.schema, idFactory)
        table.setMetadata(self.algMetadata)
        results = self.detection.run(
            table=table,
            exposure=scoreExposure,
            doSmooth=False,
        )
        # Copy the detection mask from the Score image to the difference image
        difference.mask.assign(scoreExposure.mask, scoreExposure.getBBox())

        return self.processResults(science, matchedTemplate, difference, results.sources, table,
                                   positiveFootprints=results.positive, negativeFootprints=results.negative)
