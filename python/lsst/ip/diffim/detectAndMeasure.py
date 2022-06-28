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

import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
from lsst.meas.algorithms import SkyObjectsTask, SourceDetectionTask
from lsst.meas.base import ForcedMeasurementTask, ApplyApCorrTask
import lsst.meas.extensions.trailedSources  # noqa: F401
from lsst.obs.base import ExposureIdInfo
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.utils
from lsst.utils.timer import timeMethod

from . import DipoleFitTask

__all__ = ["DetectAndMeasureConfig", "DetectAndMeasureTask"]


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
    template = pipeBase.connectionTypes.Input(
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
    selectSources = pipeBase.connectionTypes.Input(
        doc="Sources measured on the science exposure; will be matched to the "
            "detected sources on the difference image.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="SourceCatalog",
        name="{fakesType}src",
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

    def setDefaults(self):
        # DiaSource Detection
        self.detection.thresholdPolarity = "both"
        self.detection.thresholdValue = 5.0
        self.detection.reEstimateBackground = False
        self.detection.thresholdType = "pixel_stdev"

        # Add filtered flux measurement, the correct measurement for pre-convolved images.
        self.measurement.algorithms.names.add('base_PeakLikelihoodFlux')
        self.measurement.plugins.names |= ['ext_trailedSources_Naive',
                                           'base_LocalPhotoCalib',
                                           'base_LocalWcs']

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
        expId, expBits = butlerQC.quantum.dataId.pack("visit_detector",
                                                      returnMaxBits=True)
        idFactory = self.makeIdFactory(expId=expId, expBits=expBits)

        outputs = self.run(inputs['science'],
                           inputs['template'],
                           inputs['difference'],
                           inputs['selectSources'],
                           idFactory=idFactory)
        butlerQC.put(outputs, outputRefs)

    @timeMethod
    def run(self, science, template, difference, selectSources,
            idFactory=None):
        """Detect and measure sources on a difference image.

        Parameters
        ----------
        science : `lsst.afw.image.ExposureF`
            Science exposure that the template was subtracted from.
        template : `lsst.afw.image.ExposureF`
            Warped and PSF-matched template that was used produce the
            difference image.
        difference : `lsst.afw.image.ExposureF`
            Result of subtracting template from the science image.
        selectSources : `lsst.afw.table.SourceCatalog`, optional
            Identified sources on the science exposure.
        idFactory : `lsst.afw.table.IdFactory`
            Generator object to assign ids to detected sources in the difference image.

        Returns
        -------
        results : `lsst.pipe.base.Struct`

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

        if self.config.doMerge:
            fpSet = results.fpSets.positive
            fpSet.merge(results.fpSets.negative, self.config.growFootprint,
                        self.config.growFootprint, False)
            diaSources = afwTable.SourceCatalog(table)
            fpSet.makeSources(diaSources)
            self.log.info("Merging detections into %d sources", len(diaSources))
        else:
            diaSources = results.sources

        if self.config.doSkySources:
            self.addSkySources(diaSources, difference.mask, difference.info.id)

        self.measureDiaSources(diaSources, science, difference, template)

        if self.config.doForcedMeasurement:
            self.measureForcedSources(diaSources, science, difference.getWcs())

        return pipeBase.Struct(
            subtractedMeasuredExposure=difference,
            diaSources=diaSources,
        )

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
            self.applyApCorr.run(
                catalog=diaSources,
                apCorrMap=difference.getInfo().getApCorrMap()
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