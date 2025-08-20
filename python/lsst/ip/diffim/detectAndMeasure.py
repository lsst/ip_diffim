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
import requests
import os

import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
import lsst.geom
from lsst.ip.diffim.utils import (evaluateMaskFraction, computeDifferenceImageMetrics,
                                  populate_sattle_visit_cache)
from lsst.meas.algorithms import SkyObjectsTask, SourceDetectionTask, SetPrimaryFlagsTask, MaskStreaksTask
from lsst.meas.algorithms import FindGlintTrailsTask
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


class BadSubtractionError(pipeBase.AlgorithmError):
    """Raised when the residuals in footprints of stars used to compute the
    psf-matching kernel exceeds the configured maximum.
    """
    def __init__(self, *, ratio, threshold):
        msg = ("The ratio of residual power in source footprints on the"
               " difference image to the power in the footprints on the"
               f" science image was {ratio}, which exceeds the maximum"
               f" threshold of {threshold}")
        super().__init__(msg)
        self.ratio = ratio
        self.threshold = threshold

    @property
    def metadata(self):
        return {"ratio": self.ratio,
                "threshold": self.threshold
                }


class NoDiaSourcesError(pipeBase.AlgorithmError):
    """Raised when there are no diaSources detected on an image difference.
    """
    def __init__(self):
        msg = ("No diaSources detected!")
        super().__init__(msg)

    @property
    def metadata(self):
        return {}


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
    kernelSources = pipeBase.connectionTypes.Input(
        doc="Final selection of sources used for psf matching.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="SourceCatalog",
        name="{fakesType}{coaddName}Diff_psfMatchSources"
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
    differenceBackground = pipeBase.connectionTypes.Output(
        doc="Background model that was subtracted from the difference image.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="Background",
        name="difference_background",
    )
    maskedStreaks = pipeBase.connectionTypes.Output(
        doc='Catalog of streak fit parameters for the difference image.',
        storageClass="ArrowNumpyDict",
        dimensions=("instrument", "visit", "detector"),
        name="{fakesType}{coaddName}Diff_streaks",
    )
    glintTrailInfo = pipeBase.connectionTypes.Output(
        doc='Dict of fit parameters for glint trails in the catalog.',
        storageClass="ArrowNumpyDict",
        dimensions=("instrument", "visit", "detector"),
        name="trailed_glints",
    )

    def __init__(self, *, config):
        super().__init__(config=config)
        if not (self.config.writeStreakInfo and self.config.doMaskStreaks):
            self.outputs.remove("maskedStreaks")
        if not (self.config.doSubtractBackground and self.config.doWriteBackground):
            self.outputs.remove("differenceBackground")
        if not (self.config.writeGlintInfo):
            self.outputs.remove("glintTrailInfo")


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
        doc="Force photometer diaSource locations on PVI?"
    )
    doAddMetrics = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Add columns to the source table to hold analysis metrics?"
    )
    doSubtractBackground = pexConfig.Field(
        dtype=bool,
        doc="Subtract a background model from the image before detection?",
        default=True,
    )
    doWriteBackground = pexConfig.Field(
        dtype=bool,
        doc="Persist the fitted background model?",
        default=False,
    )
    subtractInitialBackground = pexConfig.ConfigurableField(
        target=lsst.meas.algorithms.SubtractBackgroundTask,
        doc="Task to perform intial background subtraction, before first detection pass.",
    )
    subtractFinalBackground = pexConfig.ConfigurableField(
        target=lsst.meas.algorithms.SubtractBackgroundTask,
        doc="Task to perform final background subtraction, after first detection pass.",
    )
    detection = pexConfig.ConfigurableField(
        target=SourceDetectionTask,
        doc="Final source detection for diaSource measurement",
    )
    streakDetection = pexConfig.ConfigurableField(
        target=SourceDetectionTask,
        doc="Separate source detection used only for streak masking",
    )
    doDeblend = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Deblend DIASources after detection?"
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
    doMaskStreaks = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Turn on streak masking",
    )
    maskStreaks = pexConfig.ConfigurableField(
        target=MaskStreaksTask,
        doc="Subtask for masking streaks. Only used if doMaskStreaks is True. "
            "Adds a mask plane to an exposure, with the mask plane name set by streakMaskName.",
    )
    streakBinFactor = pexConfig.Field(
        dtype=int,
        default=4,
        doc="Bin scale factor to use when rerunning detection for masking streaks. "
            "Only used if doMaskStreaks is True.",
    )
    writeStreakInfo = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Record the parameters of any detected streaks. For LSST, this should be turned off except for "
            "development work."
    )
    findGlints = pexConfig.ConfigurableField(
        target=FindGlintTrailsTask,
        doc="Subtask for finding glint trails, usually caused by satellites or debris."
    )
    writeGlintInfo = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Record the parameters of any detected glint trails."
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
                 "base_PixelFlags_flag_nodataCenterAll",
                 "base_PixelFlags_flag_saturatedCenterAll",
                 ),
    )
    clearMaskPlanes = lsst.pex.config.ListField(
        dtype=str,
        doc="Mask planes to clear before running detection.",
        default=("DETECTED", "DETECTED_NEGATIVE", "NOT_DEBLENDED", "STREAK"),
    )
    raiseOnBadSubtractionRatio = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Raise an error if the ratio of power in detected footprints"
            " on the difference image to the power in footprints on the science"
            " image exceeds ``badSubtractionRatioThreshold``",
    )
    badSubtractionRatioThreshold = pexConfig.Field(
        dtype=float,
        default=0.2,
        doc="Maximum ratio of power in footprints on the difference image to"
            " the same footprints on the science image."
            "Only used if ``raiseOnBadSubtractionRatio`` is set",
    )
    badSubtractionVariationThreshold = pexConfig.Field(
        dtype=float,
        default=0.4,
        doc="Maximum standard deviation of the ratio of power in footprints on"
            " the difference image to the same footprints on the science image."
            "Only used if ``raiseOnBadSubtractionRatio`` is set",
    )
    raiseOnNoDiaSources = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Raise an algorithm error if no diaSources are detected.",
    )
    run_sattle = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="If true, dia source bounding boxes will be sent for verification"
            "to the sattle service."
    )
    sattle_historical = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="If re-running a pipeline that requires sattle, this should be set "
            "to True. This will populate sattle's cache with the historic data "
            "closest in time to the exposure."
    )
    idGenerator = DetectorVisitIdGeneratorConfig.make_field()

    def setDefaults(self):
        # Background subtraction
        # Use a small binsize for the first pass to reduce detections on glints
        #  and extended structures. Should not affect the detectability of
        #  faint diaSources
        self.subtractInitialBackground.binSize = 8
        self.subtractInitialBackground.useApprox = False
        self.subtractInitialBackground.statisticsProperty = "MEDIAN"
        self.subtractInitialBackground.doFilterSuperPixels = True
        self.subtractInitialBackground.ignoredPixelMask = ["BAD",
                                                           "EDGE",
                                                           "DETECTED",
                                                           "DETECTED_NEGATIVE",
                                                           "NO_DATA",
                                                           ]
        # Use a larger binsize for the final background subtraction, to reduce
        #  over-subtraction of bright objects.
        self.subtractFinalBackground.binSize = 40
        self.subtractFinalBackground.useApprox = False
        self.subtractFinalBackground.statisticsProperty = "MEDIAN"
        self.subtractFinalBackground.doFilterSuperPixels = True
        self.subtractFinalBackground.ignoredPixelMask = ["BAD",
                                                         "EDGE",
                                                         "DETECTED",
                                                         "DETECTED_NEGATIVE",
                                                         "NO_DATA",
                                                         ]
        # DiaSource Detection
        self.detection.thresholdPolarity = "both"
        self.detection.thresholdValue = 5.0
        self.detection.reEstimateBackground = False
        self.detection.thresholdType = "pixel_stdev"
        self.detection.excludeMaskPlanes = ["EDGE",
                                            "BAD",
                                            ]

        # Copy configs for binned streak detection from the base detection task
        self.streakDetection.thresholdType = self.detection.thresholdType
        self.streakDetection.reEstimateBackground = False
        self.streakDetection.excludeMaskPlanes = self.detection.excludeMaskPlanes
        self.streakDetection.thresholdValue = self.detection.thresholdValue
        # Only detect positive streaks
        self.streakDetection.thresholdPolarity = "positive"
        # Do not grow detected mask for streaks
        self.streakDetection.nSigmaToGrow = 0
        # Set the streak mask along the entire fit line, not only where the
        # detected mask is set.
        self.maskStreaks.onlyMaskDetected = False
        # Restrict streak masking from growing too large
        self.maskStreaks.maxStreakWidth = 100
        # Restrict the number of iterations allowed for fitting streaks
        # When the fit is good it should solve quickly, and exit a bad fit quickly
        self.maskStreaks.maxFitIter = 10
        # Only mask to 2 sigma in width
        self.maskStreaks.nSigmaMask = 2
        # Threshold for including streaks after the Hough Transform.
        # A lower value will detect more features that are less linear.
        self.maskStreaks.absMinimumKernelHeight = 2

        self.measurement.plugins.names |= ["ext_trailedSources_Naive",
                                           "base_LocalPhotoCalib",
                                           "base_LocalWcs",
                                           "ext_shapeHSM_HsmSourceMoments",
                                           "ext_shapeHSM_HsmPsfMoments",
                                           "base_ClassificationSizeExtendedness",
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

    def validate(self):
        super().validate()

        if self.run_sattle:
            if not os.getenv("SATTLE_URI_BASE"):
                raise pexConfig.FieldValidationError(DetectAndMeasureConfig.run_sattle, self,
                                                     "Sattle requested but SATTLE_URI_BASE "
                                                     "environment variable not set.")


class DetectAndMeasureTask(lsst.pipe.base.PipelineTask):
    """Detect and measure sources on a difference image.
    """
    ConfigClass = DetectAndMeasureConfig
    _DefaultName = "detectAndMeasure"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.schema = afwTable.SourceTable.makeMinimalSchema()

        self.algMetadata = dafBase.PropertyList()
        if self.config.doSubtractBackground:
            self.makeSubtask("subtractInitialBackground")
            self.makeSubtask("subtractFinalBackground")
        self.makeSubtask("detection", schema=self.schema)
        if self.config.doDeblend:
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
        # Create the sky source task for use by metrics,
        # even if sky sources are not added to the diaSource catalog
        self.makeSubtask("skySources", schema=self.schema)
        if self.config.doMaskStreaks:
            self.makeSubtask("maskStreaks")
            self.makeSubtask("streakDetection")
        self.makeSubtask("findGlints")
        self.schema.addField("glint_trail", "Flag", "DiaSource is part of a glint trail.")
        self.schema.addField("reliability", type="F", doc="Reliability score of the DiaSource")

        # To get the "merge_*" fields in the schema; have to re-initialize
        # this later, once we have a peak schema post-detection.
        lsst.afw.detection.FootprintMergeList(self.schema, ["positive", "negative"])

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
        # Specify the fields that `annotate` needs below, to ensure they
        # exist, even as None.
        measurementResults = pipeBase.Struct(
            subtractedMeasuredExposure=None,
            diaSources=None,
            maskedStreaks=None,
            differenceBackground=None,
        )
        try:
            self.run(**inputs, idFactory=idFactory, measurementResults=measurementResults)
        except pipeBase.AlgorithmError as e:
            error = pipeBase.AnnotatedPartialOutputsError.annotate(
                e,
                self,
                measurementResults.subtractedMeasuredExposure,
                measurementResults.diaSources,
                measurementResults.maskedStreaks,
                measurementResults.glintTrailInfo,
                log=self.log
            )
            butlerQC.put(measurementResults, outputRefs)
            raise error from e
        butlerQC.put(measurementResults, outputRefs)

    @timeMethod
    def run(self, science, matchedTemplate, difference, kernelSources,
            idFactory=None, measurementResults=None):
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
        kernelSources : `lsst.afw.table.SourceCatalog`
            Final selection of sources that was used for psf matching.
        idFactory : `lsst.afw.table.IdFactory`, optional
            Generator object used to assign ids to detected sources in the
            difference image. Ids from this generator are not set until after
            deblending and merging positive/negative peaks.
        measurementResults : `lsst.pipe.base.Struct`, optional
            Result struct that is modified to allow saving of partial outputs
            for some failure conditions. If the task completes successfully,
            this is also returned.

        Returns
        -------
        measurementResults : `lsst.pipe.base.Struct`

            ``subtractedMeasuredExposure`` : `lsst.afw.image.ExposureF`
                Subtracted exposure with detection mask applied.
            ``diaSources``  : `lsst.afw.table.SourceCatalog`
                The catalog of detected sources.
            ``differenceBackground`` : `lsst.afw.math.BackgroundList`
                Background that was subtracted from the difference image.
        """
        if measurementResults is None:
            measurementResults = pipeBase.Struct()
        if idFactory is None:
            idFactory = lsst.meas.base.IdGenerator().make_table_id_factory()

        if self.config.doSubtractBackground:
            # Run background subtraction before clearing the mask planes
            detectionExposure = difference.clone()
            background = self.subtractInitialBackground.run(detectionExposure).background
        else:
            detectionExposure = difference
            background = afwMath.BackgroundList()

        self._prepareInputs(detectionExposure)

        # Don't use the idFactory until after deblend+merge, so that we aren't
        # generating ids that just get thrown away (footprint merge doesn't
        # know about past ids).
        table = afwTable.SourceTable.make(self.schema)
        results = self.detection.run(
            table=table,
            exposure=detectionExposure,
            doSmooth=True,
            background=background
        )

        if self.config.doSubtractBackground:
            # Run background subtraction again after detecting peaks
            # but before measurement
            # First update the mask using the detection image
            difference.setMask(detectionExposure.mask)
            background = self.subtractFinalBackground.run(difference).background

            # Re-run detection to get final footprints
            table = afwTable.SourceTable.make(self.schema)
            results = self.detection.run(
                table=table,
                exposure=difference,
                doSmooth=True,
                background=background
            )
        measurementResults.differenceBackground = background

        if self.config.doDeblend:
            sources, positives, negatives = self._deblend(difference,
                                                          results.positive,
                                                          results.negative)

        else:
            positives = afwTable.SourceCatalog(self.schema)
            results.positive.makeSources(positives)
            negatives = afwTable.SourceCatalog(self.schema)
            results.negative.makeSources(negatives)
            sources = results.sources

        self.processResults(science, matchedTemplate, difference,
                            sources, idFactory, kernelSources,
                            positives=positives,
                            negatives=negatives,
                            measurementResults=measurementResults)
        return measurementResults

    def _prepareInputs(self, difference):
        """Ensure that we start with an empty detection and deblended mask.

        Parameters
        ----------
        difference : `lsst.afw.image.ExposureF`
            The difference image that will be used for detecting diaSources.
            The mask plane will be modified in place.

        Raises
        ------
        lsst.pipe.base.UpstreamFailureNoWorkFound
            If the PSF is not usable for measurement.
        """
        # Check that we have a valid PSF now before we do more work
        sigma = difference.psf.computeShape(difference.psf.getAveragePosition()).getDeterminantRadius()
        if np.isnan(sigma):
            raise pipeBase.UpstreamFailureNoWorkFound("Invalid PSF detected! PSF width evaluates to NaN.")
        # Ensure that we start with an empty detection and deblended mask.
        mask = difference.mask
        for mp in self.config.clearMaskPlanes:
            if mp not in mask.getMaskPlaneDict():
                mask.addMaskPlane(mp)
        mask &= ~mask.getPlaneBitMask(self.config.clearMaskPlanes)

    def processResults(self, science, matchedTemplate, difference, sources, idFactory, kernelSources,
                       positives=None, negatives=None, measurementResults=None):
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
        sources : `lsst.afw.table.SourceCatalog`
            Detected sources on the difference exposure.
        idFactory : `lsst.afw.table.IdFactory`
            Generator object used to assign ids to detected sources in the
            difference image.
        kernelSources : `lsst.afw.table.SourceCatalog`
            Final selection of sources that was used for psf matching.
        positives : `lsst.afw.table.SourceCatalog`, optional
            Positive polarity footprints.
        negatives : `lsst.afw.table.SourceCatalog`, optional
            Negative polarity footprints.
        measurementResults : `lsst.pipe.base.Struct`, optional
            Result struct that is modified to allow saving of partial outputs
            for some failure conditions. If the task completes successfully,
            this is also returned.

        Returns
        -------
        measurementResults : `lsst.pipe.base.Struct`

            ``subtractedMeasuredExposure`` : `lsst.afw.image.ExposureF`
                Subtracted exposure with detection mask applied.
            ``diaSources``  : `lsst.afw.table.SourceCatalog`
                The catalog of detected sources.
        """
        if measurementResults is None:
            measurementResults = pipeBase.Struct()
        self.metadata["nUnmergedDiaSources"] = len(sources)
        if self.config.doMerge:
            # preserve peak schema, if there are any footprints
            if len(positives) > 0:
                peakSchema = positives[0].getFootprint().peaks.schema
            elif len(negatives) > 0:
                peakSchema = negatives[0].getFootprint().peaks.schema
            else:
                peakSchema = afwDetection.PeakTable.makeMinimalSchema()
            mergeList = afwDetection.FootprintMergeList(self.schema,
                                                        ["positive", "negative"], peakSchema)
            initialDiaSources = afwTable.SourceCatalog(self.schema)
            # Start with positive, as FootprintMergeList will self-merge the
            # subsequent added catalogs, and we want to try to preserve
            # deblended positive sources.
            mergeList.addCatalog(initialDiaSources.table, positives, "positive", minNewPeakDist=0)
            mergeList.addCatalog(initialDiaSources.table, negatives, "negative", minNewPeakDist=0)
            mergeList.getFinalSources(initialDiaSources)
            # Flag as negative those sources that *only* came from the negative
            # footprint set.
            initialDiaSources["is_negative"] = initialDiaSources["merge_footprint_negative"] & \
                ~initialDiaSources["merge_footprint_positive"]
            self.log.info("Merging detections into %d sources", len(initialDiaSources))
        else:
            initialDiaSources = sources

        # Assign source ids at the end: deblend/merge mean that we don't keep
        # track of parents and children, we only care about the final ids.
        for source in initialDiaSources:
            source.setId(idFactory())
        # Ensure sources added after this get correct ids.
        initialDiaSources.getTable().setIdFactory(idFactory)
        initialDiaSources.setMetadata(self.algMetadata)

        self.metadata["nMergedDiaSources"] = len(initialDiaSources)

        if self.config.doMaskStreaks:
            streakInfo = self._runStreakMasking(difference)

        if self.config.doSkySources:
            self.addSkySources(initialDiaSources, difference.mask, difference.info.id)

        if not initialDiaSources.isContiguous():
            initialDiaSources = initialDiaSources.copy(deep=True)

        self.measureDiaSources(initialDiaSources, science, difference, matchedTemplate)

        # Remove unphysical diaSources per config.badSourceFlags
        diaSources = self._removeBadSources(initialDiaSources)

        if self.config.run_sattle:
            diaSources = self.filterSatellites(diaSources, science)

        # Flag diaSources in glint trails, but do not remove them
        diaSources, trail_parameters = self._find_glint_trails(diaSources)
        if self.config.writeGlintInfo:
            measurementResults.mergeItems(trail_parameters, 'glintTrailInfo')

        if self.config.doForcedMeasurement:
            self.measureForcedSources(diaSources, science, difference.getWcs())

        # Clear the image plane for regions with NO_DATA.
        # These regions are most often caused by insufficient template coverage.
        # Do this for the final difference image after detection and measurement
        # since the subtasks should all be configured to handle NO_DATA properly
        difference.image.array[difference.mask.array & difference.mask.getPlaneBitMask('NO_DATA') > 0] = 0

        measurementResults.subtractedMeasuredExposure = difference

        if self.config.doMaskStreaks and self.config.writeStreakInfo:
            measurementResults.mergeItems(streakInfo, 'maskedStreaks')

        self.calculateMetrics(science, difference, diaSources, kernelSources)

        if np.count_nonzero(~diaSources["sky_source"]) > 0:
            measurementResults.diaSources = diaSources
        elif self.config.raiseOnNoDiaSources:
            raise NoDiaSourcesError()
        elif len(diaSources) > 0:
            # This option allows returning sky sources,
            # even if there are no diaSources
            measurementResults.diaSources = diaSources

        return measurementResults

    def _deblend(self, difference, positiveFootprints, negativeFootprints):
        """Deblend the positive and negative footprints and return a catalog
        containing just the children, and the deblended footprints.

        Parameters
        ----------
        difference : `lsst.afw.image.Exposure`
            Result of subtracting template from the science image.
        positiveFootprints, negativeFootprints : `lsst.afw.detection.FootprintSet`
            Positive and negative polarity footprints measured on
            ``difference`` to be deblended separately.

        Returns
        -------
        sources : `lsst.afw.table.SourceCatalog`
            Positive and negative deblended children.
        positives, negatives : `lsst.afw.table.SourceCatalog`
            Deblended positive and negative polarity sources with footprints
            detected on ``difference``.
        """

        def deblend(footprints, negative=False):
            """Deblend a positive or negative footprint set,
            and return the deblended children.

            Parameters
            ----------
            footprints : `lsst.afw.detection.FootprintSet`
            negative : `bool`
                Set True if the footprints contain negative fluxes

            Returns
            -------
            sources : `lsst.afw.table.SourceCatalog`
            """
            sources = afwTable.SourceCatalog(self.schema)
            footprints.makeSources(sources)
            if negative:
                # Invert the image so the deblender can run on positive peaks
                difference_inverted = difference.clone()
                difference_inverted.image *= -1
                self.deblend.run(exposure=difference_inverted, sources=sources)
                children = sources[sources["parent"] != 0]
                # Set the heavy footprint pixel values back to reality
                for child in children:
                    footprint = child.getFootprint()
                    array = footprint.getImageArray()
                    array *= -1
            else:
                self.deblend.run(exposure=difference, sources=sources)
            self.setPrimaryFlags.run(sources)
            children = sources["detect_isDeblendedSource"] == 1
            sources = sources[children].copy(deep=True)
            # Clear parents, so that measurement plugins behave correctly.
            sources['parent'] = 0
            return sources.copy(deep=True)

        positives = deblend(positiveFootprints)
        negatives = deblend(negativeFootprints, negative=True)

        sources = afwTable.SourceCatalog(self.schema)
        sources.reserve(len(positives) + len(negatives))
        sources.extend(positives, deep=True)
        sources.extend(negatives, deep=True)
        if len(negatives) > 0:
            sources[-len(negatives):]["is_negative"] = True
        return sources, positives, negatives

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
        self.metadata["nRemovedBadFlaggedSources"] = nBadTotal
        self.log.info("Removed %d unphysical sources.", nBadTotal)
        return diaSources[selector].copy(deep=True)

    def _find_glint_trails(self, diaSources):
        """Define a new flag column for diaSources that are in a glint trail.

        Parameters
        ----------
        diaSources : `lsst.afw.table.SourceCatalog`
            The catalog of detected sources.

        Returns
        -------
        diaSources : `lsst.afw.table.SourceCatalog`
            The updated catalog of detected sources, with a new bool column
            called 'glint_trail' added.

        trail_parameters : `dict`
            Parameters of all the trails that were found.
        """
        if self.config.doSkySources:
            # Do not include sky sources in glint detection
            candidateDiaSources = diaSources[~diaSources["sky_source"]].copy(deep=True)
        else:
            candidateDiaSources = diaSources
        trailed_glints = self.findGlints.run(candidateDiaSources)
        glint_mask = [True if id in trailed_glints.trailed_ids else False for id in diaSources['id']]
        diaSources['glint_trail'] = np.array(glint_mask)

        slopes = np.array([trail.slope for trail in trailed_glints.parameters])
        intercepts = np.array([trail.intercept for trail in trailed_glints.parameters])
        stderrs = np.array([trail.stderr for trail in trailed_glints.parameters])
        lengths = np.array([trail.length for trail in trailed_glints.parameters])
        angles = np.array([trail.angle for trail in trailed_glints.parameters])
        parameters = {'slopes': slopes, 'intercepts': intercepts, 'stderrs': stderrs, 'lengths': lengths,
                      'angles': angles}

        trail_parameters = pipeBase.Struct(glintTrailInfo=parameters)

        return diaSources, trail_parameters

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
        self.metadata[f"n_{subtask.getName()}"] = len(skySourceFootprints)

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

    def calculateMetrics(self, science, difference, diaSources, kernelSources):
        """Add difference image QA metrics to the Task metadata.

        This may be used to produce corresponding metrics (see
        lsst.analysis.tools.tasks.diffimTaskDetectorVisitMetricAnalysis).

        Parameters
        ----------
        science : `lsst.afw.image.ExposureF`
            Science exposure that was subtracted.
        difference : `lsst.afw.image.Exposure`
            The target difference image to calculate metrics for.
        diaSources : `lsst.afw.table.SourceCatalog`
            The catalog of detected sources.
        kernelSources : `lsst.afw.table.SourceCatalog`
            Final selection of sources that was used for psf matching.
        """
        mask = difference.mask
        badPix = (mask.array & mask.getPlaneBitMask(self.config.detection.excludeMaskPlanes)) > 0
        self.metadata["nGoodPixels"] = np.sum(~badPix)
        self.metadata["nBadPixels"] = np.sum(badPix)
        detPosPix = (mask.array & mask.getPlaneBitMask("DETECTED")) > 0
        detNegPix = (mask.array & mask.getPlaneBitMask("DETECTED_NEGATIVE")) > 0
        self.metadata["nPixelsDetectedPositive"] = np.sum(detPosPix)
        self.metadata["nPixelsDetectedNegative"] = np.sum(detNegPix)
        detPosPix &= badPix
        detNegPix &= badPix
        self.metadata["nBadPixelsDetectedPositive"] = np.sum(detPosPix)
        self.metadata["nBadPixelsDetectedNegative"] = np.sum(detNegPix)

        metricsMaskPlanes = list(mask.getMaskPlaneDict().keys())
        for maskPlane in metricsMaskPlanes:
            try:
                self.metadata["%s_mask_fraction"%maskPlane.lower()] = evaluateMaskFraction(mask, maskPlane)
            except InvalidParameterError:
                self.metadata["%s_mask_fraction"%maskPlane.lower()] = -1
                self.log.info("Unable to calculate metrics for mask plane %s: not in image"%maskPlane)

        if self.config.doSkySources:
            skySources = diaSources[diaSources["sky_source"]]
        else:
            skySources = None
        metrics = computeDifferenceImageMetrics(science, difference, kernelSources, sky_sources=skySources)

        self.metadata["residualFootprintRatioMean"] = metrics.differenceFootprintRatioMean
        self.metadata["residualFootprintRatioStdev"] = metrics.differenceFootprintRatioStdev
        self.metadata["differenceFootprintSkyRatioMean"] = metrics.differenceFootprintSkyRatioMean
        self.metadata["differenceFootprintSkyRatioStdev"] = metrics.differenceFootprintSkyRatioStdev
        self.log.info("Mean, stdev of ratio of difference to science "
                      "pixels in star footprints: %5.4f, %5.4f",
                      self.metadata["residualFootprintRatioMean"],
                      self.metadata["residualFootprintRatioStdev"])
        if self.config.raiseOnBadSubtractionRatio:
            if metrics.differenceFootprintRatioMean > self.config.badSubtractionRatioThreshold:
                raise BadSubtractionError(ratio=metrics.differenceFootprintRatioMean,
                                          threshold=self.config.badSubtractionRatioThreshold)
            if metrics.differenceFootprintRatioStdev > self.config.badSubtractionVariationThreshold:
                raise BadSubtractionError(ratio=metrics.differenceFootprintRatioStdev,
                                          threshold=self.config.badSubtractionVariationThreshold)

    def getSattleDiaSourceAllowlist(self, diaSources, science):
        """Query the sattle service and determine which diaSources are allowed.

        Parameters
        ----------
        diaSources : `lsst.afw.table.SourceCatalog`
            The catalog of detected sources.
        science : `lsst.afw.image.ExposureF`
            Science exposure that was subtracted.

        Returns
        ----------
        allow_list : `list` of `int`
            diaSourceIds of diaSources that can be made public.

        Raises
        ------
        requests.HTTPError
            Raised if sattle call does not return success.
        """
        wcs = science.getWcs()
        visit_info = science.getInfo().getVisitInfo()
        visit_id = visit_info.getId()
        sattle_uri_base = os.getenv('SATTLE_URI_BASE')

        dia_sources_json = []
        for source in diaSources:
            source_bbox = source.getFootprint().getBBox()
            corners = wcs.pixelToSky([lsst.geom.Point2D(c) for c in source_bbox.getCorners()])
            bbox_radec = [[pt.getRa().asDegrees(), pt.getDec().asDegrees()] for pt in corners]
            dia_sources_json.append({"diasource_id": source["id"], "bbox": bbox_radec})

        payload = {"visit_id": visit_id, "detector_id": science.getDetector().getId(),
                   "diasources": dia_sources_json, "historical": self.config.sattle_historical}

        sattle_output = requests.put(f'{sattle_uri_base}/diasource_allow_list',
                                     json=payload)

        # retry once if visit cache is not populated
        if sattle_output.status_code == 404:
            self.log.warning(f'Visit {visit_id} not found in sattle cache, re-sending')
            populate_sattle_visit_cache(visit_info, historical=self.config.sattle_historical)
            sattle_output = requests.put(f'{sattle_uri_base}/diasource_allow_list', json=payload)

        sattle_output.raise_for_status()

        return sattle_output.json()['allow_list']

    def filterSatellites(self, diaSources, science):
        """Remove diaSources overlapping predicted satellite positions.

        Parameters
        ----------
        diaSources : `lsst.afw.table.SourceCatalog`
            The catalog of detected sources.
        science : `lsst.afw.image.ExposureF`
            Science exposure that was subtracted.

        Returns
        ----------
        filterdDiaSources : `lsst.afw.table.SourceCatalog`
            Filtered catalog of diaSources
        """

        allow_list = self.getSattleDiaSourceAllowlist(diaSources, science)

        if allow_list:
            allow_set = set(allow_list)
            allowed_ids = [source['id'] in allow_set for source in diaSources]
            diaSources = diaSources[np.array(allowed_ids)].copy(deep=True)
        else:
            self.log.warning('Sattle allowlist is empty, all diaSources removed')
            diaSources = diaSources[0:0].copy(deep=True)
        return diaSources

    def _runStreakMasking(self, difference):
        """Do streak masking and optionally save the resulting streak
        fit parameters in a catalog.

        Only returns non-empty streakInfo if self.config.writeStreakInfo
        is set. The difference image is binned by self.config.streakBinFactor
        (and detection is run a second time) so that regions with lower
        surface brightness streaks are more likely to fall above the
        detection threshold.

        Parameters
        ----------
        difference: `lsst.afw.image.Exposure`
            The exposure in which to search for streaks. Must have a detection
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
            ``reducedChi2`` : `np.ndarray`
                Reduced chi2 of the best-fit streak profile.
            ``modelMaximum`` : `np.ndarray`
                Peak value of the fit line profile.
        """
        maskedImage = difference.maskedImage
        # Bin the diffim to enhance low surface brightness streaks
        binnedMaskedImage = afwMath.binImage(maskedImage,
                                             self.config.streakBinFactor,
                                             self.config.streakBinFactor)
        binnedExposure = afwImage.ExposureF(binnedMaskedImage.getBBox())
        binnedExposure.setMaskedImage(binnedMaskedImage)
        # Clear the DETECTED mask plane before streak detection
        binnedExposure.mask &= ~binnedExposure.mask.getPlaneBitMask('DETECTED')
        # Rerun detection to set the DETECTED mask plane on binnedExposure
        sigma = difference.psf.computeShape(difference.psf.getAveragePosition()).getDeterminantRadius()
        _table = afwTable.SourceTable.make(afwTable.SourceTable.makeMinimalSchema())
        self.streakDetection.run(table=_table, exposure=binnedExposure, doSmooth=True,
                                 sigma=sigma/self.config.streakBinFactor)
        binnedDetectedMaskPlane = binnedExposure.mask.array & binnedExposure.mask.getPlaneBitMask('DETECTED')
        rescaledDetectedMaskPlane = binnedDetectedMaskPlane.repeat(self.config.streakBinFactor,
                                                                   axis=0).repeat(self.config.streakBinFactor,
                                                                                  axis=1)
        # Create new version of a diffim with DETECTED based on binnedExposure
        streakMaskedImage = maskedImage.clone()
        ysize, xsize = rescaledDetectedMaskPlane.shape
        streakMaskedImage.mask.array[:ysize, :xsize] |= rescaledDetectedMaskPlane
        # Detect streaks on this new version of the diffim
        streaks = self.maskStreaks.run(streakMaskedImage)
        streakMaskPlane = streakMaskedImage.mask.array & streakMaskedImage.mask.getPlaneBitMask('STREAK')
        # Apply the new STREAK mask to the original diffim
        maskedImage.mask.array |= streakMaskPlane

        if self.config.writeStreakInfo:
            rhos = np.array([line.rho for line in streaks.lines])
            thetas = np.array([line.theta for line in streaks.lines])
            sigmas = np.array([line.sigma for line in streaks.lines])
            chi2s = np.array([line.reducedChi2 for line in streaks.lines])
            modelMaximums = np.array([line.modelMaximum for line in streaks.lines])
            streakInfo = {'rho': rhos, 'theta': thetas, 'sigma': sigmas, 'reducedChi2': chi2s,
                          'modelMaximum': modelMaximums}
        else:
            streakInfo = {'rho': np.array([]), 'theta': np.array([]), 'sigma': np.array([]),
                          'reducedChi2': np.array([]), 'modelMaximum': np.array([])}
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
    def run(self, science, matchedTemplate, difference, scoreExposure, kernelSources,
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
        kernelSources : `lsst.afw.table.SourceCatalog`
            Final selection of sources that was used for psf matching.
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

        self._prepareInputs(scoreExposure)

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

        if self.config.doDeblend:
            sources, positives, negatives = self._deblend(difference,
                                                          results.positive,
                                                          results.negative)

            return self.processResults(science, matchedTemplate, difference,
                                       sources, idFactory, kernelSources,
                                       positives=positives,
                                       negatives=negatives)

        else:
            positives = afwTable.SourceCatalog(self.schema)
            results.positive.makeSources(positives)
            negatives = afwTable.SourceCatalog(self.schema)
            results.negative.makeSources(negatives)
            return self.processResults(science, matchedTemplate, difference,
                                       results.sources, idFactory, kernelSources,
                                       positives=positives,
                                       negatives=negatives)
