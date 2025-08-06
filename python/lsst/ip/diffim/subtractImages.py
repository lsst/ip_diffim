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

from astropy import units as u
import numpy as np

import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.image
from lsst.afw.math._warper import computeWarpedBBox
import lsst.afw.math
import lsst.geom
from lsst.ip.diffim.utils import (evaluateMeanPsfFwhm, getPsfFwhm, computeDifferenceImageMetrics,
                                  divideExposureByPatches,
                                  )
from lsst.meas.algorithms import ScaleVarianceTask, ScienceSourceSelectorTask

import lsst.pex.config
import lsst.pipe.base
import lsst.pex.exceptions
from lsst.pipe.base import connectionTypes
from lsst.skymap import BaseSkyMap
from . import MakeKernelTask, DecorrelateALKernelTask
from lsst.utils.timer import timeMethod

__all__ = ["AlardLuptonSubtractConfig", "AlardLuptonSubtractTask",
           "AlardLuptonPreconvolveSubtractConfig", "AlardLuptonPreconvolveSubtractTask"]

_dimensions = ("instrument", "visit", "detector")
_defaultTemplates = {"coaddName": "deep", "fakesType": ""}


class SubtractInputConnections(lsst.pipe.base.PipelineTaskConnections,
                               dimensions=_dimensions,
                               defaultTemplates=_defaultTemplates):
    template = connectionTypes.Input(
        doc="Input warped template to subtract.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Diff_templateExp"
    )
    science = connectionTypes.Input(
        doc="Input science exposure to subtract from.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}calexp"
    )
    sources = connectionTypes.Input(
        doc="Sources measured on the science exposure; "
            "used to select sources for making the matching kernel.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="SourceCatalog",
        name="{fakesType}src"
    )
    skymap = connectionTypes.Input(
        doc="Input definition of geometry/bbox and projection/wcs for template exposures",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        dimensions=("skymap", ),
        storageClass="SkyMap",
        minimum=0,
    )
    visitSummary = connectionTypes.Input(
        doc=("Per-visit catalog with final calibration objects. "
             "These catalogs use the detector id for the catalog id, "
             "sorted on id for fast lookup."),
        dimensions=("instrument", "visit"),
        storageClass="ExposureCatalog",
        name="finalVisitSummary",
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        if not config.doApplyExternalCalibrations:
            del self.visitSummary


class SubtractImageOutputConnections(lsst.pipe.base.PipelineTaskConnections,
                                     dimensions=_dimensions,
                                     defaultTemplates=_defaultTemplates):
    difference = connectionTypes.Output(
        doc="Result of subtracting convolved template from science image.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Diff_differenceTempExp",
    )
    matchedTemplate = connectionTypes.Output(
        doc="Warped and PSF-matched template used to create `subtractedExposure`.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Diff_matchedExp",
    )
    psfMatchingKernel = connectionTypes.Output(
        doc="Kernel used to PSF match the science and template images.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="MatchingKernel",
        name="{fakesType}{coaddName}Diff_psfMatchKernel",
    )
    kernelSources = connectionTypes.Output(
        doc="Final selection of sources used for psf matching.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="SourceCatalog",
        name="{fakesType}{coaddName}Diff_psfMatchSources"
    )


class SubtractScoreOutputConnections(lsst.pipe.base.PipelineTaskConnections,
                                     dimensions=_dimensions,
                                     defaultTemplates=_defaultTemplates):
    scoreExposure = connectionTypes.Output(
        doc="The maximum likelihood image, used for the detection of diaSources.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Diff_scoreExp",
    )
    psfMatchingKernel = connectionTypes.Output(
        doc="Kernel used to PSF match the science and template images.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="MatchingKernel",
        name="{fakesType}{coaddName}Diff_psfScoreMatchKernel",
    )
    kernelSources = connectionTypes.Output(
        doc="Final selection of sources used for psf matching.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="SourceCatalog",
        name="{fakesType}{coaddName}Diff_psfScoreMatchSources"
    )


class AlardLuptonSubtractConnections(SubtractInputConnections, SubtractImageOutputConnections):
    pass


class AlardLuptonSubtractBaseConfig(lsst.pex.config.Config):
    makeKernel = lsst.pex.config.ConfigurableField(
        target=MakeKernelTask,
        doc="Task to construct a matching kernel for convolution.",
    )
    doDecorrelation = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        doc="Perform diffim decorrelation to undo pixel correlation due to A&L "
        "kernel convolution? If True, also update the diffim PSF."
    )
    decorrelate = lsst.pex.config.ConfigurableField(
        target=DecorrelateALKernelTask,
        doc="Task to decorrelate the image difference.",
    )
    requiredTemplateFraction = lsst.pex.config.Field(
        dtype=float,
        default=0.1,
        doc="Raise NoWorkFound and do not attempt image subtraction if template covers less than this "
        " fraction of pixels. Setting to 0 will always attempt image subtraction."
    )
    minTemplateFractionForExpectedSuccess = lsst.pex.config.Field(
        dtype=float,
        default=0.2,
        doc="Raise NoWorkFound if PSF-matching fails and template covers less than this fraction of pixels."
        " If the fraction of pixels covered by the template is less than this value (and greater than"
        " requiredTemplateFraction) this task is attempted but failure is anticipated and tolerated."
    )
    doScaleVariance = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        doc="Scale variance of the image difference?"
    )
    scaleVariance = lsst.pex.config.ConfigurableField(
        target=ScaleVarianceTask,
        doc="Subtask to rescale the variance of the template to the statistically expected level."
    )
    doSubtractBackground = lsst.pex.config.Field(
        doc="Subtract the background fit when solving the kernel?",
        dtype=bool,
        default=False,
    )
    doApplyExternalCalibrations = lsst.pex.config.Field(
        doc=(
            "Replace science Exposure's calibration objects with those"
            " in visitSummary.  Ignored if `doApplyFinalizedPsf is True."
        ),
        dtype=bool,
        default=False,
    )
    sourceSelector = lsst.pex.config.ConfigurableField(
        target=ScienceSourceSelectorTask,
        doc="Task to select sources to be used for PSF matching.",
    )
    detectionThreshold = lsst.pex.config.Field(
        dtype=float,
        default=10,
        doc="Minimum signal to noise ratio of detected sources "
        "to use for calculating the PSF matching kernel."
    )
    detectionThresholdMax = lsst.pex.config.Field(
        dtype=float,
        default=500,
        doc="Maximum signal to noise ratio of detected sources "
        "to use for calculating the PSF matching kernel."
    )
    maxKernelSources = lsst.pex.config.Field(
        dtype=int,
        default=1000,
        doc="Maximum number of sources to use for calculating the PSF matching kernel."
        "Set to -1 to disable."
    )
    minKernelSources = lsst.pex.config.Field(
        dtype=int,
        default=3,
        doc="Minimum number of sources needed for calculating the PSF matching kernel."
    )
    excludeMaskPlanes = lsst.pex.config.ListField(
        dtype=str,
        default=("NO_DATA", "BAD", "SAT", "EDGE", "FAKE"),
        doc="Mask planes to exclude when selecting sources for PSF matching."
    )
    badMaskPlanes = lsst.pex.config.ListField(
        dtype=str,
        default=("NO_DATA", "BAD", "SAT", "EDGE"),
        doc="Mask planes to interpolate over."
    )
    preserveTemplateMask = lsst.pex.config.ListField(
        dtype=str,
        default=("NO_DATA", "BAD",),
        doc="Mask planes from the template to propagate to the image difference."
    )
    renameTemplateMask = lsst.pex.config.ListField(
        dtype=str,
        default=("SAT", "INJECTED", "INJECTED_CORE",),
        doc="Mask planes from the template to propagate to the image difference"
        "with '_TEMPLATE' appended to the name."
    )
    allowKernelSourceDetection = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc="Re-run source detection for kernel candidates if an error is"
        " encountered while calculating the matching kernel."
    )
    restrictPsfMatchByPatch = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        doc="Fit the PSF matching kernel independently in each patch of the"
        " skymap."
    )

    def setDefaults(self):
        self.makeKernel.kernel.name = "AL"
        # Always include background fitting in the kernel fit,
        # even if it is not subtracted
        self.makeKernel.kernel.active.fitForBackground = True
        self.makeKernel.kernel.active.spatialKernelOrder = 1
        self.makeKernel.kernel.active.spatialBgOrder = 2
        self.sourceSelector.doUnresolved = True  # apply star-galaxy separation
        self.sourceSelector.doIsolated = True  # apply isolated star selection
        self.sourceSelector.doRequirePrimary = True  # apply primary flag selection
        self.sourceSelector.doSkySources = False  # Do not include sky sources
        self.sourceSelector.doSignalToNoise = True  # apply signal to noise filter
        self.sourceSelector.signalToNoise.minimum = 10
        self.sourceSelector.signalToNoise.maximum = 500


class AlardLuptonSubtractConfig(AlardLuptonSubtractBaseConfig, lsst.pipe.base.PipelineTaskConfig,
                                pipelineConnections=AlardLuptonSubtractConnections):
    mode = lsst.pex.config.ChoiceField(
        dtype=str,
        default="convolveTemplate",
        allowed={"auto": "Choose which image to convolve at runtime.",
                 "convolveScience": "Only convolve the science image.",
                 "convolveTemplate": "Only convolve the template image."},
        doc="Choose which image to convolve at runtime, or require that a specific image is convolved."
    )


class AlardLuptonSubtractTask(lsst.pipe.base.PipelineTask):
    """Compute the image difference of a science and template image using
    the Alard & Lupton (1998) algorithm.
    """
    ConfigClass = AlardLuptonSubtractConfig
    _DefaultName = "alardLuptonSubtract"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask("decorrelate")
        self.makeSubtask("makeKernel")
        self.makeSubtask("sourceSelector")
        if self.config.doScaleVariance:
            self.makeSubtask("scaleVariance")

        self.convolutionControl = lsst.afw.math.ConvolutionControl()
        # Normalization is an extra, unnecessary, calculation and will result
        #  in mis-subtraction of the images if there are calibration errors.
        self.convolutionControl.setDoNormalize(False)
        self.convolutionControl.setDoCopyEdge(True)

    def _applyExternalCalibrations(self, exposure, visitSummary):
        """Replace calibrations (psf, and ApCorrMap) on this exposure with
        external ones.".

        Parameters
        ----------
        exposure : `lsst.afw.image.exposure.Exposure`
            Input exposure to adjust calibrations.
        visitSummary : `lsst.afw.table.ExposureCatalog`
            Exposure catalog with external calibrations to be applied. Catalog
            uses the detector id for the catalog id, sorted on id for fast
            lookup.

        Returns
        -------
        exposure : `lsst.afw.image.exposure.Exposure`
            Exposure with adjusted calibrations.
        """
        detectorId = exposure.info.getDetector().getId()

        row = visitSummary.find(detectorId)
        if row is None:
            self.log.warning("Detector id %s not found in external calibrations catalog; "
                             "Using original calibrations.", detectorId)
        else:
            psf = row.getPsf()
            apCorrMap = row.getApCorrMap()
            if psf is None:
                self.log.warning("Detector id %s has None for psf in "
                                 "external calibrations catalog; Using original psf and aperture correction.",
                                 detectorId)
            elif apCorrMap is None:
                self.log.warning("Detector id %s has None for apCorrMap in "
                                 "external calibrations catalog; Using original psf and aperture correction.",
                                 detectorId)
            else:
                exposure.setPsf(psf)
                exposure.info.setApCorrMap(apCorrMap)

        return exposure

    @timeMethod
    def run(self, template, science, sources, skymap=None, visitSummary=None):
        """PSF match, subtract, and decorrelate two images.

        Parameters
        ----------
        template : `lsst.afw.image.ExposureF`
            Template exposure, warped to match the science exposure.
        science : `lsst.afw.image.ExposureF`
            Science exposure to subtract from the template.
        sources : `lsst.afw.table.SourceCatalog`
            Identified sources on the science exposure. This catalog is used to
            select sources in order to perform the AL PSF matching on stamp
            images around them.
        skymap : `lsst.skymap.SkyMap`, optional
            Input definition of geometry/bbox and projection/wcs for
            template exposures.
        visitSummary : `lsst.afw.table.ExposureCatalog`, optional
            Exposure catalog with external calibrations to be applied. Catalog
            uses the detector id for the catalog id, sorted on id for fast
            lookup.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            ``difference`` : `lsst.afw.image.ExposureF`
                Result of subtracting template and science.
            ``matchedTemplate`` : `lsst.afw.image.ExposureF`
                Warped and PSF-matched template exposure.
            ``backgroundModel`` : `lsst.afw.math.Function2D`
                Background model that was fit while solving for the
                PSF-matching kernel
            ``psfMatchingKernel`` : `lsst.afw.math.Kernel`
                Kernel used to PSF-match the convolved image.

        Raises
        ------
        RuntimeError
            If an unsupported convolution mode is supplied.
        RuntimeError
            If there are too few sources to calculate the PSF matching kernel.

        lsst.pipe.base.NoWorkFound
            Raised if fraction of good pixels, defined as not having NO_DATA
            set, is less then the configured requiredTemplateFraction
        """
        self._prepareInputs(template, science, visitSummary=visitSummary)

        #  Calculate estimated image depths, i.e., limiting magnitudes
        maglim_science = self._calculateMagLim(science, fallbackPsfSize=self.sciencePsfSize)
        if np.isnan(maglim_science):
            self.log.warning("Limiting magnitude of the science image is NaN!")
        fluxlim_science = (maglim_science*u.ABmag).to_value(u.nJy)
        maglim_template = self._calculateMagLim(template, fallbackPsfSize=self.templatePsfSize)
        if np.isnan(maglim_template):
            self.log.info("Cannot evaluate template limiting mag; adopting science limiting mag for diffim")
            maglim_diffim = maglim_science
        else:
            fluxlim_template = (maglim_template*u.ABmag).to_value(u.nJy)
            maglim_diffim = (np.sqrt(fluxlim_science**2 + fluxlim_template**2)*u.nJy).to(u.ABmag).value
        self.metadata["scienceLimitingMagnitude"] = maglim_science
        self.metadata["templateLimitingMagnitude"] = maglim_template
        self.metadata["diffimLimitingMagnitude"] = maglim_diffim

        if self.config.mode == "auto":
            convolveTemplate = _shapeTest(template,
                                          science,
                                          fwhmExposureBuffer=self.config.makeKernel.fwhmExposureBuffer,
                                          fwhmExposureGrid=self.config.makeKernel.fwhmExposureGrid)
            if convolveTemplate:
                if self.sciencePsfSize < self.templatePsfSize:
                    self.log.info("Average template PSF size is greater, "
                                  "but science PSF greater in one dimension: convolving template image.")
                else:
                    self.log.info("Science PSF size is greater: convolving template image.")
            else:
                self.log.info("Template PSF size is greater: convolving science image.")
        elif self.config.mode == "convolveTemplate":
            self.log.info("`convolveTemplate` is set: convolving template image.")
            convolveTemplate = True
        elif self.config.mode == "convolveScience":
            self.log.info("`convolveScience` is set: convolving science image.")
            convolveTemplate = False
        else:
            raise RuntimeError("Cannot handle AlardLuptonSubtract mode: %s", self.config.mode)
        if convolveTemplate:
            self.metadata["convolvedExposure"] = "Template"
            convolveMethod = self.runConvolveTemplate
        else:
            self.metadata["convolvedExposure"] = "Science"
            convolveMethod = self.runConvolveScience
        sourceMask = science.mask.clone()
        sourceMask.array |= template[science.getBBox()].mask.array
        selectSources = self._sourceSelector(sources, sourceMask)

        try:
            if self.config.restrictPsfMatchByPatch:
                try:
                    subtractResults = self.matchPsfByPatch(template, science, selectSources, skymap,
                                                           convolveMethod)

                except (RuntimeError, lsst.pex.exceptions.Exception) as e:
                    self.log.warning(f"Error while running per patch: {e}")
                    subtractResults = convolveMethod(template, science, selectSources)
            else:
                subtractResults = convolveMethod(template, science, selectSources)

        except (RuntimeError, lsst.pex.exceptions.Exception) as e:
            self.log.warning("Failed to match template. Checking coverage")
            #  Raise NoWorkFound if template fraction is insufficient
            checkTemplateIsSufficient(template[science.getBBox()], science, self.log,
                                      self.config.minTemplateFractionForExpectedSuccess,
                                      exceptionMessage="Template coverage lower than expected to succeed."
                                      f" Failure is tolerable: {e}")
            #  checkTemplateIsSufficient did not raise NoWorkFound, so raise original exception
            raise e

        metrics = computeDifferenceImageMetrics(science, subtractResults.difference, sources)

        self.metadata["differenceFootprintRatioMean"] = metrics.differenceFootprintRatioMean
        self.metadata["differenceFootprintRatioStdev"] = metrics.differenceFootprintRatioStdev
        self.metadata["differenceFootprintSkyRatioMean"] = metrics.differenceFootprintSkyRatioMean
        self.metadata["differenceFootprintSkyRatioStdev"] = metrics.differenceFootprintSkyRatioStdev
        self.log.info("Mean, stdev of ratio of difference to science "
                      "pixels in star footprints: %5.4f, %5.4f",
                      self.metadata["differenceFootprintRatioMean"],
                      self.metadata["differenceFootprintRatioStdev"])

        return subtractResults

    def matchPsfByPatch(self, template, science, sources, skymap, convolveMethod, subtractResults=None):
        patches = divideExposureByPatches(template, skymap, overlapThreshold=0.05)

        matchedTemplate = np.zeros_like(science.image.array)
        difference = np.zeros_like(science.image.array)
        totalWeight = np.zeros_like(science.image.array)

        templateBorderSize = (template.getBBox().getDimensions()[0] - science.getBBox().getDimensions()[0])//2
        # Create a SpatialCellSet with the desired cell size
        # cellSize = 128  # Adjust this value based on your image scale
        # spatialCellSet = lsst.afw.math.SpatialCellSet(science.getBBox(), cellSize)

        # Use the full image to fit the background, and as a default for any
        # patches that are too small
        boxS = science.getBBox()
        subtractResults = convolveMethod(template[boxS], science, sources, usePatch=False)
        science2 = science.clone()
        # science2.maskedImage -= subtractResults.backgroundModel
        totalWeight += .1
        difference += subtractResults.difference.image.array*totalWeight

        # Loop through the overlapping patches, from largest to smallest overlap
        failedPatches = []

        for patch in patches:
            self.log.info(f"PSF matching patch {patch}")
            patchCorners = patch.wcs.pixelToSky(lsst.geom.Box2D(patch.getInnerBBox()).getCorners())
            warpedPatchBox = computeWarpedBBox(science.wcs, patch.getInnerBBox(), patch.wcs)
            warpedPatchBox.clip(science.getBBox())

            patchPolygon = afwGeom.Polygon(template.wcs.skyToPixel(patchCorners))
            inds = patchPolygon.contains(sources.getX(), sources.getY())
            selectSources = sources[inds]
            patchOuterCorners = patch.wcs.pixelToSky(lsst.geom.Box2D(patch.getOuterBBox()).getCorners())
            patchOuterPolygon = afwGeom.Polygon(template.wcs.skyToPixel(patchOuterCorners))
            patchBBox = lsst.geom.Box2I(patchOuterPolygon.getBBox())
            try:
                templatePsfSize = getPsfFwhm(template[warpedPatchBox].psf, position=patchBBox.getCenter())
                sciencePsfSize = getPsfFwhm(science[warpedPatchBox].psf, position=patchBBox.getCenter())
                self.log.info(f"Science image PSF size in patch: {sciencePsfSize}")
                self.log.info(f"Template image PSF size in patch: {templatePsfSize}")
            except (lsst.pex.exceptions.Exception):
                self.log.info(f"Could not compute PSF size in patch {patch}. ")
            boxS = science.getBBox()
            boxT = boxS.clippedTo(patchBBox)
            boxT.grow(templateBorderSize)
            boxS = boxS.clippedTo(patchBBox)
            try:
                subtractResultsPatch = convolveMethod(template[boxT], science2[boxS], selectSources,
                                                      usePatch=True)
            except (RuntimeError, lsst.pex.exceptions.Exception) as e:
                self.log.warning(f"Failed to fit patch {patch}: {e}")
                failedPatches.append(patch)
                continue
            weight = patchOuterPolygon.createImage(boxS).array
            arrayBegin = lsst.geom.Point2I(0, 0)
            arrayBegin.shift(boxS.getBegin() - science.getBBox().getBegin())
            box = lsst.geom.Box2I(arrayBegin, boxS.getDimensions())
            # weight = patchWeight[science.getBBox()].array
            # patchWeights.append(weight)
            # matchingKernels.append(subtractResults.psfMatchingKernel)
            # matchedTemplates.append(subtractResults.matchedTemplate)
            # diffims.append(subtractResults.difference)
            # if firstPatch:
            #     notInPatch = weight < 1
            #     weight[notInPatch] += 0.01
            #     subtractResults = subtractResultsPatch
            matchedTemplate[box.getSlices()] += subtractResultsPatch.matchedTemplate.image.array*weight
            difference[box.getSlices()] += subtractResultsPatch.difference.image.array*weight
            totalWeight[box.getSlices()] += weight
            # candidate = MyKernelSpatialCellCandidate(patchPolygon, subtractResults.psfMatchingKernel)
            # spatialCellSet.insertCandidate(candidate)
        # if subtractResults is None:
        #     raise RuntimeError(f"Failed to fit any patch from {patches}")

        inds = totalWeight > 0
        matchedTemplate[inds] /= totalWeight[inds]
        difference[inds] /= totalWeight[inds]
        inds = science.mask.array & science.mask.getPlaneBitMask("NO_DATA") > 0
        matchedTemplate[inds] = 0
        difference[inds] = 0

        # matchedTemplate[~inds] = subtractResults.matchedTemplate.image.array[~inds]
        # difference[~inds] = subtractResults.difference.image.array[~inds]
        subtractResults.matchedTemplate.image.array = matchedTemplate
        subtractResults.difference.image.array = difference
        # subtractResults.psfMatchingKernel = CoaddPsf(spatialCellSet)
        return subtractResults

    def runConvolveTemplate(self, template, science, selectSources, usePatch=False):
        """Convolve the template image with a PSF-matching kernel and subtract
        from the science image.

        Parameters
        ----------
        template : `lsst.afw.image.ExposureF`
            Template exposure, warped to match the science exposure.
        science : `lsst.afw.image.ExposureF`
            Science exposure to subtract from the template.
        selectSources : `lsst.afw.table.SourceCatalog`
            Identified sources on the science exposure. This catalog is used to
            select sources in order to perform the AL PSF matching on stamp
            images around them.

        Returns
        -------
        results : `lsst.pipe.base.Struct`

            ``difference`` : `lsst.afw.image.ExposureF`
                Result of subtracting template and science.
            ``matchedTemplate`` : `lsst.afw.image.ExposureF`
                Warped and PSF-matched template exposure.
            ``backgroundModel`` : `lsst.afw.math.Function2D`
                Background model that was fit while solving for the PSF-matching kernel
            ``psfMatchingKernel`` : `lsst.afw.math.Kernel`
                Kernel used to PSF-match the template to the science image.
        """
        self.metadata["convolvedExposure"] = "Template"
        try:
            kernelSources = self.makeKernel.selectKernelSources(template, science,
                                                                candidateList=selectSources,
                                                                preconvolved=False,
                                                                templateFwhmPix=self.templatePsfSize,
                                                                scienceFwhmPix=self.sciencePsfSize)
            kernelResult = self.makeKernel.run(template, science, kernelSources,
                                               preconvolved=False,
                                               templateFwhmPix=self.templatePsfSize,
                                               scienceFwhmPix=self.sciencePsfSize)
        except Exception as e:
            if self.config.allowKernelSourceDetection:
                self.log.warning("Error encountered trying to construct the matching kernel"
                                 f" Running source detection and retrying. {e}")
                kernelSize = self.makeKernel.makeKernelBasisList(
                    self.templatePsfSize, self.sciencePsfSize)[0].getWidth()
                sigmaToFwhm = 2*np.log(2*np.sqrt(2))
                candidateList = self.makeKernel.makeCandidateList(template, science, kernelSize,
                                                                  candidateList=None,
                                                                  sigma=self.sciencePsfSize/sigmaToFwhm)
                kernelSources = self.makeKernel.selectKernelSources(template, science,
                                                                    candidateList=candidateList,
                                                                    preconvolved=False,
                                                                    templateFwhmPix=self.templatePsfSize,
                                                                    scienceFwhmPix=self.sciencePsfSize)
                kernelResult = self.makeKernel.run(template, science, kernelSources,
                                                   preconvolved=False,
                                                   templateFwhmPix=self.templatePsfSize,
                                                   scienceFwhmPix=self.sciencePsfSize)
            else:
                raise e

        matchedTemplate = self._convolveExposure(template, kernelResult.psfMatchingKernel,
                                                 self.convolutionControl,
                                                 bbox=science.getBBox(),
                                                 psf=science.psf,
                                                 photoCalib=science.photoCalib)

        if self.config.doSubtractBackground and ~ usePatch:
            backgroundModel = kernelResult.backgroundModel
        else:
            backgroundModel = None
        difference = _subtractImages(science, matchedTemplate,
                                     backgroundModel=backgroundModel)

        correctedExposure = self.finalize(template, science, difference,
                                          kernelResult.psfMatchingKernel,
                                          templateMatched=True)

        return lsst.pipe.base.Struct(difference=correctedExposure,
                                     matchedTemplate=matchedTemplate,
                                     matchedScience=science,
                                     backgroundModel=backgroundModel,
                                     psfMatchingKernel=kernelResult.psfMatchingKernel,
                                     kernelSources=kernelSources)

    def runConvolveScience(self, template, science, selectSources, usePatch=False):
        """Convolve the science image with a PSF-matching kernel and subtract
        the template image.

        Parameters
        ----------
        template : `lsst.afw.image.ExposureF`
            Template exposure, warped to match the science exposure.
        science : `lsst.afw.image.ExposureF`
            Science exposure to subtract from the template.
        selectSources : `lsst.afw.table.SourceCatalog`
            Identified sources on the science exposure. This catalog is used to
            select sources in order to perform the AL PSF matching on stamp
            images around them.

        Returns
        -------
        results : `lsst.pipe.base.Struct`

            ``difference`` : `lsst.afw.image.ExposureF`
                Result of subtracting template and science.
            ``matchedTemplate`` : `lsst.afw.image.ExposureF`
                Warped template exposure. Note that in this case, the template
                is not PSF-matched to the science image.
            ``backgroundModel`` : `lsst.afw.math.Function2D`
                Background model that was fit while solving for the PSF-matching kernel
            ``psfMatchingKernel`` : `lsst.afw.math.Kernel`
               Kernel used to PSF-match the science image to the template.
        """
        self.metadata["convolvedExposure"] = "Science"
        bbox = science.getBBox()
        kernelSources = self.makeKernel.selectKernelSources(science, template,
                                                            candidateList=selectSources,
                                                            preconvolved=False,
                                                            templateFwhmPix=self.templatePsfSize,
                                                            scienceFwhmPix=self.sciencePsfSize)
        kernelResult = self.makeKernel.run(science, template, kernelSources,
                                           preconvolved=False,
                                           templateFwhmPix=self.templatePsfSize,
                                           scienceFwhmPix=self.sciencePsfSize)
        modelParams = kernelResult.backgroundModel.getParameters()
        # We must invert the background model if the matching kernel is solved for the science image.
        kernelResult.backgroundModel.setParameters([-p for p in modelParams])

        kernelImage = lsst.afw.image.ImageD(kernelResult.psfMatchingKernel.getDimensions())
        norm = kernelResult.psfMatchingKernel.computeImage(kernelImage, doNormalize=False)

        matchedScience = self._convolveExposure(science, kernelResult.psfMatchingKernel,
                                                self.convolutionControl,
                                                psf=template.psf)

        # Place back on native photometric scale
        matchedScience.maskedImage /= norm
        matchedTemplate = template.clone()[bbox]
        matchedTemplate.maskedImage /= norm
        matchedTemplate.setPhotoCalib(science.photoCalib)

        if self.config.doSubtractBackground and ~ usePatch:
            backgroundModel = kernelResult.backgroundModel
        else:
            backgroundModel = None
        difference = _subtractImages(matchedScience, matchedTemplate,
                                     backgroundModel=backgroundModel)

        correctedExposure = self.finalize(template, science, difference,
                                          kernelResult.psfMatchingKernel,
                                          templateMatched=False)

        return lsst.pipe.base.Struct(difference=correctedExposure,
                                     matchedTemplate=matchedTemplate,
                                     matchedScience=matchedScience,
                                     backgroundModel=backgroundModel,
                                     psfMatchingKernel=kernelResult.psfMatchingKernel,
                                     kernelSources=kernelSources)

    def finalize(self, template, science, difference, kernel,
                 templateMatched=True,
                 preConvMode=False,
                 preConvKernel=None,
                 spatiallyVarying=False):
        """Decorrelate the difference image to undo the noise correlations
        caused by convolution.

        Parameters
        ----------
        template : `lsst.afw.image.ExposureF`
            Template exposure, warped to match the science exposure.
        science : `lsst.afw.image.ExposureF`
            Science exposure to subtract from the template.
        difference : `lsst.afw.image.ExposureF`
            Result of subtracting template and science.
        kernel : `lsst.afw.math.Kernel`
            An (optionally spatially-varying) PSF matching kernel
        templateMatched : `bool`, optional
            Was the template PSF-matched to the science image?
        preConvMode : `bool`, optional
            Was the science image preconvolved with its own PSF
            before PSF matching the template?
        preConvKernel : `lsst.afw.detection.Psf`, optional
            If not `None`, then the science image was pre-convolved with
            (the reflection of) this kernel. Must be normalized to sum to 1.
        spatiallyVarying : `bool`, optional
            Compute the decorrelation kernel spatially varying across the image?

        Returns
        -------
        correctedExposure : `lsst.afw.image.ExposureF`
            The decorrelated image difference.
        """
        if self.config.doDecorrelation:
            self.log.info("Decorrelating image difference.")
            # We have cleared the template mask plane, so copy the mask plane of
            # the image difference so that we can calculate correct statistics
            # during decorrelation
            correctedExposure = self.decorrelate.run(science, template[science.getBBox()], difference, kernel,
                                                     templateMatched=templateMatched,
                                                     preConvMode=preConvMode,
                                                     preConvKernel=preConvKernel,
                                                     spatiallyVarying=spatiallyVarying).correctedExposure
        else:
            self.log.info("NOT decorrelating image difference.")
            correctedExposure = difference
        return correctedExposure

    def _calculateMagLim(self, exposure, nsigma=5.0, fallbackPsfSize=None):
        """Calculate an exposure's limiting magnitude.

        This method uses the photometric zeropoint together with the
        PSF size from the average position of the exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The target exposure to calculate the limiting magnitude for.
        nsigma : `float`, optional
            The detection threshold in sigma.
        fallbackPsfSize : `float`, optional
            PSF FWHM to use in the event the exposure PSF cannot be retrieved.

        Returns
        -------
        maglim : `astropy.units.Quantity`
            The limiting magnitude of the exposure, or np.nan.
        """
        if exposure.photoCalib is None:
            return np.nan
        # Set maglim to nan upfront in case on an unexpected RuntimeError
        maglim = np.nan
        try:
            psf = exposure.getPsf()
            psf_shape = psf.computeShape(psf.getAveragePosition())
        except (lsst.pex.exceptions.InvalidParameterError,
                afwDetection.InvalidPsfError,
                lsst.pex.exceptions.RangeError):
            if fallbackPsfSize is not None:
                self.log.info("Unable to evaluate PSF, using fallback FWHM %f", fallbackPsfSize)
                psf_area = np.pi*(fallbackPsfSize/2)**2
                zeropoint = exposure.photoCalib.instFluxToMagnitude(1)
                maglim = zeropoint - 2.5*np.log10(nsigma*np.sqrt(psf_area))
            else:
                self.log.info("Unable to evaluate PSF, setting maglim to nan")
                maglim = np.nan
        else:
            # Get a more accurate area than `psf_shape.getArea()` via moments
            psf_area = np.pi*np.sqrt(psf_shape.getIxx()*psf_shape.getIyy())
            zeropoint = exposure.photoCalib.instFluxToMagnitude(1)
            maglim = zeropoint - 2.5*np.log10(nsigma*np.sqrt(psf_area))
        finally:
            return maglim

    @staticmethod
    def _validateExposures(template, science):
        """Check that the WCS of the two Exposures match, the template bbox
        contains the science bbox, and that the bands match.

        Parameters
        ----------
        template : `lsst.afw.image.ExposureF`
            Template exposure, warped to match the science exposure.
        science : `lsst.afw.image.ExposureF`
            Science exposure to subtract from the template.

        Raises
        ------
        AssertionError
            Raised if the WCS of the template is not equal to the science WCS,
            if the science image is not fully contained in the template
            bounding box, or if the bands do not match.
        """
        assert template.wcs == science.wcs, \
            "Template and science exposure WCS are not identical."
        templateBBox = template.getBBox()
        scienceBBox = science.getBBox()
        assert science.filter.bandLabel == template.filter.bandLabel, \
            "Science and template exposures have different bands: %s, %s" % \
            (science.filter, template.filter)

        assert templateBBox.contains(scienceBBox), \
            "Template bbox does not contain all of the science image."

    def _convolveExposure(self, exposure, kernel, convolutionControl,
                          bbox=None,
                          psf=None,
                          photoCalib=None,
                          interpolateBadMaskPlanes=False,
                          ):
        """Convolve an exposure with the given kernel.

        Parameters
        ----------
        exposure : `lsst.afw.Exposure`
            exposure to convolve.
        kernel : `lsst.afw.math.LinearCombinationKernel`
            PSF matching kernel computed in the ``makeKernel`` subtask.
        convolutionControl : `lsst.afw.math.ConvolutionControl`
            Configuration for convolve algorithm.
        bbox : `lsst.geom.Box2I`, optional
            Bounding box to trim the convolved exposure to.
        psf : `lsst.afw.detection.Psf`, optional
            Point spread function (PSF) to set for the convolved exposure.
        photoCalib : `lsst.afw.image.PhotoCalib`, optional
            Photometric calibration of the convolved exposure.

        Returns
        -------
        convolvedExp : `lsst.afw.Exposure`
            The convolved image.
        """
        convolvedExposure = exposure.clone()
        if psf is not None:
            convolvedExposure.setPsf(psf)
        if photoCalib is not None:
            convolvedExposure.setPhotoCalib(photoCalib)
        if interpolateBadMaskPlanes and self.config.badMaskPlanes is not None:
            nInterp = _interpolateImage(convolvedExposure.maskedImage,
                                        self.config.badMaskPlanes)
            self.metadata["nInterpolated"] = nInterp
        convolvedImage = lsst.afw.image.MaskedImageF(convolvedExposure.getBBox())
        lsst.afw.math.convolve(convolvedImage, convolvedExposure.maskedImage, kernel, convolutionControl)
        convolvedExposure.setMaskedImage(convolvedImage)
        if bbox is None:
            return convolvedExposure
        else:
            return convolvedExposure[bbox]

    def _sourceSelector(self, sources, mask):
        """Select sources from a catalog that meet the selection criteria.

        Parameters
        ----------
        sources : `lsst.afw.table.SourceCatalog`
            Input source catalog to select sources from.
        mask : `lsst.afw.image.Mask`
            The image mask plane to use to reject sources
            based on their location on the ccd.

        Returns
        -------
        selectSources : `lsst.afw.table.SourceCatalog`
            The input source catalog, with flagged and low signal-to-noise
            sources removed.

        Raises
        ------
        RuntimeError
            If there are too few sources to compute the PSF matching kernel
            remaining after source selection.
        """

        selected = self.sourceSelector.selectSources(sources).selected
        nInitialSelected = np.count_nonzero(selected)
        nSelected = np.count_nonzero(selected)
        self.log.info("Rejecting %i candidate sources: an excluded template mask plane is set.",
                      nInitialSelected - nSelected)
        selectSources = sources[selected].copy(deep=True)
        # Trim selectSources if they exceed ``maxKernelSources``.
        # Keep the highest signal-to-noise sources of those selected.
        if (len(selectSources) > self.config.maxKernelSources) & (self.config.maxKernelSources > 0):
            signalToNoise = selectSources.getPsfInstFlux()/selectSources.getPsfInstFluxErr()
            indices = np.argsort(signalToNoise)
            indices = indices[-self.config.maxKernelSources:]
            selected = np.zeros(len(selectSources), dtype=bool)
            selected[indices] = True
            selectSources = selectSources[selected].copy(deep=True)

        self.log.info("%i/%i=%.1f%% of sources selected for PSF matching from the input catalog",
                      len(selectSources), len(sources), 100*len(selectSources)/len(sources))
        if len(selectSources) < self.config.minKernelSources:
            self.log.error("Too few sources to calculate the PSF matching kernel: "
                           "%i selected but %i needed for the calculation.",
                           len(selectSources), self.config.minKernelSources)
            if not self.config.allowKernelSourceDetection:
                raise RuntimeError("Cannot compute PSF matching kernel: too few sources selected.")
        self.metadata["nPsfSources"] = len(selectSources)

        return selectSources

    def _prepareInputs(self, template, science, visitSummary=None):
        """Perform preparatory calculations common to all Alard&Lupton Tasks.

        Parameters
        ----------
        template : `lsst.afw.image.ExposureF`
            Template exposure, warped to match the science exposure. The
            variance plane of the template image is modified in place.
        science : `lsst.afw.image.ExposureF`
            Science exposure to subtract from the template. The variance plane
            of the science image is modified in place.
        visitSummary : `lsst.afw.table.ExposureCatalog`, optional
            Exposure catalog with external calibrations to be applied.  Catalog
            uses the detector id for the catalog id, sorted on id for fast
            lookup.
        """
        self._validateExposures(template, science)
        if visitSummary is not None:
            self._applyExternalCalibrations(science, visitSummary=visitSummary)
        templateCoverageFraction = checkTemplateIsSufficient(
            template[science.getBBox()], science, self.log,
            requiredTemplateFraction=self.config.requiredTemplateFraction,
            exceptionMessage="Not attempting subtraction. To force subtraction,"
            " set config requiredTemplateFraction=0"
        )
        self.metadata["templateCoveragePercent"] = 100*templateCoverageFraction

        if self.config.doScaleVariance:
            # Scale the variance of the template and science images before
            # convolution, subtraction, or decorrelation so that they have the
            # correct ratio.
            templateVarFactor = self.scaleVariance.run(template.maskedImage)
            sciVarFactor = self.scaleVariance.run(science.maskedImage)
            self.log.info("Template variance scaling factor: %.2f", templateVarFactor)
            self.metadata["scaleTemplateVarianceFactor"] = templateVarFactor
            self.log.info("Science variance scaling factor: %.2f", sciVarFactor)
            self.metadata["scaleScienceVarianceFactor"] = sciVarFactor

        # Erase existing detection mask planes.
        #  We don't want the detection mask from the science image
        self.updateMasks(template, science)

        # Calling getPsfFwhm on template.psf fails on some rare occasions when
        # the template has no input exposures at the average position of the
        # stars. So we try getPsfFwhm first on template, and if that fails we
        # evaluate the PSF on a grid specified by fwhmExposure* fields.
        # To keep consistent definitions for PSF size on the template and
        # science images, we use the same method for both.
        # In the try block below, we catch two exceptions:
        # 1. InvalidParameterError, in case the point where we are evaluating
        #    the PSF lands in a gap in the template.
        # 2. RangeError, in case the template coverage is so poor that we end
        #    up near a region with no data.
        try:
            self.templatePsfSize = getPsfFwhm(template.psf)
            self.sciencePsfSize = getPsfFwhm(science.psf)
        except lsst.pex.exceptions.Exception:
            # Catch a broad range of exceptions, since some are C++ only
            # Catching:
            #  - lsst::geom::SingularTransformException
            #  - lsst.pex.exceptions.InvalidParameterError
            #  - lsst.pex.exceptions.RangeError
            self.log.info("Unable to evaluate PSF at the average position. "
                          "Evaluting PSF on a grid of points."
                          )
            self.templatePsfSize = evaluateMeanPsfFwhm(
                template,
                fwhmExposureBuffer=self.config.makeKernel.fwhmExposureBuffer,
                fwhmExposureGrid=self.config.makeKernel.fwhmExposureGrid
            )
            self.sciencePsfSize = evaluateMeanPsfFwhm(
                science,
                fwhmExposureBuffer=self.config.makeKernel.fwhmExposureBuffer,
                fwhmExposureGrid=self.config.makeKernel.fwhmExposureGrid
            )
        self.log.info("Science PSF FWHM: %f pixels", self.sciencePsfSize)
        self.log.info("Template PSF FWHM: %f pixels", self.templatePsfSize)
        self.metadata["sciencePsfSize"] = self.sciencePsfSize
        self.metadata["templatePsfSize"] = self.templatePsfSize

    def updateMasks(self, template, science):
        """Update the science and template mask planes before differencing.

        Parameters
        ----------
        template : `lsst.afw.image.Exposure`
            Template exposure, warped to match the science exposure.
            The template mask planes will be erased, except for a few specified
            in the task config.
        science : `lsst.afw.image.Exposure`
            Science exposure to subtract from the template.
            The DETECTED and DETECTED_NEGATIVE mask planes of the science image
            will be erased.
        """
        self._clearMask(science.mask, clearMaskPlanes=["DETECTED", "DETECTED_NEGATIVE"])

        # We will clear ALL template mask planes, except for those specified
        # via the `preserveTemplateMask` config. Mask planes specified via
        # the `renameTemplateMask` config will be copied to new planes with
        # "_TEMPLATE" appended to their names, and the original mask plane will
        # be cleared.
        clearMaskPlanes = [mp for mp in template.mask.getMaskPlaneDict().keys()
                           if mp not in self.config.preserveTemplateMask]
        renameMaskPlanes = [mp for mp in self.config.renameTemplateMask
                            if mp in template.mask.getMaskPlaneDict().keys()]

        # propagate the mask plane related to Fake source injection
        # NOTE: the fake source injection sets FAKE plane, but it should be INJECTED
        # NOTE: This can be removed in DM-40796
        if "FAKE" in science.mask.getMaskPlaneDict().keys():
            self.log.info("Adding injected mask plane to science image")
            self._renameMaskPlanes(science.mask, "FAKE", "INJECTED")
        if "FAKE" in template.mask.getMaskPlaneDict().keys():
            self.log.info("Adding injected mask plane to template image")
            self._renameMaskPlanes(template.mask, "FAKE", "INJECTED_TEMPLATE")
            if "INJECTED" in renameMaskPlanes:
                renameMaskPlanes.remove("INJECTED")
            if "INJECTED_TEMPLATE" in clearMaskPlanes:
                clearMaskPlanes.remove("INJECTED_TEMPLATE")

        for maskPlane in renameMaskPlanes:
            self._renameMaskPlanes(template.mask, maskPlane, maskPlane + "_TEMPLATE")
        self._clearMask(template.mask, clearMaskPlanes=clearMaskPlanes)

    @staticmethod
    def _renameMaskPlanes(mask, maskPlane, newMaskPlane):
        """Rename a mask plane by adding the new name and copying the data.

        Parameters
        ----------
        mask : `lsst.afw.image.Mask`
            The mask image to update in place.
        maskPlane : `str`
            The name of the existing mask plane to copy.
        newMaskPlane : `str`
            The new name of the mask plane that will be added.
            If the mask plane already exists, it will be updated in place.
        """
        mask.addMaskPlane(newMaskPlane)
        originBitMask = mask.getPlaneBitMask(maskPlane)
        destinationBitMask = mask.getPlaneBitMask(newMaskPlane)
        mask.array |= ((mask.array & originBitMask) > 0)*destinationBitMask

    def _clearMask(self, mask, clearMaskPlanes=None):
        """Clear the mask plane of an exposure.

        Parameters
        ----------
        mask : `lsst.afw.image.Mask`
            The mask plane to erase, which will be modified in place.
        clearMaskPlanes : `list` of `str`, optional
            Erase the specified mask planes.
            If not supplied, the entire mask will be erased.
        """
        if clearMaskPlanes is None:
            clearMaskPlanes = list(mask.getMaskPlaneDict().keys())

        bitMaskToClear = mask.getPlaneBitMask(clearMaskPlanes)
        mask &= ~bitMaskToClear


class AlardLuptonPreconvolveSubtractConnections(SubtractInputConnections,
                                                SubtractScoreOutputConnections):
    pass


class AlardLuptonPreconvolveSubtractConfig(AlardLuptonSubtractBaseConfig, lsst.pipe.base.PipelineTaskConfig,
                                           pipelineConnections=AlardLuptonPreconvolveSubtractConnections):
    pass


class AlardLuptonPreconvolveSubtractTask(AlardLuptonSubtractTask):
    """Subtract a template from a science image, convolving the science image
    before computing the kernel, and also convolving the template before
    subtraction.
    """
    ConfigClass = AlardLuptonPreconvolveSubtractConfig
    _DefaultName = "alardLuptonPreconvolveSubtract"

    def run(self, template, science, sources, skymap=None, visitSummary=None):
        """Preconvolve the science image with its own PSF,
        convolve the template image with a PSF-matching kernel and subtract
        from the preconvolved science image.

        Parameters
        ----------
        template : `lsst.afw.image.ExposureF`
            The template image, which has previously been warped to the science
            image. The template bbox will be padded by a few pixels compared to
            the science bbox.
        science : `lsst.afw.image.ExposureF`
            The science exposure.
        sources : `lsst.afw.table.SourceCatalog`
            Identified sources on the science exposure. This catalog is used to
            select sources in order to perform the AL PSF matching on stamp
            images around them.
        skymap : `lsst.skymap.SkyMap`, optional
            Input definition of geometry/bbox and projection/wcs for
            template exposures.
        visitSummary : `lsst.afw.table.ExposureCatalog`, optional
            Exposure catalog with complete external calibrations. Catalog uses
            the detector id for the catalog id, sorted on id for fast lookup.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            ``scoreExposure`` : `lsst.afw.image.ExposureF`
                Result of subtracting the convolved template and science
                images. Attached PSF is that of the original science image.
            ``matchedTemplate`` : `lsst.afw.image.ExposureF`
                Warped and PSF-matched template exposure. Attached PSF is that
                of the original science image.
            ``matchedScience`` : `lsst.afw.image.ExposureF`
                The science exposure after convolving with its own PSF.
                Attached PSF is that of the original science image.
            ``backgroundModel`` : `lsst.afw.math.Function2D`
                Background model that was fit while solving for the
                PSF-matching kernel
            ``psfMatchingKernel`` : `lsst.afw.math.Kernel`
                Final kernel used to PSF-match the template to the science
                image.
        """
        self._prepareInputs(template, science, visitSummary=visitSummary)

        # TODO: DM-37212 we need to mirror the kernel in order to get correct cross correlation
        scienceKernel = science.psf.getKernel()
        matchedScience = self._convolveExposure(science, scienceKernel, self.convolutionControl,
                                                interpolateBadMaskPlanes=True)
        self.metadata["convolvedExposure"] = "Preconvolution"
        try:
            selectSources = self._sourceSelector(sources, matchedScience.mask)
            subtractResults = self.runPreconvolve(template, science, matchedScience,
                                                  selectSources, scienceKernel)

        except (RuntimeError, lsst.pex.exceptions.Exception) as e:
            self.log.warning("Failed to match template. Checking coverage")
            #  Raise NoWorkFound if template fraction is insufficient
            checkTemplateIsSufficient(template[science.getBBox()], science, self.log,
                                      self.config.minTemplateFractionForExpectedSuccess,
                                      exceptionMessage="Template coverage lower than expected to succeed."
                                      f" Failure is tolerable: {e}")
            #  checkTemplateIsSufficient did not raise NoWorkFound, so raise original exception
            raise e

        return subtractResults

    def runPreconvolve(self, template, science, matchedScience, selectSources, preConvKernel):
        """Convolve the science image with its own PSF, then convolve the
        template with a matching kernel and subtract to form the Score
        exposure.

        Parameters
        ----------
        template : `lsst.afw.image.ExposureF`
            Template exposure, warped to match the science exposure.
        science : `lsst.afw.image.ExposureF`
            Science exposure to subtract from the template.
        matchedScience : `lsst.afw.image.ExposureF`
            The science exposure, convolved with the reflection of its own PSF.
        selectSources : `lsst.afw.table.SourceCatalog`
            Identified sources on the science exposure. This catalog is used to
            select sources in order to perform the AL PSF matching on stamp
            images around them.
        preConvKernel : `lsst.afw.math.Kernel`
            The reflection of the kernel that was used to preconvolve the
            `science` exposure. Must be normalized to sum to 1.

        Returns
        -------
        results : `lsst.pipe.base.Struct`

            ``scoreExposure`` : `lsst.afw.image.ExposureF`
                Result of subtracting the convolved template and science
                images. Attached PSF is that of the original science image.
            ``matchedTemplate`` : `lsst.afw.image.ExposureF`
                Warped and PSF-matched template exposure. Attached PSF is that
                of the original science image.
            ``matchedScience`` : `lsst.afw.image.ExposureF`
                The science exposure after convolving with its own PSF.
                Attached PSF is that of the original science image.
            ``backgroundModel`` : `lsst.afw.math.Function2D`
                Background model that was fit while solving for the
                PSF-matching kernel
            ``psfMatchingKernel`` : `lsst.afw.math.Kernel`
                Final kernel used to PSF-match the template to the science
                image.
        """
        bbox = science.getBBox()
        innerBBox = preConvKernel.shrinkBBox(bbox)

        kernelSources = self.makeKernel.selectKernelSources(template[innerBBox], matchedScience[innerBBox],
                                                            candidateList=selectSources,
                                                            preconvolved=True,
                                                            templateFwhmPix=self.templatePsfSize,
                                                            scienceFwhmPix=self.sciencePsfSize)
        kernelResult = self.makeKernel.run(template[innerBBox], matchedScience[innerBBox], kernelSources,
                                           preconvolved=True,
                                           templateFwhmPix=self.templatePsfSize,
                                           scienceFwhmPix=self.sciencePsfSize)

        matchedTemplate = self._convolveExposure(template, kernelResult.psfMatchingKernel,
                                                 self.convolutionControl,
                                                 bbox=bbox,
                                                 psf=science.psf,
                                                 interpolateBadMaskPlanes=True,
                                                 photoCalib=science.photoCalib)
        score = _subtractImages(matchedScience, matchedTemplate,
                                backgroundModel=(kernelResult.backgroundModel
                                                 if self.config.doSubtractBackground else None))
        correctedScore = self.finalize(template[bbox], science, score,
                                       kernelResult.psfMatchingKernel,
                                       templateMatched=True, preConvMode=True,
                                       preConvKernel=preConvKernel)

        return lsst.pipe.base.Struct(scoreExposure=correctedScore,
                                     matchedTemplate=matchedTemplate,
                                     matchedScience=matchedScience,
                                     backgroundModel=kernelResult.backgroundModel,
                                     psfMatchingKernel=kernelResult.psfMatchingKernel,
                                     kernelSources=kernelSources)


def checkTemplateIsSufficient(templateExposure, scienceExposure, logger, requiredTemplateFraction=0.,
                              exceptionMessage=""):
    """Raise NoWorkFound if template coverage < requiredTemplateFraction

    Parameters
    ----------
    templateExposure : `lsst.afw.image.ExposureF`
        The template exposure to check
    logger : `logging.Logger`
        Logger for printing output.
    requiredTemplateFraction : `float`, optional
        Fraction of pixels of the science image required to have coverage
        in the template.
    exceptionMessage : `str`, optional
        Message to include in the exception raised if the template coverage
        is insufficient.

    Returns
    -------
    templateCoverageFraction: `float`
        Fraction of pixels in the template with data.

    Raises
    ------
    lsst.pipe.base.NoWorkFound
        Raised if fraction of good pixels, defined as not having NO_DATA
        set, is less than the requiredTemplateFraction
    """
    # Count the number of pixels with the NO_DATA mask bit set
    # counting NaN pixels is insufficient because pixels without data are often intepolated over)
    noTemplate = templateExposure.mask.array & templateExposure.mask.getPlaneBitMask('NO_DATA')
    # Also need to account for missing data in the science image,
    # because template coverage there doesn't help
    noScience = scienceExposure.mask.array & scienceExposure.mask.getPlaneBitMask('NO_DATA')
    pixNoData = np.count_nonzero(noTemplate | noScience)
    pixGood = templateExposure.getBBox().getArea() - pixNoData
    templateCoverageFraction = pixGood/templateExposure.getBBox().getArea()
    logger.info("template has %d good pixels (%.1f%%)", pixGood, 100*templateCoverageFraction)

    if templateCoverageFraction < requiredTemplateFraction:
        message = ("Insufficient Template Coverage. (%.1f%% < %.1f%%)" % (
                   100*templateCoverageFraction,
                   100*requiredTemplateFraction))
        raise lsst.pipe.base.NoWorkFound(message + " " + exceptionMessage)
    return templateCoverageFraction


def _subtractImages(science, template, backgroundModel=None):
    """Subtract template from science, propagating relevant metadata.

    Parameters
    ----------
    science : `lsst.afw.Exposure`
        The input science image.
    template : `lsst.afw.Exposure`
        The template to subtract from the science image.
    backgroundModel : `lsst.afw.MaskedImage`, optional
        Differential background model

    Returns
    -------
    difference : `lsst.afw.Exposure`
        The subtracted image.
    """
    difference = science.clone()
    if backgroundModel is not None:
        difference.maskedImage -= backgroundModel
    difference.maskedImage -= template.maskedImage
    return difference


def _shapeTest(exp1, exp2, fwhmExposureBuffer, fwhmExposureGrid):
    """Determine that the PSF of ``exp1`` is not wider than that of ``exp2``.

    Parameters
    ----------
    exp1 : `~lsst.afw.image.Exposure`
        Exposure with the reference point spread function (PSF) to evaluate.
    exp2 : `~lsst.afw.image.Exposure`
        Exposure with a candidate point spread function (PSF) to evaluate.
    fwhmExposureBuffer : `float`
        Fractional buffer margin to be left out of all sides of the image
        during the construction of the grid to compute mean PSF FWHM in an
        exposure, if the PSF is not available at its average position.
    fwhmExposureGrid : `int`
        Grid size to compute the mean FWHM in an exposure, if the PSF is not
        available at its average position.
    Returns
    -------
    result : `bool`
        True if ``exp1`` has a PSF that is not wider than that of ``exp2`` in
        either dimension.
    """
    try:
        shape1 = getPsfFwhm(exp1.psf, average=False)
        shape2 = getPsfFwhm(exp2.psf, average=False)
    except (lsst.pex.exceptions.InvalidParameterError, lsst.pex.exceptions.RangeError):
        shape1 = evaluateMeanPsfFwhm(exp1,
                                     fwhmExposureBuffer=fwhmExposureBuffer,
                                     fwhmExposureGrid=fwhmExposureGrid
                                     )
        shape2 = evaluateMeanPsfFwhm(exp2,
                                     fwhmExposureBuffer=fwhmExposureBuffer,
                                     fwhmExposureGrid=fwhmExposureGrid
                                     )
        return shape1 <= shape2

    # Results from getPsfFwhm is a tuple of two values, one for each dimension.
    xTest = shape1[0] <= shape2[0]
    yTest = shape1[1] <= shape2[1]
    return xTest | yTest


def _interpolateImage(maskedImage, badMaskPlanes, fallbackValue=None):
    """Replace masked image pixels with interpolated values.

    Parameters
    ----------
    maskedImage : `lsst.afw.image.MaskedImage`
        Image on which to perform interpolation.
    badMaskPlanes : `list` of `str`
        List of mask planes to interpolate over.
    fallbackValue : `float`, optional
        Value to set when interpolation fails.

    Returns
    -------
    result: `float`
        The number of masked pixels that were replaced.
    """
    imgBadMaskPlanes = [
        maskPlane for maskPlane in badMaskPlanes if maskPlane in maskedImage.mask.getMaskPlaneDict()
    ]

    image = maskedImage.image.array
    badPixels = (maskedImage.mask.array & maskedImage.mask.getPlaneBitMask(imgBadMaskPlanes)) > 0
    image[badPixels] = np.nan
    if fallbackValue is None:
        fallbackValue = np.nanmedian(image)
    # For this initial implementation, skip the interpolation and just fill with
    # the median value.
    image[badPixels] = fallbackValue
    return np.sum(badPixels)
