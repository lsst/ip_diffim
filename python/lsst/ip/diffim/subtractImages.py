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

import lsst.afw.image
import lsst.afw.math
import lsst.geom
from lsst.ip.diffim.utils import getPsfFwhm
from lsst.meas.algorithms import ScaleVarianceTask
import lsst.pex.config
import lsst.pipe.base
from lsst.pipe.base import connectionTypes
from . import MakeKernelTask, DecorrelateALKernelTask

__all__ = ["AlardLuptonSubtractConfig", "AlardLuptonSubtractTask"]

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
    finalizedPsfApCorrCatalog = connectionTypes.Input(
        doc=("Per-visit finalized psf models and aperture correction maps. "
             "These catalogs use the detector id for the catalog id, "
             "sorted on id for fast lookup."),
        dimensions=("instrument", "visit"),
        storageClass="ExposureCatalog",
        name="finalized_psf_ap_corr_catalog",
    )


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


class AlardLuptonSubtractConnections(SubtractInputConnections, SubtractImageOutputConnections):

    def __init__(self, *, config=None):
        super().__init__(config=config)
        if not config.doApplyFinalizedPsf:
            self.inputs.remove("finalizedPsfApCorrCatalog")


class AlardLuptonSubtractConfig(lsst.pipe.base.PipelineTaskConfig,
                                pipelineConnections=AlardLuptonSubtractConnections):
    mode = lsst.pex.config.ChoiceField(
        dtype=str,
        default="auto",
        allowed={"auto": "Choose which image to convolve at runtime.",
                 "convolveScience": "Only convolve the science image.",
                 "convolveTemplate": "Only convolve the template image."},
        doc="Choose which image to convolve at runtime, or require that a specific image is convolved."
    )
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
        doc="Abort task if template covers less than this fraction of pixels."
        " Setting to 0 will always attempt image subtraction."
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
        default=True,
    )
    doApplyFinalizedPsf = lsst.pex.config.Field(
        doc="Replace science Exposure's psf and aperture correction map"
        " with those in finalizedPsfApCorrCatalog.",
        dtype=bool,
        default=False,
    )

    forceCompatibility = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc="Set up and run diffim using settings that ensure the results"
        "are compatible with the old version in pipe_tasks.",
        deprecated="This option is only for backwards compatibility purposes"
        " and will be removed after v24.",
    )

    def setDefaults(self):
        self.makeKernel.kernel.name = "AL"
        self.makeKernel.kernel.active.fitForBackground = self.doSubtractBackground
        self.makeKernel.kernel.active.spatialKernelOrder = 1
        self.makeKernel.kernel.active.spatialBgOrder = 2

    def validate(self):
        if self.forceCompatibility:
            self.mode = "convolveTemplate"


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
        if self.config.doScaleVariance:
            self.makeSubtask("scaleVariance")

        self.convolutionControl = lsst.afw.math.ConvolutionControl()
        # Normalization is an extra, unnecessary, calculation and will result
        #  in mis-subtraction of the images if there are calibration errors.
        self.convolutionControl.setDoNormalize(False)
        self.convolutionControl.setDoCopyEdge(True)

    def _applyExternalCalibrations(self, exposure, finalizedPsfApCorrCatalog):
        """Replace calibrations (psf, and ApCorrMap) on this exposure with external ones.".

        Parameters
        ----------
        exposure : `lsst.afw.image.exposure.Exposure`
            Input exposure to adjust calibrations.
        finalizedPsfApCorrCatalog : `lsst.afw.table.ExposureCatalog`
            Exposure catalog with finalized psf models and aperture correction
            maps to be applied if config.doApplyFinalizedPsf=True.  Catalog uses
            the detector id for the catalog id, sorted on id for fast lookup.

        Returns
        -------
        exposure : `lsst.afw.image.exposure.Exposure`
            Exposure with adjusted calibrations.
        """
        detectorId = exposure.info.getDetector().getId()

        row = finalizedPsfApCorrCatalog.find(detectorId)
        if row is None:
            self.log.warning("Detector id %s not found in finalizedPsfApCorrCatalog; "
                             "Using original psf.", detectorId)
        else:
            psf = row.getPsf()
            apCorrMap = row.getApCorrMap()
            if psf is None:
                self.log.warning("Detector id %s has None for psf in "
                                 "finalizedPsfApCorrCatalog; Using original psf and aperture correction.",
                                 detectorId)
            elif apCorrMap is None:
                self.log.warning("Detector id %s has None for apCorrMap in "
                                 "finalizedPsfApCorrCatalog; Using original psf and aperture correction.",
                                 detectorId)
            else:
                exposure.setPsf(psf)
                exposure.info.setApCorrMap(apCorrMap)

        return exposure

    def run(self, template, science, sources, finalizedPsfApCorrCatalog=None):
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
        finalizedPsfApCorrCatalog : `lsst.afw.table.ExposureCatalog`, optional
            Exposure catalog with finalized psf models and aperture correction
            maps to be applied if config.doApplyFinalizedPsf=True.  Catalog uses
            the detector id for the catalog id, sorted on id for fast lookup.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            ``difference`` : `lsst.afw.image.ExposureF`
                Result of subtracting template and science.
            ``matchedTemplate`` : `lsst.afw.image.ExposureF`
                Warped and PSF-matched template exposure.
            ``backgroundModel`` : `lsst.afw.math.Chebyshev1Function2D`
                Background model that was fit while solving for the PSF-matching kernel
            ``psfMatchingKernel`` : `lsst.afw.math.Kernel`
                Kernel used to PSF-match the convolved image.

        Raises
        ------
        RuntimeError
            If an unsupported convolution mode is supplied.
        lsst.pipe.base.NoWorkFound
            Raised if fraction of good pixels, defined as not having NO_DATA
            set, is less then the configured requiredTemplateFraction
        """
        self._validateExposures(template, science)
        if self.config.doApplyFinalizedPsf:
            self._applyExternalCalibrations(science,
                                            finalizedPsfApCorrCatalog=finalizedPsfApCorrCatalog)
        checkTemplateIsSufficient(template, self.log,
                                  requiredTemplateFraction=self.config.requiredTemplateFraction)
        if self.config.forceCompatibility:
            # Compatibility option to maintain old functionality
            # This should be removed in the future!
            self.log.warning("Running with `config.forceCompatibility=True`")
            sources = None
        sciencePsfSize = getPsfFwhm(science.psf)
        templatePsfSize = getPsfFwhm(template.psf)
        self.log.info("Science PSF size: %f", sciencePsfSize)
        self.log.info("Template PSF size: %f", templatePsfSize)
        if self.config.mode == "auto":
            if sciencePsfSize < templatePsfSize:
                self.log.info("Template PSF size is greater: convolving science image.")
                convolveTemplate = False
            else:
                self.log.info("Science PSF size is greater: convolving template image.")
                convolveTemplate = True
        elif self.config.mode == "convolveTemplate":
            self.log.info("`convolveTemplate` is set: convolving template image.")
            convolveTemplate = True
        elif self.config.mode == "convolveScience":
            self.log.info("`convolveScience` is set: convolving science image.")
            convolveTemplate = False
        else:
            raise RuntimeError("Cannot handle AlardLuptonSubtract mode: %s", self.config.mode)

        if self.config.doScaleVariance and ~self.config.forceCompatibility:
            # Scale the variance of the template and science images before
            # convolution, subtraction, or decorrelation so that they have the
            # correct ratio.
            templateVarFactor = self.scaleVariance.run(template.maskedImage)
            sciVarFactor = self.scaleVariance.run(science.maskedImage)
            self.log.info("Template variance scaling factor: %.2f", templateVarFactor)
            self.metadata.add("scaleTemplateVarianceFactor", templateVarFactor)
            self.log.info("Science variance scaling factor: %.2f", sciVarFactor)
            self.metadata.add("scaleScienceVarianceFactor", sciVarFactor)

        kernelSources = self.makeKernel.selectKernelSources(template, science,
                                                            candidateList=sources,
                                                            preconvolved=False)
        if convolveTemplate:
            subtractResults = self.runConvolveTemplate(template, science, kernelSources)
        else:
            subtractResults = self.runConvolveScience(template, science, kernelSources)

        if self.config.doScaleVariance and self.config.forceCompatibility:
            # The old behavior scaled the variance of the final image difference.
            diffimVarFactor = self.scaleVariance.run(subtractResults.difference.maskedImage)
            self.log.info("Diffim variance scaling factor: %.2f", diffimVarFactor)
            self.metadata.add("scaleDiffimVarianceFactor", diffimVarFactor)

        return subtractResults

    def runConvolveTemplate(self, template, science, sources):
        """Convolve the template image with a PSF-matching kernel and subtract
        from the science image.

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

        Returns
        -------
        results : `lsst.pipe.base.Struct`

            ``difference`` : `lsst.afw.image.ExposureF`
                Result of subtracting template and science.
            ``matchedTemplate`` : `lsst.afw.image.ExposureF`
                Warped and PSF-matched template exposure.
            ``backgroundModel`` : `lsst.afw.math.Chebyshev1Function2D`
                Background model that was fit while solving for the PSF-matching kernel
            ``psfMatchingKernel`` : `lsst.afw.math.Kernel`
                Kernel used to PSF-match the template to the science image.
        """
        if self.config.forceCompatibility:
            # Compatibility option to maintain old behavior
            # This should be removed in the future!
            template = template[science.getBBox()]
        kernelResult = self.makeKernel.run(template, science, sources, preconvolved=False)

        matchedTemplate = self._convolveExposure(template, kernelResult.psfMatchingKernel,
                                                 self.convolutionControl,
                                                 bbox=science.getBBox(),
                                                 psf=science.psf,
                                                 photoCalib=science.getPhotoCalib())
        difference = _subtractImages(science, matchedTemplate,
                                     backgroundModel=(kernelResult.backgroundModel
                                                      if self.config.doSubtractBackground else None))
        correctedExposure = self.finalize(template, science, difference, kernelResult.psfMatchingKernel,
                                          templateMatched=True)

        return lsst.pipe.base.Struct(difference=correctedExposure,
                                     matchedTemplate=matchedTemplate,
                                     matchedScience=science,
                                     backgroundModel=kernelResult.backgroundModel,
                                     psfMatchingKernel=kernelResult.psfMatchingKernel)

    def runConvolveScience(self, template, science, sources):
        """Convolve the science image with a PSF-matching kernel and subtract the template image.

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

        Returns
        -------
        results : `lsst.pipe.base.Struct`

            ``difference`` : `lsst.afw.image.ExposureF`
                Result of subtracting template and science.
            ``matchedTemplate`` : `lsst.afw.image.ExposureF`
                Warped template exposure. Note that in this case, the template
                is not PSF-matched to the science image.
            ``backgroundModel`` : `lsst.afw.math.Chebyshev1Function2D`
                Background model that was fit while solving for the PSF-matching kernel
            ``psfMatchingKernel`` : `lsst.afw.math.Kernel`
               Kernel used to PSF-match the science image to the template.
        """
        if self.config.forceCompatibility:
            # Compatibility option to maintain old behavior
            # This should be removed in the future!
            template = template[science.getBBox()]
        kernelResult = self.makeKernel.run(science, template, sources, preconvolved=False)
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
        matchedTemplate = template.clone()[science.getBBox()]
        matchedTemplate.maskedImage /= norm
        matchedTemplate.setPhotoCalib(science.getPhotoCalib())

        difference = _subtractImages(matchedScience, matchedTemplate,
                                     backgroundModel=(kernelResult.backgroundModel
                                                      if self.config.doSubtractBackground else None))

        correctedExposure = self.finalize(template, science, difference, kernelResult.psfMatchingKernel,
                                          templateMatched=False)

        return lsst.pipe.base.Struct(difference=correctedExposure,
                                     matchedTemplate=matchedTemplate,
                                     matchedScience=matchedScience,
                                     backgroundModel=kernelResult.backgroundModel,
                                     psfMatchingKernel=kernelResult.psfMatchingKernel,)

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
        # Erase existing detection mask planes.
        #  We don't want the detection mask from the science image
        mask = difference.mask
        mask &= ~(mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE"))

        if self.config.doDecorrelation:
            self.log.info("Decorrelating image difference.")
            correctedExposure = self.decorrelate.run(science, template[science.getBBox()], difference, kernel,
                                                     templateMatched=templateMatched,
                                                     preConvMode=preConvMode,
                                                     preConvKernel=preConvKernel,
                                                     spatiallyVarying=spatiallyVarying).correctedExposure
        else:
            self.log.info("NOT decorrelating image difference.")
            correctedExposure = difference
        return correctedExposure

    @staticmethod
    def _validateExposures(template, science):
        """Check that the WCS of the two Exposures match, and the template bbox
        contains the science bbox.

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
            or if the science image is not fully contained in the template
            bounding box.
        """
        assert template.wcs == science.wcs,\
            "Template and science exposure WCS are not identical."
        templateBBox = template.getBBox()
        scienceBBox = science.getBBox()

        assert templateBBox.contains(scienceBBox),\
            "Template bbox does not contain all of the science image."

    @staticmethod
    def _convolveExposure(exposure, kernel, convolutionControl,
                          bbox=None,
                          psf=None,
                          photoCalib=None):
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
        convolvedImage = lsst.afw.image.MaskedImageF(exposure.getBBox())
        lsst.afw.math.convolve(convolvedImage, exposure.maskedImage, kernel, convolutionControl)
        convolvedExposure.setMaskedImage(convolvedImage)
        if bbox is None:
            return convolvedExposure
        else:
            return convolvedExposure[bbox]


def checkTemplateIsSufficient(templateExposure, logger, requiredTemplateFraction=0.):
    """Raise NoWorkFound if template coverage < requiredTemplateFraction

    Parameters
    ----------
    templateExposure : `lsst.afw.image.ExposureF`
        The template exposure to check
    logger : `lsst.log.Log`
        Logger for printing output.
    requiredTemplateFraction : `float`, optional
        Fraction of pixels of the science image required to have coverage
        in the template.

    Raises
    ------
    lsst.pipe.base.NoWorkFound
        Raised if fraction of good pixels, defined as not having NO_DATA
        set, is less then the configured requiredTemplateFraction
    """
    # Count the number of pixels with the NO_DATA mask bit set
    # counting NaN pixels is insufficient because pixels without data are often intepolated over)
    pixNoData = np.count_nonzero(templateExposure.mask.array
                                 & templateExposure.mask.getPlaneBitMask('NO_DATA'))
    pixGood = templateExposure.getBBox().getArea() - pixNoData
    logger.info("template has %d good pixels (%.1f%%)", pixGood,
                100*pixGood/templateExposure.getBBox().getArea())

    if pixGood/templateExposure.getBBox().getArea() < requiredTemplateFraction:
        message = ("Insufficient Template Coverage. (%.1f%% < %.1f%%) Not attempting subtraction. "
                   "To force subtraction, set config requiredTemplateFraction=0." % (
                       100*pixGood/templateExposure.getBBox().getArea(),
                       100*requiredTemplateFraction))
        raise lsst.pipe.base.NoWorkFound(message)


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
