#
# LSST Data Management System
# Copyright 2012 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
from __future__ import absolute_import, division, print_function
from builtins import zip
import math
import random
import numpy

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
from lsst.meas.astrom import AstrometryTask
from lsst.meas.extensions.astrometryNet import LoadAstrometryNetObjectsTask
from lsst.meas.algorithms import SingleGaussianPsf, \
    ObjectSizeStarSelectorTask
from . import ImagePsfMatchTask, ZogyImagePsfMatchTask, makeKernelBasisList, \
    KernelCandidateQa, DiaCatalogSourceSelectorTask, DiaCatalogSourceSelectorConfig, \
    RegisterTask, GetCoaddAsTemplateTask, GetCalexpAsTemplateTask, \
    DecorrelateALKernelSpatialTask

__all__ = ["MakeDiffimTask", "MakeDiffimConfig"]

FwhmPerSigma = 2 * math.sqrt(2 * math.log(2))
IqrToSigma = 0.741

subtractAlgorithmRegistry = pexConfig.makeRegistry(
    doc="A registry of subtraction algorithms for use as a subtask in imageDifference",
)

subtractAlgorithmRegistry.register('al', ImagePsfMatchTask)
subtractAlgorithmRegistry.register('zogy', ZogyImagePsfMatchTask)


class MakeDiffimConfig(pexConfig.Config):
    """Config for MakeDiffimTask
    """
    doAddCalexpBackground = pexConfig.Field(dtype=bool, default=True,
                                            doc="Add background to calexp before processing it. "
                                                "Useful as ipDiffim does background matching.")
    doUseRegister = pexConfig.Field(dtype=bool, default=True,
                                    doc="Use image-to-image registration to align template with "
                                        "science image")
    doDebugRegister = pexConfig.Field(dtype=bool, default=False,
                                      doc="Writing debugging data for doUseRegister")
    doSelectSources = pexConfig.Field(dtype=bool, default=True,
                                      doc="Select stars to use for kernel fitting")
    doSelectDcrCatalog = pexConfig.Field(dtype=bool, default=False,
                                         doc="Select stars of extreme color as part of the control sample")
    doSelectVariableCatalog = pexConfig.Field(dtype=bool, default=False,
                                              doc="Select stars that are variable to be part "
                                                  "of the control sample")
    doPreConvolve = pexConfig.Field(dtype=bool, default=True,
                                    doc="Convolve science image by its PSF before PSF-matching?")
    useGaussianForPreConvolution = pexConfig.Field(dtype=bool, default=True,
                                                   doc="Use a simple gaussian PSF model for pre-convolution "
                                                       "(else use fit PSF)? Ignored if doPreConvolve false.")
    doDecorrelation = pexConfig.Field(dtype=bool, default=False,
                                      doc="Perform diffim decorrelation to undo pixel correlation due to A&L "
                                          "kernel convolution? If True, also update the diffim PSF.")
    doWriteSubtractedExp = pexConfig.Field(dtype=bool, default=True, doc="Write difference exposure?")
    doWriteMatchedExp = pexConfig.Field(dtype=bool, default=False,
                                        doc="Write warped and PSF-matched template coadd exposure?")
    doAddMetrics = pexConfig.Field(dtype=bool, default=True,
                                   doc="Add columns to the source table to hold analysis metrics?")

    coaddName = pexConfig.Field(
        doc="coadd name: typically one of deep or goodSeeing",
        dtype=str,
        default="deep",
    )
    convolveTemplate = pexConfig.Field(
        doc="Which image gets convolved (default = template)",
        dtype=bool,
        default=True
    )
    refObjLoader = pexConfig.ConfigurableField(
        target=LoadAstrometryNetObjectsTask,
        doc="reference object loader",
    )
    astrometer = pexConfig.ConfigurableField(
        target=AstrometryTask,
        doc="astrometry task; used to match sources to reference objects, but not to fit a WCS",
    )
    sourceSelector = pexConfig.ConfigurableField(
        target=ObjectSizeStarSelectorTask,
        doc="Source selection algorithm",
    )
    subtract = subtractAlgorithmRegistry.makeField("Subtraction Algorithm", default="al")
    decorrelate = pexConfig.ConfigurableField(
        target=DecorrelateALKernelSpatialTask,
        doc="Decorrelate effects of A&L kernel convolution on image difference, only if doSubtract is True. "
        "If this option is enabled, then detection.thresholdValue should be set to 5.0 (rather than the "
        "default of 5.5).",
    )
    doSpatiallyVarying = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="If using Zogy or A&L decorrelation, perform these on a grid across the "
        "image in order to allow for spatial variations"
    )
    getTemplate = pexConfig.ConfigurableField(
        target=GetCoaddAsTemplateTask,
        doc="Subtask to retrieve template exposure and sources",
    )
    controlStepSize = pexConfig.Field(
        doc="What step size (every Nth one) to select a control sample from the kernelSources",
        dtype=int,
        default=5
    )
    controlRandomSeed = pexConfig.Field(
        doc = "Random seed for shuffing the control sample",
        dtype = int,
        default = 10
    )
    register = pexConfig.ConfigurableField(
        target=RegisterTask,
        doc="Task to enable image-to-image image registration (warping)",
    )
    kernelSourcesFromRef = pexConfig.Field(
        doc="Select sources to measure kernel from reference catalog if True, template if false",
        dtype=bool,
        default=False
    )
    templateSipOrder = pexConfig.Field(dtype=int, default=2,
                                       doc="Sip Order for fitting the Template Wcs "
                                           "(default is too high, overfitting)")

    def setDefaults(self):
        # defaults are OK for catalog and diacatalog

        self.subtract['al'].kernel.name = "AL"
        self.subtract['al'].kernel.active.fitForBackground = True
        self.subtract['al'].kernel.active.spatialKernelOrder = 1
        self.subtract['al'].kernel.active.spatialBgOrder = 0
        self.doPreConvolve = False
        self.doAddMetrics = False
        self.doUseRegister = False

        # For shuffling the control sample
        random.seed(self.controlRandomSeed)

    def validate(self):
        pexConfig.Config.validate(self)
        if self.doUseRegister and not self.doSelectSources:
            raise ValueError("doUseRegister=True and doSelectSources=False. " +
                             "Cannot run RegisterTask without selecting sources.")


class MakeDiffimTaskRunner(pipeBase.ButlerInitializedTaskRunner):

    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return pipeBase.TaskRunner.getTargetList(parsedCmd, templateIdList=parsedCmd.templateId.idList,
                                                 **kwargs)


class MakeDiffimTask(pipeBase.CmdLineTask):
    """Subtract an image from a template and measure the result
    """
    ConfigClass = MakeDiffimConfig
    RunnerClass = MakeDiffimTaskRunner
    _DefaultName = "makeDiffim"

    def __init__(self, butler=None, **kwargs):
        """!Construct an MakeDiffim Task

        @param[in] butler  Butler object to use in constructing reference object loaders
        """
        pipeBase.CmdLineTask.__init__(self, **kwargs)
        self.makeSubtask("getTemplate")

        self.makeSubtask("subtract")

        if self.config.subtract.name == 'al' and self.config.doDecorrelation:
            self.makeSubtask("decorrelate")

        if self.config.doUseRegister:
            self.makeSubtask("register")
        self.schema = afwTable.SourceTable.makeMinimalSchema()

        if self.config.doSelectSources:
            self.makeSubtask("sourceSelector", schema=self.schema)
            self.makeSubtask('refObjLoader', butler=butler)
            self.makeSubtask("astrometer", refObjLoader=self.refObjLoader)

        self.algMetadata = dafBase.PropertyList()

    @pipeBase.timeMethod
    def run(self, sensorRef, templateIdList=None):
        """Subtract an image from a template coadd

        Steps include:
        - warp template coadd to match WCS of image
        - PSF match image to warped template
        - subtract image from PSF-matched, warped template
        - persist difference image

        @param sensorRef: sensor-level butler data reference, used for the following data products:
        Input only:
        - calexp
        - psf
        - ccdExposureId
        - ccdExposureId_bits
        - self.config.coaddName + "Coadd_skyMap"
        - self.config.coaddName + "Coadd"
        Input or output, depending on config:
        - self.config.coaddName + "Diff_differenceExp"
        Output, depending on config:
        - self.config.coaddName + "Diff_matchedExp"

        @return pipe_base Struct containing these fields:
        - subtractedExposure: exposure after subtracting template;
            the unpersisted version if subtraction not run but detection run
            None if neither subtraction nor detection run (i.e. nothing useful done)
        - subtractRes: results of subtraction task; None if subtraction not run
        """
        self.log.info("Processing %s" % (sensorRef.dataId))

        # We make one IdFactory that will be used by both icSrc and src datasets;
        # I don't know if this is the way we ultimately want to do things, but at least
        # this ensures the source IDs are fully unique.
        expBits = sensorRef.get("ccdExposureId_bits")
        expId = int(sensorRef.get("ccdExposureId"))
        idFactory = afwTable.IdFactory.makeSource(expId, 64 - expBits)

        # Retrieve the science image we wish to analyze
        exposure = sensorRef.get("calexp", immediate=True)
        template = self.getTemplate.run(exposure, sensorRef, templateIdList=templateIdList)

        result = self.doMakeDiffim(template, exposure, idFactory=idFactory,
                                   sensorRef=sensorRef)

        return result

    @pipeBase.timeMethod
    def doMakeDiffim(self, template, exposure, idFactory=None, sensorRef=None):
        templateExposure = template.exposure   # Stitched coadd exposure
        templateSources = template.sources    # Sources on the template image

        if not exposure.hasPsf():
            raise pipeBase.TaskError("Exposure has no psf")
        sciencePsf = exposure.getPsf()

        # compute scienceSigmaOrig: sigma of PSF of science image
        scienceSigmaOrig = sciencePsf.computeShape().getDeterminantRadius()
        scienceSigmaPost = scienceSigmaOrig
        # compute scienceSigmaPost: sigma of PSF of science image after pre-convolution (A&L)
        if self.config.doPreConvolve:
            scienceSigmaPost = scienceSigmaOrig * math.sqrt(2)

        # If requested, find sources in the image
        selectSourceRes = None
        selectSources = None
        kernelSources = None
        matches = None
        if self.config.doSelectSources:
            selectSourcesRes = self.runSelectSources(exposure, templateExposure,
                                                     templateSources, scienceSigmaPost,
                                                     idFactory, sensorRef=sensorRef)
            selectSources = selectSourcesRes.selectSources
            kernelSources = selectSourcesRes.kernelSources
            matches = selectSourcesRes.matches

        allresids = {}
        if self.config.doUseRegister:
            res = self.runRegister(exposure, templateExposure, templateSources,
                                   selectSources, idFactory=None)
            wcsResults = res.wcsResults
            templateExposure = res.templateExposure
            # Create debugging outputs on the astrometric
            # residuals as a function of position.  Persistence
            # not yet implemented; expected on (I believe) #2636.
            if self.config.doDebugRegister:
                allresids = self.runDebugRegister(wcsResults, matches)

        if self.config.subtract.name == 'zogy':
            spatiallyVarying = self.config.doSpatiallyVarying
            subtractRes = self.subtract.subtractExposures(templateExposure, exposure,
                                                          doWarping=True,
                                                          spatiallyVarying=spatiallyVarying,
                                                          doPreConvolve=self.config.doPreConvolve)
            subtractedExposure = subtractRes.subtractedExposure

        elif self.config.subtract.name == 'al':
            # Useful since AL can match backgrounds (Zogy cannot)
            if self.config.doAddCalexpBackground:
                mi = exposure.getMaskedImage()
                mi += sensorRef.get("calexpBackground").getImage()

            # if requested, convolve the science exposure with its PSF
            # (properly, this should be a cross-correlation, but our code does not yet support that)
            # compute scienceSigmaPost: sigma of science exposure with pre-convolution, if done,
            # else sigma of original science exposure
            preConvKernel = None
            if self.config.doPreConvolve:
                convControl = afwMath.ConvolutionControl()
                # cannot convolve in place, so make a new MI to receive convolved image
                srcMI = exposure.getMaskedImage()
                destMI = srcMI.Factory(srcMI.getDimensions())
                srcPsf = sciencePsf
                if self.config.useGaussianForPreConvolution:
                    # convolve with a simplified PSF model: a double Gaussian
                    kWidth, kHeight = sciencePsf.getLocalKernel().getDimensions()
                    preConvPsf = SingleGaussianPsf(kWidth, kHeight, scienceSigmaOrig)
                else:
                    # convolve with science exposure's PSF model
                    preConvPsf = srcPsf
                afwMath.convolve(destMI, srcMI, preConvPsf.getLocalKernel(), convControl)
                exposure.setMaskedImage(destMI)

            # warp template exposure to match exposure,
            # PSF match template exposure to exposure,
            # then return the difference

            # Return warped template...
            self.log.info("Subtracting images")
            subtractRes = self.subtract.subtractExposures(
                templateExposure=templateExposure,
                scienceExposure=exposure,
                candidateList=kernelSources,
                convolveTemplate=self.config.convolveTemplate,
                doWarping=not self.config.doUseRegister
            )
            subtractedExposure = subtractRes.subtractedExposure

            # Add the selectSources result to the final result Struct.
            if selectSourceRes is not None:
                subtractRes.selectSourceResult = selectSourceRes
            subtractRes.allresids = allresids

            self.log.info("Computing diffim PSF")

            # Get Psf from the appropriate input image if it doesn't exist
            if not subtractedExposure.hasPsf():
                if self.config.convolveTemplate:
                    subtractedExposure.setPsf(exposure.getPsf())
                else:
                    subtractedExposure.setPsf(template.exposure.getPsf())

            # Perform diffim decorrelation
            if self.config.doDecorrelation:
                spatiallyVarying = self.config.doSpatiallyVarying
                preConvKernel = None
                if preConvPsf is not None:
                    preConvKernel = preConvPsf.getLocalKernel()
                decorrResult = self.decorrelate.run(exposure, templateExposure,
                                                    subtractedExposure,
                                                    subtractRes.psfMatchingKernel,
                                                    spatiallyVarying=spatiallyVarying,
                                                    preConvKernel=preConvKernel)
                subtractedExposure = decorrResult.correctedExposure
                subtractRes.subtractedExposure = subtractedExposure

        if sensorRef is not None:
            if self.config.doWriteMatchedExp and 'matchedImage' in subtractRes.getDict():
                self.log.info("Writing matched exposure, %s" %
                              (self.config.coaddName + "Diff_matchedExp"))
                sensorRef.put(subtractRes.matchedImage, self.config.coaddName + "Diff_matchedExp")

            if self.config.doWriteSubtractedExp and 'subtractedExposure' in subtractRes.getDict():
                self.log.info("Writing subtracted exposure, %s" %
                              (self.config.coaddName + "Diff_differenceExp"))
                subtractedExposureName = self.config.coaddName + "Diff_differenceExp"
                sensorRef.put(subtractRes.subtractedExposure, subtractedExposureName)

        return subtractRes

    def runSelectSources(self, exposure, templateExposure, templateSources,
                         scienceSigmaPost, idFactory=None, sensorRef=None):
        if sensorRef is not None and sensorRef.datasetExists("src"):
            self.log.info("Source selection via src product")
            selectSources = sensorRef.get("src")
        else:
            self.log.warn("Src product does not exist; running detection, measurement, selection")
            # Run own detection and measurement; necessary in nightly processing
            selectSources = self.subtract.getSelectSources(
                exposure,
                sigma=scienceSigmaPost,
                doSmooth=not self.config.doPreConvolve,
                idFactory=idFactory,
            )

        kcQa = nparam = None
        if self.config.doAddMetrics:
            # sigma of PSF of template image before warping
            templateSigma = templateExposure.getPsf().computeShape().getDeterminantRadius()
            # Number of basis functions
            nparam = len(makeKernelBasisList(self.subtract.config.kernel.active,
                                             referenceFwhmPix=scienceSigmaPost * FwhmPerSigma,
                                             targetFwhmPix=templateSigma * FwhmPerSigma))
            # Modify the schema of all Sources
            kcQa = KernelCandidateQa(nparam)
            selectSources = kcQa.addToSchema(selectSources)

        if self.config.kernelSourcesFromRef:
            # match exposure sources to reference catalog
            astromRet = self.astrometer.loadAndMatch(exposure=exposure, sourceCat=selectSources)
            matches = astromRet.matches
        elif templateSources is not None:
            # match exposure sources to template sources
            mc = afwTable.MatchControl()
            mc.findOnlyClosest = False
            matches = afwTable.matchRaDec(templateSources, selectSources, 1.0*afwGeom.arcseconds,
                                          mc)
        else:
            raise RuntimeError("doSelectSources=True and kernelSourcesFromRef=False," +
                               "but template sources not available. Cannot match science " +
                               "sources with template sources. Run process* on data from " +
                               "which templates are built.")

        kernelSources = self.sourceSelector.selectStars(exposure, selectSources,
                                                        matches=matches).starCat

        random.shuffle(kernelSources, random.random)
        controlSources = kernelSources[::self.config.controlStepSize]
        kernelSources = [k for i, k in enumerate(kernelSources) if i % self.config.controlStepSize]

        if self.config.doSelectDcrCatalog:
            redSelector = DiaCatalogSourceSelectorTask(
                DiaCatalogSourceSelectorConfig(grMin=self.sourceSelector.config.grMax, grMax=99.999))
            redSources = redSelector.selectStars(exposure, selectSources, matches=matches).starCat
            controlSources.extend(redSources)

            blueSelector = DiaCatalogSourceSelectorTask(
                DiaCatalogSourceSelectorConfig(grMin=-99.999, grMax=self.sourceSelector.config.grMin))
            blueSources = blueSelector.selectStars(exposure, selectSources, matches=matches).starCat
            controlSources.extend(blueSources)

        if self.config.doSelectVariableCatalog:
            varSelector = DiaCatalogSourceSelectorTask(
                DiaCatalogSourceSelectorConfig(includeVariable=True))
            varSources = varSelector.selectStars(exposure, selectSources, matches=matches).starCat
            controlSources.extend(varSources)

        self.log.info("Selected %d / %d sources for Psf matching (%d for control sample)"
                      % (len(kernelSources), len(selectSources), len(controlSources)))

        return pipeBase.Struct(selectSources=selectSources,
                               kernelSources=kernelSources,
                               controlSources=controlSources,
                               matches=matches,
                               kcQa=kcQa, nparam=nparam)

    def runRegister(self, exposure, templateExposure, templateSources, selectSources, idFactory=None):
        self.log.info("Registering images")

        if templateSources is None:
            # sigma of PSF of template image before warping
            templateSigma = templateExposure.getPsf().computeShape().getDeterminantRadius()
            # Run detection on the template, which is
            # temporarily background-subtracted
            templateSources = self.subtract.getSelectSources(
                templateExposure,
                sigma=templateSigma,
                doSmooth=True,
                idFactory=idFactory
            )

        # Third step: we need to fit the relative astrometry.
        wcsResults = self.fitAstrometry(templateSources, templateExposure, selectSources)
        warpedExp = self.register.warpExposure(templateExposure, wcsResults.wcs,
                                               exposure.getWcs(), exposure.getBBox())
        # Psf seems to get removed by the registration, so re-add it.
        templPsf = templateExposure.getPsf()
        templateExposure = warpedExp
        templateExposure.setPsf(templPsf)
        return pipeBase.Struct(wcsResults=wcsResults,
                               templateExposure=templateExposure)

    def runDebugRegister(self, wcsResults, matches):
        # Create debugging outputs on the astrometric
        # residuals as a function of position.  Persistence
        # not yet implemented; expected on (I believe) #2636.

        # Grab matches to reference catalog
        srcToMatch = {x.second.getId(): x.first for x in matches}

        refCoordKey = wcsResults.matches[0].first.getTable().getCoordKey()
        inCentroidKey = wcsResults.matches[0].second.getTable().getCentroidKey()
        sids = [m.first.getId() for m in wcsResults.matches]
        positions = [m.first.get(refCoordKey) for m in wcsResults.matches]
        residuals = [m.first.get(refCoordKey).getOffsetFrom(wcsResults.wcs.pixelToSky(
            m.second.get(inCentroidKey))) for m in wcsResults.matches]
        allresids = dict(zip(sids, zip(positions, residuals)))

        cresiduals = [m.first.get(refCoordKey).getTangentPlaneOffset(
            wcsResults.wcs.pixelToSky(
                m.second.get(inCentroidKey))) for m in wcsResults.matches]
        colors = numpy.array([-2.5*numpy.log10(srcToMatch[x].get("g")) +
                              2.5*numpy.log10(srcToMatch[x].get("r"))
                              for x in sids if x in srcToMatch.keys()])
        dlong = numpy.array([r[0].asArcseconds() for s, r in zip(sids, cresiduals)
                             if s in srcToMatch.keys()])
        dlat = numpy.array([r[1].asArcseconds() for s, r in zip(sids, cresiduals)
                            if s in srcToMatch.keys()])
        idx1 = numpy.where(colors < self.sourceSelector.config.grMin)
        idx2 = numpy.where((colors >= self.sourceSelector.config.grMin) &
                           (colors <= self.sourceSelector.config.grMax))
        idx3 = numpy.where(colors > self.sourceSelector.config.grMax)
        rms1Long = IqrToSigma * \
            (numpy.percentile(dlong[idx1], 75)-numpy.percentile(dlong[idx1], 25))
        rms1Lat = IqrToSigma*(numpy.percentile(dlat[idx1], 75)-numpy.percentile(dlat[idx1], 25))
        rms2Long = IqrToSigma * \
            (numpy.percentile(dlong[idx2], 75)-numpy.percentile(dlong[idx2], 25))
        rms2Lat = IqrToSigma*(numpy.percentile(dlat[idx2], 75)-numpy.percentile(dlat[idx2], 25))
        rms3Long = IqrToSigma * \
            (numpy.percentile(dlong[idx3], 75)-numpy.percentile(dlong[idx3], 25))
        rms3Lat = IqrToSigma*(numpy.percentile(dlat[idx3], 75)-numpy.percentile(dlat[idx3], 25))
        self.log.info("Blue star offsets'': %.3f %.3f, %.3f %.3f" % (numpy.median(dlong[idx1]),
                                                                     rms1Long,
                                                                     numpy.median(dlat[idx1]),
                                                                     rms1Lat))
        self.log.info("Green star offsets'': %.3f %.3f, %.3f %.3f" % (numpy.median(dlong[idx2]),
                                                                      rms2Long,
                                                                      numpy.median(dlat[idx2]),
                                                                      rms2Lat))
        self.log.info("Red star offsets'': %.3f %.3f, %.3f %.3f" % (numpy.median(dlong[idx3]),
                                                                    rms3Long,
                                                                    numpy.median(dlat[idx3]),
                                                                    rms3Lat))

        self.metadata.add("RegisterBlueLongOffsetMedian", numpy.median(dlong[idx1]))
        self.metadata.add("RegisterGreenLongOffsetMedian", numpy.median(dlong[idx2]))
        self.metadata.add("RegisterRedLongOffsetMedian", numpy.median(dlong[idx3]))
        self.metadata.add("RegisterBlueLongOffsetStd", rms1Long)
        self.metadata.add("RegisterGreenLongOffsetStd", rms2Long)
        self.metadata.add("RegisterRedLongOffsetStd", rms3Long)

        self.metadata.add("RegisterBlueLatOffsetMedian", numpy.median(dlat[idx1]))
        self.metadata.add("RegisterGreenLatOffsetMedian", numpy.median(dlat[idx2]))
        self.metadata.add("RegisterRedLatOffsetMedian", numpy.median(dlat[idx3]))
        self.metadata.add("RegisterBlueLatOffsetStd", rms1Lat)
        self.metadata.add("RegisterGreenLatOffsetStd", rms2Lat)
        self.metadata.add("RegisterRedLatOffsetStd", rms3Lat)
        return allresids

    def fitAstrometry(self, templateSources, templateExposure, selectSources):
        """Fit the relative astrometry between templateSources and selectSources

        @todo remove this method. It originally fit a new WCS to the template before calling register.run
        because our TAN-SIP fitter behaved badly for points far from CRPIX, but that's been fixed.
        It remains because a subtask overrides it.
        """
        results = self.register.run(templateSources, templateExposure.getWcs(),
                                    templateExposure.getBBox(), selectSources)
        return results

    def _getConfigName(self):
        """Return the name of the config dataset
        """
        return "%sDiff_config" % (self.config.coaddName,)

    def _getMetadataName(self):
        """Return the name of the metadata dataset
        """
        return "%sDiff_metadata" % (self.config.coaddName,)

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "calexp", help="data ID, e.g. --id visit=12345 ccd=1,2")
        parser.add_id_argument("--templateId", "calexp", doMakeDataRefList=True,
                               help="Optional template data ID (visit only), e.g. --templateId visit=6789")
        return parser


class Winter2013MakeDiffimConfig(MakeDiffimConfig):
    winter2013WcsShift = pexConfig.Field(dtype=float, default=0.0,
                                         doc="Shift stars going into RegisterTask by this amount")
    winter2013WcsRms = pexConfig.Field(dtype=float, default=0.0,
                                       doc="Perturb stars going into RegisterTask by this amount")

    def setDefaults(self):
        MakeDiffimConfig.setDefaults(self)
        self.getTemplate.retarget(GetCalexpAsTemplateTask)


class Winter2013MakeDiffimTask(MakeDiffimTask):
    """!Image difference Task used in the Winter 2013 data challege.
    Enables testing the effects of registration shifts and scatter.

    For use with winter 2013 simulated images:
    Use --templateId visit=88868666 for sparse data
        --templateId visit=22222200 for dense data (g)
        --templateId visit=11111100 for dense data (i)
    """
    ConfigClass = Winter2013MakeDiffimConfig
    _DefaultName = "winter2013MakeDiffim"

    def __init__(self, **kwargs):
        MakeDiffimTask.__init__(self, **kwargs)

    def fitAstrometry(self, templateSources, templateExposure, selectSources):
        """Fit the relative astrometry between templateSources and selectSources"""
        if self.config.winter2013WcsShift > 0.0:
            offset = afwGeom.Extent2D(self.config.winter2013WcsShift,
                                      self.config.winter2013WcsShift)
            cKey = templateSources[0].getTable().getCentroidKey()
            for source in templateSources:
                centroid = source.get(cKey)
                source.set(cKey, centroid+offset)
        elif self.config.winter2013WcsRms > 0.0:
            cKey = templateSources[0].getTable().getCentroidKey()
            for source in templateSources:
                offset = afwGeom.Extent2D(self.config.winter2013WcsRms*numpy.random.normal(),
                                          self.config.winter2013WcsRms*numpy.random.normal())
                centroid = source.get(cKey)
                source.set(cKey, centroid+offset)

        results = self.register.run(templateSources, templateExposure.getWcs(),
                                    templateExposure.getBBox(), selectSources)
        return results
