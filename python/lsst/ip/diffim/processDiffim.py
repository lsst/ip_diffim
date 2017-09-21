from __future__ import absolute_import, division, print_function
#
# LSST Data Management System
# Copyright 2016 AURA/LSST.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import random
import math

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.daf.base as dafBase
import lsst.afw.table as afwTable
from lsst.meas.extensions.astrometryNet import LoadAstrometryNetObjectsTask
from lsst.meas.astrom import AstrometryTask
from lsst.meas.algorithms import SourceDetectionTask, LoadIndexedReferenceObjectsTask
from . import MakeDiffimTask, GetCoaddAsTemplateTask, DipoleFitTask

__all__ = ["ProcessDiffimTask", "ProcessDiffimConfig"]


class ProcessDiffimConfig(pexConfig.Config):
    """Config for ProcessDiffimTask
    """
    doAddCalexpBackground = pexConfig.Field(dtype=bool, default=True,
                                            doc="Add background to calexp before processing it.  "
                                                "Useful as ipDiffim does background matching.")
    doMerge = pexConfig.Field(dtype=bool, default=True,
                              doc="Merge positive and negative diaSources with grow radius "
                                  "set by growFootprint")
    doMatchSources = pexConfig.Field(dtype=bool, default=True,
                                     doc="Match diaSources with input calexp sources and ref catalog sources")
    doMeasurement = pexConfig.Field(dtype=bool, default=True, doc="Measure diaSources?")
    doDipoleFitting = pexConfig.Field(dtype=bool, default=True, doc="Measure dipoles using new algorithm?")
    doWriteSources = pexConfig.Field(dtype=bool, default=True, doc="Write sources?")

    coaddName = pexConfig.Field(
        doc="coadd name: typically one of deep or goodSeeing",
        dtype=str,
        default="deep",
    )
    detection = pexConfig.ConfigurableField(
        target=SourceDetectionTask,
        doc="Low-threshold detection for final measurement",
    )
    measurement = pexConfig.ConfigurableField(
        target=DipoleFitTask,
        doc="Enable updated dipole fitting method",
    )
    refObjLoader = pexConfig.ConfigurableField(
        target=LoadAstrometryNetObjectsTask,
        doc="reference object loader",
    )
    astrometer = pexConfig.ConfigurableField(
        target=AstrometryTask,
        doc="astrometry task; used to match sources to reference objects, but not to fit a WCS",
    )
    controlRandomSeed = pexConfig.Field(
        doc = "Random seed for shuffing the control sample",
        dtype = int,
        default = 10
    )

    growFootprint = pexConfig.Field(
        doc="Grow positive and negative footprints by this amount before merging",
        dtype=int,
        default=2
    )

    diaSourceMatchRadius = pexConfig.Field(
        doc="Match radius (in arcseconds) for DiaSource to Source association",
        dtype=float,
        default=0.5
    )
    isPreConvolved = pexConfig.Field(
        doc="Diffim is pre-convolved (doPreConvolve=True in makeDiffim)",
        dtype=bool,
        default=False
    )
    isDecorrelated = pexConfig.Field(
        doc="Diffim is decorrelated (doDecorrelation=True in makeDiffim",
        dtype=bool,
        default=True
    )

    def setDefaults(self):
        # defaults are OK for catalog and diacatalog
        self.doMatchSources = False

        # DiaSource Detection
        self.detection.thresholdPolarity = "both"
        self.detection.thresholdValue = 5.5
        self.detection.reEstimateBackground = False
        self.detection.thresholdType = "pixel_stdev"

        # Add filtered flux measurement, the correct measurement for pre-convolved images.
        # Enable all measurements, regardless of doPreConvolved, as it makes data harvesting easier.
        # To change that you must modify algorithms.names in the task's applyOverrides method,
        # after the user has set doPreConvolved.
        self.measurement.algorithms.names.add('base_PeakLikelihoodFlux')

        # For shuffling the control sample
        random.seed(self.controlRandomSeed)

    def validate(self):
        pexConfig.Config.validate(self)


class ProcessDiffimTaskRunner(pipeBase.ButlerInitializedTaskRunner):

    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return pipeBase.TaskRunner.getTargetList(parsedCmd, **kwargs)

## \addtogroup LSST_task_documentation
## \{
## \page ProcessDiffimTask
## \ref ProcessDiffimTask_ "ProcessDiffimTask"
## \copybrief ProcessDiffimTask
## \}


class ProcessDiffimTask(pipeBase.CmdLineTask):
    """!
\anchor ProcessDiffimTask_

\brief Subtract an image from a template and save the results

\section ip_diffim_processDiffim_Contents Contents

 - \ref ip_diffim_processDiffim_Purpose
 - \ref ip_diffim_processDiffim_Initialize
 - \ref ip_diffim_processDiffim_IO
 - \ref ip_diffim_processDiffim_Config
 - \ref ip_diffim_processDiffim_Example

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_processDiffim_Purpose   Description

This task serves to process the image difference (diffim) produced by the
lsst.ip.diffim.makeDiffim.MakeDiffimTask. This includes detection and measurement
of sources in the diffim, and optionally saving the resulting catalog to disk.

The detection and measurement algorithms are provided by the `detection` and `measurement`
subtasks, respectively.

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_processDiffim_Initialize    Task initialization

\copydoc \_\_init\_\_

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_processDiffim_IO        Invoking the Task

This task may be invoked in two ways:

1. its run() method which expects a Butler sensorRef for loading all relevant data
(template, science image, and diffim).
2. the doProcessDiffim() method which is called from run() but may also be called separately
and expects the template and new exposures and optionally an idFactory and sensorRef for
loading additional data.

See each method's returned lsst.pipe.base.Struct for more details.

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_processDiffim_Config       Configuration parameters

See \ref ProcessDiffimConfig

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

\section ip_diffim_processDiffim_Example   A complete example of using ProcessDiffimTask

This task is called from lsst.ip.diffim.imageDifference.ImageDifferenceTask.run(), which
may be used as an example. This task itself may be invoked via the command line, for
example (assuming obs_decam has been set up and decamDiffimRepo was passed to
makeDiffim.py as the output repo):

\code
processDiffim.py decamDiffimRepo --output decamDiaSrcRepo --id visit=289820 ccdnum=11
\endcode

This will add a new diaSrc record in the decamDiaSrcRepo.
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    """
    ConfigClass = ProcessDiffimConfig
    RunnerClass = ProcessDiffimTaskRunner
    _DefaultName = "processDiffim"

    def __init__(self, schema=None, butler=None, **kwargs):
        """Construct a ProcessDiffim Task

        Parameters
        ----------
        butler : Butler object to use in constructing reference object loaders
        """
        pipeBase.CmdLineTask.__init__(self, **kwargs)

        self.schema = schema
        if schema is None:
            self.schema = afwTable.SourceTable.makeMinimalSchema()

        self.algMetadata = dafBase.PropertyList()
        self.makeSubtask("detection", schema=self.schema)
        if self.config.doMeasurement:
            self.makeSubtask("measurement", schema=self.schema,
                             algMetadata=self.algMetadata)
        if self.config.doMatchSources:
            self.makeSubtask('refObjLoader', butler=butler)
            self.makeSubtask("astrometer", refObjLoader=self.refObjLoader)
            self.schema.addField("refMatchId", "L", "unique id of reference catalog match")
            self.schema.addField("srcMatchId", "L", "unique id of source match")

    @pipeBase.timeMethod
    def run(self, sensorRef):
        """Measure the result of an image subtraction by calling doProcessDiffim

        Parameters
        ----------
        sensorRef : sensor-level butler data reference, used for the following data products:
            Input only:
            - calexp
            - psf
            - ccdExposureId
            - ccdExposureId_bits
            - self.config.coaddName + "Coadd_skyMap"
            - self.config.coaddName + "Coadd"
            - self.config.coaddName + "Diff_differenceExp"
            - self.config.coaddName + "Diff_matchedExp"
            Output, depending on config:
            - self.config.coaddName + "Diff_src"

        Returns
        -------
        pipe_base Struct containing these fields:
            - sources: detected and possibly measured sources
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
        if self.config.doAddCalexpBackground:
            mi = exposure.getMaskedImage()
            mi += sensorRef.get("calexpBackground").getImage()
        if not exposure.hasPsf():
            raise pipeBase.TaskError("Exposure has no psf")

        # Get the subtracted and matched exposures, and `selectSources` for metrics evaluation
        subtractedExposure = sensorRef.get(self.config.coaddName + "Diff_differenceExp")
        matchedExposure = sensorRef.get(self.config.coaddName + "Diff_matchedExp")

        result = self.doProcessDiffim(subtractedExposure, exposure, matchedExposure,
                                      isPreConvolved=self.config.isPreConvolved,
                                      isDecorrelated=self.config.isDecorrelated,
                                      selectSources=None, idFactory=idFactory,
                                      sensorRef=sensorRef)

        return result

    @pipeBase.timeMethod
    def doProcessDiffim(self, subtractedExposure, scienceExposure, matchedTemplateExposure,
                        isPreConvolved=False, isDecorrelated=True, selectSources=None,
                        idFactory=None, sensorRef=None):
        """Measure the result of an image subtraction

        Steps include:
        - detect sources
        - measure sources

        Parameters
        ----------
        subtractedExposure : lsst.afw.image.Exposure
            Subtracted exposure
        scienceExposure : lsst.afw.image.Exposure
            'Science', or new exposure
        templateExposure : lsst.afw.image.Exposure
            Template exposure
        isPreConvolved : bool
            Was the subtractedExposure pre-convolved? (see `makeDiffim.MakeDiffimConfig`)
            This is used to select the type of detection and measurement performed
        isDecorrelated : bool
            Was the subtracted exposure decorrelated? (see `makeDiffim.MakeDiffimConfig`)
            This is used to select the type of dipole measurement performed
        selectSources : lsst.afw.table.SourceCatalog
            Source catalog detected on scienceExposure and selected for A&L PSF
            matching (see `makeDiffim.MakeDiffimTask`). If None, it is read from the
            butler.
        idFactory : lsst.afw.table.Factory
             Factory for the generation of Source ids
        sensorRef :
             Sensor-level butler data reference, used for the following data products:
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

        Returns
        -------
        pipe_base Struct containing these fields:
            - sources: detected and possibly measured sources
        """
        self.log.info("Running diaSource detection")
        # Erase existing detection mask planes
        mask = subtractedExposure.getMaskedImage().getMask()
        mask &= ~(mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE"))

        table = afwTable.SourceTable.make(self.schema, idFactory)
        table.setMetadata(self.algMetadata)
        results = self.detection.makeSourceCatalog(
            table=table,
            exposure=subtractedExposure,
            doSmooth=not isPreConvolved,
        )

        if self.config.doMerge:
            fpSet = results.fpSets.positive
            fpSet.merge(results.fpSets.negative, self.config.growFootprint,
                        self.config.growFootprint, False)
            diaSources = afwTable.SourceCatalog(table)
            fpSet.makeSources(diaSources)
            self.log.info("Merging detections into %d sources" % (len(diaSources)))
        else:
            diaSources = results.sources

        if self.config.doMeasurement:
            newDipoleFitting = self.config.doDipoleFitting
            self.log.info("Running diaSource measurement: newDipoleFitting=%r", newDipoleFitting)
            if not newDipoleFitting:
                # Just fit dipole in diffim
                self.measurement.run(diaSources, subtractedExposure)
            else:
                # Use science image (if avail.) to constrain dipole fitting; matchedTemplate is
                # computed in measurement.run() from scienceExposure - subtractedExposure.
                if matchedTemplateExposure is not None and not isDecorrelated:
                    self.measurement.run(diaSources, subtractedExposure, scienceExposure,
                                         matchedTemplateExposure)
                else:
                    # If decorrelation is turned on, we don't want to use the matchedTemplate, since it
                    # is not decorrelated. If we don't pass the matchedTemplate, then the dipoleFitter
                    # automatically creates from scienceExposure - subtractedExposure.
                    self.measurement.run(diaSources, subtractedExposure, scienceExposure)

        if self.config.doMatchSources and sensorRef is not None:
            # Match with the calexp sources if possible
            selectSources = MakeDiffimTask.getSelectSources(scienceExposure, self.config.isPreConvolved,
                                                            self.log, idFactory, sensorRef)
            refObjLoader = LoadIndexedReferenceObjectsTask(butler=sensorRef.getButler())
            diaSources = self.runMatchSources(diaSources, selectSources, scienceExposure,
                                              refObjLoader)

        if self.config.doWriteSources and diaSources is not None and sensorRef is not None:
            self.log.info("Writing diaSources, %s" %
                          (self.config.coaddName + "Diff_diaSrc"))
            sensorRef.put(diaSources, self.config.coaddName + "Diff_diaSrc")

        return pipeBase.Struct(
            sources=diaSources
        )

    def runMatchSources(self, diaSources, selectSources, exposure, refObjLoader):
        """Match diaSources and selectSources catalogs by RA/Dec

        Parameters
        ----------
        diaSources : lsst.afw.table.SourceCatalog
           Catalog of sources detected on diffim
        selectSources : lsst.afw.table.SourceCatalog
           Catalog of sources detected on science image and used for PSF matching

        """
        matchRadAsec = self.config.diaSourceMatchRadius
        if selectSources is not None:
            # Create key,val pair where key=diaSourceId and val=sourceId
            matchRadPixel = matchRadAsec / exposure.getWcs().pixelScale().asArcseconds()

            srcMatches = afwTable.matchXy(selectSources, diaSources, matchRadPixel)
            srcMatchDict = dict([(srcMatch.second.getId(), srcMatch.first.getId()) for
                                 srcMatch in srcMatches])
            self.log.info("Matched %d / %d diaSources to sources" % (len(srcMatchDict),
                                                                     len(diaSources)))
        else:
            self.log.warn("Src product does not exist; cannot match with diaSources")
            srcMatchDict = {}

        # Create key,val pair where key=diaSourceId and val=refId
        self.astrometer.matcher.maxMatchDistArcSec = matchRadAsec
        astromRet = self.astrometer.loadAndMatch(exposure=exposure, sourceCat=selectSources)
        refMatches = astromRet.matches
        if refMatches is None:
            self.log.warn("No diaSource matches with reference catalog")
            refMatchDict = {}
        else:
            self.log.info("Matched %d / %d diaSources to reference catalog" % (len(refMatches),
                                                                               len(diaSources)))
            refMatchDict = dict([(refMatch.second.getId(), refMatch.first.getId()) for
                                 refMatch in refMatches])

        # Assign source Ids
        for diaSource in diaSources:
            sid = diaSource.getId()
            if sid in srcMatchDict:
                diaSource.set("srcMatchId", srcMatchDict[sid])
            if sid in refMatchDict:
                diaSource.set("refMatchId", refMatchDict[sid])

        return diaSources

    def _getConfigName(self):
        """Return the name of the config dataset
        """
        return "%sDiff_config" % (self.config.coaddName,)

    def _getMetadataName(self):
        """Return the name of the metadata dataset
        """
        return "%sDiff_metadata" % (self.config.coaddName,)

    def getSchemaCatalogs(self):
        """Return a dict of empty catalogs for each catalog dataset produced by this task."""
        diaSrc = afwTable.SourceCatalog(self.schema)
        diaSrc.getTable().setMetadata(self.algMetadata)
        return {self.config.coaddName + "Diff_diaSrc": diaSrc}

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "calexp", help="data ID, e.g. --id visit=12345 ccd=1,2")
        return parser
