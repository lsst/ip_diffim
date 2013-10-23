#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import numpy as np
import lsst.daf.base as dafBase
import lsst.pex.logging as pexLog
import lsst.pex.config as pexConfig
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.detection as afwDetect
import lsst.pipe.base as pipeBase
from lsst.meas.algorithms import (SourceDetectionTask, SourceMeasurementTask,
                                  getBackground, BackgroundConfig)
from .makeKernelBasisList import makeKernelBasisList
from .psfMatch import PsfMatch, PsfMatchConfigDF, PsfMatchConfigAL
from . import utils as diUtils 
from . import diffimLib
from . import diffimTools
import lsst.afw.display.ds9 as ds9

sigma2fwhm = 2. * np.sqrt(2. * np.log(2.))

class ImagePsfMatchConfig(pexConfig.Config):
    kernel = pexConfig.ConfigChoiceField(
        doc="kernel type",
        typemap=dict(
            AL=PsfMatchConfigAL,
            DF=PsfMatchConfigDF
        ),
        default="AL",
    )
    selectDetection = pexConfig.ConfigurableField(
        target=SourceDetectionTask,
        doc="Initial detections used to feed stars to kernel fitting",
    )
    selectMeasurement = pexConfig.ConfigurableField(
        target=SourceMeasurementTask,
        doc="Initial measurements used to feed stars to kernel fitting",
    )

    def setDefaults(self):
        # High sigma detections only
        self.selectDetection.reEstimateBackground = False
        self.selectDetection.thresholdValue = 10.0

        # Minimal set of measurments for star selection
        self.selectMeasurement.algorithms.names.clear()
        self.selectMeasurement.algorithms.names = ('flux.psf', 'flags.pixel', 'shape.sdss',
                                                   'flux.gaussian', 'skycoord')
        self.selectMeasurement.slots.modelFlux = None
        self.selectMeasurement.slots.apFlux = None

class ImagePsfMatchTask(PsfMatch):
    """PSF-match images to reference images

    Fits the following model:
    image to not convolve = (image to convolve convolved with PSF matching kernel) + background model
    """
    ConfigClass = ImagePsfMatchConfig

    def __init__(self, *args, **kwargs):
        """Create a PsfMatchToImage

        @param config: see lsst.ip.diffim.PsfMatchConfig
        @param logName: name by which messages are logged
        """
        PsfMatch.__init__(self, *args, **kwargs)
        self.kconfig = self.config.kernel.active
        self._warper = afwMath.Warper.fromConfig(self.kconfig.warpingConfig)
        self.selectSchema = afwTable.SourceTable.makeMinimalSchema()
        self.selectAlgMetadata = dafBase.PropertyList()
        self.makeSubtask("selectDetection", schema=self.selectSchema)
        self.makeSubtask("selectMeasurement", schema=self.selectSchema, algMetadata=self.selectAlgMetadata)

    @pipeBase.timeMethod
    def matchExposures(self, templateExposure, scienceExposure,
                       templateFwhmPix=None, scienceFwhmPix=None,
                       candidateList=None, doWarping=True, convolveTemplate=True):
        """Warp and PSF-match an exposure to the reference

        Do the following, in order:
        - Warp templateExposure to match scienceExposure,
            if doWarping True and their WCSs do not already match
        - Determine a PSF matching kernel and differential background model
            that matches templateExposure to scienceExposure
        - Convolve templateExposure by PSF matching kernel

        @param templateExposure: Exposure to warp and PSF-match to the reference masked image
        @param scienceExposure: Exposure whose WCS and PSF are to be matched to
        @param templateFwhmPix: FWHM (in pixels) of the Psf in the template image (image to convolve)
        @param scienceFwhmPix: FWHM (in pixels) of the Psf in the science image
        @param candidateList: a list of footprints/maskedImages for kernel candidates; 
                              if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF
        @param doWarping: what to do if templateExposure's and scienceExposure's WCSs do not match:
            - if True then warp templateExposure to match scienceExposure
            - if False then raise an Exception
        @param convolveTemplate: convolve the template image or the science image
            - if True, templateExposure is warped if doWarping, templateExposure is convolved
            - if False, templateExposure is warped if doWarping, scienceExposure is convolved

        @return a pipeBase.Struct containing these fields:
        - matchedImage: the PSF-matched exposure =
            warped templateExposure convolved by psfMatchingKernel. This has:
            - the same parent bbox, Wcs and Calib as scienceExposure
            - the same filter as templateExposure
            - no Psf (because the PSF-matching process does not compute one)
        - psfMatchingKernel: the PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to solve for the PSF matching kernel

        @raise RuntimeError if doWarping is False and templateExposure's and scienceExposure's
            WCSs do not match
        """
        if not self._validateWcs(templateExposure, scienceExposure):
            if doWarping:
                self.log.info("Astrometrically registering template to science image")
                templateExposure = self._warper.warpExposure(scienceExposure.getWcs(),
                    templateExposure, destBBox=scienceExposure.getBBox(afwImage.PARENT))
            else:
                pexLog.Trace(self.log.getName(), 1, "ERROR: Input images not registered")
                raise RuntimeError("Input images not registered")
        if templateFwhmPix is None:
            if not templateExposure.hasPsf():
                self.log.warn("No estimate of Psf FWHM for template image")
            else:
                psf = templateExposure.getPsf()
                templateSigPix = psf.computeShape().getDeterminantRadius()
                templateFwhmPix = templateSigPix * sigma2fwhm
        if scienceFwhmPix is None:
            if not scienceExposure.hasPsf():
                self.log.warn("No estimate of Psf FWHM for science image")
            else:
                psf = scienceExposure.getPsf()
                scienceSigPix = psf.computeShape().getDeterminantRadius()
                scienceFwhmPix = scienceSigPix * sigma2fwhm

        kernelSize = makeKernelBasisList(self.kconfig, templateFwhmPix, scienceFwhmPix)[0].getWidth()
        candidateList = self.makeCandidateList(templateExposure, scienceExposure, kernelSize, candidateList)

        if convolveTemplate:
            results = self.matchMaskedImages(
                templateExposure.getMaskedImage(), scienceExposure.getMaskedImage(), candidateList,
                templateFwhmPix=templateFwhmPix, scienceFwhmPix=scienceFwhmPix)
        else:
            results = self.matchMaskedImages(
                scienceExposure.getMaskedImage(), templateExposure.getMaskedImage(), candidateList,
                templateFwhmPix=scienceFwhmPix, scienceFwhmPix=templateFwhmPix)

        psfMatchedExposure = afwImage.makeExposure(results.matchedImage, scienceExposure.getWcs())
        psfMatchedExposure.setFilter(templateExposure.getFilter())
        psfMatchedExposure.setCalib(scienceExposure.getCalib())
        results.warpedExposure  = templateExposure
        results.matchedExposure = psfMatchedExposure
        return results

    @pipeBase.timeMethod
    def matchMaskedImages(self, templateMaskedImage, scienceMaskedImage, candidateList,
                          templateFwhmPix=None, scienceFwhmPix=None):
        """PSF-match a MaskedImage (templateMaskedImage) to a reference MaskedImage (scienceMaskedImage)

        Do the following, in order:
        - Determine a PSF matching kernel and differential background model
            that matches templateMaskedImage to scienceMaskedImage
        - Convolve templateMaskedImage by the PSF matching kernel

        @param templateMaskedImage: masked image to PSF-match to the reference masked image;
            must be warped to match the reference masked image
        @param scienceMaskedImage: maskedImage whose PSF is to be matched to
        @param templateFwhmPix: FWHM (in pixels) of the Psf in the template image (image to convolve)
        @param scienceFwhmPix: FWHM (in pixels) of the Psf in the science image
        @param candidateList: a list of footprints/maskedImages for kernel candidates; 
                              if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF

        @return a pipeBase.Struct containing these fields:
        - psfMatchedMaskedImage: the PSF-matched masked image =
            templateMaskedImage convolved with psfMatchingKernel.
            This has the same xy0, dimensions and wcs as scienceMaskedImage.
        - psfMatchingKernel: the PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to solve for the PSF matching kernel

        @raise RuntimeError if input images have different dimensions
        """

        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayTemplate = lsstDebug.Info(__name__).displayTemplate
        displaySciIm = lsstDebug.Info(__name__).displaySciIm
        displaySpatialCells = lsstDebug.Info(__name__).displaySpatialCells
        maskTransparency = lsstDebug.Info(__name__).maskTransparency
        if not maskTransparency:
            maskTransparency = 0
        ds9.setMaskTransparency(maskTransparency)

        if not candidateList:
            raise RuntimeError("Candidate list must be populated by makeCandidateList")

        if not self._validateSize(templateMaskedImage, scienceMaskedImage):
            pexLog.Trace(self.log.getName(), 1, "ERROR: Input images different size")
            raise RuntimeError("Input images different size")

        if display and displayTemplate:
            ds9.mtv(templateMaskedImage, frame=lsstDebug.frame, title="Image to convolve")
            lsstDebug.frame += 1

        if display and  displaySciIm:
            ds9.mtv(scienceMaskedImage, frame=lsstDebug.frame, title="Image to not convolve")
            lsstDebug.frame += 1

        kernelCellSet = self._buildCellSet(templateMaskedImage,
                                           scienceMaskedImage,
                                           candidateList)

        if display and displaySpatialCells:
            diUtils.showKernelSpatialCells(scienceMaskedImage, kernelCellSet,
                                           symb="o", ctype=ds9.CYAN, ctypeUnused=ds9.YELLOW, ctypeBad=ds9.RED,
                                           size=4, frame=lsstDebug.frame, title="Image to not convolve")
            lsstDebug.frame += 1

        if templateFwhmPix and scienceFwhmPix:
            self.log.info("Matching Psf FWHM %.2f -> %.2f pix" % (templateFwhmPix, scienceFwhmPix))

        if self.kconfig.useBicForKernelBasis:
            tmpKernelCellSet = self._buildCellSet(templateMaskedImage,
                                                  scienceMaskedImage,
                                                  candidateList)
            nbe = diffimTools.NbasisEvaluator(self.kconfig, templateFwhmPix, scienceFwhmPix)
            bicDegrees = nbe(tmpKernelCellSet, self.log)
            basisList = makeKernelBasisList(self.kconfig, templateFwhmPix, scienceFwhmPix,
                                            alardDegGauss=bicDegrees[0], metadata=self.metadata)
            del tmpKernelCellSet
        else:
            basisList = makeKernelBasisList(self.kconfig, templateFwhmPix, scienceFwhmPix,
                                            metadata=self.metadata)

        spatialSolution, psfMatchingKernel, backgroundModel = self._solve(kernelCellSet, basisList)




        psfMatchedMaskedImage = afwImage.MaskedImageF(templateMaskedImage.getBBox(afwImage.PARENT))
        doNormalize = False
        afwMath.convolve(psfMatchedMaskedImage, templateMaskedImage, psfMatchingKernel, doNormalize)
        return pipeBase.Struct(
            matchedImage=psfMatchedMaskedImage,
            psfMatchingKernel=psfMatchingKernel,
            backgroundModel=backgroundModel,
            kernelCellSet=kernelCellSet,
        )

    @pipeBase.timeMethod
    def subtractExposures(self, templateExposure, scienceExposure,
                          templateFwhmPix=None, scienceFwhmPix=None,
                          candidateList=None, doWarping=True, convolveTemplate=True):
        """Subtract two Exposures

        Do the following, in order:
        - Warp templateExposure to match scienceExposure, if their WCSs do not already match
        - Determine a PSF matching kernel and differential background model
            that matches templateExposure to scienceExposure
        - PSF-match templateExposure to scienceExposure
        - Compute subtracted exposure (see return values for equation).

        @param templateExposure: exposure to PSF-match to scienceExposure
        @param scienceExposure: reference Exposure
        @param templateFwhmPix: FWHM (in pixels) of the Psf in the template image (image to convolve)
        @param scienceFwhmPix: FWHM (in pixels) of the Psf in the science image
        @param candidateList: a list of footprints/maskedImages for kernel candidates;
                              if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF
        @param doWarping: what to do if templateExposure's and scienceExposure's WCSs do not match:
            - if True then warp templateExposure to match scienceExposure
            - if False then raise an Exception
        @param convolveTemplate: convolve the template image or the science image
            - if True, templateExposure is warped if doWarping, templateExposure is convolved
            - if False, templateExposure is warped if doWarping, scienceExposure is convolved

        @return a pipeBase.Struct containing these fields:
        - subtractedExposure: subtracted Exposure = scienceExposure - (matchedImage + backgroundModel)
        - matchedImage: templateExposure after warping to match templateExposure (if doWarping true),
            and convolving with psfMatchingKernel
        - psfMatchingKernel: PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to determine PSF matching kernel
        """
        results = self.matchExposures(
            templateExposure=templateExposure,
            scienceExposure=scienceExposure,
            templateFwhmPix=templateFwhmPix,
            scienceFwhmPix=scienceFwhmPix,
            candidateList=candidateList,
            doWarping=doWarping,
            convolveTemplate=convolveTemplate
        )

        subtractedExposure = afwImage.ExposureF(scienceExposure, True)
        if convolveTemplate:
            subtractedMaskedImage  = subtractedExposure.getMaskedImage()
            subtractedMaskedImage -= results.matchedExposure.getMaskedImage()
            subtractedMaskedImage -= results.backgroundModel
        else:
            subtractedExposure.setMaskedImage(results.warpedExposure.getMaskedImage())
            subtractedMaskedImage  = subtractedExposure.getMaskedImage()
            subtractedMaskedImage -= results.matchedExposure.getMaskedImage()
            subtractedMaskedImage -= results.backgroundModel

            # Preserve polarity of differences
            subtractedMaskedImage *= -1

            # Place back on native photometric scale
            subtractedMaskedImage /= results.psfMatchingKernel.computeImage(
                afwImage.ImageD(results.psfMatchingKernel.getDimensions()), False)

        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayDiffIm = lsstDebug.Info(__name__).displayDiffIm
        maskTransparency = lsstDebug.Info(__name__).maskTransparency
        if not maskTransparency:
            maskTransparency = 0
        ds9.setMaskTransparency(maskTransparency)
        if display and displayDiffIm:
            ds9.mtv(templateExposure, frame=lsstDebug.frame, title="Template")
            lsstDebug.frame += 1
            ds9.mtv(results.matchedExposure, frame=lsstDebug.frame, title="Matched template")
            lsstDebug.frame += 1
            ds9.mtv(scienceExposure, frame=lsstDebug.frame, title="Science Image")
            lsstDebug.frame += 1
            ds9.mtv(subtractedExposure, frame=lsstDebug.frame, title="Difference Image")
            lsstDebug.frame += 1

        results.subtractedExposure = subtractedExposure
        return results

    @pipeBase.timeMethod
    def subtractMaskedImages(self, templateMaskedImage, scienceMaskedImage, candidateList,
            templateFwhmPix=None, scienceFwhmPix=None):
        """Subtract two MaskedImages

        Do the following, in order:
        - PSF-match templateMaskedImage to scienceMaskedImage
        - Determine the differential background
        - Return the difference: scienceMaskedImage -
            ((warped templateMaskedImage convolved with psfMatchingKernel) + backgroundModel)

        @param templateMaskedImage: MaskedImage to PSF-match to scienceMaskedImage
        @param scienceMaskedImage: reference MaskedImage
        @param templateFwhmPix: FWHM (in pixels) of the Psf in the template image (image to convolve)
        @param scienceFwhmPix: FWHM (in pixels) of the Psf in the science image
        @param candidateList: a list of footprints/maskedImages for kernel candidates;
                              if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF

        @return a pipeBase.Struct containing these fields:
        - subtractedMaskedImage = scienceMaskedImage - (matchedImage + backgroundModel)
        - matchedImage: templateMaskedImage convolved with psfMatchingKernel
        - psfMatchingKernel: PSF matching kernel
        - backgroundModel: differential background model
        - kernelCellSet: SpatialCellSet used to determine PSF matching kernel
        """
        if not candidateList:
            raise RuntimeError("Candidate list must be populated by makeCandidateList")

        results = self.matchMaskedImages(
            templateMaskedImage=templateMaskedImage,
            scienceMaskedImage=scienceMaskedImage,
            candidateList=candidateList,
            templateFwhmPix=templateFwhmPix,
            scienceFwhmPix=scienceFwhmPix,
            )

        subtractedMaskedImage  = afwImage.MaskedImageF(scienceMaskedImage, True)
        subtractedMaskedImage -= results.matchedImage
        subtractedMaskedImage -= results.backgroundModel
        results.subtractedMaskedImage = subtractedMaskedImage

        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayDiffIm = lsstDebug.Info(__name__).displayDiffIm
        maskTransparency = lsstDebug.Info(__name__).maskTransparency
        if not maskTransparency:
            maskTransparency = 0
        ds9.setMaskTransparency(maskTransparency)
        if display and displayDiffIm:
            ds9.mtv(subtractedMaskedImage, frame=lsstDebug.frame)
            lsstDebug.frame += 1

        return results

    def getSelectSources(self, exposure, sigma=None, doSmooth=True, idFactory=None):
        """Get sources to use for Psf-matching

        This method runs detection and measurement on an exposure.
        The returned set of sources will be used as candidates for
        Psf-matching.

        @param exposure: Exposure on which to run detection/measurement
        @param sigma: Detection threshold
        @param doSmooth: Whether or not to smooth the Exposure with Psf before detection
        @param idFactory: Factory for the generation of Source ids
        @param binSize: Binsize for background subtraction of Exposure

        @return source catalog containing candidates for the Psf-matching
        """

        if idFactory:
            table = afwTable.SourceTable.make(self.selectSchema, idFactory)
        else:
            table = afwTable.SourceTable.make(self.selectSchema)
        mi = exposure.getMaskedImage()
        # If binsize is not set, fall back to simple median estimation
        imArr = mi.getImage().getArray()
        maskArr = mi.getMask().getArray()
        miArr = np.ma.masked_array(imArr, mask=maskArr)
        if not binSize:
            bkgd = np.ma.extras.median(miArr)
        else:
            try:
                bkgd = getBackground(mi, self.afwBackgroundConfig).getImageF()
            except:
                self.log.warn("Failed to get background model.  Falling back to median background estimation")
                bkgd = np.ma.extras.median(miArr)


        #Take off background for detection
        mi -= bkgd
        table.setMetadata(self.selectAlgMetadata) 
        detRet = self.selectDetection.makeSourceCatalog(
            table=table,
            exposure=exposure,
            sigma=sigma,
            doSmooth=doSmooth
        )
        selectSources = detRet.sources
        self.selectMeasurement.measure(exposure, selectSources)
        #Put back on the background in case it is needed down stream
        mi += bkgd
        del bkgd
        return selectSources

    def makeCandidateList(self, templateExposure, scienceExposure, kernelSize, candidateList=None):
        """Accept or generate a list of candidate sources for
        Psf-matching, and examine the Mask planes in both of the
        images for indications of bad pixels

        @param templateExposure: Exposure that will be convolved
        @param scienceExposure: Exposure that will be matched-to
        @param kernelSize: Dimensions of the Psf-matching Kernel, used to grow detection footprints
        @param candidateList: List of Sources to examine

        @return a list of dicts having a "source" and "footprint"
        field for the Sources deemed to be appropriate for Psf
        matching
        """
        if candidateList is None:
            candidateList = self.getSelectSources(scienceExposure)

        listTypes = set(type(x) for x in candidateList)
        if (not len(listTypes) == 1) or (listTypes.pop() == afwTable.SourceRecord):
            raise RuntimeError("Can only make candidate list from set of SourceRecords.  Got %s instead." \
                                   % (type(candidateList[0])))
        candidateList = diffimTools.sourceToFootprintList(candidateList,
                                                          templateExposure, scienceExposure,
                                                          kernelSize,
                                                          self.kconfig.detectionConfig,
                                                          self.log)
        if len(candidateList) == 0:
            raise RuntimeError("Cannot find any objects suitable for KernelCandidacy")

        return candidateList

    def _adaptCellSize(self, candidateList):
        """ NOT IMPLEMENTED YET"""
        nCand = len(candidateList)
        return self.kconfig.sizeCellX, self.kconfig.sizeCellY

    def _buildCellSet(self, templateMaskedImage, scienceMaskedImage, candidateList):
        """Build a SpatialCellSet for use with the solve method

        @param templateMaskedImage: MaskedImage to PSF-matched to scienceMaskedImage
        @param scienceMaskedImage: reference MaskedImage
        @param candidateList: a list of footprints/maskedImages for kernel candidates;
                              if None then source detection is run.
            - Currently supported: list of Footprints or measAlg.PsfCandidateF

        @return kernelCellSet: a SpatialCellSet for use with self._solve
        """
        if not candidateList:
            raise RuntimeError("Candidate list must be populated by makeCandidateList")

        sizeCellX, sizeCellY = self._adaptCellSize(candidateList)

        # Object to store the KernelCandidates for spatial modeling
        kernelCellSet = afwMath.SpatialCellSet(templateMaskedImage.getBBox(afwImage.PARENT),
                                               sizeCellX, sizeCellY)

        policy = pexConfig.makePolicy(self.kconfig)
        # Place candidates within the spatial grid
        for cand in candidateList:
            bbox = cand['footprint'].getBBox()

            tmi  = afwImage.MaskedImageF(templateMaskedImage, bbox, afwImage.PARENT)
            smi  = afwImage.MaskedImageF(scienceMaskedImage, bbox, afwImage.PARENT)
            cand = diffimLib.makeKernelCandidate(cand['source'], tmi, smi, policy)

            self.log.logdebug("Candidate %d at %f, %f" % (cand.getId(), cand.getXCenter(), cand.getYCenter()))
            kernelCellSet.insertCandidate(cand)

        return kernelCellSet

    def _validateSize(self, templateMaskedImage, scienceMaskedImage):
        """Return True if two image-like objects are the same size
        """
        return templateMaskedImage.getDimensions() == scienceMaskedImage.getDimensions()

    def _validateWcs(self, templateExposure, scienceExposure):
        """Return True if the WCS of the two Exposures have the same origin and extent
        """
        templateWcs    = templateExposure.getWcs() 
        scienceWcs     = scienceExposure.getWcs()
        templateBBox   = templateExposure.getBBox(afwImage.PARENT)
        scienceBBox    = scienceExposure.getBBox(afwImage.PARENT)

        # LLC
        templateOrigin = templateWcs.pixelToSky(afwGeom.Point2D(templateBBox.getBegin()))
        scienceOrigin  = scienceWcs.pixelToSky(afwGeom.Point2D(scienceBBox.getBegin()))

        # URC
        templateLimit  = templateWcs.pixelToSky(afwGeom.Point2D(templateBBox.getEnd()))
        scienceLimit   = scienceWcs.pixelToSky(afwGeom.Point2D(scienceBBox.getEnd()))

        self.log.info("Template Wcs : %f,%f -> %f,%f" %
                      (templateOrigin[0], templateOrigin[1],
                       templateLimit[0], templateLimit[1]))
        self.log.info("Science Wcs : %f,%f -> %f,%f" %
                      (scienceOrigin[0], scienceOrigin[1],
                       scienceLimit[0], scienceLimit[1]))

        templateBBox = afwGeom.Box2D(templateOrigin.getPosition(), templateLimit.getPosition())
        scienceBBox  = afwGeom.Box2D(scienceOrigin.getPosition(), scienceLimit.getPosition())
        if not (templateBBox.overlaps(scienceBBox)):
            raise RuntimeError("Input images do not overlap at all")

        if ( (templateOrigin.getPosition() != scienceOrigin.getPosition()) or
             (templateLimit.getPosition()  != scienceLimit.getPosition())  or
             (templateExposure.getDimensions() != scienceExposure.getDimensions())):
            return False
        return True
