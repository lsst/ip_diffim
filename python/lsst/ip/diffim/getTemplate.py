#
# LSST Data Management System
# Copyright 2016 LSST Corporation.
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

import numpy as np

import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.sphgeom as sphgeom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.afw.math as afwMath
from lsst.ip.diffim.dcrModel import DcrModel
from lsst.meas.algorithms import CoaddPsf, WarpedPsf, CoaddPsfConfig

__all__ = ["GetCoaddAsTemplateTask", "GetCoaddAsTemplateConfig",
           "GetCalexpAsTemplateTask", "GetCalexpAsTemplateConfig"]


class GetCoaddAsTemplateConfig(pexConfig.Config):
    templateBorderSize = pexConfig.Field(
        dtype=int,
        default=10,
        doc="Number of pixels to grow the requested template image to account for warping"
    )
    coaddName = pexConfig.Field(
        doc="coadd name: typically one of 'deep', 'goodSeeing', or 'dcr'",
        dtype=str,
        default="deep",
    )
    numSubfilters = pexConfig.Field(
        doc="Number of subfilters in the DcrCoadd. Used only if ``coaddName``='dcr'",
        dtype=int,
        default=3,
    )
    effectiveWavelength = pexConfig.Field(
        doc="Effective wavelength of the filter. Used only if ``coaddName``='dcr'",
        optional=True,
        dtype=float,
    )
    bandwidth = pexConfig.Field(
        doc="Bandwidth of the physical filter. Used only if ``coaddName``='dcr'",
        optional=True,
        dtype=float,
    )
    warpType = pexConfig.Field(
        doc="Warp type of the coadd template: one of 'direct' or 'psfMatched'",
        dtype=str,
        default="direct",
    )
    coaddPsf = pexConfig.ConfigField(
        doc="Configuration for CoaddPsf",
        dtype=CoaddPsfConfig,
    )

    def validate(self):
        if self.coaddName == 'dcr':
            if self.effectiveWavelength is None or self.bandwidth is None:
                raise ValueError("The effective wavelength and bandwidth of the physical filter "
                                 "must be set in the getTemplate config for DCR coadds. "
                                 "Required until transmission curves are used in DM-13668.")


class GetCoaddAsTemplateTask(pipeBase.Task):
    """Subtask to retrieve coadd for use as an image difference template.

    This is the default getTemplate Task to be run as a subtask by
    ``pipe.tasks.ImageDifferenceTask``. The main methods are ``run()`` and
    ``runGen3()``.

    Notes
    -----
    From the given skymap, the closest tract is selected;  multiple tracts  are
    not supported. The assembled template inherits the WCS of the selected
    skymap tract and the resolution of the template exposures. Overlapping box
    regions of the input template patches are pixel by pixel copied into the
    assembled template image. There is no warping or pixel resampling.

    Pixels with no overlap of any available input patches are set to ``nan`` value
    and ``NO_DATA`` flagged.
    """

    ConfigClass = GetCoaddAsTemplateConfig
    _DefaultName = "GetCoaddAsTemplateTask"

    def runDataRef(self, exposure, sensorRef, templateIdList=None):
        """Gen2 task entry point. Retrieve and mosaic a template coadd exposure
        that overlaps the science exposure.

        Parameters
        ----------
        exposure: `lsst.afw.image.Exposure`
            an exposure for which to generate an overlapping template
        sensorRef : TYPE
            a Butler data reference that can be used to obtain coadd data
        templateIdList : TYPE, optional
            list of data ids, unused here, in the case of coadd template

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            - ``exposure`` : `lsst.afw.image.ExposureF`
                a template coadd exposure assembled out of patches
            - ``sources`` :  None for this subtask
        """
        skyMap = sensorRef.get(datasetType=self.config.coaddName + "Coadd_skyMap")
        tractInfo, patchList, skyCorners = self.getOverlapPatchList(exposure, skyMap)

        availableCoaddRefs = dict()
        for patchInfo in patchList:
            patchNumber = tractInfo.getSequentialPatchIndex(patchInfo)
            patchArgDict = dict(
                datasetType=self.getCoaddDatasetName() + "_sub",
                bbox=patchInfo.getOuterBBox(),
                tract=tractInfo.getId(),
                patch="%s,%s" % (patchInfo.getIndex()[0], patchInfo.getIndex()[1]),
                subfilter=0,
                numSubfilters=self.config.numSubfilters,
            )

            if sensorRef.datasetExists(**patchArgDict):
                self.log.info("Reading patch %s" % patchArgDict)
                availableCoaddRefs[patchNumber] = patchArgDict

        templateExposure = self.run(
            tractInfo, patchList, skyCorners, availableCoaddRefs,
            sensorRef=sensorRef, visitInfo=exposure.getInfo().getVisitInfo()
        )
        return pipeBase.Struct(exposure=templateExposure, sources=None)

    def runQuantum(self, exposure, butlerQC, skyMapRef, coaddExposureRefs):
        """Gen3 task entry point. Retrieve and mosaic a template coadd exposure
        that overlaps the science exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The science exposure to define the sky region of the template coadd.
        butlerQC : `lsst.pipe.base.ButlerQuantumContext`
            Butler like object that supports getting data by DatasetRef.
        skyMapRef : `lsst.daf.butler.DatasetRef`
            Reference to SkyMap object that corresponds to the template coadd.
        coaddExposureRefs : iterable of `lsst.daf.butler.DeferredDatasetRef`
            Iterable of references to the available template coadd patches.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            - ``exposure`` : `lsst.afw.image.ExposureF`
                a template coadd exposure assembled out of patches
            - ``sources`` :  `None` for this subtask
        """
        skyMap = butlerQC.get(skyMapRef)
        coaddExposureRefs = butlerQC.get(coaddExposureRefs)
        tracts = [ref.dataId['tract'] for ref in coaddExposureRefs]
        if tracts.count(tracts[0]) == len(tracts):
            tractInfo = skyMap[tracts[0]]
        else:
            raise RuntimeError("Templates constructed from multiple Tracts not yet supported")

        detectorBBox = exposure.getBBox()
        detectorWcs = exposure.getWcs()
        detectorCorners = detectorWcs.pixelToSky(geom.Box2D(detectorBBox).getCorners())

        availableCoaddRefs = dict()
        for coaddRef in coaddExposureRefs:
            dataId = coaddRef.dataId
            patchWcs = skyMap[dataId['tract']].getWcs()
            patchBBox = skyMap[dataId['tract']][dataId['patch']].getOuterBBox()
            if self.bboxIntersectsCorners(patchBBox, patchWcs, detectorCorners):
                if self.config.coaddName == 'dcr':
                    self.log.info("Using template input tract=%s, patch=%s, subfilter=%s" %
                                  (dataId['tract'], dataId['patch'], dataId['subfilter']))
                    if dataId['patch'] in availableCoaddRefs:
                        availableCoaddRefs[dataId['patch']].append(coaddRef)
                    else:
                        availableCoaddRefs[dataId['patch']] = [coaddRef, ]
                else:
                    self.log.info("Using template input tract=%s, patch=%s" %
                                  (dataId['tract'], dataId['patch']))
                    availableCoaddRefs[dataId['patch']] = coaddRef

        patchList = [tractInfo[patch] for patch in availableCoaddRefs.keys()]
        #templateExposure = self.run(tractInfo, patchList, detectorCorners, availableCoaddRefs,
        #                            visitInfo=exposure.getInfo().getVisitInfo(), wcs=detectorWcs)
        templateExposure = self.runBob(exposure, availableCoaddRefs, skyMap, templateIdList=None).exposure
        return pipeBase.Struct(exposure=templateExposure, sources=None)

    def getOverlapPatchList(self, exposure, skyMap):
        """Select the relevant tract and its patches that overlap with the science exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The science exposure to define the sky region of the template coadd.

        skyMap : `lsst.skymap.BaseSkyMap`
            SkyMap object that corresponds to the template coadd.

        Returns
        -------
        result : `tuple` of
         - ``tractInfo`` : `lsst.skymap.TractInfo`
             The selected tract.
         - ``patchList`` : `list` of `lsst.skymap.PatchInfo`
             List of all overlap patches of the selected tract.
         - ``skyCorners`` : `list` of `lsst.geom.SpherePoint`
             Corners of the exposure in the sky in the order given by `lsst.geom.Box2D.getCorners`.
        """
        expWcs = exposure.getWcs()
        expBoxD = geom.Box2D(exposure.getBBox())
        expBoxD.grow(self.config.templateBorderSize)
        ctrSkyPos = expWcs.pixelToSky(expBoxD.getCenter())
        tractInfo = skyMap.findTract(ctrSkyPos)
        self.log.info("Using skyMap tract %s" % (tractInfo.getId(),))
        skyCorners = [expWcs.pixelToSky(pixPos) for pixPos in expBoxD.getCorners()]
        patchList = tractInfo.findPatchList(skyCorners)

        if not patchList:
            raise RuntimeError("No suitable tract found")

        self.log.info("Assembling %s coadd patches" % (len(patchList),))
        self.log.info("exposure dimensions=%s" % exposure.getDimensions())

        return (tractInfo, patchList, skyCorners)

    def run(self, tractInfo, patchList, skyCorners, availableCoaddRefs,
            sensorRef=None, visitInfo=None, wcs=None):
        """Gen2 and gen3 shared code: determination of exposure dimensions and
        copying of pixels from overlapping patch regions.

        Parameters
        ----------
        skyMap : `lsst.skymap.BaseSkyMap`
            SkyMap object that corresponds to the template coadd.
        tractInfo : `lsst.skymap.TractInfo`
            The selected tract.
        patchList : iterable of `lsst.skymap.patchInfo.PatchInfo`
            Patches to consider for making the template exposure.
        skyCorners : list of `lsst.geom.SpherePoint`
            Sky corner coordinates to be covered by the template exposure.
        availableCoaddRefs : `dict` [`int`]
            Dictionary of spatially relevant retrieved coadd patches,
            indexed by their sequential patch number. In Gen3 mode, values are
            `lsst.daf.butler.DeferredDatasetHandle` and ``.get()`` is called,
            in Gen2 mode, ``sensorRef.get(**coaddef)`` is called to retrieve the coadd.
        sensorRef : `lsst.daf.persistence.ButlerDataRef`, Gen2 only
            Butler data reference to get coadd data.
            Must be `None` for Gen3.
        visitInfo : `lsst.afw.image.VisitInfo`, Gen2 only
            VisitInfo to make dcr model.

        Returns
        -------
        templateExposure: `lsst.afw.image.ExposureF`
            The created template exposure.
        """
        coaddWcs = tractInfo.getWcs()

        # compute coadd bbox
        coaddBBox = geom.Box2D()
        for skyPos in skyCorners:
            coaddBBox.include(coaddWcs.skyToPixel(skyPos))
        coaddBBox = geom.Box2I(coaddBBox)
        self.log.info("coadd dimensions=%s" % coaddBBox.getDimensions())

        # Table for CoaddPSF
        tractsSchema = afwTable.ExposureTable.makeMinimalSchema()
        tractKey = tractsSchema.addField('tract', type=np.int32, doc='Which tract')
        patchKey = tractsSchema.addField('patch', type=np.int32, doc='Which patch')
        weightKey = tractsSchema.addField('weight', type=float, doc='Weight for each tract, should be 1')
        tractsCatalog = afwTable.ExposureCatalog(tractsSchema)

        coaddExposure = afwImage.ExposureF(coaddBBox, coaddWcs)
        coaddExposure.maskedImage.set(np.nan, afwImage.Mask.getPlaneBitMask("NO_DATA"), np.nan)
        nPatchesFound = 0
        coaddFilterLabel = None
        coaddPsf = None
        coaddPhotoCalib = None
        for patchInfo in patchList:
            patchNumber = tractInfo.getSequentialPatchIndex(patchInfo)
            patchSubBBox = patchInfo.getInnerBBox()
            patchSubBBox.clip(coaddBBox)
            if patchNumber not in availableCoaddRefs:
                self.log.warn(f"skip patch={patchNumber}; patch does not exist for this coadd")
                continue
            if patchSubBBox.isEmpty():
                self.log.info(f"skip tract={availableCoaddRefs[patchNumber]['tract']}, "
                              f"patch={patchNumber}; no overlapping pixels")
                continue

            if self.config.coaddName == 'dcr':
                patchInnerBBox = patchInfo.getInnerBBox()
                patchInnerBBox.clip(coaddBBox)
                if np.min(patchInnerBBox.getDimensions()) <= 2*self.config.templateBorderSize:
                    self.log.info("skip tract=%(tract)s, patch=%(patch)s; too few pixels."
                                  % availableCoaddRefs[patchNumber])
                    continue
                self.log.info("Constructing DCR-matched template for patch %s"
                              % availableCoaddRefs[patchNumber])

                if sensorRef:
                    dcrModel = DcrModel.fromDataRef(sensorRef,
                                                    self.config.effectiveWavelength,
                                                    self.config.bandwidth,
                                                    **availableCoaddRefs[patchNumber])
                else:
                    dcrModel = DcrModel.fromQuantum(availableCoaddRefs[patchNumber],
                                                    self.config.effectiveWavelength,
                                                    self.config.bandwidth)
                # The edge pixels of the DcrCoadd may contain artifacts due to missing data.
                # Each patch has significant overlap, and the contaminated edge pixels in
                # a new patch will overwrite good pixels in the overlap region from
                # previous patches.
                # Shrink the BBox to remove the contaminated pixels,
                # but make sure it is only the overlap region that is reduced.
                dcrBBox = geom.Box2I(patchSubBBox)
                dcrBBox.grow(-self.config.templateBorderSize)
                dcrBBox.include(patchInnerBBox)
                coaddPatch = dcrModel.buildMatchedExposure(bbox=dcrBBox,
                                                           wcs=coaddWcs,
                                                           visitInfo=visitInfo)
            else:
                if sensorRef is None:
                    # Gen3
                    coaddPatch = availableCoaddRefs[patchNumber].get()
                else:
                    # Gen2
                    coaddPatch = sensorRef.get(**availableCoaddRefs[patchNumber])
            nPatchesFound += 1

            # Gen2 get() seems to clip based on bbox kwarg but we removed bbox
            # calculation from caller code. Gen3 also does not do this.
            overlapBox = coaddPatch.getBBox()
            overlapBox.clip(coaddBBox)
            coaddExposure.maskedImage.assign(coaddPatch.maskedImage[overlapBox], overlapBox)

            if coaddFilterLabel is None:
                coaddFilterLabel = coaddPatch.getFilterLabel()

            # Retrieve the PSF for this coadd tract, if not already retrieved
            # if coaddPsf is None and coaddPatch.hasPsf():
            record = tractsCatalog.addNew()
            record.setPsf(coaddPatch.getPsf())
            record.setWcs(coaddPatch.getWcs())
            record.setPhotoCalib(coaddPatch.getPhotoCalib())
            record.setBBox(patchSubBBox)
            record.set(tractKey, tractInfo.getId())
            record.set(patchKey, int(patchNumber))
            record.set(weightKey, 1.)

            # Retrieve the calibration for this coadd tract, if not already retrieved
            if coaddPhotoCalib is None:
                coaddPhotoCalib = coaddPatch.getPhotoCalib()

        import pdb; pdb.set_trace()
        if coaddPhotoCalib is None:
            raise RuntimeError("No coadd PhotoCalib found!")
        if nPatchesFound == 0:
            raise RuntimeError("No patches found!")

        coaddExposure.setPhotoCalib(coaddPhotoCalib)
        coaddPsf = CoaddPsf(tractsCatalog, coaddWcs, self.config.coaddPsf.makeControl())
        coaddExposure.setPsf(coaddPsf)
        coaddExposure.setFilterLabel(coaddFilterLabel)

        return coaddExposure

    def runBob(self, exposure, availableCoaddRefs, skyMap, templateIdList=None):
        """Retrieve and mosaic a template coadd that overlaps the exposure where
        the template spans multiple tracts.
        The resulting template image will be an average of all the input templates from
        the separate tracts.
        The PSF on the template is created by combining the CoaddPsf on each template image
        into a meta-CoaddPsf.
        Parameters
        ----------
        exposure: `lsst.afw.image.Exposure`
            an exposure for which to generate an overlapping template
        sensorRef : TYPE
            a Butler data reference that can be used to obtain coadd data
        templateIdList : TYPE, optional
            list of data ids (unused)
        Returns
        -------
        result : `struct`
            return a pipeBase.Struct:
            - ``exposure`` : a template coadd exposure assembled out of patches
            - ``sources`` :  None for this subtask
        """

        # Table for CoaddPSF
        tractsSchema = afwTable.ExposureTable.makeMinimalSchema()
        tractKey = tractsSchema.addField('tract', type=np.int32, doc='Which tract')
        patchKey = tractsSchema.addField('patch', type=np.int32, doc='Which patch')
        weightKey = tractsSchema.addField('weight', type=float, doc='Weight for each tract, should be 1')
        tractsCatalog = afwTable.ExposureCatalog(tractsSchema)

        expWcs = exposure.getWcs()
        expBoxD = geom.Box2D(exposure.getBBox())
        expBoxD.grow(self.config.templateBorderSize)
        ctrSkyPos = expWcs.pixelToSky(expBoxD.getCenter())

        centralTractInfo = skyMap.findTract(ctrSkyPos)
        if not centralTractInfo:
            raise RuntimeError("No suitable tract found for central point")

        self.log.info("Central skyMap tract %s" % (centralTractInfo.getId(),))

        skyCorners = [expWcs.pixelToSky(pixPos) for pixPos in expBoxD.getCorners()]
        tractPatchList = skyMap.findTractPatchList(skyCorners)
        if not tractPatchList:
            raise RuntimeError("No suitable tract found")

        self.log.info("All overlapping skyMap tracts %s" % ([a[0].getId() for a in tractPatchList]))

        # Move central tract to front of the list and use as the reference
        tracts = [tract[0].getId() for tract in tractPatchList]
        centralIndex = tracts.index(centralTractInfo.getId())
        tracts.insert(0, tracts.pop(centralIndex))
        tractPatchList.insert(0, tractPatchList.pop(centralIndex))

        coaddPsf = None
        coaddFilter = None
        nPatchesFound = 0

        maskedImageList = []
        weightList = []

        for itract, tract in enumerate([tracts[0],]):
            tractInfo = tractPatchList[itract][0]

            coaddWcs = tractInfo.getWcs()
            coaddBBox = geom.Box2D()
            for skyPos in skyCorners:
                coaddBBox.include(coaddWcs.skyToPixel(skyPos))
            coaddBBox = geom.Box2I(coaddBBox)

            if itract == 0:
                # Define final wcs and bounding box from the reference tract
                finalWcs = coaddWcs
                finalBBox = coaddBBox

            patchList = tractPatchList[itract][1]
            #import pdb; pdb.set_trace()
            for patchInfo in patchList:
                self.log.info('Adding patch %s from tract %s' % (patchInfo.getIndex(),tract))

                # Local patch information
                patchSubBBox = geom.Box2I(patchInfo.getInnerBBox())
                patchSubBBox.clip(coaddBBox)
                patchNumber = tractInfo.getSequentialPatchIndex(patchInfo)
                patchInt = int(f"{patchInfo.getIndex()[0]}{patchInfo.getIndex()[1]}")
                innerBBox = geom.Box2I(tractInfo._minimumBoundingBox(finalWcs))

                if itract == 0:
                    # clip to image and tract boundaries
                    patchSubBBox.clip(finalBBox)
                    patchSubBBox.clip(innerBBox)
                    if patchSubBBox.getArea() == 0:
                        self.log.debug("No ovlerap for patch %s" % patchInfo)
                        continue

                    coaddPatch = availableCoaddRefs[patchNumber].get()
                    if coaddFilter is None:
                        coaddFilter = coaddPatch.getFilter()

                    # create full image from final bounding box
                    exp = afwImage.ExposureF(finalBBox, finalWcs)
                    exp.maskedImage.set(np.nan, afwImage.Mask.getPlaneBitMask("NO_DATA"), np.nan)
                    overlapBox = coaddPatch.getBBox()
                    overlapBox.clip(coaddBBox)
                    exp.maskedImage.assign(coaddPatch.maskedImage[overlapBox], overlapBox)
                    maskedImageList.append(exp.maskedImage)
                    weightList.append(1)

                    record = tractsCatalog.addNew()
                    record.setPsf(coaddPatch.getPsf())
                    record.setWcs(coaddPatch.getWcs())
                    record.setPhotoCalib(coaddPatch.getPhotoCalib())
                    record.setBBox(patchSubBBox)
                    record.set(tractKey, tract)
                    record.set(patchKey, patchInt)
                    record.set(weightKey, 1.)
                    nPatchesFound += 1
                else:
                    # compute the exposure bounding box in a tract that is not the reference tract
                    localBox = geom.Box2I()
                    for skyPos in skyCorners:
                        localBox.include(geom.Point2I(tractInfo.getWcs().skyToPixel(skyPos)))

                    # clip to patch bounding box
                    localBox.clip(patchSubBBox)

                    # grow border to deal with warping at edges
                    localBox.grow(self.config.templateBorderSize)

                    # clip to tract inner bounding box
                    localInnerBBox = geom.Box2I(tractInfo._minimumBoundingBox(tractInfo.getWcs()))
                    localBox.clip(localInnerBBox)

                    coaddPatch = availableCoaddRefs[patchNumber].get()

                    # warp to reference tract wcs
                    xyTransform = geom.makeWcsPairTransform(coaddPatch.getWcs(), finalWcs)
                    psfWarped = WarpedPsf(coaddPatch.getPsf(), xyTransform)
                    warped = self.warper.warpExposure(finalWcs, coaddPatch, maxBBox=finalBBox)

                    # check if warpped image is viable
                    if warped.getBBox().getArea() == 0:
                        self.log.info("No ovlerap for warped patch %s. Skipping" % patchInfo)
                        continue

                    warped.setPsf(psfWarped)

                    exp = afwImage.ExposureF(finalBBox, finalWcs)
                    exp.maskedImage.set(np.nan, afwImage.Mask.getPlaneBitMask("NO_DATA"), np.nan)
                    exp.maskedImage.assign(warped.maskedImage, warped.getBBox())

                    maskedImageList.append(exp.maskedImage)
                    weightList.append(1)
                    record = tractsCatalog.addNew()
                    record.setPsf(psfWarped)
                    record.setWcs(finalWcs)
                    record.setPhotoCalib(coaddPatch.getPhotoCalib())
                    record.setBBox(warped.getBBox())
                    record.set(tractKey, tract)
                    record.set(patchKey, patchInt)
                    record.set(weightKey, 1.)
                    nPatchesFound += 1

        if nPatchesFound == 0:
            raise RuntimeError("No patches found!")

        # Combine images from individual patches together
        #  import pdb; pdb.set_trace()
        # Do not mask any values
        statsFlags = afwMath.stringToStatisticsProperty('MEAN')
        maskMap = []
        statsCtrl = afwMath.StatisticsControl()
        statsCtrl.setNanSafe(True)
        statsCtrl.setWeighted(True)
        statsCtrl.setCalcErrorFromInputVariance(True)

        coaddExposure = afwImage.ExposureF(finalBBox, finalWcs)
        coaddExposure.maskedImage.set(np.nan, afwImage.Mask.getPlaneBitMask("NO_DATA"), np.nan)
        xy0 = coaddExposure.getXY0()
        coaddExposure.maskedImage = afwMath.statisticsStack(maskedImageList,
                                                            statsFlags, statsCtrl, weightList, 0, maskMap)
        coaddExposure.maskedImage.setXY0(xy0)
        center = exposure.getBBox().getCenter()
        centerCoord = finalWcs.skyToPixel(expWcs.pixelToSky(center))
        ctrl = self.config.coaddPsf.makeControl()

        coaddPsf = CoaddPsf(tractsCatalog, finalWcs, centerCoord, ctrl.warpingKernelName, ctrl.cacheSize)
        if coaddPsf is None:
            raise RuntimeError("No coadd Psf found!")

        coaddExposure.setPsf(coaddPsf)
        coaddExposure.setFilter(coaddFilter)
        return pipeBase.Struct(exposure=coaddExposure, sources=None)


    def getCoaddDatasetName(self):
        """Return coadd name for given task config

        Returns
        -------
        CoaddDatasetName : `string`

        TODO: This nearly duplicates a method in CoaddBaseTask (DM-11985)
        """
        warpType = self.config.warpType
        suffix = "" if warpType == "direct" else warpType[0].upper() + warpType[1:]
        return self.config.coaddName + "Coadd" + suffix

    def bboxIntersectsCorners(self, bbox, wcs, otherCorners):
        """Returns true if the bbox with wcs intersects otherCorners

        Parameters
        ----------
        bbox : `lsst.geom.Box2I`
            specifying the bounding box of test region
        wcs : `lsst.afw.geom.SkyWcs`
            specifying the WCS of test region
        otherCorners : `list` of `lsst.geom.SpherePoint`
            ICRS coordinates specifying boundary of the other sky region

        Returns
        -------
        result: `bool`
           Does bbox/wcs intersect other corners?
        """
        bboxCorners = wcs.pixelToSky(geom.Box2D(bbox).getCorners())
        bboxPolygon = sphgeom.ConvexPolygon.convexHull([coord.getVector() for coord in bboxCorners])
        otherPolygon = sphgeom.ConvexPolygon.convexHull([coord.getVector() for coord in otherCorners])
        return otherPolygon.intersects(bboxPolygon)


class GetCalexpAsTemplateConfig(pexConfig.Config):
    doAddCalexpBackground = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Add background to calexp before processing it."
    )


class GetCalexpAsTemplateTask(pipeBase.Task):
    """Subtask to retrieve calexp of the same ccd number as the science image SensorRef
    for use as an image difference template. Only gen2 supported.

    To be run as a subtask by pipe.tasks.ImageDifferenceTask.
    Intended for use with simulations and surveys that repeatedly visit the same pointing.
    This code was originally part of Winter2013ImageDifferenceTask.
    """

    ConfigClass = GetCalexpAsTemplateConfig
    _DefaultName = "GetCalexpAsTemplateTask"

    def run(self, exposure, sensorRef, templateIdList):
        """Return a calexp exposure with based on input sensorRef.

        Construct a dataId based on the sensorRef.dataId combined
        with the specifications from the first dataId in templateIdList

        Parameters
        ----------
        exposure :  `lsst.afw.image.Exposure`
            exposure (unused)
        sensorRef : `list` of `lsst.daf.persistence.ButlerDataRef`
            Data reference of the calexp(s) to subtract from.
        templateIdList : `list` of `lsst.daf.persistence.ButlerDataRef`
            Data reference of the template calexp to be subtraced.
            Can be incomplete, fields are initialized from `sensorRef`.
            If there are multiple items, only the first one is used.

        Returns
        -------
        result : `struct`

            return a pipeBase.Struct:

                - ``exposure`` : a template calexp
                - ``sources`` : source catalog measured on the template
        """

        if len(templateIdList) == 0:
            raise RuntimeError("No template data reference supplied.")
        if len(templateIdList) > 1:
            self.log.warn("Multiple template data references supplied. Using the first one only.")

        templateId = sensorRef.dataId.copy()
        templateId.update(templateIdList[0])

        self.log.info("Fetching calexp (%s) as template." % (templateId))

        butler = sensorRef.getButler()
        template = butler.get(datasetType="calexp", dataId=templateId)
        if self.config.doAddCalexpBackground:
            templateBg = butler.get(datasetType="calexpBackground", dataId=templateId)
            mi = template.getMaskedImage()
            mi += templateBg.getImage()

        if not template.hasPsf():
            raise pipeBase.TaskError("Template has no psf")

        templateSources = butler.get(datasetType="src", dataId=templateId)
        return pipeBase.Struct(exposure=template,
                               sources=templateSources)

    def runDataRef(self, *args, **kwargs):
        return self.run(*args, **kwargs)

    def runQuantum(self, **kwargs):
        raise NotImplementedError("Calexp template is not supported with gen3 middleware")
