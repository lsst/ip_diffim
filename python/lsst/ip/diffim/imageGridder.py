from __future__ import absolute_import, division, print_function
from future import standard_library
standard_library.install_aliases()
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

import numpy as np

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

__all__ = ("ImageGridderTask", "ImageGridderConfig",
           "ImageGridSubtask", "ImageGridSubtaskConfig")


class ImageGridSubtaskConfig(pexConfig.Config):
    """!
    \anchor ImageGridSubtaskConfig_

    \brief Configuration parameters for the ImageGridSubtask
    """
    pass


class ImageGridSubtask(pipeBase.Task):
    """!
    \anchor ImageGridSubtask_

    \brief Simple example and base class for any task that is to be
    used as `ImageGridderConfig.gridSubtask`

    """
    ConfigClass = ImageGridSubtaskConfig
    _DefaultName = "ip_diffim_ImageGridSubtask"

    def __init__(self, *args, **kwargs):
        """! Create the image gridding subTask
        @param *args arguments to be passed to
        lsst.pipe.base.task.Task.__init__
        @param **kwargs keyword arguments to be passed to
        lsst.pipe.base.task.Task.__init__
        """
        pipeBase.Task.__init__(self, *args, **kwargs)

    def run(self, subExp, expandedSubExp, fullBBox, **kwargs):
        """! Perform operation on given sub-xposure

        @param[in] subExp the sub-exposure upon which to
        operate
        @param[in] expandedSubExp the expanded sub-exposure
        upon which to operate
        @param[in] fullBBox the bounding box of the original exposure
        @return anything, including an a `afw.Exposure`
        """
        pass


class ImageGridderConfig(pexConfig.Config):
    """!
    \anchor ImageGridderConfig_

    \brief Configuration parameters for the ImageGridderTask
    """

    gridSubtask = pexConfig.ConfigurableField(
        doc="Subtask to run on each subimage",
        target=ImageGridSubtask,
    )

    # Separate gridCentroidsX and gridCentroidsY since pexConfig.ListField accepts limited dtypes
    #  (i.e., no Point2D)
    gridCentroidsX = pexConfig.ListField(
        dtype=float,
        doc="Input X centroids around which to place subimages. If None, use grid config options below.",
        optional=True,
        default=None
    )

    gridCentroidsY = pexConfig.ListField(
        dtype=float,
        doc="Input Y centroids around which to place subimages. If None, use grid config options below.",
        optional=True,
        default=None
    )

    gridSizeX = pexConfig.Field(
        dtype=float,
        doc="""Dimensions of each grid cell in x direction""",
        default=10.,
        check=lambda x: x > 0.
    )

    gridSizeY = pexConfig.Field(
        dtype=float,
        doc="""Dimensions of each grid cell in y direction""",
        default=10.,
        check=lambda x: x > 0.
    )

    gridStepX = pexConfig.Field(
        dtype=float,
        doc="""Spacing between subsequent grid cells in x direction. If equal to
               gridSizeX, then there is no overlap in the x direction.""",
        default=10.,
        check=lambda x: x > 0.
    )

    gridStepY = pexConfig.Field(
        dtype=float,
        doc="""Spacing between subsequent grid cells in y direction. If equal to
               gridSizeY, then there is no overlap in the y direction.""",
        default=10.,
        check=lambda x: x > 0.
    )

    borderSizeX = pexConfig.Field(
        dtype=float,
        doc="""Dimensions of cell border in +/- x direction""",
        default=5.,
        check=lambda x: x > 0.
    )

    borderSizeY = pexConfig.Field(
        dtype=float,
        doc="""Dimensions of cell border in +/- y direction""",
        default=5.,
        check=lambda x: x > 0.
    )

    adjustGridOption = pexConfig.ChoiceField(
        dtype=str,
        doc="""Whether and how to adjust grid to fit evenly within image""",
        default="spacing",
        allowed={
            "spacing": "adjust spacing between centers of grid cells (allowing overlaps)",
            "size": "adjust the sizes of the grid cells (disallowing overlaps)",
            "none": "do not adjust the grid sizes or spacing"
        }
    )

    scaleByFwhm = pexConfig.Field(
        dtype=bool,
        doc="""Scale gridSize/gridStep/borderSize/overlapSize by PSF FWHM rather
        than pixels?""",
        default=True
    )

    reduceOperation = pexConfig.ChoiceField(
        dtype=str,
        doc="""Operation to use for combining subimages into new image.""",
        default="copy",
        allowed={
            "none": """simply return a list of values and not re-mapping results into
                       a new image""",
            "copy": "copy pixels directly from subimage (non-deterministic for overlaps)",
            "sum": "add pixels from overlaps (probably never wanted; used for testing)",
            "average": "average pixels from overlaps"
        }
    )

    ignoreMaskPlanes = pexConfig.ListField(
        dtype=str,
        doc="""Mask planes to ignore for sigma-clipped statistics""",
        default=("INTRP", "EDGE", "DETECTED", "SAT", "CR", "BAD", "NO_DATA", "DETECTED_NEGATIVE")
    )


## \addtogroup LSST_task_documentation
## \{
## \page ImageGridderTask
## \ref ImageGridderTask_ "ImageGridderTask"
##      Task for performing operations on an image over a regular-spaced grid
## \}


class ImageGridderTask(pipeBase.Task):
    """!
    \anchor ImageGridderTask_

    \brief Break an image in to subimages on a grid and perform
    the same operation on each.

    \section ip_diffim_imageGridder_ImageGridderTask_Contents Contents

      - \ref ip_diffim_imageGridder_ImageGridderTask_Purpose
      - \ref ip_diffim_imageGridder_ImageGridderTask_Config
      - \ref ip_diffim_imageGridder_ImageGridderTask_Run
      - \ref ip_diffim_imageGridder_ImageGridderTask_Debug
      - \ref ip_diffim_imageGridder_ImageGridderTask_Example

    \section ip_diffim_imageGridder_ImageGridderTask_Purpose	Description

    Task that can perform 'simple' operations on a gridded set of
    subimages of a larger image, and then have those subimages
    stitched back together into a new, full-sized modified image.

    The actual operation is performed by a subTask passed to the
    config. The input exposure will be pre-subimaged, but the `run`
    method of the subtask will also have access to the entire
    (original) image.

    \section ip_diffim_imageGridder_ImageGridderTask_Initialize       Task initialization

    \copydoc \_\_init\_\_

    \section ip_diffim_imageGridder_ImageGridderTask_Run       Invoking the Task

    \copydoc run

    \section ip_diffim_imageGridder_ImageGridderTask_Config       Configuration parameters

    See \ref ImageGridderConfig

    \section ip_diffim_imageGridder_ImageGridderTask_Debug		Debug variables

    This task has no debug variables

    \section ip_diffim_imageGridder_ImageGridderTask_Example	Example of using ImageGridderTask

    This task has an example simple implementation of the use of
    ImageGridSubtask/Config in its unit test
    \link tests/testImageGridder testImageGridder\endlink.
    """
    ConfigClass = ImageGridderConfig
    _DefaultName = "ip_diffim_imageGridder"

    def __init__(self, *args, **kwargs):
        """! Create the image gridding Task
        @param *args arguments to be passed to
        lsst.pipe.base.task.Task.__init__
        @param **kwargs keyword arguments to be passed to
        lsst.pipe.base.task.Task.__init__
        """
        pipeBase.Task.__init__(self, *args, **kwargs)

        self.makeSubtask("gridSubtask")
        self.statsControl = afwMath.StatisticsControl()
        self.statsControl.setNumSigmaClip(3.)
        self.statsControl.setNumIter(3)
        self.statsControl.setAndMask(afwImage.MaskU.getPlaneBitMask(self.config.ignoreMaskPlanes))

    @pipeBase.timeMethod
    def run(self, exposure, returnPatches=False, **kwargs):
        """! Perform an operation on the given exposure.

        Break the exposure into sub-exps on a grid (parameters given
        by `ImageGridderConfig`) and perform `runSubtask` on
        each. Stitch together the resulting sub-exps generated by (or
        modified by) `runSubtask` into a final exposure of the same
        dimensions as the input `exposure`.

        @param[in] exposure the full exposure to process
        @return anything, including an `afw.Exposure`

        """
        self.log.info("Processing.")

        patches = self._makePatches(exposure, **kwargs)
        if returnPatches:
            return patches

        newMI = self._reduceImage(patches, exposure, **kwargs)
        return newMI

    def runSubtask(self, subExp, expandedSubExp, fullBBox, **kwargs):
        """! Run subtask to perform an operation on the given subExposure.

        Return an exposure of the same dimensions as `subExp`.
        Can use `expandedSubExp`, an expanded region of the original exposure,
        to perform computations. `subExp` and `expandedSubExp` should be considered
        *read-only* -- if you want to make changes to `subExp` it is recommended that
        you clone it first, and return the cloned, modified sub-exposure. This
        will be done in `gridSubtask.run()`.

        @param[in] subExp the sub-exposure of `exposure` upon which to operate
        @param[in] expandedSubExp the expanded sub-exposure of `exposure` upon which to operate
        @param[in] fullBBox the bounding box of the original exposure
        @return a `afw.Image` or `afw.Exposure`
        """
        subExp = self.gridSubtask.run(subExp, expandedSubExp, fullBBox, **kwargs)
        return subExp

    def _makePatches(self, exposure, doClone=False, **kwargs):
        """! Perform `runSubtask` on each patch

        Perform `runSubtask` on each patch across a grid on `exposure` generated by
        `generateGrid`.

        @param[in] exposure the original exposure which is used as the template.
        @param[in] doClone if True, clone the subimages before passing to subtask
        In this case, the subimages do not have to be considered as read-only.
        @return a list of `afwExposure`s
        """
        boxes0, boxes1 = self._generateGrid(exposure)
        if len(boxes0) != len(boxes1):
            raise Exception('Uh oh!')   # TBD: define a specific exception to raise

        self.log.info("Processing %d sub-exposures" % len(boxes0))
        patches = []
        for i in range(len(boxes0)):
            subExp = afwImage.ExposureF(exposure, boxes0[i])
            expandedSubExp = afwImage.ExposureF(exposure, boxes1[i])
            if doClone:
                subExp = subExp.clone()
                expandedSubExp = expandedSubExp.clone()
            result = self.runSubtask(subExp, expandedSubExp, exposure.getBBox(), **kwargs)
            patches.append(result)

        return patches

    def _reduceImage(self, patches, exposure, **kwargs):
        """! Reduce a set of sub-exposure patches into a final exposure

        Return an exposure of the same dimensions as `exposure`.
        `patches` is expected to have been produced by `makePatches`.

        @param[in] patches list of subExposures of `exposure` to merge
        @param[in] exposure the original exposure which is used as the template
        @return a `afw.Exposure`

        @notes Notes and currently known issues:
        1. This currently *should* correctly handle overlapping patches.
           For overlapping patches, use `config.reduceOperation='average'.
        2. This currently does not correctly handle varying PSFs (in fact,
           it just copies over the PSF from the original exposure)
        3. This logic currently makes *two* copies of the original exposure
           (one here and one in `gridSubtask.run()`).
        """
        if self.config.reduceOperation == 'none':
            return patches  # this is likely a non-image-type, such as a list of floats.

        newExp = afwImage.ExposureF(exposure).clone()
        newMI = newExp.getMaskedImage()
        newMI.getImage()[:, :] = 0.
        newMI.getVariance()[:, :] = 0.

        reduceOp = self.config.reduceOperation
        if reduceOp == 'average':  # make an array to keep track of weights
            weights = afwImage.ImageF(newMI.getBBox())  # Needs to be a float to divide later.

        for patch in patches:
            subExp = afwImage.ExposureF(newExp, patch.getBBox())
            subMI = subExp.getMaskedImage()
            patchMI = patch.getMaskedImage()
            if reduceOp == 'copy':
                subMI[:, :] = patchMI
            elif reduceOp == 'sum':
                subMI += patchMI
            elif reduceOp == 'average':
                subMI += patchMI
                wsubim = afwImage.ImageF(weights, patch.getBBox())
                wsubim += 1.

        if reduceOp == 'average':
            newMI /= weights
        return newExp

    def _generateGrid(self, exposure, forceEvenSized=False):
        """! Generate two lists of bounding boxes that evenly grid `exposure`

        Grid (subimage) centers will be spaced by gridStepX/Y. Then
        the grid will be adjusted as little as possible to evenly
        cover the input exposure (if adjustGridOption is not
        'none').  Then the bounding boxes will be expanded by
        borderSizeX/Y. The expanded bounding boxes will be adjusted to
        ensure that they intersect the exposure's bounding box.  The
        resulting lists of bounding boxes and corresponding expanded
        bounding boxes will be returned, and also set to
        `self.boxes0`, `self.boxes1`.

        @param[in] exposure an `afwImage.Exposure` whose full bounding
        box is to be evenly gridded.
        @param[in] forceEvenSized force grid elements to have even x- and
        y- dimensions (for ZOGY)

        @return tuple containing two lists of `afwGeom.BoundingBox`es
        """

        # Extract the config parameters for conciseness.
        gridSizeX = self.config.gridSizeX
        gridSizeY = self.config.gridSizeY
        gridStepX = self.config.gridStepX
        gridStepY = self.config.gridStepY
        borderSizeX = self.config.borderSizeX
        borderSizeY = self.config.borderSizeY
        adjustGridOption = self.config.adjustGridOption
        scaleByFwhm = self.config.scaleByFwhm
        bbox = exposure.getBBox()

        psfFwhm = (exposure.getPsf().computeShape().getDeterminantRadius() *
                   2. * np.sqrt(2. * np.log(2.)))

        def rescaleValue(val):
            if scaleByFwhm:
                return np.rint(val * psfFwhm).astype(int)
            else:
                return np.rint(val).astype(int)

        gridSizeX = rescaleValue(gridSizeX)
        if gridSizeX > bbox.getWidth():
            gridSizeX = bbox.getWidth()
        gridSizeY = rescaleValue(gridSizeY)
        if gridSizeY > bbox.getHeight():
            gridSizeY = bbox.getHeight()
        gridStepX = rescaleValue(gridStepX)
        if gridStepX > bbox.getWidth():
            gridStepX = bbox.getWidth()
        gridStepY = rescaleValue(gridStepY)
        if gridStepY > bbox.getHeight():
            gridStepY = bbox.getHeight()
        borderSizeX = rescaleValue(borderSizeX)
        borderSizeY = rescaleValue(borderSizeY)

        nGridX = bbox.getWidth() // gridStepX
        nGridY = bbox.getHeight() // gridStepY

        if adjustGridOption == 'spacing':
            nGridX = bbox.getWidth() / gridStepX
            # Readjust gridStepX so that it fits perfectly in the image.
            gridStepX = float(bbox.getWidth() - gridSizeX) / float(nGridX)
            if gridStepX <= 0:
                gridStepX = 1

        if adjustGridOption == 'spacing':
            nGridY = bbox.getWidth() / gridStepY
            # Readjust gridStepY so that it fits perfectly in the image.
            gridStepY = float(bbox.getHeight() - gridSizeY) / float(nGridY)
            if gridStepY <= 0:
                gridStepY = 1

        # first "main" box at 0,0
        bbox0 = afwGeom.Box2I(afwGeom.Point2I(bbox.getBegin()), afwGeom.Extent2I(gridSizeX, gridSizeY))
        # first expanded box
        bbox1 = afwGeom.Box2I(bbox0)
        bbox1.grow(afwGeom.Extent2I(borderSizeX, borderSizeY))

        self.boxes0 = []  # "main" boxes; store in task so can be extracted if needed
        self.boxes1 = []  # "expanded" boxes

        # use given centroids as centers for bounding boxes
        if self.config.gridCentroidsX is not None and len(self.config.gridCentroidsX) > 0:
            for i, centroidX in enumerate(self.config.gridCentroidsX):
                centroidY = self.config.gridCentroidsY[i]
                centroid = afwGeom.Point2D(centroidX, centroidY)
                bb0 = afwGeom.Box2I(bbox0)
                xoff = int(np.floor(centroid.getX())) - bb0.getWidth()//2
                yoff = int(np.floor(centroid.getY())) - bb0.getHeight()//2
                bb0.shift(afwGeom.Extent2I(xoff, yoff))
                bb0.clip(bbox)
                self.boxes0.append(bb0)
                bb1 = afwGeom.Box2I(bbox1)
                bb1.shift(afwGeom.Extent2I(xoff, yoff))
                bb1.clip(bbox)
                self.boxes1.append(bb1)
            return self.boxes0, self.boxes1

        # Offset the "main" (bbox0) and "expanded" (bbox1) bboxes by xoff, yoff. Clip them by the
        # exposure's bbox.
        def offsetAndClipBoxes(bbox0, bbox1, xoff, yoff, bbox):
            xoff = int(np.floor(xoff))
            yoff = int(np.floor(yoff))
            bb0 = afwGeom.Box2I(bbox0)
            bb0.shift(afwGeom.Extent2I(xoff, yoff))
            bb0.clip(bbox)
            bb1 = afwGeom.Box2I(bbox1)
            bb1.shift(afwGeom.Extent2I(xoff, yoff))
            bb1.clip(bbox)
            if forceEvenSized:
                bb0.grow(afwGeom.Extent2I(bb0.getWidth() % 2, bb0.getHeight() % 2))
                bb1.grow(afwGeom.Extent2I(bb1.getWidth() % 2, bb1.getHeight() % 2))
            return bb0, bb1

        xoff = 0
        while(xoff <= bbox.getWidth()):
            yoff = 0
            while(yoff <= bbox.getHeight()):
                bb0, bb1 = offsetAndClipBoxes(bbox0, bbox1, xoff, yoff, bbox)
                yoff += gridStepY
                if bb0.getArea() > 0 and bb1.getArea() > 0:
                    self.boxes0.append(bb0)
                    self.boxes1.append(bb1)
            xoff += gridStepX

        return self.boxes0, self.boxes1

    def _plotBoxGrid(self, boxes, bbox, **kwargs):
        """! Plot a grid of boxes using matplotlib.

        @param[in] boxes a list of `afwGeom.BoundginBox`es
        @param[in] bbox an overall bounding box
        """
        import matplotlib.pyplot as plt

        def plotBox(box):
            corners = np.array([np.array([pt.getX(), pt.getY()]) for pt in box.getCorners()])
            corners = np.vstack([corners, corners[0, :]])
            plt.plot(corners[:, 0], corners[:, 1], **kwargs)

        for b in boxes:
            plotBox(b)
        plt.xlim(bbox.getBeginX(), bbox.getEndX())
        plt.ylim(bbox.getBeginY(), bbox.getEndY())

    def _plotBoxes(self, exposure):
        """! Plot both grids of boxes using matplotlib.

        `self.boxes0` and `self.boxes1` must have been set.
        @param[in] bbox an overall bounding box
        """
        import matplotlib.pyplot as plt

        bbox = exposure.getBBox()
        boxes0, boxes1 = self._generateGrid(exposure)
        self._plotBoxGrid(boxes0[::3], bbox, ls='--')
        # reset the color cycle -- see
        # http://stackoverflow.com/questions/24193174/reset-color-cycle-in-matplotlib
        plt.gca().set_prop_cycle(None)
        self._plotBoxGrid(boxes1[::3], bbox, ls=':')

    def _computeVarianceMean(self, subImage):
        """! Utility function: compute mean of variance plane of subimage

        @param[in] subImage the sub-image of `exposure` upon which to operate
        @return float clipped mean of masked variance plane of subImage
        """
        statObj = afwMath.makeStatistics(subImage.getMaskedImage().getVariance(),
                                         subImage.getMaskedImage().getMask(),
                                         afwMath.MEANCLIP, self.statsControl)
        var = statObj.getValue(afwMath.MEANCLIP)
        return var

    def _computePsf(self, subImage):
        """! Utility function: compute Psf at center of subImage.

        TBD: is this computing the Psf at the center of the subimage
        (i.e. center of its bounding box)?

        @param[in] subImage the sub-image of `exposure` upon which to operate
        @return 2d numpy.array of Psf for calculations.

        """
        return subImage.getPsf().computeImage().getArray()
