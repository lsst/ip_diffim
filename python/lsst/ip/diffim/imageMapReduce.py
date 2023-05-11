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
import abc

import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.geom as geom
import lsst.meas.algorithms as measAlg
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod

__all__ = ("ImageMapReduceTask", "ImageMapReduceConfig",
           "ImageMapper", "ImageMapperConfig",
           "ImageReducer", "ImageReducerConfig")


"""Tasks for processing an exposure via processing on
multiple sub-exposures and then collecting the results
to either re-stitch the sub-exposures back into a new
exposure, or return summary results for each sub-exposure.

This provides a framework for arbitrary mapper-reducer
operations on an exposure by implementing simple operations in
subTasks. It currently is not parallelized, although it could be in
the future. It does enable operations such as spatially-mapped
processing on a grid across an image, processing regions surrounding
centroids (such as for PSF processing), etc.

It is implemented as primary Task, `ImageMapReduceTask` which contains
two subtasks, `ImageMapper` and `ImageReducer`.
`ImageMapReduceTask` configures the centroids and sub-exposure
dimensions to be processed, and then calls the `run` methods of the
`ImageMapper` and `ImageReducer` on those sub-exposures.
`ImageMapReduceTask` may be configured with a list of sub-exposure
centroids (`config.cellCentroidsX` and `config.cellCentroidsY`) and a
single pair of bounding boxes defining their dimensions, or a set of
parameters defining a regular grid of centroids (`config.gridStepX`
and `config.gridStepY`).

`ImageMapper` is an abstract class and must be subclassed with
an implemented `run` method to provide the desired operation for
processing individual sub-exposures. It is called from
`ImageMapReduceTask.run`, and may return a new, processed sub-exposure
which is to be "stitched" back into a new resulting larger exposure
(depending on the configured `ImageMapReduceTask.mapper`);
otherwise if it does not return an lsst.afw.image.Exposure, then the results are
passed back directly to the caller.

`ImageReducer` will either stitch the `mapperResults` list of
results generated by the `ImageMapper` together into a new
Exposure (by default) or pass it through to the
caller. `ImageReducer` has an implemented `run` method for
basic reducing operations (`reduceOperation`) such as `average` (which
will average all overlapping pixels from sub-exposures produced by the
`ImageMapper` into the new exposure). Another notable
implemented `reduceOperation` is 'none', in which case the
`mapperResults` list is simply returned directly.
"""


class ImageMapperConfig(pexConfig.Config):
    """Configuration parameters for ImageMapper
    """
    pass


class ImageMapper(pipeBase.Task, metaclass=abc.ABCMeta):
    """Abstract base class for any task that is to be
    used as `ImageMapReduceConfig.mapper`.

    Notes
    -----
    An `ImageMapper` is responsible for processing individual
    sub-exposures in its `run` method, which is called from
    `ImageMapReduceTask.run`. `run` may return a processed new
    sub-exposure which can be be "stitched" back into a new resulting
    larger exposure (depending on the configured
    `ImageReducer`); otherwise if it does not return an
    lsst.afw.image.Exposure, then the
    `ImageReducer.config.reducer.reduceOperation`
    should be set to 'none' and the result will be propagated
    as-is.
    """
    ConfigClass = ImageMapperConfig
    _DefaultName = "ip_diffim_ImageMapper"

    @abc.abstractmethod
    def run(self, subExposure, expandedSubExposure, fullBBox, **kwargs):
        """Perform operation on `subExposure`.

        To be implemented by subclasses. See class docstring for more
        details. This method is given the `subExposure` which
        is to be operated upon, and an `expandedSubExposure` which
        will contain `subExposure` with additional surrounding
        pixels. This allows for, for example, convolutions (which
        should be performed on `expandedSubExposure`), to prevent the
        returned sub-exposure from containing invalid pixels.

        This method may return a new, processed sub-exposure which can
        be be "stitched" back into a new resulting larger exposure
        (depending on the paired, configured `ImageReducer`);
        otherwise if it does not return an lsst.afw.image.Exposure, then the
        `ImageReducer.config.mapper.reduceOperation`
        should be set to 'none' and the result will be propagated
        as-is.

        Parameters
        ----------
        subExposure : `lsst.afw.image.Exposure`
            the sub-exposure upon which to operate
        expandedSubExposure : `lsst.afw.image.Exposure`
            the expanded sub-exposure upon which to operate
        fullBBox : `lsst.geom.Box2I`
            the bounding box of the original exposure
        kwargs :
            additional keyword arguments propagated from
            `ImageMapReduceTask.run`.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            A structure containing the result of the `subExposure` processing,
            which may itself be of any type. See above for details. If it is an
            `lsst.afw.image.Exposure` (processed sub-exposure), then the name in
            the Struct should be 'subExposure'. This is implemented here as a
            pass-through example only.
        """
        return pipeBase.Struct(subExposure=subExposure)


class ImageReducerConfig(pexConfig.Config):
    """Configuration parameters for the ImageReducer
    """
    reduceOperation = pexConfig.ChoiceField(
        dtype=str,
        doc="""Operation to use for reducing subimages into new image.""",
        default="average",
        allowed={
            "none": """simply return a list of values and don't re-map results into
                       a new image (noop operation)""",
            "copy": """copy pixels directly from subimage into correct location in
                       new exposure (potentially non-deterministic for overlaps)""",
            "sum": """add pixels from overlaps (probably never wanted; used for testing)
                       into correct location in new exposure""",
            "average": """same as copy, but also average pixels from overlapped regions
                       (NaNs ignored)""",
            "coaddPsf": """Instead of constructing an Exposure, take a list of returned
                       PSFs and use CoaddPsf to construct a single PSF that covers the
                       entire input exposure""",
        }
    )
    badMaskPlanes = pexConfig.ListField(
        dtype=str,
        doc="""Mask planes to set for invalid pixels""",
        default=('INVALID_MAPREDUCE', 'BAD', 'NO_DATA')
    )


class ImageReducer(pipeBase.Task):
    """Base class for any 'reduce' task that is to be
    used as `ImageMapReduceConfig.reducer`.

    Basic reduce operations are provided by the `run` method
    of this class, to be selected by its config.
    """
    ConfigClass = ImageReducerConfig
    _DefaultName = "ip_diffim_ImageReducer"

    def run(self, mapperResults, exposure, **kwargs):
        """Reduce a list of items produced by `ImageMapper`.

        Either stitch the passed `mapperResults` list
        together into a new Exposure (default) or pass it through
        (if `self.config.reduceOperation` is 'none').

        If `self.config.reduceOperation` is not 'none', then expect
        that the `pipeBase.Struct`s in the `mapperResults` list
        contain sub-exposures named 'subExposure', to be stitched back
        into a single Exposure with the same dimensions, PSF, and mask
        as the input `exposure`. Otherwise, the `mapperResults` list
        is simply returned directly.

        Parameters
        ----------
        mapperResults : `list`
            list of `lsst.pipe.base.Struct` returned by `ImageMapper.run`.
        exposure : `lsst.afw.image.Exposure`
            the original exposure which is cloned to use as the
            basis for the resulting exposure (if
            ``self.config.mapper.reduceOperation`` is not 'None')
        kwargs :
            additional keyword arguments propagated from
            `ImageMapReduceTask.run`.

        Returns
        -------
        A `lsst.pipe.base.Struct` containing either an `lsst.afw.image.Exposure`
        (named 'exposure') or a list (named 'result'),
        depending on `config.reduceOperation`.

        Notes
        -----
        1. This currently correctly handles overlapping sub-exposures.
           For overlapping sub-exposures, use `config.reduceOperation='average'`.
        2. This correctly handles varying PSFs, constructing the resulting
           exposure's PSF via CoaddPsf (DM-9629).

        Known issues

        1. To be done: correct handling of masks (nearly there)
        2. This logic currently makes *two* copies of the original exposure
           (one here and one in `mapper.run()`). Possibly of concern
           for large images on memory-constrained systems.
        """
        # No-op; simply pass mapperResults directly to ImageMapReduceTask.run
        if self.config.reduceOperation == 'none':
            return pipeBase.Struct(result=mapperResults)

        if self.config.reduceOperation == 'coaddPsf':
            # Each element of `mapperResults` should contain 'psf' and 'bbox'
            coaddPsf = self._constructPsf(mapperResults, exposure)
            return pipeBase.Struct(result=coaddPsf)

        newExp = exposure.clone()
        newMI = newExp.getMaskedImage()

        reduceOp = self.config.reduceOperation
        if reduceOp == 'copy':
            weights = None
            newMI.image.array[:, :] = np.nan
            newMI.variance.array[:, :] = np.nan
        else:
            newMI.image.array[:, :] = 0.
            newMI.variance.array[:, :] = 0.
            if reduceOp == 'average':  # make an array to keep track of weights
                weights = afwImage.ImageI(newMI.getBBox())

        for item in mapperResults:
            item = item.subExposure  # Expected named value in the pipeBase.Struct
            if not (isinstance(item, afwImage.ExposureF) or isinstance(item, afwImage.ExposureI)
                    or isinstance(item, afwImage.ExposureU) or isinstance(item, afwImage.ExposureD)):
                raise TypeError("""Expecting an Exposure type, got %s.
                                   Consider using `reduceOperation="none".""" % str(type(item)))
            subExp = newExp.Factory(newExp, item.getBBox())
            subMI = subExp.getMaskedImage()
            patchMI = item.getMaskedImage()
            isValid = ~np.isnan(patchMI.image.array * patchMI.variance.array)

            if reduceOp == 'copy':
                subMI.image.array[isValid] = patchMI.image.array[isValid]
                subMI.variance.array[isValid] = patchMI.variance.array[isValid]
                subMI.mask.array[:, :] |= patchMI.mask.array

            if reduceOp == 'sum' or reduceOp == 'average':  # much of these two options is the same
                subMI.image.array[isValid] += patchMI.image.array[isValid]
                subMI.variance.array[isValid] += patchMI.variance.array[isValid]
                subMI.mask.array[:, :] |= patchMI.mask.array
                if reduceOp == 'average':
                    # wtsView is a view into the `weights` Image
                    wtsView = afwImage.ImageI(weights, item.getBBox())
                    wtsView.array[isValid] += 1

        # New mask plane - for debugging map-reduced images
        mask = newMI.mask
        for m in self.config.badMaskPlanes:
            mask.addMaskPlane(m)
        bad = mask.getPlaneBitMask(self.config.badMaskPlanes)

        isNan = np.where(np.isnan(newMI.image.array * newMI.variance.array))
        if len(isNan[0]) > 0:
            # set mask to INVALID for pixels where produced exposure is NaN
            mask.array[isNan[0], isNan[1]] |= bad

        if reduceOp == 'average':
            wts = weights.array.astype(float)
            self.log.info('AVERAGE: Maximum overlap: %f', np.nanmax(wts))
            self.log.info('AVERAGE: Average overlap: %f', np.nanmean(wts))
            self.log.info('AVERAGE: Minimum overlap: %f', np.nanmin(wts))
            wtsZero = np.equal(wts, 0.)
            wtsZeroInds = np.where(wtsZero)
            wtsZeroSum = len(wtsZeroInds[0])
            self.log.info('AVERAGE: Number of zero pixels: %f (%f%%)', wtsZeroSum,
                          wtsZeroSum * 100. / wtsZero.size)
            notWtsZero = ~wtsZero
            tmp = newMI.image.array
            np.divide(tmp, wts, out=tmp, where=notWtsZero)
            tmp = newMI.variance.array
            np.divide(tmp, wts, out=tmp, where=notWtsZero)
            if len(wtsZeroInds[0]) > 0:
                newMI.image.array[wtsZeroInds] = np.nan
                newMI.variance.array[wtsZeroInds] = np.nan
                # set mask to something for pixels where wts == 0.
                # happens sometimes if operation failed on a certain subexposure
                mask.array[wtsZeroInds] |= bad

        # Not sure how to construct a PSF when reduceOp=='copy'...
        if reduceOp == 'sum' or reduceOp == 'average':
            psf = self._constructPsf(mapperResults, exposure)
            newExp.setPsf(psf)

        return pipeBase.Struct(exposure=newExp)

    def _constructPsf(self, mapperResults, exposure):
        """Construct a CoaddPsf based on PSFs from individual subExposures

        Currently uses (and returns) a CoaddPsf. TBD if we want to
        create a custom subclass of CoaddPsf to differentiate it.

        Parameters
        ----------
        mapperResults : `list`
            list of `pipeBase.Struct` returned by `ImageMapper.run`.
            For this to work, each element of `mapperResults` must contain
            a `subExposure` element, from which the component Psfs are
            extracted (thus the reducerTask cannot have
            `reduceOperation = 'none'`.
        exposure : `lsst.afw.image.Exposure`
            the original exposure which is used here solely for its
            bounding-box and WCS.

        Returns
        -------
        psf : `lsst.meas.algorithms.CoaddPsf`
            A psf constructed from the PSFs of the individual subExposures.
        """
        schema = afwTable.ExposureTable.makeMinimalSchema()
        schema.addField("weight", type="D", doc="Coadd weight")
        mycatalog = afwTable.ExposureCatalog(schema)

        # We're just using the exposure's WCS (assuming that the subExposures'
        # WCSs are the same, which they better be!).
        wcsref = exposure.getWcs()
        for i, res in enumerate(mapperResults):
            record = mycatalog.getTable().makeRecord()
            if 'subExposure' in res.getDict():
                subExp = res.subExposure
                if subExp.getWcs() != wcsref:
                    raise ValueError('Wcs of subExposure is different from exposure')
                record.setPsf(subExp.getPsf())
                record.setWcs(subExp.getWcs())
                record.setBBox(subExp.getBBox())
            elif 'psf' in res.getDict():
                record.setPsf(res.psf)
                record.setWcs(wcsref)
                record.setBBox(res.bbox)
            record['weight'] = 1.0
            record['id'] = i
            mycatalog.append(record)

        # create the coaddpsf
        psf = measAlg.CoaddPsf(mycatalog, wcsref, 'weight')
        return psf


class ImageMapReduceConfig(pexConfig.Config):
    """Configuration parameters for the ImageMapReduceTask
    """
    mapper = pexConfig.ConfigurableField(
        doc="Task to run on each subimage",
        target=ImageMapper,
    )

    reducer = pexConfig.ConfigurableField(
        doc="Task to combine results of mapper task",
        target=ImageReducer,
    )

    # Separate cellCentroidsX and cellCentroidsY since pexConfig.ListField accepts limited dtypes
    #  (i.e., no Point2D). The resulting set of centroids is the "vertical stack" of
    #  `cellCentroidsX` and `cellCentroidsY`, i.e. for (1,2), (3,4) respectively, the
    #   resulting centroids are ((1,3), (2,4)).
    cellCentroidsX = pexConfig.ListField(
        dtype=float,
        doc="""Input X centroids around which to place subimages.
               If None, use grid config options below.""",
        optional=True,
        default=None
    )

    cellCentroidsY = pexConfig.ListField(
        dtype=float,
        doc="""Input Y centroids around which to place subimages.
               If None, use grid config options below.""",
        optional=True,
        default=None
    )

    cellSizeX = pexConfig.Field(
        dtype=float,
        doc="""Dimensions of each grid cell in x direction""",
        default=10.,
        check=lambda x: x > 0.
    )

    cellSizeY = pexConfig.Field(
        dtype=float,
        doc="""Dimensions of each grid cell in y direction""",
        default=10.,
        check=lambda x: x > 0.
    )

    gridStepX = pexConfig.Field(
        dtype=float,
        doc="""Spacing between subsequent grid cells in x direction. If equal to
               cellSizeX, then there is no overlap in the x direction.""",
        default=10.,
        check=lambda x: x > 0.
    )

    gridStepY = pexConfig.Field(
        dtype=float,
        doc="""Spacing between subsequent grid cells in y direction. If equal to
               cellSizeY, then there is no overlap in the y direction.""",
        default=10.,
        check=lambda x: x > 0.
    )

    borderSizeX = pexConfig.Field(
        dtype=float,
        doc="""Dimensions of grid cell border in +/- x direction, to be used
               for generating `expandedSubExposure`.""",
        default=5.,
        check=lambda x: x > 0.
    )

    borderSizeY = pexConfig.Field(
        dtype=float,
        doc="""Dimensions of grid cell border in +/- y direction, to be used
               for generating `expandedSubExposure`.""",
        default=5.,
        check=lambda x: x > 0.
    )

    adjustGridOption = pexConfig.ChoiceField(
        dtype=str,
        doc="""Whether and how to adjust grid to fit evenly within, and cover entire
               image""",
        default="spacing",
        allowed={
            "spacing": "adjust spacing between centers of grid cells (allowing overlaps)",
            "size": "adjust the sizes of the grid cells (disallowing overlaps)",
            "none": "do not adjust the grid sizes or spacing"
        }
    )

    scaleByFwhm = pexConfig.Field(
        dtype=bool,
        doc="""Scale cellSize/gridStep/borderSize/overlapSize by PSF FWHM rather
               than pixels?""",
        default=True
    )

    returnSubImages = pexConfig.Field(
        dtype=bool,
        doc="""Return the input subExposures alongside the processed ones (for debugging)""",
        default=False
    )

    ignoreMaskPlanes = pexConfig.ListField(
        dtype=str,
        doc="""Mask planes to ignore for sigma-clipped statistics""",
        default=("INTRP", "EDGE", "DETECTED", "SAT", "CR", "BAD", "NO_DATA", "DETECTED_NEGATIVE")
    )


class ImageMapReduceTask(pipeBase.Task):
    """Split an Exposure into subExposures (optionally on a grid) and
    perform the same operation on each.

    Perform 'simple' operations on a gridded set of subExposures of a
    larger Exposure, and then (by default) have those subExposures
    stitched back together into a new, full-sized image.

    Contrary to the expectation given by its name, this task does not
    perform these operations in parallel, although it could be updatd
    to provide such functionality.

    The actual operations are performed by two subTasks passed to the
    config. The exposure passed to this task's `run` method will be
    divided, and those subExposures will be passed to the subTasks,
    along with the original exposure. The reducing operation is
    performed by the second subtask.
    """
    ConfigClass = ImageMapReduceConfig
    _DefaultName = "ip_diffim_imageMapReduce"

    def __init__(self, *args, **kwargs):
        """Create the image map-reduce task

        Parameters
        ----------
        args :
            arguments to be passed to
            `lsst.pipe.base.task.Task.__init__`
        kwargs :
            additional keyword arguments to be passed to
            `lsst.pipe.base.task.Task.__init__`
        """
        pipeBase.Task.__init__(self, *args, **kwargs)

        self.boxes0 = self.boxes1 = None
        self.makeSubtask("mapper")
        self.makeSubtask("reducer")

    @timeMethod
    def run(self, exposure, **kwargs):
        """Perform a map-reduce operation on the given exposure.

        Split the exposure into sub-expposures on a grid (parameters
        given by `ImageMapReduceConfig`) and perform
        `config.mapper.run()` on each. Reduce the resulting
        sub-exposures by running `config.reducer.run()`.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            the full exposure to process
        kwargs :
            additional keyword arguments to be passed to
            subtask `run` methods

        Returns
        -------
        output of `reducer.run()`

        """
        self.log.info("Mapper sub-task: %s", self.mapper._DefaultName)
        mapperResults = self._runMapper(exposure, **kwargs)
        self.log.info("Reducer sub-task: %s", self.reducer._DefaultName)
        result = self._reduceImage(mapperResults, exposure, **kwargs)
        return result

    def _runMapper(self, exposure, doClone=False, **kwargs):
        """Perform `mapper.run` on each sub-exposure

        Perform `mapper.run` on each sub-exposure across a
        grid on `exposure` generated by `_generateGrid`. Also pass to
        `mapper.run` an 'expanded sub-exposure' containing the
        same region as the sub-exposure but with an expanded bounding box.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            the original exposure which is used as the template
        doClone : `bool`
            if True, clone the subimages before passing to subtask;
            in that case, the sub-exps do not have to be considered as read-only
        kwargs :
            additional keyword arguments to be passed to
            `mapper.run` and `self._generateGrid`, including `forceEvenSized`.

        Returns
        -------
        a list of `pipeBase.Struct`s as returned by `mapper.run`.
        """
        if self.boxes0 is None:
            self._generateGrid(exposure, **kwargs)  # possibly pass `forceEvenSized`
        if len(self.boxes0) != len(self.boxes1):
            raise ValueError('Bounding boxes list and expanded bounding boxes list are of different lengths')

        self.log.info("Processing %d sub-exposures", len(self.boxes0))
        mapperResults = []
        for box0, box1 in zip(self.boxes0, self.boxes1):
            subExp = exposure.Factory(exposure, box0)
            expandedSubExp = exposure.Factory(exposure, box1)
            if doClone:
                subExp = subExp.clone()
                expandedSubExp = expandedSubExp.clone()
            result = self.mapper.run(subExp, expandedSubExp, exposure.getBBox(), **kwargs)
            if self.config.returnSubImages:
                toAdd = pipeBase.Struct(inputSubExposure=subExp,
                                        inputExpandedSubExposure=expandedSubExp)
                result.mergeItems(toAdd, 'inputSubExposure', 'inputExpandedSubExposure')
            mapperResults.append(result)

        return mapperResults

    def _reduceImage(self, mapperResults, exposure, **kwargs):
        """Reduce/merge a set of sub-exposures into a final result

        Return an exposure of the same dimensions as `exposure`.
        `mapperResults` is expected to have been produced by `runMapper`.

        Parameters
        ----------
        mapperResults : `list`
            `list` of `lsst.pipe.base.Struct`, each of which was produced by
            `config.mapper`
        exposure : `lsst.afw.image.Exposure`
            the original exposure
        **kwargs
            additional keyword arguments

        Returns
        -------
        Output of `reducer.run` which is a `pipeBase.Struct`.
        """
        result = self.reducer.run(mapperResults, exposure, **kwargs)
        return result

    def _generateGrid(self, exposure, forceEvenSized=False, **kwargs):
        """Generate two lists of bounding boxes that evenly grid `exposure`

        Unless the config was provided with `cellCentroidsX` and
        `cellCentroidsY`, grid (subimage) centers are spaced evenly
        by gridStepX/Y. Then the grid is adjusted as little as
        possible to evenly cover the input exposure (if
        adjustGridOption is not 'none'). Then the second set of
        bounding boxes is expanded by borderSizeX/Y. The expanded
        bounding boxes are adjusted to ensure that they intersect the
        exposure's bounding box. The resulting lists of bounding boxes
        and corresponding expanded bounding boxes are set to
        `self.boxes0`, `self.boxes1`.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            input exposure whose full bounding box is to be evenly gridded.
        forceEvenSized : `bool`
            force grid elements to have even-valued x- and y- dimensions?
            (Potentially useful if doing Fourier transform of subExposures.)
        """
        # kwargs are ignored, but necessary to enable optional passing of
        # `forceEvenSized` from `_runMapper`.
        bbox = exposure.getBBox()

        # Extract the config parameters for conciseness.
        cellCentroidsX = self.config.cellCentroidsX
        cellCentroidsY = self.config.cellCentroidsY
        cellSizeX = self.config.cellSizeX
        cellSizeY = self.config.cellSizeY
        gridStepX = self.config.gridStepX
        gridStepY = self.config.gridStepY
        borderSizeX = self.config.borderSizeX
        borderSizeY = self.config.borderSizeY
        adjustGridOption = self.config.adjustGridOption
        scaleByFwhm = self.config.scaleByFwhm

        if cellCentroidsX is None or len(cellCentroidsX) <= 0:
            # Not given centroids; construct them from cellSize/gridStep
            psf = exposure.getPsf()
            psfFwhm = (psf.computeShape(psf.getAveragePosition()).getDeterminantRadius()
                       * 2.*np.sqrt(2.*np.log(2.)))
            if scaleByFwhm:
                self.log.info("Scaling grid parameters by %f", psfFwhm)

            def rescaleValue(val):
                if scaleByFwhm:
                    return np.rint(val*psfFwhm).astype(int)
                else:
                    return np.rint(val).astype(int)

            cellSizeX = rescaleValue(cellSizeX)
            cellSizeY = rescaleValue(cellSizeY)
            gridStepX = rescaleValue(gridStepX)
            gridStepY = rescaleValue(gridStepY)
            borderSizeX = rescaleValue(borderSizeX)
            borderSizeY = rescaleValue(borderSizeY)

            nGridX = bbox.getWidth()//gridStepX
            nGridY = bbox.getHeight()//gridStepY

            if adjustGridOption == 'spacing':
                # Readjust spacings so that they fit perfectly in the image.
                nGridX = bbox.getWidth()//cellSizeX + 1
                nGridY = bbox.getHeight()//cellSizeY + 1
                xLinSpace = np.linspace(cellSizeX//2, bbox.getWidth() - cellSizeX//2, nGridX)
                yLinSpace = np.linspace(cellSizeY//2, bbox.getHeight() - cellSizeY//2, nGridY)

            elif adjustGridOption == 'size':
                cellSizeX = gridStepX
                cellSizeY = gridStepY
                xLinSpace = np.arange(cellSizeX//2, bbox.getWidth() + cellSizeX//2, cellSizeX)
                yLinSpace = np.arange(cellSizeY//2, bbox.getHeight() + cellSizeY//2, cellSizeY)
                cellSizeX += 1  # add 1 to make sure there are no gaps
                cellSizeY += 1

            else:
                xLinSpace = np.arange(cellSizeX//2, bbox.getWidth() + cellSizeX//2, gridStepX)
                yLinSpace = np.arange(cellSizeY//2, bbox.getHeight() + cellSizeY//2, gridStepY)

            cellCentroids = [(x, y) for x in xLinSpace for y in yLinSpace]

        else:
            # in py3 zip returns an iterator, but want to test length below, so use this instead:
            cellCentroids = [(cellCentroidsX[i], cellCentroidsY[i]) for i in range(len(cellCentroidsX))]

        # first "main" box at 0,0
        bbox0 = geom.Box2I(geom.Point2I(bbox.getBegin()), geom.Extent2I(cellSizeX, cellSizeY))
        # first expanded box
        bbox1 = geom.Box2I(bbox0)
        bbox1.grow(geom.Extent2I(borderSizeX, borderSizeY))

        self.boxes0 = []  # "main" boxes; store in task so can be extracted if needed
        self.boxes1 = []  # "expanded" boxes

        def _makeBoxEvenSized(bb):
            """Force a bounding-box to have dimensions that are modulo 2."""

            if bb.getWidth() % 2 == 1:  # grow to the right
                bb.include(geom.Point2I(bb.getMaxX()+1, bb.getMaxY()))  # Expand by 1 pixel!
                bb.clip(bbox)
                if bb.getWidth() % 2 == 1:  # clipped at right -- so grow to the left
                    bb.include(geom.Point2I(bb.getMinX()-1, bb.getMaxY()))
                    bb.clip(bbox)
            if bb.getHeight() % 2 == 1:  # grow upwards
                bb.include(geom.Point2I(bb.getMaxX(), bb.getMaxY()+1))  # Expand by 1 pixel!
                bb.clip(bbox)
                if bb.getHeight() % 2 == 1:  # clipped upwards -- so grow down
                    bb.include(geom.Point2I(bb.getMaxX(), bb.getMinY()-1))
                    bb.clip(bbox)
            if bb.getWidth() % 2 == 1 or bb.getHeight() % 2 == 1:  # Box is probably too big
                raise RuntimeError('Cannot make bounding box even-sized. Probably too big.')

            return bb

        # Use given or grid-parameterized centroids as centers for bounding boxes
        if cellCentroids is not None and len(cellCentroids) > 0:
            for x, y in cellCentroids:
                centroid = geom.Point2D(x, y)
                bb0 = geom.Box2I(bbox0)
                xoff = int(np.floor(centroid.getX())) - bb0.getWidth()//2
                yoff = int(np.floor(centroid.getY())) - bb0.getHeight()//2
                bb0.shift(geom.Extent2I(xoff, yoff))
                bb0.clip(bbox)
                if forceEvenSized:
                    bb0 = _makeBoxEvenSized(bb0)
                bb1 = geom.Box2I(bbox1)
                bb1.shift(geom.Extent2I(xoff, yoff))
                bb1.clip(bbox)
                if forceEvenSized:
                    bb1 = _makeBoxEvenSized(bb1)

                if bb0.getArea() > 1 and bb1.getArea() > 1:
                    self.boxes0.append(bb0)
                    self.boxes1.append(bb1)

        return self.boxes0, self.boxes1

    def plotBoxes(self, fullBBox, skip=3):
        """Plot both grids of boxes using matplotlib.

        Will compute the grid via `_generateGrid` if
        `self.boxes0` and `self.boxes1` have not already been set.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure whose bounding box is gridded by this task.
        skip : `int`
            Plot every skip-ped box (help make plots less confusing)
        """
        import matplotlib.pyplot as plt

        if self.boxes0 is None:
            raise RuntimeError('Cannot plot boxes. Run _generateGrid first.')
        self._plotBoxGrid(self.boxes0[::skip], fullBBox, ls='--')
        # reset the color cycle -- see
        # http://stackoverflow.com/questions/24193174/reset-color-cycle-in-matplotlib
        plt.gca().set_prop_cycle(None)
        self._plotBoxGrid(self.boxes1[::skip], fullBBox, ls=':')

    def _plotBoxGrid(self, boxes, bbox, **kwargs):
        """Plot a grid of boxes using matplotlib.

        Parameters
        ----------
        boxes : `list` of `lsst.geom.Box2I`
            a list of bounding boxes.
        bbox : `lsst.geom.Box2I`
            an overall bounding box
        **kwargs
            additional keyword arguments for matplotlib
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
