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
import lsst.geom as geom
import lsst.log
import lsst.meas.algorithms as measAlg
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase


from .imageMapReduce import (ImageMapReduceConfig, ImageMapReduceTask,
                             ImageMapper)

__all__ = ("DecorrelateALKernelTask", "DecorrelateALKernelConfig",
           "DecorrelateALKernelMapper", "DecorrelateALKernelMapReduceConfig",
           "DecorrelateALKernelSpatialConfig", "DecorrelateALKernelSpatialTask")


class DecorrelateALKernelConfig(pexConfig.Config):
    """Configuration parameters for the DecorrelateALKernelTask
    """

    ignoreMaskPlanes = pexConfig.ListField(
        dtype=str,
        doc="""Mask planes to ignore for sigma-clipped statistics""",
        default=("INTRP", "EDGE", "DETECTED", "SAT", "CR", "BAD", "NO_DATA", "DETECTED_NEGATIVE")
    )


class DecorrelateALKernelTask(pipeBase.Task):
    """Decorrelate the effect of convolution by Alard-Lupton matching kernel in image difference

    Notes
    -----

    Pipe-task that removes the neighboring-pixel covariance in an
    image difference that are added when the template image is
    convolved with the Alard-Lupton PSF matching kernel.

    The image differencing pipeline task @link
    ip.diffim.psfMatch.PsfMatchTask PSFMatchTask@endlink and @link
    ip.diffim.psfMatch.PsfMatchConfigAL PSFMatchConfigAL@endlink uses
    the Alard and Lupton (1998) method for matching the PSFs of the
    template and science exposures prior to subtraction. The
    Alard-Lupton method identifies a matching kernel, which is then
    (typically) convolved with the template image to perform PSF
    matching. This convolution has the effect of adding covariance
    between neighboring pixels in the template image, which is then
    added to the image difference by subtraction.

    The pixel covariance may be corrected by whitening the noise of
    the image difference. This task performs such a decorrelation by
    computing a decorrelation kernel (based upon the A&L matching
    kernel and variances in the template and science images) and
    convolving the image difference with it. This process is described
    in detail in [DMTN-021](http://dmtn-021.lsst.io).

    This task has no standalone example, however it is applied as a
    subtask of pipe.tasks.imageDifference.ImageDifferenceTask.
    """
    ConfigClass = DecorrelateALKernelConfig
    _DefaultName = "ip_diffim_decorrelateALKernel"

    def __init__(self, *args, **kwargs):
        """Create the image decorrelation Task

        Parameters
        ----------
        args :
            arguments to be passed to ``lsst.pipe.base.task.Task.__init__``
        kwargs :
            keyword arguments to be passed to ``lsst.pipe.base.task.Task.__init__``
        """
        pipeBase.Task.__init__(self, *args, **kwargs)

        self.statsControl = afwMath.StatisticsControl()
        self.statsControl.setNumSigmaClip(3.)
        self.statsControl.setNumIter(3)
        self.statsControl.setAndMask(afwImage.Mask.getPlaneBitMask(self.config.ignoreMaskPlanes))

    def computeVarianceMean(self, exposure):
        statObj = afwMath.makeStatistics(exposure.getMaskedImage().getVariance(),
                                         exposure.getMaskedImage().getMask(),
                                         afwMath.MEANCLIP, self.statsControl)
        var = statObj.getValue(afwMath.MEANCLIP)
        return var

    @pipeBase.timeMethod
    def run(self, exposure, templateExposure, subtractedExposure, psfMatchingKernel,
            preConvKernel=None, xcen=None, ycen=None, svar=None, tvar=None):
        """Perform decorrelation of an image difference exposure.

        Decorrelates the diffim due to the convolution of the templateExposure with the
        A&L PSF matching kernel. Currently can accept a spatially varying matching kernel but in
        this case it simply uses a static kernel from the center of the exposure. The decorrelation
        is described in [DMTN-021, Equation 1](http://dmtn-021.lsst.io/#equation-1), where
        `exposure` is I_1; templateExposure is I_2; `subtractedExposure` is D(k);
        `psfMatchingKernel` is kappa; and svar and tvar are their respective
        variances (see below).

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The science afwImage.Exposure used for PSF matching
        templateExposure : `lsst.afw.image.Exposure`
            The template exposure used for PSF matching
        subtractedExposure :
            the subtracted exposure produced by
            `ip_diffim.ImagePsfMatchTask.subtractExposures()`
        psfMatchingKernel :
            An (optionally spatially-varying) PSF matching kernel produced
            by `ip_diffim.ImagePsfMatchTask.subtractExposures()`
        preConvKernel :
            if not None, then the `exposure` was pre-convolved with this kernel
        xcen : `float`, optional
            X-pixel coordinate to use for computing constant matching kernel to use
            If `None` (default), then use the center of the image.
        ycen : `float`, optional
            Y-pixel coordinate to use for computing constant matching kernel to use
            If `None` (default), then use the center of the image.
        svar : `float`, optional
            image variance for science image
            If `None` (default) then compute the variance over the entire input science image.
        tvar : `float`, optional
            Image variance for template image
            If `None` (default) then compute the variance over the entire input template image.

        Returns
        -------
        result : `Struct`
            a `lsst.pipe.base.Struct` containing:

            - ``correctedExposure`` : the decorrelated diffim
            - ``correctionKernel`` : the decorrelation correction kernel (which may be ignored)

        Notes
        -----
        The `subtractedExposure` is NOT updated

        The returned `correctedExposure` has an updated PSF as well.

        Here we currently convert a spatially-varying matching kernel into a constant kernel,
        just by computing it at the center of the image (tickets DM-6243, DM-6244).

        We are also using a constant accross-the-image measure of sigma (sqrt(variance)) to compute
        the decorrelation kernel.

        Still TBD (ticket DM-6580): understand whether the convolution is correctly modifying
        the variance plane of the new subtractedExposure.
        """
        spatialKernel = psfMatchingKernel
        kimg = afwImage.ImageD(spatialKernel.getDimensions())
        bbox = subtractedExposure.getBBox()
        if xcen is None:
            xcen = (bbox.getBeginX() + bbox.getEndX()) / 2.
        if ycen is None:
            ycen = (bbox.getBeginY() + bbox.getEndY()) / 2.
        self.log.info("Using matching kernel computed at (%d, %d)", xcen, ycen)
        spatialKernel.computeImage(kimg, True, xcen, ycen)

        if svar is None:
            svar = self.computeVarianceMean(exposure)
        if tvar is None:
            tvar = self.computeVarianceMean(templateExposure)
        self.log.info("Variance (science, template): (%f, %f)", svar, tvar)

        # Should not happen unless entire image has been masked, which could happen
        # if this is a small subimage of the main exposure. In this case, just return a full NaN
        # exposure
        if np.isnan(svar) or np.isnan(tvar):
            # Double check that one of the exposures is all NaNs
            if (np.all(np.isnan(exposure.getMaskedImage().getImage().getArray())) or
                    np.all(np.isnan(templateExposure.getMaskedImage().getImage().getArray()))):
                self.log.warn('Template or science image is entirely NaNs: skipping decorrelation.')
                outExposure = subtractedExposure.clone()
                return pipeBase.Struct(correctedExposure=outExposure, correctionKernel=None)

        tOverSVar = tvar/svar
        if tOverSVar > 1e8:
            self.log.warn("Science image variance is much smaller than template"
                          f", tvar/svar:{tOverSVar:.2e}")

        var = self.computeVarianceMean(subtractedExposure)
        self.log.info("Variance (uncorrected diffim): %f", var)

        if preConvKernel is not None:
            self.log.info('Using a pre-convolution kernel as part of decorrelation correction.')
            kimg2 = afwImage.ImageD(preConvKernel.getDimensions())
            preConvKernel.computeImage(kimg2, False)
            pckArr = kimg2.getArray()

        kArr = kimg.getArray()
        diffExpArr = subtractedExposure.getMaskedImage().getImage().getArray()
        psfArr = subtractedExposure.getPsf().computeKernelImage(geom.Point2D(xcen, ycen)).getArray()

        # Determine the common shape
        if preConvKernel is None:
            self.computeCommonShape(kArr.shape, psfArr.shape, diffExpArr.shape)
            corrft = self.computeCorrection(kArr, svar, tvar)
        else:
            self.computeCommonShape(pckArr.shape, kArr.shape,
                                    psfArr.shape, diffExpArr.shape)
            corrft = self.computeCorrection(kArr, svar, tvar, preConvArr=pckArr)

        diffExpArr = self.computeCorrectedImage(corrft, diffExpArr)
        # The whitening should remove the correlation and scale the diffexp variance
        # to svar + tvar on average
        psfArr = self.computeCorrectedDiffimPsf(corrft, psfArr)

        psfcI = afwImage.ImageD(psfArr.shape[0], psfArr.shape[1])
        psfcI.getArray()[...] = psfArr
        psfcK = afwMath.FixedKernel(psfcI)
        psfNew = measAlg.KernelPsf(psfcK)

        correctedExposure = subtractedExposure.clone()
        correctedExposure.getMaskedImage().getImage().getArray()[...] = diffExpArr
        var = correctedExposure.getMaskedImage().getVariance()
        var.assign(exposure.getMaskedImage().getVariance())
        var += templateExposure.getMaskedImage().getVariance()
        correctedExposure.setPsf(psfNew)

        var = self.computeVarianceMean(correctedExposure)
        self.log.info(f"Variance (corrected diffim): {var:.2e}")

        return pipeBase.Struct(correctedExposure=correctedExposure, )

    def computeCommonShape(self, *shapes):
        """Calculate and sets internally the common shape for FFT operations.

        Parameters
        ----------
        shapes : one or more`tuple` of `int`
            Shapes of the arrays. All must have the same dimensionality.
            At least one shape must be provided.

        Returns
        -------
        None

        Notes
        -----
        For each dimension, gets the smallest even number greater than or equal to
        `N1+N2-1` where `N1` and `N2` are the two largest values.
        In case of only one shape given, rounds up to even each dimension value.
        """
        S = np.array(shapes, dtype=int)
        if len(shapes) > 2:
            S.sort(axis=0)
            S = S[-2:]
        if len(shapes) > 1:
            commonShape = np.sum(S, axis=0) - 1
        else:
            commonShape = S[0]
        commonShape[commonShape % 2 != 0] += 1
        self.freqSpaceShape = tuple(commonShape)
        self.log.info(f"Common frequency space shape {self.freqSpaceShape}")

    @staticmethod
    def padCenterOriginArray(A, newShape: tuple, onwardOp=True):
        """Zero pad an image where the origin is at the center and replace the
        origin to the corner as required by the periodic input of FFT. Implement also
        the inverse operation, crop the padding and re-center data.

        Parameters
        ----------
        A : `numpy.ndarray`
            An array to copy from.
        newShape : `tuple` of `int`
            The dimensions of the resulting array. For padding, the resulting array
            must be larger than A in each dimension. For the inverse operation this
            must be the original, before padding size of the array.
        onwardOp : bool, optional
            Selector of the padding (True) or its inverse (False) operation.

        Returns
        -------
        R : `numpy.ndarray`
            The padded or unpadded array with shape of `newShape` and the same dtype as A.

        Notes
        -----
        Supports n-dimension arrays. For odd dimensions, the splitting is rounded to
        put the center element into the new origin (eg. the center pixel of an odd sized
        kernel will be located at (0,0) ready for FFT).

        """

        # The onward and inverse operations should round odd dimension halves at the opposite
        # sides to get the pixels back to their original positions.
        if onwardOp:
            firstHalves = [x//2 for x in A.shape]
            secondHalves = [x-y for x, y in zip(A.shape, firstHalves)]
        else:
            secondHalves = [x//2 for x in newShape]
            firstHalves = [x-y for x, y in zip(newShape, secondHalves)]

        R = np.zeros_like(A, shape=newShape)
        R[-firstHalves[0]:, -firstHalves[1]:] = A[:firstHalves[0], :firstHalves[1]]
        R[:secondHalves[0], -firstHalves[1]:] = A[-secondHalves[0]:, :firstHalves[1]]
        R[:secondHalves[0], :secondHalves[1]] = A[-secondHalves[0]:, -secondHalves[1]:]
        R[-firstHalves[0]:, :secondHalves[1]] = A[:firstHalves[0], -secondHalves[1]:]
        return R

    def computeCorrection(self, kappa, svar, tvar, preConvArr=None):
        """Compute the Lupton decorrelation post-conv. kernel for decorrelating an
        image difference, based on the PSF-matching kernel.

        Parameters
        ----------
        kappa : `numpy.ndarray`
            A matching kernel 2-d numpy.array derived from Alard & Lupton PSF matching
        svar : `float`, optional
            Average variance of science image used for PSF matching
        tvar : `float`, optional
            Average variance of the template (matched) image used for PSF matching
        preConvKernel If not None, then pre-filtering was applied
            to science exposure, and this is the pre-convolution kernel.

        Returns
        -------
        corrft : `numpy.ndarray` dtype complex but contains real numbers
            TBD

        Notes
        -----
        kappa, and preConvKernel must have shape of self.freqSpaceShape.
        The maximum correction factor converges to sqrt(tvar/svar) towards high frequencies.
        This should be a plausible value.
        """
        kappa = self.padCenterOriginArray(kappa, self.freqSpaceShape)
        kft = np.fft.fft2(kappa)
        kft2 = np.conj(kft) * kft
        if preConvArr is None:
            denom = svar + tvar * kft2
        else:
            preConvArr = self.padCenterOriginArray(preConvArr, self.freqSpaceShape)
            mk = np.fft.fft2(preConvArr)
            mk2 = np.conj(mk) * mk
            denom = svar * mk2 + tvar * kft2
        kft = np.sqrt((svar + tvar) / denom)
        return kft

    def computeCorrectedDiffimPsf(self, corrft, psfArr):
        """Compute the (decorrelated) difference image's new PSF.
        new_psf = psf(k) * sqrt((svar + tvar) / (svar + tvar * kappa_ft(k)**2))

        Parameters
        ----------
        TBD

        Returns
        -------
        TBD
        """
        psfShape = psfArr.shape
        psfArr = self.padCenterOriginArray(psfArr, self.freqSpaceShape)
        psfArr = np.fft.fft2(psfArr)
        psfArr *= corrft
        psfArr = np.fft.ifft2(psfArr)
        psfArr = psfArr.real
        psfArr = self.padCenterOriginArray(psfArr, psfShape, onwardOp=False)
        psfArr = psfArr/psfArr.sum()
        return psfArr

    def computeCorrectedImage(self, corrft, expArr):
        """Convolve an Exposure with a decorrelation convolution kernel.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Input exposure to be convolved.
        kernel : `numpy.array`
            Input 2-d numpy.array to convolve the image with

        Returns
        -------
        out :
            TBD

        Notes
        -----

        """
        expShape = expArr.shape
        filtInf = np.isinf(expArr)
        filtNan = np.isnan(expArr)
        expArr[filtInf] = np.nan
        expArr[filtInf | filtNan] = np.nanmean(expArr)
        expArr = self.padCenterOriginArray(expArr, self.freqSpaceShape)
        expArr = np.fft.fft2(expArr)
        expArr *= corrft
        expArr = np.fft.ifft2(expArr)
        expArr = expArr.real
        expArr = self.padCenterOriginArray(expArr, expShape, onwardOp=False)
        expArr[filtNan] = np.nan
        expArr[filtInf] = np.inf
        return expArr


class DecorrelateALKernelMapper(DecorrelateALKernelTask, ImageMapper):
    """Task to be used as an ImageMapper for performing
    A&L decorrelation on subimages on a grid across a A&L difference image.

    This task subclasses DecorrelateALKernelTask in order to implement
    all of that task's configuration parameters, as well as its `run` method.
    """

    ConfigClass = DecorrelateALKernelConfig
    _DefaultName = 'ip_diffim_decorrelateALKernelMapper'

    def __init__(self, *args, **kwargs):
        DecorrelateALKernelTask.__init__(self, *args, **kwargs)

    def run(self, subExposure, expandedSubExposure, fullBBox,
            template, science, alTaskResult=None, psfMatchingKernel=None,
            preConvKernel=None, **kwargs):
        """Perform decorrelation operation on `subExposure`, using
        `expandedSubExposure` to allow for invalid edge pixels arising from
        convolutions.

        This method performs A&L decorrelation on `subExposure` using
        local measures for image variances and PSF. `subExposure` is a
        sub-exposure of the non-decorrelated A&L diffim. It also
        requires the corresponding sub-exposures of the template
        (`template`) and science (`science`) exposures.

        Parameters
        ----------
        subExposure : `lsst.afw.image.Exposure`
            the sub-exposure of the diffim
        expandedSubExposure : `lsst.afw.image.Exposure`
            the expanded sub-exposure upon which to operate
        fullBBox : `lsst.geom.Box2I`
            the bounding box of the original exposure
        template : `lsst.afw.image.Exposure`
            the corresponding sub-exposure of the template exposure
        science : `lsst.afw.image.Exposure`
            the corresponding sub-exposure of the science exposure
        alTaskResult : `lsst.pipe.base.Struct`
            the result of A&L image differencing on `science` and
            `template`, importantly containing the resulting
            `psfMatchingKernel`. Can be `None`, only if
            `psfMatchingKernel` is not `None`.
        psfMatchingKernel : Alternative parameter for passing the
            A&L `psfMatchingKernel` directly.
        preConvKernel : If not None, then pre-filtering was applied
            to science exposure, and this is the pre-convolution
            kernel.
        kwargs :
            additional keyword arguments propagated from
            `ImageMapReduceTask.run`.

        Returns
        -------
        A `pipeBase.Struct` containing:

            - ``subExposure`` : the result of the `subExposure` processing.
            - ``decorrelationKernel`` : the decorrelation kernel, currently
                not used.

        Notes
        -----
        This `run` method accepts parameters identical to those of
        `ImageMapper.run`, since it is called from the
        `ImageMapperTask`. See that class for more information.
        """
        templateExposure = template  # input template
        scienceExposure = science  # input science image
        if alTaskResult is None and psfMatchingKernel is None:
            raise RuntimeError('Both alTaskResult and psfMatchingKernel cannot be None')
        psfMatchingKernel = alTaskResult.psfMatchingKernel if alTaskResult is not None else psfMatchingKernel

        # subExp and expandedSubExp are subimages of the (un-decorrelated) diffim!
        # So here we compute corresponding subimages of templateExposure and scienceExposure
        subExp2 = scienceExposure.Factory(scienceExposure, expandedSubExposure.getBBox())
        subExp1 = templateExposure.Factory(templateExposure, expandedSubExposure.getBBox())

        # Prevent too much log INFO verbosity from DecorrelateALKernelTask.run
        logLevel = self.log.getLevel()
        self.log.setLevel(lsst.log.WARN)
        res = DecorrelateALKernelTask.run(self, subExp2, subExp1, expandedSubExposure,
                                          psfMatchingKernel, preConvKernel)
        self.log.setLevel(logLevel)  # reset the log level

        diffim = res.correctedExposure.Factory(res.correctedExposure, subExposure.getBBox())
        out = pipeBase.Struct(subExposure=diffim, )
        return out


class DecorrelateALKernelMapReduceConfig(ImageMapReduceConfig):
    """Configuration parameters for the ImageMapReduceTask to direct it to use
       DecorrelateALKernelMapper as its mapper for A&L decorrelation.
    """
    mapper = pexConfig.ConfigurableField(
        doc='A&L decorrelation task to run on each sub-image',
        target=DecorrelateALKernelMapper
    )


class DecorrelateALKernelSpatialConfig(pexConfig.Config):
    """Configuration parameters for the DecorrelateALKernelSpatialTask.
    """
    decorrelateConfig = pexConfig.ConfigField(
        dtype=DecorrelateALKernelConfig,
        doc='DecorrelateALKernel config to use when running on complete exposure (non spatially-varying)',
    )

    decorrelateMapReduceConfig = pexConfig.ConfigField(
        dtype=DecorrelateALKernelMapReduceConfig,
        doc='DecorrelateALKernelMapReduce config to use when running on each sub-image (spatially-varying)',
    )

    ignoreMaskPlanes = pexConfig.ListField(
        dtype=str,
        doc="""Mask planes to ignore for sigma-clipped statistics""",
        default=("INTRP", "EDGE", "DETECTED", "SAT", "CR", "BAD", "NO_DATA", "DETECTED_NEGATIVE")
    )

    def setDefaults(self):
        self.decorrelateMapReduceConfig.gridStepX = self.decorrelateMapReduceConfig.gridStepY = 40
        self.decorrelateMapReduceConfig.cellSizeX = self.decorrelateMapReduceConfig.cellSizeY = 41
        self.decorrelateMapReduceConfig.borderSizeX = self.decorrelateMapReduceConfig.borderSizeY = 8
        self.decorrelateMapReduceConfig.reducer.reduceOperation = 'average'


class DecorrelateALKernelSpatialTask(pipeBase.Task):
    """Decorrelate the effect of convolution by Alard-Lupton matching kernel in image difference

    Notes
    -----

    Pipe-task that removes the neighboring-pixel covariance in an
    image difference that are added when the template image is
    convolved with the Alard-Lupton PSF matching kernel.

    This task is a simple wrapper around @ref DecorrelateALKernelTask,
    which takes a `spatiallyVarying` parameter in its `run` method. If
    it is `False`, then it simply calls the `run` method of @ref
    DecorrelateALKernelTask. If it is True, then it uses the @ref
    ImageMapReduceTask framework to break the exposures into
    subExposures on a grid, and performs the `run` method of @ref
    DecorrelateALKernelTask on each subExposure. This enables it to
    account for spatially-varying PSFs and noise in the exposures when
    performing the decorrelation.

    This task has no standalone example, however it is applied as a
    subtask of pipe.tasks.imageDifference.ImageDifferenceTask.
    There is also an example of its use in `tests/testImageDecorrelation.py`.
    """
    ConfigClass = DecorrelateALKernelSpatialConfig
    _DefaultName = "ip_diffim_decorrelateALKernelSpatial"

    def __init__(self, *args, **kwargs):
        """Create the image decorrelation Task

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

        self.statsControl = afwMath.StatisticsControl()
        self.statsControl.setNumSigmaClip(3.)
        self.statsControl.setNumIter(3)
        self.statsControl.setAndMask(afwImage.Mask.getPlaneBitMask(self.config.ignoreMaskPlanes))

    def computeVarianceMean(self, exposure):
        """Compute the mean of the variance plane of `exposure`.
        """
        statObj = afwMath.makeStatistics(exposure.getMaskedImage().getVariance(),
                                         exposure.getMaskedImage().getMask(),
                                         afwMath.MEANCLIP, self.statsControl)
        var = statObj.getValue(afwMath.MEANCLIP)
        return var

    def run(self, scienceExposure, templateExposure, subtractedExposure, psfMatchingKernel,
            spatiallyVarying=True, preConvKernel=None):
        """Perform decorrelation of an image difference exposure.

        Decorrelates the diffim due to the convolution of the
        templateExposure with the A&L psfMatchingKernel. If
        `spatiallyVarying` is True, it utilizes the spatially varying
        matching kernel via the `imageMapReduce` framework to perform
        spatially-varying decorrelation on a grid of subExposures.

        Parameters
        ----------
        scienceExposure : `lsst.afw.image.Exposure`
           the science Exposure used for PSF matching
        templateExposure : `lsst.afw.image.Exposure`
           the template Exposure used for PSF matching
        subtractedExposure : `lsst.afw.image.Exposure`
           the subtracted Exposure produced by `ip_diffim.ImagePsfMatchTask.subtractExposures()`
        psfMatchingKernel :
           an (optionally spatially-varying) PSF matching kernel produced
           by `ip_diffim.ImagePsfMatchTask.subtractExposures()`
        spatiallyVarying : `bool`
           if True, perform the spatially-varying operation
        preConvKernel : `lsst.meas.algorithms.Psf`
           if not none, the scienceExposure has been pre-filtered with this kernel. (Currently
           this option is experimental.)

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            a structure containing:

            - ``correctedExposure`` : the decorrelated diffim

        """
        self.log.info('Running A&L decorrelation: spatiallyVarying=%r' % spatiallyVarying)

        svar = self.computeVarianceMean(scienceExposure)
        tvar = self.computeVarianceMean(templateExposure)
        if np.isnan(svar) or np.isnan(tvar):  # Should not happen unless entire image has been masked.
            # Double check that one of the exposures is all NaNs
            if (np.all(np.isnan(scienceExposure.getMaskedImage().getImage().getArray())) or
                    np.all(np.isnan(templateExposure.getMaskedImage().getImage().getArray()))):
                self.log.warn('Template or science image is entirely NaNs: skipping decorrelation.')
                if np.isnan(svar):
                    svar = 1e-9
                if np.isnan(tvar):
                    tvar = 1e-9

        var = self.computeVarianceMean(subtractedExposure)

        if spatiallyVarying:
            self.log.info("Variance (science, template): (%f, %f)", svar, tvar)
            self.log.info("Variance (uncorrected diffim): %f", var)
            config = self.config.decorrelateMapReduceConfig
            task = ImageMapReduceTask(config=config)
            results = task.run(subtractedExposure, science=scienceExposure,
                               template=templateExposure, psfMatchingKernel=psfMatchingKernel,
                               preConvKernel=preConvKernel, forceEvenSized=True)
            results.correctedExposure = results.exposure

            # Make sure masks of input image are propagated to diffim
            def gm(exp):
                return exp.getMaskedImage().getMask()
            gm(results.correctedExposure)[:, :] = gm(subtractedExposure)

            var = self.computeVarianceMean(results.correctedExposure)
            self.log.info("Variance (corrected diffim): %f", var)

        else:
            config = self.config.decorrelateConfig
            task = DecorrelateALKernelTask(config=config)
            results = task.run(scienceExposure, templateExposure,
                               subtractedExposure, psfMatchingKernel, preConvKernel=preConvKernel)

        return results
