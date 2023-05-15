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
import lsst.meas.algorithms as measAlg
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod

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
    completeVarPlanePropagation = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Compute the full effect of the decorrelated matching kernel on the variance plane."
            " Otherwise use a model weighed sum of the input variances."
    )


class DecorrelateALKernelTask(pipeBase.Task):
    """Decorrelate the effect of convolution by Alard-Lupton matching kernel in image difference

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
        statObj = afwMath.makeStatistics(exposure.variance,
                                         exposure.mask,
                                         afwMath.MEANCLIP, self.statsControl)
        var = statObj.getValue(afwMath.MEANCLIP)
        return var

    @timeMethod
    def run(self, scienceExposure, templateExposure, subtractedExposure, psfMatchingKernel,
            preConvKernel=None, xcen=None, ycen=None, svar=None, tvar=None,
            templateMatched=True, preConvMode=False, **kwargs):
        """Perform decorrelation of an image difference or of a score difference exposure.

        Corrects the difference or score image due to the convolution of the
        templateExposure with the A&L PSF matching kernel.
        See [DMTN-021, Equation 1](http://dmtn-021.lsst.io/#equation-1) and
        [DMTN-179](http://dmtn-179.lsst.io/) for details.

        Parameters
        ----------
        scienceExposure : `lsst.afw.image.Exposure`
            The original science exposure (before pre-convolution, if ``preConvMode==True``).
        templateExposure : `lsst.afw.image.Exposure`
            The original template exposure warped, but not psf-matched, to the science exposure.
        subtractedExposure : `lsst.afw.image.Exposure`
            the subtracted exposure produced by
            `ip_diffim.ImagePsfMatchTask.subtractExposures()`. The `subtractedExposure` must
            inherit its PSF from `exposure`, see notes below.
        psfMatchingKernel : `lsst.afw.detection.Psf`
            An (optionally spatially-varying) PSF matching kernel produced
            by `ip_diffim.ImagePsfMatchTask.subtractExposures()`.
        preConvKernel : `lsst.afw.math.Kernel`, optional
            If not `None`, then the `scienceExposure` was pre-convolved with (the reflection of)
            this kernel. Must be normalized to sum to 1.
            Allowed only if ``templateMatched==True`` and ``preConvMode==True``.
            Defaults to the PSF of the science exposure at the image center.
        xcen : `float`, optional
            X-pixel coordinate to use for computing constant matching kernel to use
            If `None` (default), then use the center of the image.
        ycen : `float`, optional
            Y-pixel coordinate to use for computing constant matching kernel to use
            If `None` (default), then use the center of the image.
        svar : `float`, optional
            Image variance for science image
            If `None` (default) then compute the variance over the entire input science image.
        tvar : `float`, optional
            Image variance for template image
            If `None` (default) then compute the variance over the entire input template image.
        templateMatched : `bool`, optional
            If True, the template exposure was matched (convolved) to the science exposure.
            See also notes below.
        preConvMode : `bool`, optional
            If True, ``subtractedExposure`` is assumed to be a likelihood difference image
            and will be noise corrected as a likelihood image.
        **kwargs
            Additional keyword arguments propagated from DecorrelateALKernelSpatialTask.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            - ``correctedExposure`` : the decorrelated diffim

        Notes
        -----
        If ``preConvMode==True``, ``subtractedExposure`` is assumed to be a
        score image and the noise correction for likelihood images
        is applied. The resulting image is an optimal detection likelihood image
        when the templateExposure has noise. (See DMTN-179) If ``preConvKernel`` is
        not specified, the PSF of ``scienceExposure`` is assumed as pre-convolution kernel.

        The ``subtractedExposure`` is NOT updated. The returned ``correctedExposure``
        has an updated but spatially fixed PSF. It is calculated as the center of
        image PSF corrected by the center of image matching kernel.

        If ``templateMatched==True``, the templateExposure was matched (convolved)
        to the ``scienceExposure`` by ``psfMatchingKernel`` during image differencing.
        Otherwise the ``scienceExposure`` was matched (convolved) by ``psfMatchingKernel``.
        In either case, note that the original template and science images are required,
        not the psf-matched version.

        This task discards the variance plane of ``subtractedExposure`` and re-computes
        it from the variance planes of ``scienceExposure`` and ``templateExposure``.
        The image plane of ``subtractedExposure`` must be at the photometric level
        set by the AL PSF matching in `ImagePsfMatchTask.subtractExposures`.
        The assumptions about the photometric level are controlled by the
        `templateMatched` option in this task.

        Here we currently convert a spatially-varying matching kernel into a constant kernel,
        just by computing it at the center of the image (tickets DM-6243, DM-6244).

        We are also using a constant accross-the-image measure of sigma (sqrt(variance)) to compute
        the decorrelation kernel.

        TODO DM-23857 As part of the spatially varying correction implementation
        consider whether returning a Struct is still necessary.
        """
        if preConvKernel is not None and not (templateMatched and preConvMode):
            raise ValueError("Pre-convolution kernel is allowed only if "
                             "preConvMode==True and templateMatched==True.")

        spatialKernel = psfMatchingKernel
        kimg = afwImage.ImageD(spatialKernel.getDimensions())
        bbox = subtractedExposure.getBBox()
        if xcen is None:
            xcen = (bbox.getBeginX() + bbox.getEndX()) / 2.
        if ycen is None:
            ycen = (bbox.getBeginY() + bbox.getEndY()) / 2.
        self.log.info("Using matching kernel computed at (%d, %d)", xcen, ycen)
        spatialKernel.computeImage(kimg, False, xcen, ycen)

        preConvImg = None
        if preConvMode:
            if preConvKernel is None:
                pos = scienceExposure.getPsf().getAveragePosition()
                preConvKernel = scienceExposure.getPsf().getLocalKernel(pos)
            preConvImg = afwImage.ImageD(preConvKernel.getDimensions())
            preConvKernel.computeImage(preConvImg, True)

        if svar is None:
            svar = self.computeVarianceMean(scienceExposure)
        if tvar is None:
            tvar = self.computeVarianceMean(templateExposure)
        self.log.info("Original variance plane means. Science:%.5e, warped template:%.5e)",
                      svar, tvar)

        # Should not happen unless entire image has been masked, which could happen
        # if this is a small subimage of the main exposure. In this case, just return a full NaN
        # exposure
        if np.isnan(svar) or np.isnan(tvar):
            # Double check that one of the exposures is all NaNs
            if (np.all(np.isnan(scienceExposure.image.array))
                    or np.all(np.isnan(templateExposure.image.array))):
                self.log.warning('Template or science image is entirely NaNs: skipping decorrelation.')
                outExposure = subtractedExposure.clone()
                return pipeBase.Struct(correctedExposure=outExposure, )

        if templateMatched:
            # Regular subtraction, we convolved the template
            self.log.info("Decorrelation after template image convolution")
            varianceMean = svar
            targetVarianceMean = tvar
            # variance plane of the image that is not convolved
            variance = scienceExposure.variance.array
            # Variance plane of the convolved image, before convolution.
            targetVariance = templateExposure.variance.array
        else:
            # We convolved the science image
            self.log.info("Decorrelation after science image convolution")
            varianceMean = tvar
            targetVarianceMean = svar
            # variance plane of the image that is not convolved
            variance = templateExposure.variance.array
            # Variance plane of the convolved image, before convolution.
            targetVariance = scienceExposure.variance.array

        # The maximal correction value converges to sqrt(targetVarianceMean/varianceMean).
        # Correction divergence warning if the correction exceeds 4 orders of magnitude.
        mOverExpVar = targetVarianceMean/varianceMean
        if mOverExpVar > 1e8:
            self.log.warning("Diverging correction: matched image variance is "
                             " much larger than the unconvolved one's"
                             ", targetVarianceMean/varianceMean:%.2e", mOverExpVar)

        oldVarMean = self.computeVarianceMean(subtractedExposure)
        self.log.info("Variance plane mean of uncorrected diffim: %f", oldVarMean)

        kArr = kimg.array
        diffimShape = subtractedExposure.image.array.shape
        psfImg = subtractedExposure.getPsf().computeKernelImage(geom.Point2D(xcen, ycen))
        psfShape = psfImg.array.shape

        if preConvMode:
            self.log.info("Decorrelation of likelihood image")
            self.computeCommonShape(preConvImg.array.shape, kArr.shape,
                                    psfShape, diffimShape)
            corr = self.computeScoreCorrection(kArr, varianceMean, targetVarianceMean, preConvImg.array)
        else:
            self.log.info("Decorrelation of difference image")
            self.computeCommonShape(kArr.shape, psfShape, diffimShape)
            corr = self.computeDiffimCorrection(kArr, varianceMean, targetVarianceMean)

        correctedImage = self.computeCorrectedImage(corr.corrft, subtractedExposure.image.array)
        correctedPsf = self.computeCorrectedDiffimPsf(corr.corrft, psfImg.array)

        # The subtracted exposure variance plane is already correlated, we cannot propagate
        # it through another convolution; instead we need to use the uncorrelated originals
        # The whitening should scale it to varianceMean + targetVarianceMean on average
        if self.config.completeVarPlanePropagation:
            self.log.debug("Using full variance plane calculation in decorrelation")
            correctedVariance = self.calculateVariancePlane(
                variance, targetVariance,
                varianceMean, targetVarianceMean, corr.cnft, corr.crft)
        else:
            self.log.debug("Using estimated variance plane calculation in decorrelation")
            correctedVariance = self.estimateVariancePlane(
                variance, targetVariance,
                corr.cnft, corr.crft)

        # Determine the common shape
        kSum = np.sum(kArr)
        kSumSq = kSum*kSum
        self.log.debug("Matching kernel sum: %.3e", kSum)
        if not templateMatched:
            # ImagePsfMatch.subtractExposures re-scales the difference in
            # the science image convolution mode
            correctedVariance /= kSumSq
        subtractedExposure.image.array[...] = correctedImage  # Allow for numpy type casting
        subtractedExposure.variance.array[...] = correctedVariance
        subtractedExposure.setPsf(correctedPsf)

        newVarMean = self.computeVarianceMean(subtractedExposure)
        self.log.info("Variance plane mean of corrected diffim: %.5e", newVarMean)

        # TODO DM-23857 As part of the spatially varying correction implementation
        # consider whether returning a Struct is still necessary.
        return pipeBase.Struct(correctedExposure=subtractedExposure, )

    def computeCommonShape(self, *shapes):
        """Calculate the common shape for FFT operations. Set `self.freqSpaceShape`
        internally.

        Parameters
        ----------
        shapes : one or more `tuple` of `int`
            Shapes of the arrays. All must have the same dimensionality.
            At least one shape must be provided.

        Returns
        -------
        None.

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
        self.log.info("Common frequency space shape %s", self.freqSpaceShape)

    @staticmethod
    def padCenterOriginArray(A, newShape: tuple, useInverse=False):
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
        useInverse : bool, optional
            Selector of forward, add padding, operation (False)
            or its inverse, crop padding, operation (True).

        Returns
        -------
        R : `numpy.ndarray`
            The padded or unpadded array with shape of `newShape` and the same dtype as A.

        Notes
        -----
        For odd dimensions, the splitting is rounded to
        put the center pixel into the new corner origin (0,0). This is to be consistent
        e.g. for a dirac delta kernel that is originally located at the center pixel.
        """

        # The forward and inverse operations should round odd dimension halves at the opposite
        # sides to get the pixels back to their original positions.
        if not useInverse:
            # Forward operation: First and second halves with respect to the axes of A.
            firstHalves = [x//2 for x in A.shape]
            secondHalves = [x-y for x, y in zip(A.shape, firstHalves)]
        else:
            # Inverse operation: Opposite rounding
            secondHalves = [x//2 for x in newShape]
            firstHalves = [x-y for x, y in zip(newShape, secondHalves)]

        R = np.zeros_like(A, shape=newShape)
        R[-firstHalves[0]:, -firstHalves[1]:] = A[:firstHalves[0], :firstHalves[1]]
        R[:secondHalves[0], -firstHalves[1]:] = A[-secondHalves[0]:, :firstHalves[1]]
        R[:secondHalves[0], :secondHalves[1]] = A[-secondHalves[0]:, -secondHalves[1]:]
        R[-firstHalves[0]:, :secondHalves[1]] = A[:firstHalves[0], -secondHalves[1]:]
        return R

    def computeDiffimCorrection(self, kappa, svar, tvar):
        """Compute the Lupton decorrelation post-convolution kernel for decorrelating an
        image difference, based on the PSF-matching kernel.

        Parameters
        ----------
        kappa : `numpy.ndarray` of `float`
            A matching kernel 2-d numpy.array derived from Alard & Lupton PSF matching.
        svar : `float` > 0.
            Average variance of science image used for PSF matching.
        tvar : `float` > 0.
            Average variance of the template (matched) image used for PSF matching.

        Returns
        -------
        corrft : `numpy.ndarray` of `float`
            The frequency space representation of the correction. The array is real (dtype float).
            Shape is `self.freqSpaceShape`.

        cnft, crft : `numpy.ndarray` of `complex`
            The overall convolution (pre-conv, PSF matching, noise correction) kernel
            for the science and template images, respectively for the variance plane
            calculations. These are intermediate results in frequency space.

        Notes
        -----
        The maximum correction factor converges to `sqrt(tvar/svar)` towards high frequencies.
        This should be a plausible value.
        """
        kSum = np.sum(kappa)  # We scale the decorrelation to preserve fluxes
        kappa = self.padCenterOriginArray(kappa, self.freqSpaceShape)
        kft = np.fft.fft2(kappa)
        kftAbsSq = np.real(np.conj(kft) * kft)

        denom = svar + tvar * kftAbsSq
        corrft = np.sqrt((svar + tvar * kSum*kSum) / denom)
        cnft = corrft
        crft = kft*corrft
        return pipeBase.Struct(corrft=corrft, cnft=cnft, crft=crft)

    def computeScoreCorrection(self, kappa, svar, tvar, preConvArr):
        """Compute the correction kernel for a score image.

        Parameters
        ----------
        kappa : `numpy.ndarray`
            A matching kernel 2-d numpy.array derived from Alard & Lupton PSF matching.
        svar : `float`
            Average variance of science image used for PSF matching (before pre-convolution).
        tvar : `float`
            Average variance of the template (matched) image used for PSF matching.
        preConvArr : `numpy.ndarray`
            The pre-convolution kernel of the science image. It should be the PSF
            of the science image or an approximation of it. It must be normed to sum 1.

        Returns
        -------
        corrft : `numpy.ndarray` of `float`
            The frequency space representation of the correction. The array is real (dtype float).
            Shape is `self.freqSpaceShape`.
        cnft, crft : `numpy.ndarray` of `complex`
            The overall convolution (pre-conv, PSF matching, noise correction) kernel
            for the science and template images, respectively for the variance plane
            calculations. These are intermediate results in frequency space.

        Notes
        -----
        To be precise, the science image should be _correlated_ by ``preConvArray`` but this
        does not matter for this calculation.

        ``cnft``, ``crft`` contain the scaling factor as well.

        """
        kSum = np.sum(kappa)
        kappa = self.padCenterOriginArray(kappa, self.freqSpaceShape)
        kft = np.fft.fft2(kappa)
        preConvArr = self.padCenterOriginArray(preConvArr, self.freqSpaceShape)
        preFt = np.fft.fft2(preConvArr)
        preFtAbsSq = np.real(np.conj(preFt) * preFt)
        kftAbsSq = np.real(np.conj(kft) * kft)
        # Avoid zero division, though we don't normally approach `tiny`.
        # We have numerical noise instead.
        tiny = np.finfo(preFtAbsSq.dtype).tiny * 1000.
        flt = preFtAbsSq < tiny
        # If we pre-convolve to avoid deconvolution in AL, then kftAbsSq / preFtAbsSq
        # theoretically expected to diverge to +inf. But we don't care about the convergence
        # properties here, S goes to 0 at these frequencies anyway.
        preFtAbsSq[flt] = tiny
        denom = svar + tvar * kftAbsSq / preFtAbsSq
        corrft = (svar + tvar * kSum*kSum) / denom
        cnft = np.conj(preFt)*corrft
        crft = kft*corrft
        return pipeBase.Struct(corrft=corrft, cnft=cnft, crft=crft)

    @staticmethod
    def estimateVariancePlane(vplane1, vplane2, c1ft, c2ft):
        """Estimate the variance planes.

        The estimation assumes that around each pixel the surrounding
        pixels' sigmas within the convolution kernel are the same.

        Parameters
        ----------
        vplane1, vplane2 : `numpy.ndarray` of `float`
            Variance planes of the original (before pre-convolution or matching)
            exposures.
        c1ft, c2ft : `numpy.ndarray` of `complex`
            The overall convolution that includes the matching and the
            afterburner in frequency space. The result of either
            ``computeScoreCorrection`` or ``computeDiffimCorrection``.

        Returns
        -------
        vplaneD : `numpy.ndarray` of `float`
          The estimated variance plane of the difference/score image
          as a weighted sum of the input variances.

        Notes
        ------
        See DMTN-179 Section 5 about the variance plane calculations.
        """
        w1 = np.sum(np.real(np.conj(c1ft)*c1ft)) / c1ft.size
        w2 = np.sum(np.real(np.conj(c2ft)*c2ft)) / c2ft.size
        # w1, w2: the frequency space sum of abs(c1)^2 is the same as in image
        # space.
        return vplane1*w1 + vplane2*w2

    def calculateVariancePlane(self, vplane1, vplane2, varMean1, varMean2, c1ft, c2ft):
        """Full propagation of the variance planes of the original exposures.

        The original variance planes of independent pixels are convolved with the
        image space square of the overall kernels.

        Parameters
        ----------
        vplane1, vplane2 : `numpy.ndarray` of `float`
            Variance planes of the original (before pre-convolution or matching)
            exposures.
        varMean1, varMean2 : `float`
            Replacement average values for non-finite ``vplane1`` and ``vplane2`` values respectively.

        c1ft, c2ft : `numpy.ndarray` of `complex`
            The overall convolution that includes the matching and the
            afterburner in frequency space. The result of either
            ``computeScoreCorrection`` or ``computeDiffimCorrection``.

        Returns
        -------
        vplaneD : `numpy.ndarray` of `float`
          The variance plane of the difference/score images.

        Notes
        ------
        See DMTN-179 Section 5 about the variance plane calculations.

        Infs and NaNs are allowed and kept in the returned array.
        """
        D = np.real(np.fft.ifft2(c1ft))
        c1SqFt = np.fft.fft2(D*D)

        v1shape = vplane1.shape
        filtInf = np.isinf(vplane1)
        filtNan = np.isnan(vplane1)
        # This copy could be eliminated if inf/nan handling were go into padCenterOriginArray
        vplane1 = np.copy(vplane1)
        vplane1[filtInf | filtNan] = varMean1
        D = self.padCenterOriginArray(vplane1, self.freqSpaceShape)
        v1 = np.real(np.fft.ifft2(np.fft.fft2(D) * c1SqFt))
        v1 = self.padCenterOriginArray(v1, v1shape, useInverse=True)
        v1[filtNan] = np.nan
        v1[filtInf] = np.inf

        D = np.real(np.fft.ifft2(c2ft))
        c2SqFt = np.fft.fft2(D*D)

        v2shape = vplane2.shape
        filtInf = np.isinf(vplane2)
        filtNan = np.isnan(vplane2)
        vplane2 = np.copy(vplane2)
        vplane2[filtInf | filtNan] = varMean2
        D = self.padCenterOriginArray(vplane2, self.freqSpaceShape)
        v2 = np.real(np.fft.ifft2(np.fft.fft2(D) * c2SqFt))
        v2 = self.padCenterOriginArray(v2, v2shape, useInverse=True)
        v2[filtNan] = np.nan
        v2[filtInf] = np.inf

        return v1 + v2

    def computeCorrectedDiffimPsf(self, corrft, psfOld):
        """Compute the (decorrelated) difference image's new PSF.

        Parameters
        ----------
        corrft : `numpy.ndarray`
            The frequency space representation of the correction calculated by
            `computeCorrection`. Shape must be `self.freqSpaceShape`.
        psfOld : `numpy.ndarray`
            The psf of the difference image to be corrected.

        Returns
        -------
        correctedPsf : `lsst.meas.algorithms.KernelPsf`
            The corrected psf, same shape as `psfOld`, sum normed to 1.

        Notes
        -----
        There is no algorithmic guarantee that the corrected psf can
        meaningfully fit to the same size as the original one.
        """
        psfShape = psfOld.shape
        psfNew = self.padCenterOriginArray(psfOld, self.freqSpaceShape)
        psfNew = np.fft.fft2(psfNew)
        psfNew *= corrft
        psfNew = np.fft.ifft2(psfNew)
        psfNew = psfNew.real
        psfNew = self.padCenterOriginArray(psfNew, psfShape, useInverse=True)
        psfNew = psfNew/psfNew.sum()

        psfcI = afwImage.ImageD(geom.Extent2I(psfShape[1], psfShape[0]))
        psfcI.array = psfNew
        psfcK = afwMath.FixedKernel(psfcI)
        correctedPsf = measAlg.KernelPsf(psfcK)
        return correctedPsf

    def computeCorrectedImage(self, corrft, imgOld):
        """Compute the decorrelated difference image.

        Parameters
        ----------
        corrft : `numpy.ndarray`
            The frequency space representation of the correction calculated by
            `computeCorrection`. Shape must be `self.freqSpaceShape`.
        imgOld : `numpy.ndarray`
            The difference image to be corrected.

        Returns
        -------
        imgNew : `numpy.ndarray`
            The corrected image, same size as the input.
        """
        expShape = imgOld.shape
        imgNew = np.copy(imgOld)
        filtInf = np.isinf(imgNew)
        filtNan = np.isnan(imgNew)
        imgNew[filtInf] = np.nan
        imgNew[filtInf | filtNan] = np.nanmean(imgNew)
        imgNew = self.padCenterOriginArray(imgNew, self.freqSpaceShape)
        imgNew = np.fft.fft2(imgNew)
        imgNew *= corrft
        imgNew = np.fft.ifft2(imgNew)
        imgNew = imgNew.real
        imgNew = self.padCenterOriginArray(imgNew, expShape, useInverse=True)
        imgNew[filtNan] = np.nan
        imgNew[filtInf] = np.inf
        return imgNew


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
        logLevel = self.log.level
        self.log.setLevel(self.log.WARNING)
        res = DecorrelateALKernelTask.run(self, subExp2, subExp1, expandedSubExposure,
                                          psfMatchingKernel, preConvKernel, **kwargs)
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
        statObj = afwMath.makeStatistics(exposure.variance,
                                         exposure.mask,
                                         afwMath.MEANCLIP, self.statsControl)
        var = statObj.getValue(afwMath.MEANCLIP)
        return var

    def run(self, scienceExposure, templateExposure, subtractedExposure, psfMatchingKernel,
            spatiallyVarying=True, preConvKernel=None, templateMatched=True, preConvMode=False):
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
        psfMatchingKernel : an (optionally spatially-varying) PSF matching kernel produced
           by `ip_diffim.ImagePsfMatchTask.subtractExposures()`
        spatiallyVarying : `bool`
           if True, perform the spatially-varying operation
        preConvKernel : `lsst.meas.algorithms.Psf`
           if not none, the scienceExposure has been pre-filtered with this kernel. (Currently
           this option is experimental.)
        templateMatched : `bool`, optional
           If True, the template exposure was matched (convolved) to the science exposure.
        preConvMode : `bool`, optional
            If True, ``subtractedExposure`` is assumed to be a likelihood difference image
            and will be noise corrected as a likelihood image.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            a structure containing:
            - ``correctedExposure`` : the decorrelated diffim
        """
        self.log.info('Running A&L decorrelation: spatiallyVarying=%r', spatiallyVarying)

        svar = self.computeVarianceMean(scienceExposure)
        tvar = self.computeVarianceMean(templateExposure)
        if np.isnan(svar) or np.isnan(tvar):  # Should not happen unless entire image has been masked.
            # Double check that one of the exposures is all NaNs
            if (np.all(np.isnan(scienceExposure.image.array))
                    or np.all(np.isnan(templateExposure.image.array))):
                self.log.warning('Template or science image is entirely NaNs: skipping decorrelation.')
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
                               preConvKernel=preConvKernel, forceEvenSized=True,
                               templateMatched=templateMatched, preConvMode=preConvMode)
            results.correctedExposure = results.exposure

            # Make sure masks of input image are propagated to diffim
            def gm(exp):
                return exp.mask
            gm(results.correctedExposure)[:, :] = gm(subtractedExposure)

            var = self.computeVarianceMean(results.correctedExposure)
            self.log.info("Variance (corrected diffim): %f", var)

        else:
            config = self.config.decorrelateConfig
            task = DecorrelateALKernelTask(config=config)
            results = task.run(scienceExposure, templateExposure,
                               subtractedExposure, psfMatchingKernel, preConvKernel=preConvKernel,
                               templateMatched=templateMatched, preConvMode=preConvMode)

        return results
