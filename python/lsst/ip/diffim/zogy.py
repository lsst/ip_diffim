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

import lsst.geom as geom
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.meas.algorithms as measAlg
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase


from .imagePsfMatch import (ImagePsfMatchTask, ImagePsfMatchConfig,
                            subtractAlgorithmRegistry)

__all__ = ["ZogyTask", "ZogyConfig",
           "ZogyImagePsfMatchConfig", "ZogyImagePsfMatchTask"]


"""Tasks for performing the "Proper image subtraction" algorithm of
Zackay, et al. (2016), hereafter simply referred to as 'ZOGY (2016)'.

`ZogyTask` contains methods to perform the basic estimation of the
ZOGY diffim `D`, its updated PSF, and the variance-normalized
likelihood image `S_corr`. We have implemented ZOGY using the
proscribed methodology, computing all convolutions in Fourier space,
and also variants in which the convolutions are performed in real
(image) space. The former is faster and results in fewer artifacts
when the PSFs are noisy (i.e., measured, for example, via
`PsfEx`). The latter is presumed to be preferred as it can account for
masks correctly with fewer "ringing" artifacts from edge effects or
saturated stars, but noisy PSFs result in their own smaller
artifacts. Removal of these artifacts is a subject of continuing
research. Currently, we "pad" the PSFs when performing the
subtractions in real space, which reduces, but does not entirely
eliminate these artifacts.

All methods in `ZogyTask` assume template and science images are
already accurately photometrically and astrometrically registered.

`ZogyMapper` is a wrapper which runs `ZogyTask` in the
`ImageMapReduce` framework, computing of ZOGY diffim's on small,
overlapping sub-images, thereby enabling complete ZOGY diffim's which
account for spatially-varying noise and PSFs across the two input
exposures. An example of the use of this task is in the `testZogy.py`
unit test.
"""


class ZogyConfig(pexConfig.Config):
    """Configuration parameters for the ZogyTask
    """

    templateFluxScaling = pexConfig.Field(
        dtype=float,
        default=1.,
        doc="Template flux scaling factor (Fr in ZOGY paper)"
    )

    scienceFluxScaling = pexConfig.Field(
        dtype=float,
        default=1.,
        doc="Science flux scaling factor (Fn in ZOGY paper)"
    )

    scaleByCalibration = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Compute the flux normalization scaling based on the image calibration."
        "This overrides 'templateFluxScaling' and 'scienceFluxScaling'."
    )

    correctBackground = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Subtract exposure background mean to have zero expectation value."
    )

    ignoreMaskPlanes = pexConfig.ListField(
        dtype=str,
        default=("INTRP", "EDGE", "DETECTED", "SAT", "CR", "BAD", "NO_DATA", "DETECTED_NEGATIVE"),
        doc="Mask planes to ignore for statistics"
    )
    maxPsfCentroidDist = pexConfig.Field(
        dtype=float,
        default=0.2,
        doc="Maximum centroid difference allowed between the two exposure PSFs (pixels)."
    )


class ZogyTask(pipeBase.Task):
    """Task to perform ZOGY proper image subtraction. See module-level documentation for
    additional details.

    """
    ConfigClass = ZogyConfig
    _DefaultName = "ip_diffim_Zogy"

    def _computeVarianceMean(self, exposure):
        """Compute the sigma-clipped mean of the variance image of `exposure`.
        """
        statObj = afwMath.makeStatistics(exposure.getMaskedImage().getVariance(),
                                         exposure.getMaskedImage().getMask(),
                                         afwMath.MEANCLIP, self.statsControl)
        var = statObj.getValue(afwMath.MEANCLIP)
        return var

    @staticmethod
    def padCenterOriginArray(A, newShape, useInverse=False, dtype=None):
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
        dtype: `numpy.dtype`, optional
            Dtype of output array. Values must be implicitly castable to this type.
            Use to get expected result type, e.g. single float (nympy.float32).
            If not specified, dtype is inherited from ``A``.

        Returns
        -------
        R : `numpy.ndarray`
            The padded or unpadded array with shape of `newShape` and dtype of ``dtype``.

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
        if dtype is None:
            dtype = A.dtype

        R = np.zeros(newShape, dtype=dtype)
        R[-firstHalves[0]:, -firstHalves[1]:] = A[:firstHalves[0], :firstHalves[1]]
        R[:secondHalves[0], -firstHalves[1]:] = A[-secondHalves[0]:, :firstHalves[1]]
        R[:secondHalves[0], :secondHalves[1]] = A[-secondHalves[0]:, -secondHalves[1]:]
        R[-firstHalves[0]:, :secondHalves[1]] = A[:firstHalves[0], -secondHalves[1]:]
        return R

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

    def padAndFftImage(self, imgArr):
        """Prepare and forward FFT an image array.

        Parameters
        ----------
        imgArr : `numpy.ndarray` of `float`
            Original array. In-place modified as `numpy.nan` and `numpy.inf` are replaced by
            array mean.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            - ``imFft`` : `numpy.ndarray` of `numpy.complex`.
                FFT of image.
            - ``filtInf``, ``filtNaN`` : `numpy.ndarray` of `bool`

        Notes
        -----
        Save location of non-finite values for restoration, and replace them
        with image mean values. Re-center and zero pad array by `padCenterOriginArray`.
        """
        filtInf = np.isinf(imgArr)
        filtNaN = np.isnan(imgArr)
        imgArr[filtInf] = np.nan
        imgArr[filtInf | filtNaN] = np.nanmean(imgArr)
        self.log.debug(f"Replacing {np.sum(filtInf)} Inf and {np.sum(filtNaN)} NaN values.")
        imgArr = self.padCenterOriginArray(imgArr, self.freqSpaceShape)
        imgArr = np.fft.fft2(imgArr)
        return pipeBase.Struct(imFft=imgArr, filtInf=filtInf, filtNaN=filtNaN)

    def inverseFftAndCropImage(self, imgArr, origSize, filtInf=None, filtNaN=None, dtype=None):
        """Inverse FFT and crop padding from image array.

        Parameters
        ----------
        imgArr : `numpy.ndarray` of `numpy.complex`
            Fourier space array representing a real image.

        origSize : `tuple` of `int`
            Original unpadded shape tuple of the image to be cropped to.

        filtInf, filtNan : `numpy.ndarray` of indices, optional
            If specified, they are used as index arrays for ``result`` to set values to
            `numpy.inf` and `numpy.nan` respectively at these positions.

        dtype : `numpy.dtype`, optional
            Dtype of result array to cast return values to implicitly. This is to
            spare one array copy operation at reducing double precision to single.
            If `None` result inherits dtype of `imgArr`.

        Returns
        -------
        result : `numpy.ndarray` of `dtype`
        """
        imgNew = np.fft.ifft2(imgArr)
        imgNew = imgNew.real
        imgNew = self.padCenterOriginArray(imgNew, origSize, useInverse=True, dtype=dtype)
        if filtInf is not None:
            imgNew[filtInf] = np.inf
        if filtNaN is not None:
            imgNew[filtNaN] = np.nan
        return imgNew

    def computePsfAtCenter(self, exposure, assumeGaussianPsf=False):
        """Computes the PSF image at the bbox center point. This may be a fractional pixel position.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure with psf.

        Returns
        -------
        psfImg : `lsst.afw.image.Image`
            Calculated psf image.
        """

        bbox = exposure.getBBox()
        xcen = (bbox.getBeginX() + bbox.getEndX()) / 2.
        ycen = (bbox.getBeginY() + bbox.getEndY()) / 2.
        psf = exposure.getPsf()
        psfImg = psf.computeKernelImage(geom.Point2D(xcen, ycen))  # Centered and normed
        if assumeGaussianPsf:
            psf_sig = psf.computeShape().getDeterminantRadius()
            kWidth = (int(psf_sig * 30 + 0.5)//2)*2 + 1  # make sure it is odd
            psf_bbox = psfImg.getBBox()
            self.log.debug(f"PSF bbox {psf_bbox}.")
            self.log.debug(f"Using Gaussian PSF width {psf_sig:.3f} and {kWidth} size.")
            gaussFunc = afwMath.GaussianFunction1D(psf_sig)
            gaussKernel = afwMath.SeparableKernel(kWidth, kWidth, gaussFunc, gaussFunc)
            psfImg = afwImage.Image(geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(kWidth, kWidth)),
                                    dtype=psfImg.dtype)
            gaussKernel.computeImage(psfImg, True)
            self.log.debug(f"psfImg min {np.amin(psfImg.array):.3e}")
        return psfImg

    @staticmethod
    def subtractImageMean(image, mask, statsControl):
        """In-place subtraction of sigma-clipped mean of the image.

        Parameters
        ----------
        image : `lsst.afw.image.Image`
            Image to manipulate. Its sigma clipped mean is in-place subtracted.

        mask : `lsst.afw.image.Mask`
            Mask to use for ignoring pixels.

        statsControl : `lsst.afw.math.StatisticsControl`
            Config of sigma clipped mean statistics calculation.

        Returns
        -------
        None
        """
        statObj = afwMath.makeStatistics(image, mask,
                                         afwMath.MEANCLIP, statsControl)
        mean = statObj.getValue(afwMath.MEANCLIP)
        if not np.isnan(mean):
            image -= mean

    def prepareFullExposure(self, exposure1, exposure2, correctBackground=False):
        """Performs calculations that apply to the full exposures once only in the psf matching.

        Parameters
        ----------

        correctBackground : `bool`, optional
            If True, subtracts sigma-clipped mean of exposures. The algorithm
            assumes zero expectation value at background pixels.

        Returns
        -------
        None

        Notes
        -----
        Set a number of instance fields with pre-calculated values.

        Raises
        ------
        ValueError : If photometric calibrations are not available.
        """

        self.statsControl = afwMath.StatisticsControl()
        self.statsControl.setNumSigmaClip(3.)
        self.statsControl.setNumIter(3)
        self.statsControl.setAndMask(afwImage.Mask.getPlaneBitMask(
            self.config.ignoreMaskPlanes))

        exposure1 = exposure1.clone()
        exposure2 = exposure2.clone()
        # If 'scaleByCalibration' is True then these norms are overwritten
        if self.config.scaleByCalibration:
            calibObj1 = exposure1.getPhotoCalib()
            calibObj2 = exposure2.getPhotoCalib()
            if calibObj1 is None or calibObj2 is None:
                raise ValueError("Photometric calibrations are not available for both exposures.")
            mImg1 = calibObj1.calibrateImage(exposure1.maskedImage)
            mImg2 = calibObj2.calibrateImage(exposure2.maskedImage)
            self.F1 = 1.
            self.F2 = 1.
        else:
            self.F1 = self.config.templateFluxScaling  # default is 1
            self.F2 = self.config.scienceFluxScaling  # default is 1
            mImg1 = exposure1.maskedImage
            mImg2 = exposure2.maskedImage

        # mImgs can be in-place modified
        if correctBackground:
            self.subtractImageMean(mImg1.image, mImg1.mask, self.statsControl)
            self.subtractImageMean(mImg2.image, mImg2.mask, self.statsControl)

        psfBBox1 = exposure1.getPsf().computeBBox()
        psfBBox2 = exposure2.getPsf().computeBBox()
        self.psfShape1 = (psfBBox1.getHeight(), psfBBox1.getWidth())
        self.psfShape2 = (psfBBox2.getHeight(), psfBBox2.getWidth())
        self.imgShape = (mImg1.getHeight(), mImg1.getWidth())
        # We need the calibrated, full size original
        # MaskedImages for the variance plane calculations
        exposure1.maskedImage = mImg1
        exposure2.maskedImage = mImg2
        # TODO DM-25174 : Here we need actually not psfShape but an
        # estimation of the size of Pd and Ps
        # worst case scenario a padding of imgShape (? TBC)
        self.computeCommonShape(self.imgShape, self.psfShape1, self.psfShape2)

        self.fullExp1 = exposure1
        self.fullExp2 = exposure2

        self.fftFullIm1 = self.padAndFftImage(mImg1.image.array)
        self.fftVarPl1 = self.padAndFftImage(mImg1.variance.array)
        self.fftFullIm2 = self.padAndFftImage(mImg2.image.array)
        self.fftVarPl2 = self.padAndFftImage(mImg2.variance.array)

    def prepareSubExposure(self, bbox1=None, bbox2=None, psf1=None, psf2=None, sig1=None, sig2=None):
        """Perform per-sub exposure preparations.

        Parameters
        ----------
        sig1, sig2 : `float`, optional
            For debug purposes only, copnsider that the image
            may already be rescaled by the photometric calibration.
        bbox1, bbox2 : `lsst.geom.Box2I`, optional
            If specified, the region of the full exposure to use.

        psf1, psf2 : `lsst.afw.detection.Psf`, optional
            If specified, use given psf as the sub exposure psf. For debug purposes.

        sig1, sig2 : `float`, optional
            If specified, use value as the sub-exposures' background noise sigma value.

        Returns
        -------
        None

        Notes
        -----
        TODO DM-23855: Performing ZOGY on a grid is not yet implemented.

        Raises
        ------
        ValueError: If sub-exposure dimensions mismatch.
        """
        if bbox1 is None:
            subExposure1 = self.fullExp1
        else:
            subExposure1 = self.fullExp1.Factory(self.exposure1, bbox1)
        if bbox2 is None:
            subExposure2 = self.fullExp2
        else:
            subExposure2 = self.fullExp2.Factory(self.exposure2, bbox2)

        if subExposure1.getDimensions() != subExposure2.getDimensions():
            raise ValueError("Subexposure dimensions do not match.")

        if psf1 is None:
            self.subExpPsf1 = self.computePsfAtCenter(subExposure1)
        if psf2 is None:
            self.subExpPsf2 = self.computePsfAtCenter(subExposure2)
        self.checkCentroids(self.subExpPsf1.array, self.subExpPsf2.array)
        # sig1 and sig2  should not be set externally, just for debug purpose
        if sig1 is None:
            sig1 = np.sqrt(self._computeVarianceMean(subExposure1))
        self.subExpVar1 = sig1*sig1
        if sig2 is None:
            sig2 = np.sqrt(self._computeVarianceMean(subExposure2))
        self.subExpVar2 = sig2*sig2

        D = self.padCenterOriginArray(self.subExpPsf1.array, self.freqSpaceShape)
        self.psfFft1 = np.fft.fft2(D)
        D = self.padCenterOriginArray(self.subExpPsf2.array, self.freqSpaceShape)
        self.psfFft2 = np.fft.fft2(D)

        self.subExposure1 = subExposure1
        self.subExposure2 = subExposure2

    @staticmethod
    def realSq(D):
        """Square the argument in pixel space.

        Parameters
        ----------
        D : 2D `numpy.ndarray` of `numpy.complex`
            Fourier transform of a real valued array.

        Returns
        -------
        R : `numpy.ndarray` of `numpy.complex`

        Notes
        -----
        ``D`` is to be inverse Fourier transformed, squared and then
        forward Fourier transformed again, i.e. an autoconvolution in Fourier space.
        This operation is not distributive over multiplication.
        ``realSq(A*B) != realSq(A)*realSq(B)``
        """
        D = np.real(np.fft.ifft2(D))
        D *= D
        D = np.fft.fft2(D)
        return D

    @staticmethod
    def getCentroid(A):
        """Calculate the centroid coordinates of a 2D array.

        Parameters
        ----------
        A : 2D `numpy.ndarray` of `float`
            The input array. Must not be all exact zero.

        Notes
        -----
        Calculates the centroid as if the array represented a 2D geometrical shape with
        weights per cell, allowing for "negative" weights. If sum equals to exact (float) zero,
        calculates centroid of absolute value array.

        The geometrical center is defined as (0,0), independently of the array shape.
        For an odd dimension, this is the center of the center pixel,
        for an even dimension, this is between the two center pixels.

        Returns
        -------
        ycen, xcen : `tuple` of `float`

        """
        s = np.sum(A)
        if s == 0.:
            A = np.fabs(A)
            s = np.sum(A)
        w = np.arange(A.shape[0], dtype=float) - (A.shape[0] - 1.)/2
        ycen = np.sum(w[:, np.newaxis]*A)/s
        w = np.arange(A.shape[1], dtype=float) - (A.shape[1] - 1.)/2
        xcen = np.sum(w[np.newaxis, :]*A)/s

        return ycen, xcen

    def checkCentroids(self, psfArr1, psfArr2):
        """Check whether two PSF array centroids' distance is within tolerance.

        Parameters
        ----------
        psfArr1, psfArr2 : `numpy.ndarray` of `float`
            Input PSF arrays to check.

        Returns
        -------
        None

        Raises
        ------
        ValueError:
            Centroid distance exceeds `config.maxPsfCentroidDist` pixels.
        """
        yc1, xc1 = self.getCentroid(psfArr1)
        yc2, xc2 = self.getCentroid(psfArr2)
        dy = yc2 - yc1
        dx = xc2 - xc1
        if dy*dy + dx*dx > self.config.maxPsfCentroidDist*self.config.maxPsfCentroidDist:
            raise ValueError(
                f"PSF centroids are offset by more than {self.config.maxPsfCentroidDist:.2f} pixels.")

    def calculateFourierDiffim(self, psf1, im1, varPlane1, F1, varMean1,
                               psf2, im2, varPlane2, F2, varMean2, calculateS):
        """Calculates the difference image, ``im1-im2``.

        Parameters
        ----------
        psf1, psf2, im1, im2, varPlane1, varPlane2 : `numpy.ndarray` of `numpy.complex`,
        shape ``self.freqSpaceShape``
            Psf, image and variance plane arrays respectively.
            All arrays must be already in Fourier space.

        varMean1, varMean2: `numpy.float` > 0.
            Average per-pixel noise variance in im1, im2 respectively. Used as weighing
            of input images. Must be greater than zero.

        F1, F2 : `numpy.float` > 0.
            Photometric scaling of the images. See eqs. (5)--(9)

        calculateS : `bool`
            If True, calculate and return the detection significance (score) image.

        Returns
        -------
        result : `pipe.base.Struct`
            All arrays are in Fourier space and have shape ``self.freqSpaceShape``.
            - ``Fd`` : `float`
                Photometric level of ``D``.
            - ``D`` : `numpy.ndarray` of `numpy.complex`
                The difference image.
            - ``varplaneD`` : `numpy.ndarray` of `numpy.complex`
                Variance plane of ``D``.
            - ``Pd`` : `numpy.ndarray` of `numpy.complex`
                PSF of ``D``.
            - ``S`` : `numpy.ndarray` of `numpy.complex` or `None`
                Significance (score) image.
            - ``varplaneS`` : `numpy.ndarray` of `numpy.complex` or `None`
                Variance plane of ``S``.
            - ``Ps`` : `numpy.ndarray` of `numpy.complex`
                PSF of ``S``.

        Notes
        -----
        All array inputs and outputs are Fourier-space images with size of
        `self.freqSpaceShape` in this method.

        ``varMean1``, ``varMean2`` quantities are part of the noise model and not to be confused
        with the variance of image frequency components or with ``varPlane1``, ``varPlane2`` that
        are the Fourier transform of the variance planes.
        """
        var1F2Sq = varMean1*F2*F2
        var2F1Sq = varMean2*F1*F1
        # We need reals for comparison, also real operations are usually faster
        psfAbsSq1 = np.real(np.conj(psf1)*psf1)
        psfAbsSq2 = np.real(np.conj(psf2)*psf2)
        FdDenom = np.sqrt(var1F2Sq + var2F1Sq)  # one number

        # PSF smoothing
        tiny = np.finfo(psf1.dtype).tiny * 100
        epsSq = np.finfo(psf1.dtype).eps

        self.log.debug(f"fPsf1 min {np.amin(np.abs(psf1))}")
        fltZero1 = psfAbsSq1 < epsSq
        self.log.debug(f"fPsf1 cut to zero {np.sum(fltZero1)}.")
        psf1[fltZero1] = 0.
        psfAbsSq1[fltZero1] = 0.

        self.log.debug(f"fPsf2 min {np.amin(np.abs(psf2))}")
        fltZero2 = psfAbsSq2 < epsSq
        self.log.debug(f"fPsf2 cut to zero {np.sum(fltZero2)}.")
        psf2[fltZero2] = 0.
        psfAbsSq2[fltZero2] = 0.

        fltZero = np.logical_and(fltZero1, fltZero2)
        nZero = np.sum(fltZero)
        fltZero1 = None
        fltZero2 = None

        # Secure positive limit to avoid floating point operations resulting in exact zero

        sDenom = var1F2Sq*psfAbsSq2 + var2F1Sq*psfAbsSq1  # array, eq. (12)
        # Frequencies where both psfs are too close to zero.
        # We expect this only in cases when psf1, psf2 are identical,
        # and either having very well sampled Gaussian tails
        # or having "edges" such that some sinc-like zero crossings are found at symmetry points
        #
        # if sDenom < tiny then it can be == 0. -> `denom` = 0. and 0/0 occur at `c1` , `c2`
        # if we keep SDenom = tiny, denom ~ O(sqrt(tiny)), Pd ~ O(sqrt(tiny)), S ~ O(sqrt(tiny)*tiny) == 0
        # Where S = 0 then Pd = 0 and D should still yield the same variance ~ O(1)
        # For safety, we set S = 0 explicitly, too, though it should be unnecessary.
        self.log.debug(f"Handling {nZero} both PSFs are zero points.")
        if nZero > 0:
            fltZero = np.nonzero(fltZero)  # We expect only a small fraction of such frequencies
            sDenom[fltZero] = tiny  # Avoid division problem but overwrite result anyway
        denom = np.sqrt(sDenom)  # array, eq. (13)

        c1 = F2*psf2/denom
        c2 = F1*psf1/denom
        if nZero > 0:
            c1[fltZero] = F2/FdDenom
            c2[fltZero] = F1/FdDenom
        D = c1*im1 - c2*im2  # Difference image eq. (13)
        varPlaneD = self.realSq(c1)*varPlane1 + self.realSq(c2)*varPlane2  # eq. (26)

        Pd = FdDenom*psf1*psf2/denom  # Psf of D eq. (14)
        if nZero > 0:
            Pd[fltZero] = 0.

        Fd = F1*F2/FdDenom  # Flux scaling of D eq. (15)
        if calculateS:
            c1 = F1*F2*F2*np.conj(psf1)*psfAbsSq2/sDenom
            c2 = F2*F1*F1*np.conj(psf2)*psfAbsSq1/sDenom
            if nZero > 0:
                c1[fltZero] = 0.
                c2[fltZero] = 0.
            S = c1*im1 - c2*im2  # eq. (12)
            varPlaneS = self.realSq(c1)*varPlane1 + self.realSq(c2)*varPlane2
            Ps = np.conj(Pd)*Pd  # eq. (17) Source detection expects a PSF
        else:
            S = None
            Ps = None
            varPlaneS = None
        return pipeBase.Struct(D=D, Pd=Pd, varPlaneD=varPlaneD, Fd=Fd,
                               S=S, Ps=Ps, varPlaneS=varPlaneS)

    @staticmethod
    def calculateMaskPlane(mask1, mask2, effPsf1=None, effPsf2=None):
        """Calculate the mask plane of the difference image.

        Parameters
        ----------
        mask1, maks2 : `lsst.afw.image.Mask`
            Mask planes of the two exposures.


        Returns
        -------
        diffmask : `lsst.afw.image.Mask`
            Mask plane for the subtraction result.

        Notes
        -----
        TODO DM-25174 : Specification of effPsf1, effPsf2 are not yet supported.
        """

        # mask1 x effPsf2 | mask2 x effPsf1
        if effPsf1 is not None or effPsf2 is not None:
            # TODO: DM-25174 effPsf1, effPsf2: the effective psf for cross-blurring.
            # We need a "size" approximation of the c1 and c2 coefficients to make effPsfs
            # Also convolution not yet supports mask-only operation
            raise NotImplementedError("Mask plane only 'convolution' operation is not yet supported")
            # if effPsf1 is not None:
            #     mask1 = mask1.clone()
            #     afwMath.convolve(mask1, effPsf2)
            # if effPsf2 is not None:
            #     mask2 = mask2.clone()
            #     afwMath.convolve(mask2, effPsf1)
        R = mask1.clone()
        R |= mask2
        return R

    def makeDiffimSubExposure(self, ftDiff):
        """Wrap array results into Exposure objects.

        Parameters
        ----------
        ftDiff : `lsst.pipe.base.Struct`
            Result struct by `calculateFourierDiffim`.

        Returns
        -------
        resultName : `lsst.pipe.base.Struct`
            - ``diffSubExp`` : `lsst.afw.image.Exposure`
                The difference (sub)exposure. The exposure is calibrated
                in its pixel values, and has a constant `PhotoCalib` object of 1.
            - ``scoreSubExp`` : `lsst.afw.image.Exposure` or `None`
                The score (sub)exposure if it was calculated.
        """
        D = self.inverseFftAndCropImage(
            ftDiff.D, self.imgShape, np.logical_or(self.fftFullIm1.filtInf, self.fftFullIm2.filtInf),
            np.logical_or(self.fftFullIm1.filtNaN, self.fftFullIm2.filtNaN),
            dtype=self.subExposure1.image.dtype)
        varPlaneD = self.inverseFftAndCropImage(
            ftDiff.varPlaneD, self.imgShape, np.logical_or(self.fftVarPl1.filtInf, self.fftVarPl2.filtInf),
            np.logical_or(self.fftVarPl1.filtNaN, self.fftVarPl2.filtNaN),
            dtype=self.subExposure1.variance.dtype)
        Pd = self.inverseFftAndCropImage(
            ftDiff.Pd, self.psfShape1, dtype=self.subExpPsf1.dtype)
        sumPd = np.sum(Pd)
        # If this is smaller than 1. it is an indicator that it does not fit its original dimensions
        self.log.info(f"Pd sum before normalization: {sumPd:.3f}")
        Pd /= sumPd

        diffSubExposure = self.subExposure1.clone()
        # Indices of the subexposure bbox in the full image array
        bbox = self.subExposure1.getBBox()
        arrIndex = bbox.getMin() - self.fullExp1.getXY0()
        diffSubExposure.image.array = D[
            arrIndex.getY():arrIndex.getY() + bbox.getHeight(),
            arrIndex.getX():arrIndex.getX() + bbox.getWidth()]
        diffSubExposure.variance.array = varPlaneD[
            arrIndex.getY():arrIndex.getY() + bbox.getHeight(),
            arrIndex.getX():arrIndex.getX() + bbox.getWidth()]
        diffSubExposure.mask = self.calculateMaskPlane(self.subExposure1.mask, self.subExposure2.mask)

        # PhotoCalib does not support ImageD.
        # calib = afwImage.PhotoCalib(1./ftDiff.Fd)
        # calibImg = calib.calibrateImage(
        #     afwImage.MaskedImage(diffSubExposure.maskedImage, deep=True, dtype=np.float32))
        # diffSubExposure.maskedImage = calibImg
        diffSubExposure.maskedImage /= ftDiff.Fd

        # Now the subExposure calibration is 1. everywhere
        calibOne = afwImage.PhotoCalib(1.)
        diffSubExposure.setPhotoCalib(calibOne)

        # Set the PSF of this subExposure
        self.log.debug(f"PSF dim {self.subExpPsf1.getDimensions()}")
        self.log.debug(f"Pd shape {Pd.shape}")

        psfImg = self.subExpPsf1.Factory(geom.Extent2I(Pd.shape[1], Pd.shape[0]))
        psfImg.array = Pd
        psfNew = measAlg.KernelPsf(afwMath.FixedKernel(psfImg))
        diffSubExposure.setPsf(psfNew)

        if ftDiff.S is not None:
            S = self.inverseFftAndCropImage(
                ftDiff.S, self.imgShape, np.logical_or(self.fftFullIm1.filtInf, self.fftFullIm2.filtInf),
                np.logical_or(self.fftFullIm1.filtNaN, self.fftFullIm2.filtNaN),
                dtype=self.subExposure1.image.dtype)
            varPlaneS = self.inverseFftAndCropImage(
                ftDiff.varPlaneS, self.imgShape,
                np.logical_or(self.fftVarPl1.filtInf, self.fftVarPl2.filtInf),
                np.logical_or(self.fftVarPl1.filtNaN, self.fftVarPl2.filtNaN),
                dtype=self.subExposure1.variance.dtype)

            S = S[arrIndex.getY():arrIndex.getY() + bbox.getHeight(),
                  arrIndex.getX():arrIndex.getX() + bbox.getWidth()]
            varPlaneS = varPlaneS[arrIndex.getY():arrIndex.getY() + bbox.getHeight(),
                                  arrIndex.getX():arrIndex.getX() + bbox.getWidth()]
            # PSF of S
            Ps = self.inverseFftAndCropImage(ftDiff.Ps, self.psfShape1, dtype=self.subExpPsf1.dtype)
            sumPs = np.sum(Ps)
            self.log.info(f"Ps sum before normalization: {sumPs:.3f}")
            Ps /= sumPs

            # Ensure that no 0/0 occur in S/var(S).
            tiny = np.finfo(varPlaneS.dtype).tiny * 10
            fltZero = np.nonzero(varPlaneS < tiny)
            varPlaneS[fltZero] = tiny
            S[fltZero] = 0

            # TODO DM-23855 : Scorr corrections may be done here

            scoreSubExposure = self.subExposure1.clone()
            scoreSubExposure.image.array = S
            scoreSubExposure.variance.array = varPlaneS
            scoreSubExposure.mask = diffSubExposure.mask
            scoreSubExposure.setPhotoCalib(None)

            # Set the PSF of this subExposure
            psfSImg = self.subExpPsf1.Factory(geom.Extent2I(Ps.shape[1], Ps.shape[0]))
            psfSImg.array = Ps
            psfSNew = measAlg.KernelPsf(afwMath.FixedKernel(psfSImg))
            scoreSubExposure.setPsf(psfSNew)
        else:
            scoreSubExposure = None

        return pipeBase.Struct(diffSubExp=diffSubExposure, scoreSubExp=scoreSubExposure)

    def run(self, exposure1, exposure2, calculateS=True):
        """Perform zogy subtraction of exposure1 - exposure2. Task entry point.

        Parameters
        ----------
        exposure1, exposure2 : `lsst.afw.image.Exposure`
            Two exposures warped and matched into matching pixel dimensions.

        Returns
        -------
        resultName : `lsst.pipe.base.Struct`
            - ``diffExp`` : `lsst.afw.image.Exposure`
                The Zogy difference exposure (``exposure1`` - ``exposure2``).
            - ``scoreExp`` : `lsst.afw.image.Exposure`
                The Zogy score exposure.
            - ``ftDiff`` : `lsst.pipe.base.Struct`
                Lower level return struct by `calculateFourierDiffim` w/ added fields from the task instance.
                For debug purposes.

        Notes
        -----
        TODO DM-23855 : spatially varying solution on a grid is not yet implemented
        """
        # We use the dimensions of the 1st image only in the code
        if exposure1.getDimensions() != exposure2.getDimensions():
            raise ValueError("Exposure dimensions do not match.")

        self.prepareFullExposure(exposure1, exposure2, correctBackground=self.config.correctBackground)

        # TODO DM-23855: Add grid splitting support here for spatially varying PSF support
        # Passing exposure1,2 won't be ok here: they're not photometrically scaled.
        # Use the modified full maskedImages here
        self.prepareSubExposure()
        ftDiff = self.calculateFourierDiffim(
            self.psfFft1, self.fftFullIm1.imFft, self.fftVarPl1.imFft, self.F1, self.subExpVar1,
            self.psfFft2, self.fftFullIm2.imFft, self.fftVarPl2.imFft, self.F2, self.subExpVar2,
            calculateS=calculateS)
        diffExp = self.makeDiffimSubExposure(ftDiff)
        # Add debug info from the task instance
        ftDiff.freqSpaceShape = self.freqSpaceShape
        ftDiff.imgShape = self.imgShape
        ftDiff.psfShape1 = self.psfShape1
        ftDiff.psfShape2 = self.psfShape2
        return pipeBase.Struct(diffExp=diffExp.diffSubExp,
                               scoreExp=diffExp.scoreSubExp,
                               ftDiff=ftDiff)


class ZogyImagePsfMatchConfig(ImagePsfMatchConfig):
    """Config for the ZogyImagePsfMatchTask"""

    zogyConfig = pexConfig.ConfigField(
        dtype=ZogyConfig,
        doc='ZogyTask config to use when running on complete exposure (non spatially-varying)',
    )


class ZogyImagePsfMatchTask(ImagePsfMatchTask):
    """Task to perform Zogy PSF matching and image subtraction.

    This class inherits from ImagePsfMatchTask to contain the _warper
    subtask and related methods.
    """

    ConfigClass = ZogyImagePsfMatchConfig

    def __init__(self, *args, **kwargs):
        ImagePsfMatchTask.__init__(self, *args, **kwargs)

    def run(self, scienceExposure, templateExposure, doWarping=True, spatiallyVarying=False):
        """Register, PSF-match, and subtract two Exposures, ``scienceExposure - templateExposure``
        using the ZOGY algorithm.

        Parameters
        ----------
        templateExposure : `lsst.afw.image.Exposure`
            exposure to be warped to scienceExposure.
        scienceExposure : `lsst.afw.image.Exposure`
            reference Exposure.
        doWarping : `bool`
            what to do if templateExposure's and scienceExposure's WCSs do not match:
            - if True then warp templateExposure to match scienceExposure
            - if False then raise an Exception
        spatiallyVarying : `bool`
            If True, perform the operation over a grid of patches across the two exposures

        Notes
        -----
        Do the following, in order:
            - Warp templateExposure to match scienceExposure, if their WCSs do not already match
            - Compute subtracted exposure ZOGY image subtraction algorithm on the two exposures

        This is the new entry point of the task as of DM-25115.


        Returns
        -------
        results : `lsst.pipe.base.Struct` containing these fields:
            - subtractedExposure: `lsst.afw.image.Exposure`
                The subtraction result.
            - warpedExposure: `lsst.afw.image.Exposure` or `None`
                templateExposure after warping to match scienceExposure
        """

        if spatiallyVarying:
            raise NotImplementedError(
                "DM-25115 Spatially varying zogy subtraction is not implemented.")

        if not self._validateWcs(scienceExposure, templateExposure):
            if doWarping:
                self.log.info("Warping templateExposure to scienceExposure")
                xyTransform = afwGeom.makeWcsPairTransform(templateExposure.getWcs(),
                                                           scienceExposure.getWcs())
                psfWarped = measAlg.WarpedPsf(templateExposure.getPsf(), xyTransform)
                templateExposure = self._warper.warpExposure(
                    scienceExposure.getWcs(), templateExposure, destBBox=scienceExposure.getBBox())
                templateExposure.setPsf(psfWarped)
            else:
                self.log.error("ERROR: Input images not registered")
                raise RuntimeError("Input images not registered")

        config = self.config.zogyConfig
        task = ZogyTask(config=config)
        results = task.run(scienceExposure, templateExposure)
        results.warpedExposure = templateExposure
        return results

    def subtractExposures(self, templateExposure, scienceExposure,
                          doWarping=True, spatiallyVarying=True, inImageSpace=False,
                          doPreConvolve=False):
        raise NotImplementedError

    def subtractMaskedImages(self, templateExposure, scienceExposure,
                             doWarping=True, spatiallyVarying=True, inImageSpace=False,
                             doPreConvolve=False):
        raise NotImplementedError


subtractAlgorithmRegistry.register('zogy', ZogyImagePsfMatchTask)
