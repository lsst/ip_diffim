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

from .imageMapReduce import (ImageMapReduceConfig, ImageMapper,
                             ImageMapReduceTask)
from .imagePsfMatch import (ImagePsfMatchTask, ImagePsfMatchConfig,
                            subtractAlgorithmRegistry)

__all__ = ["ZogyTask", "ZogyConfig",
           "ZogyMapper", "ZogyMapReduceConfig",
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
    inImageSpace = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Perform all convolutions in real (image) space rather than Fourier space. "
        "Currently if True, this results in artifacts when using real (noisy) PSFs."
    )

    padSize = pexConfig.Field(
        dtype=int,
        default=7,
        doc="Number of pixels to pad PSFs to avoid artifacts (when inImageSpace is True)"
    )

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

    doTrimKernels = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Trim kernels for image-space ZOGY. Speeds up convolutions and shrinks artifacts. "
        "Subject of future research."
    )

    doFilterPsfs = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Filter PSFs for image-space ZOGY. Aids in reducing artifacts. "
        "Subject of future research."
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


MIN_KERNEL = 1.0e-4


class ZogyTask(pipeBase.Task):
    """Task to perform ZOGY proper image subtraction. See module-level documentation for
    additional details.

    As of DM-25115, there are two entry points, `run()` and `subtractExposures()` to keep backward
    compatibility with the previous implementation.
    """
    ConfigClass = ZogyConfig
    _DefaultName = "ip_diffim_Zogy"

    def __init__(self, templateExposure=None, scienceExposure=None, sig1=None, sig2=None,
                 psf1=None, psf2=None, *args, **kwargs):
        """Create the ZOGY task.

        Parameters
        ----------
        templateExposure : `lsst.afw.image.Exposure`
            Template exposure ("Reference image" in ZOGY (2016)).
        scienceExposure : `lsst.afw.image.Exposure`
            Science exposure ("New image" in ZOGY (2016)). Must have already been
            registered and photmetrically matched to template.
        sig1 : `float`
            (Optional) sqrt(variance) of `templateExposure`. If `None`, it is
            computed from the sqrt(mean) of the `templateExposure` variance image.
        sig2 : `float`
            (Optional) sqrt(variance) of `scienceExposure`. If `None`, it is
            computed from the sqrt(mean) of the `scienceExposure` variance image.
        psf1 : 2D `numpy.array`
            (Optional) 2D array containing the PSF image for the template. If
            `None`, it is extracted from the PSF taken at the center of `templateExposure`.
        psf2 : 2D `numpy.array`
            (Optional) 2D array containing the PSF image for the science img. If
            `None`, it is extracted from the PSF taken at the center of `scienceExposure`.
        *args
            additional arguments to be passed to
            `lsst.pipe.base.Task`
        **kwargs
            additional keyword arguments to be passed to
            `lsst.pipe.base.Task`
        """
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.template = self.science = None
        self.setup(templateExposure=templateExposure, scienceExposure=scienceExposure,
                   sig1=sig1, sig2=sig2, psf1=psf1, psf2=psf2, *args, **kwargs)

    def setup(self, templateExposure=None, scienceExposure=None, sig1=None, sig2=None,
              psf1=None, psf2=None, correctBackground=False, *args, **kwargs):
        """Set up the ZOGY task.

        Parameters
        ----------
        templateExposure : `lsst.afw.image.Exposure`
            Template exposure ("Reference image" in ZOGY (2016)).
        scienceExposure : `lsst.afw.image.Exposure`
            Science exposure ("New image" in ZOGY (2016)). Must have already been
            registered and photometrically matched to template.
        sig1 : `float`
            (Optional) sqrt(variance) of `templateExposure`. If `None`, it is
            computed from the sqrt(mean) of the `templateExposure` variance image.
        sig2 : `float`
            (Optional) sqrt(variance) of `scienceExposure`. If `None`, it is
            computed from the sqrt(mean) of the `scienceExposure` variance image.
        psf1 : 2D `numpy.array`
            (Optional) 2D array containing the PSF image for the template. If
            `None`, it is extracted from the PSF taken at the center of `templateExposure`.
        psf2 : 2D `numpy.array`
            (Optional) 2D array containing the PSF image for the science img. If
            `None`, it is extracted from the PSF taken at the center of `scienceExposure`.
        correctBackground : `bool`
            (Optional) subtract sigma-clipped mean of exposures. Zogy doesn't correct
            nonzero backgrounds (unlike AL) so subtract them here.
        *args
            additional arguments to be passed to
            `lsst.pipe.base.Task`
        **kwargs
            additional keyword arguments to be passed to
            `lsst.pipe.base.Task`
        """
        if self.template is None and templateExposure is None:
            return
        if self.science is None and scienceExposure is None:
            return

        self.template = templateExposure
        self.science = scienceExposure

        self.statsControl = afwMath.StatisticsControl()
        self.statsControl.setNumSigmaClip(3.)
        self.statsControl.setNumIter(3)
        self.statsControl.setAndMask(afwImage.Mask.getPlaneBitMask(
            self.config.ignoreMaskPlanes))

        self.im1 = self.template.getMaskedImage().getImage().getArray()
        self.im2 = self.science.getMaskedImage().getImage().getArray()
        self.im1_var = self.template.getMaskedImage().getVariance().getArray()
        self.im2_var = self.science.getMaskedImage().getVariance().getArray()

        def selectPsf(psf, exposure):
            if psf is not None:
                return psf
            else:
                bbox1 = self.template.getBBox()
                xcen = (bbox1.getBeginX() + bbox1.getEndX()) / 2.
                ycen = (bbox1.getBeginY() + bbox1.getEndY()) / 2.
                return exposure.getPsf().computeKernelImage(geom.Point2D(xcen, ycen)).getArray()

        self.im1_psf = selectPsf(psf1, self.template)
        self.im2_psf = selectPsf(psf2, self.science)

        # Make sure PSFs are the same size. Messy, but should work for all cases.
        psf1 = self.im1_psf
        psf2 = self.im2_psf
        pShape1 = psf1.shape
        pShape2 = psf2.shape
        if (pShape1[0] < pShape2[0]):
            psf1 = np.pad(psf1, ((0, pShape2[0] - pShape1[0]), (0, 0)), mode='constant', constant_values=0.)
        elif (pShape2[0] < pShape1[0]):
            psf2 = np.pad(psf2, ((0, pShape1[0] - pShape2[0]), (0, 0)), mode='constant', constant_values=0.)
        if (pShape1[1] < pShape2[1]):
            psf1 = np.pad(psf1, ((0, 0), (0, pShape2[1] - pShape1[1])), mode='constant', constant_values=0.)
        elif (pShape2[1] < pShape1[1]):
            psf2 = np.pad(psf2, ((0, 0), (0, pShape1[1] - pShape2[1])), mode='constant', constant_values=0.)

        # PSFs' centers may be offset relative to each other; now fix that!
        maxLoc1 = np.unravel_index(np.argmax(psf1), psf1.shape)
        maxLoc2 = np.unravel_index(np.argmax(psf2), psf2.shape)
        # *Very* rarely happens but if they're off by >1 pixel, do it more than once.
        while (maxLoc1[0] != maxLoc2[0]) or (maxLoc1[1] != maxLoc2[1]):
            if maxLoc1[0] > maxLoc2[0]:
                psf2[1:, :] = psf2[:-1, :]
            elif maxLoc1[0] < maxLoc2[0]:
                psf1[1:, :] = psf1[:-1, :]
            if maxLoc1[1] > maxLoc2[1]:
                psf2[:, 1:] = psf2[:, :-1]
            elif maxLoc1[1] < maxLoc2[1]:
                psf1[:, 1:] = psf1[:, :-1]
            maxLoc1 = np.unravel_index(np.argmax(psf1), psf1.shape)
            maxLoc2 = np.unravel_index(np.argmax(psf2), psf2.shape)

        # Make sure there are no div-by-zeros
        psf1[psf1 < MIN_KERNEL] = MIN_KERNEL
        psf2[psf2 < MIN_KERNEL] = MIN_KERNEL

        self.im1_psf = psf1
        self.im2_psf = psf2

        self.sig1 = np.sqrt(self._computeVarianceMean(self.template)) if sig1 is None else sig1
        self.sig2 = np.sqrt(self._computeVarianceMean(self.science)) if sig2 is None else sig2
        # if sig1 or sig2 are NaN, then the entire region being Zogy-ed is masked.
        # Don't worry about it - the result will be masked but avoid warning messages.
        if np.isnan(self.sig1) or self.sig1 == 0:
            self.sig1 = 1.
        if np.isnan(self.sig2) or self.sig2 == 0:
            self.sig2 = 1.

        # Zogy doesn't correct nonzero backgrounds (unlike AL) so subtract them here.
        if correctBackground:
            def _subtractImageMean(exposure):
                """Compute the sigma-clipped mean of the image of `exposure`."""
                mi = exposure.getMaskedImage()
                statObj = afwMath.makeStatistics(mi.getImage(), mi.getMask(),
                                                 afwMath.MEANCLIP, self.statsControl)
                mean = statObj.getValue(afwMath.MEANCLIP)
                if not np.isnan(mean):
                    mi -= mean

            _subtractImageMean(self.template)
            _subtractImageMean(self.science)

        # Define the normalization of each image from the config
        self.Fr = self.config.templateFluxScaling  # default is 1
        self.Fn = self.config.scienceFluxScaling  # default is 1
        # If 'scaleByCalibration' is True then these norms are overwritten
        if self.config.scaleByCalibration:
            calib_template = self.template.getPhotoCalib()
            calib_science = self.science.getPhotoCalib()
            if calib_template is None:
                self.log.warning("No calibration information available for template image.")
            if calib_science is None:
                self.log.warning("No calibration information available for science image.")
            if calib_template is None or calib_science is None:
                self.log.warning("Due to lack of calibration information, "
                                 "reverting to templateFluxScaling and scienceFluxScaling.")
            else:
                self.Fr = 1/calib_template.getCalibrationMean()
                self.Fn = 1/calib_science.getCalibrationMean()

        self.log.info("Setting template image scaling to Fr=%f" % self.Fr)
        self.log.info("Setting science  image scaling to Fn=%f" % self.Fn)

        self.padSize = self.config.padSize  # default is 7

    def _computeVarianceMean(self, exposure):
        """Compute the sigma-clipped mean of the variance image of `exposure`.
        """
        statObj = afwMath.makeStatistics(exposure.getMaskedImage().getVariance(),
                                         exposure.getMaskedImage().getMask(),
                                         afwMath.MEANCLIP, self.statsControl)
        var = statObj.getValue(afwMath.MEANCLIP)
        return var

    @staticmethod
    def _padPsfToSize(psf, size):
        """Zero-pad `psf` to the dimensions given by `size`.

        Parameters
        ----------
        psf : 2D `numpy.array`
            Input psf to be padded
        size : `list`
            Two element list containing the dimensions to pad the `psf` to

        Returns
        -------
        psf : 2D `numpy.array`
            The padded copy of the input `psf`.
        """
        newArr = np.zeros(size)
        # The center of the PSF sould be placed in the center-right.
        offset = [size[0]//2 - psf.shape[0]//2, size[1]//2 - psf.shape[1]//2]
        tmp = newArr[offset[0]:(psf.shape[0] + offset[0]), offset[1]:(psf.shape[1] + offset[1])]
        tmp[:, :] = psf
        return newArr

    def computePrereqs(self, psf1=None, psf2=None, padSize=0):
        """Compute standard ZOGY quantities used by (nearly) all methods.

        Many of the ZOGY calculations require similar quantities, including
        FFTs of the PSFs, and the "denominator" term (e.g. in eq. 13 of
        ZOGY manuscript (2016). This function consolidates many of those
        operations.

        Parameters
        ----------
        psf1 : 2D `numpy.array`
            (Optional) Input psf of template, override if already padded
        psf2 : 2D `numpy.array`
            (Optional) Input psf of science image, override if already padded
        padSize : `int`, optional
            Number of pixels to pad the image on each side with zeroes.

        Returns
        -------
        A `lsst.pipe.base.Struct` containing:
        - Pr : 2D `numpy.array`, the (possibly zero-padded) template PSF
        - Pn : 2D `numpy.array`, the (possibly zero-padded) science PSF
        - Pr_hat : 2D `numpy.array`, the FFT of `Pr`
        - Pn_hat : 2D `numpy.array`, the FFT of `Pn`
        - denom : 2D `numpy.array`, the denominator of equation (13) in ZOGY (2016) manuscript
        - Fd : `float`, the relative flux scaling factor between science and template
        """
        psf1 = self.im1_psf if psf1 is None else psf1
        psf2 = self.im2_psf if psf2 is None else psf2
        padSize = self.padSize if padSize is None else padSize
        Pr, Pn = psf1, psf2
        if padSize > 0:
            Pr = ZogyTask._padPsfToSize(psf1, (psf1.shape[0] + padSize, psf1.shape[1] + padSize))
            Pn = ZogyTask._padPsfToSize(psf2, (psf2.shape[0] + padSize, psf2.shape[1] + padSize))
        # Make sure there are no div-by-zeros
        psf1[np.abs(psf1) <= MIN_KERNEL] = MIN_KERNEL
        psf2[np.abs(psf2) <= MIN_KERNEL] = MIN_KERNEL

        sigR, sigN = self.sig1, self.sig2
        Pr_hat = np.fft.fft2(Pr)
        Pr_hat2 = np.conj(Pr_hat) * Pr_hat
        Pn_hat = np.fft.fft2(Pn)
        Pn_hat2 = np.conj(Pn_hat) * Pn_hat
        denom = np.sqrt((sigN**2 * self.Fr**2 * Pr_hat2) + (sigR**2 * self.Fn**2 * Pn_hat2))
        Fd = self.Fr * self.Fn / np.sqrt(sigN**2 * self.Fr**2 + sigR**2 * self.Fn**2)

        res = pipeBase.Struct(
            Pr=Pr, Pn=Pn, Pr_hat=Pr_hat, Pn_hat=Pn_hat, denom=denom, Fd=Fd
        )
        return res

    def computeDiffimFourierSpace(self, debug=False, returnMatchedTemplate=False, **kwargs):
        r"""Compute ZOGY diffim `D` as proscribed in ZOGY (2016) manuscript

        Parameters
        ----------
        debug : `bool`, optional
            If set to True, filter the kernels by setting the edges to zero.
        returnMatchedTemplate : `bool`, optional
            Calculate the template image.
            If not set, the returned template will be None.

        Notes
        -----
        In all functions, im1 is R (reference, or template) and im2 is N (new, or science)
        Compute the ZOGY eqn. (13):

        .. math::

            \widehat{D} = \frac{Fr\widehat{Pr}\widehat{N} -
            F_n\widehat{Pn}\widehat{R}}{\sqrt{\sigma_n^2 Fr^2
            \|\widehat{Pr}\|^2 + \sigma_r^2 F_n^2 \|\widehat{Pn}\|^2}}

        where :math:`D` is the optimal difference image, :math:`R` and :math:`N` are the
        reference and "new" image, respectively, :math:`Pr` and :math:`P_n` are their
        PSFs, :math:`Fr` and :math:`Fn` are their flux-based zero-points (which we
        will set to one here), :math:`\sigma_r^2` and :math:`\sigma_n^2` are their
        variance, and :math:`\widehat{D}` denotes the FT of :math:`D`.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            - ``D`` : 2D `numpy.array`, the proper image difference
            - ``D_var`` : 2D `numpy.array`, the variance image for `D`
        """
        # Do all in fourier space (needs image-sized PSFs)
        psf1 = ZogyTask._padPsfToSize(self.im1_psf, self.im1.shape)
        psf2 = ZogyTask._padPsfToSize(self.im2_psf, self.im2.shape)

        preqs = self.computePrereqs(psf1, psf2, padSize=0)  # already padded the PSFs

        def _filterKernel(K, trim_amount):
            # Filter the wings of Kn, Kr, set to zero
            ps = trim_amount
            K[:ps, :] = K[-ps:, :] = 0
            K[:, :ps] = K[:, -ps:] = 0
            return K

        Kr_hat = self.Fr * preqs.Pr_hat / preqs.denom
        Kn_hat = self.Fn * preqs.Pn_hat / preqs.denom
        if debug and self.config.doTrimKernels:  # default False
            # Suggestion from Barak to trim Kr and Kn to remove artifacts
            # Here we just filter them (in image space) to keep them the same size
            ps = (Kn_hat.shape[1] - 80)//2
            Kn = _filterKernel(np.fft.ifft2(Kn_hat), ps)
            Kn_hat = np.fft.fft2(Kn)
            Kr = _filterKernel(np.fft.ifft2(Kr_hat), ps)
            Kr_hat = np.fft.fft2(Kr)

        def processImages(im1, im2, doAdd=False):
            # Some masked regions are NaN or infinite!, and FFTs no likey.
            im1[np.isinf(im1)] = np.nan
            im1[np.isnan(im1)] = np.nanmean(im1)
            im2[np.isinf(im2)] = np.nan
            im2[np.isnan(im2)] = np.nanmean(im2)

            R_hat = np.fft.fft2(im1)
            N_hat = np.fft.fft2(im2)

            D_hat = Kr_hat * N_hat
            D_hat_R = Kn_hat * R_hat
            if not doAdd:
                D_hat -= D_hat_R
            else:
                D_hat += D_hat_R

            D = np.fft.ifft2(D_hat)
            D = np.fft.ifftshift(D.real) / preqs.Fd

            R = None
            if returnMatchedTemplate:
                R = np.fft.ifft2(D_hat_R)
                R = np.fft.ifftshift(R.real) / preqs.Fd

            return D, R

        # First do the image
        D, R = processImages(self.im1, self.im2, doAdd=False)
        # Do the exact same thing to the var images, except add them
        D_var, R_var = processImages(self.im1_var, self.im2_var, doAdd=True)

        return pipeBase.Struct(D=D, D_var=D_var, R=R, R_var=R_var)

    def _doConvolve(self, exposure, kernel, recenterKernel=False):
        """Convolve an Exposure with a decorrelation convolution kernel.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Input exposure to be convolved.
        kernel : `numpy.array`
            2D `numpy.array` to convolve the image with
        recenterKernel : `bool`, optional
            Force the kernel center to the pixel with the maximum value.

        Returns
        -------
        A new `lsst.afw.image.Exposure` with the convolved pixels and the (possibly
        re-centered) kernel.

        Notes
        -----
        - We optionally re-center the kernel if necessary and return the possibly
          re-centered kernel
        """
        kernelImg = afwImage.ImageD(kernel.shape[1], kernel.shape[0])
        kernelImg.getArray()[:, :] = kernel
        kern = afwMath.FixedKernel(kernelImg)
        if recenterKernel:
            maxloc = np.unravel_index(np.argmax(kernel), kernel.shape)
            kern.setCtr(geom.Point2I(maxloc))
        outExp = exposure.clone()  # Do this to keep WCS, PSF, masks, etc.
        convCntrl = afwMath.ConvolutionControl(doNormalize=False, doCopyEdge=False,
                                               maxInterpolationDistance=0)
        try:
            afwMath.convolve(outExp.getMaskedImage(), exposure.getMaskedImage(), kern, convCntrl)
        except AttributeError:
            # Allow exposure to actually be an image/maskedImage
            # (getMaskedImage will throw AttributeError in that case)
            afwMath.convolve(outExp, exposure, kern, convCntrl)

        return outExp, kern

    def computeDiffimImageSpace(self, padSize=None, debug=False, **kwargs):
        """Compute ZOGY diffim `D` using image-space convlutions

        This method is still being debugged as it results in artifacts
        when the PSFs are noisy (see module-level docstring). Thus
        there are several options still enabled by the `debug` flag,
        which are disabled by defult.

        Parameters
        ----------
        padSize : `int`
           The amount to pad the PSFs by
        debug : `bool`
           Flag to enable debugging tests and options

        Returns
        -------
        D : `lsst.afw.Exposure`
           the proper image difference, including correct variance,
           masks, and PSF
        """
        preqs = self.computePrereqs(padSize=padSize)

        delta = 0.
        if debug:
            delta = 1.  # Regularize the ratio, a possible option to remove artifacts
        Kr_hat = (preqs.Pr_hat + delta) / (preqs.denom + delta)
        Kn_hat = (preqs.Pn_hat + delta) / (preqs.denom + delta)
        Kr = np.fft.ifft2(Kr_hat).real
        Kr = np.roll(np.roll(Kr, -1, 0), -1, 1)
        Kn = np.fft.ifft2(Kn_hat).real
        Kn = np.roll(np.roll(Kn, -1, 0), -1, 1)

        def _trimKernel(self, K, trim_amount):
            # Trim out the wings of Kn, Kr (see notebook #15)
            # only necessary if it's from a measured psf and PsfEx seems to always make PSFs of size 41x41
            ps = trim_amount
            K = K[ps:-ps, ps:-ps]
            return K

        padSize = self.padSize if padSize is None else padSize
        # Enabling this block (debug=True) makes it slightly faster, but ~25% worse artifacts:
        if debug and self.config.doTrimKernels:  # default False
            # Filtering also makes it slightly faster (zeros are ignored in convolution)
            # but a bit worse. Filter the wings of Kn, Kr (see notebook #15)
            Kn = _trimKernel(Kn, padSize)
            Kr = _trimKernel(Kr, padSize)

        # Note these are reverse-labelled, this is CORRECT!
        exp1, _ = self._doConvolve(self.template, Kn)
        exp2, _ = self._doConvolve(self.science, Kr)
        D = exp2
        tmp = D.getMaskedImage()
        tmp -= exp1.getMaskedImage()
        tmp /= preqs.Fd
        return pipeBase.Struct(D=D, R=exp1)

    def _setNewPsf(self, exposure, psfArr):
        """Utility method to set an exposure's PSF when provided as a 2-d numpy.array
        """
        bbox = exposure.getBBox()
        center = ((bbox.getBeginX() + bbox.getEndX()) // 2., (bbox.getBeginY() + bbox.getEndY()) // 2.)
        center = geom.Point2D(center[0], center[1])
        psfI = afwImage.ImageD(psfArr.shape[1], psfArr.shape[0])
        psfI.getArray()[:, :] = psfArr
        psfK = afwMath.FixedKernel(psfI)
        psfNew = measAlg.KernelPsf(psfK, center)
        exposure.setPsf(psfNew)
        return exposure

    def computeDiffim(self, inImageSpace=None, padSize=None,
                      returnMatchedTemplate=False, **kwargs):
        """Wrapper method to compute ZOGY proper diffim

        This method should be used as the public interface for
        computing the ZOGY diffim.

        Parameters
        ----------
        inImageSpace : `bool`
           Override config `inImageSpace` parameter
        padSize : `int`
           Override config `padSize` parameter
        returnMatchedTemplate : `bool`
           Include the PSF-matched template in the results Struct
        **kwargs
            additional keyword arguments to be passed to
            `computeDiffimFourierSpace` or `computeDiffimImageSpace`.

        Returns
        -------
        An lsst.pipe.base.Struct containing:
           - D : `lsst.afw.Exposure`
               the proper image difference, including correct variance,
               masks, and PSF
           - R : `lsst.afw.Exposure`
               If `returnMatchedTemplate` is True, the PSF-matched template
               exposure
        """
        R = None
        inImageSpace = self.config.inImageSpace if inImageSpace is None else inImageSpace
        if inImageSpace:
            padSize = self.padSize if padSize is None else padSize
            res = self.computeDiffimImageSpace(padSize=padSize, **kwargs)
            D = res.D
            if returnMatchedTemplate:
                R = res.R
        else:
            res = self.computeDiffimFourierSpace(**kwargs)
            D = self.science.clone()
            D.getMaskedImage().getImage().getArray()[:, :] = res.D
            D.getMaskedImage().getVariance().getArray()[:, :] = res.D_var
            if returnMatchedTemplate:
                R = self.science.clone()
                R.getMaskedImage().getImage().getArray()[:, :] = res.R
                R.getMaskedImage().getVariance().getArray()[:, :] = res.R_var

        psf = self.computeDiffimPsf()
        D = self._setNewPsf(D, psf)
        return pipeBase.Struct(D=D, R=R)

    def computeDiffimPsf(self, padSize=0, keepFourier=False, psf1=None, psf2=None):
        """Compute the ZOGY diffim PSF (ZOGY manuscript eq. 14)

        Parameters
        ----------
        padSize : `int`
           Override config `padSize` parameter
        keepFourier : `bool`
           Return the FFT of the diffim PSF (do not inverse-FFT it)
        psf1 : 2D `numpy.array`
            (Optional) Input psf of template, override if already padded
        psf2 : 2D `numpy.array`
            (Optional) Input psf of science image, override if already padded

        Returns
        -------
        Pd : 2D `numpy.array`
            The diffim PSF (or FFT of PSF if `keepFourier=True`)
        """
        preqs = self.computePrereqs(psf1=psf1, psf2=psf2, padSize=padSize)

        Pd_hat_numerator = (self.Fr * self.Fn * preqs.Pr_hat * preqs.Pn_hat)
        Pd_hat = Pd_hat_numerator / (preqs.Fd * preqs.denom)

        if keepFourier:
            return Pd_hat

        Pd = np.fft.ifft2(Pd_hat)
        Pd = np.fft.ifftshift(Pd).real

        return Pd

    def _computeVarAstGradients(self, xVarAst=0., yVarAst=0., inImageSpace=False,
                                R_hat=None, Kr_hat=None, Kr=None,
                                N_hat=None, Kn_hat=None, Kn=None):
        """Compute the astrometric noise correction terms

        Compute the correction for estimated astrometric noise as
        proscribed in ZOGY (2016), section 3.3. All convolutions
        performed either in real (image) or Fourier space.

        Parameters
        ----------
        xVarAst, yVarAst : `float`
           estimated astrometric noise (variance of astrometric registration errors)
        inImageSpace : `bool`
           Perform all convolutions in real (image) space rather than Fourier space
        R_hat : 2-D `numpy.array`
           (Optional) FFT of template image, only required if `inImageSpace=False`
        Kr_hat : 2-D `numpy.array`
           FFT of Kr kernel (eq. 28 of ZOGY (2016)), only required if `inImageSpace=False`
        Kr : 2-D `numpy.array`
           Kr kernel (eq. 28 of ZOGY (2016)), only required if `inImageSpace=True`.
           Kr is associated with the template (reference).
        N_hat : 2-D `numpy.array`
           FFT of science image, only required if `inImageSpace=False`
        Kn_hat : 2-D `numpy.array`
           FFT of Kn kernel (eq. 29 of ZOGY (2016)), only required if `inImageSpace=False`
        Kn : 2-D `numpy.array`
           Kn kernel (eq. 29 of ZOGY (2016)), only required if `inImageSpace=True`.
           Kn is associated with the science (new) image.

        Returns
        -------
        VastSR, VastSN : 2-D `numpy.array`
           Arrays containing the values in eqs. 30 and 32 of ZOGY (2016).
        """
        VastSR = VastSN = 0.
        if xVarAst + yVarAst > 0:  # Do the astrometric variance correction
            if inImageSpace:
                S_R, _ = self._doConvolve(self.template, Kr)
                S_R = S_R.getMaskedImage().getImage().getArray()
            else:
                S_R = np.fft.ifft2(R_hat * Kr_hat)
            gradRx, gradRy = np.gradient(S_R)
            VastSR = xVarAst * gradRx**2. + yVarAst * gradRy**2.

            if inImageSpace:
                S_N, _ = self._doConvolve(self.science, Kn)
                S_N = S_N.getMaskedImage().getImage().getArray()
            else:
                S_N = np.fft.ifft2(N_hat * Kn_hat)
            gradNx, gradNy = np.gradient(S_N)
            VastSN = xVarAst * gradNx**2. + yVarAst * gradNy**2.

        return VastSR, VastSN

    def computeScorrFourierSpace(self, xVarAst=0., yVarAst=0., **kwargs):
        """Compute corrected likelihood image, optimal for source detection

        Compute ZOGY S_corr image. This image can be thresholded for
        detection without optimal filtering, and the variance image is
        corrected to account for astrometric noise (errors in
        astrometric registration whether systematic or due to effects
        such as DCR). The calculations here are all performed in
        Fourier space, as proscribed in ZOGY (2016).

        Parameters
        ----------
        xVarAst, yVarAst : `float`
           estimated astrometric noise (variance of astrometric registration errors)

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            - ``S`` : `numpy.array`, the likelihood image S (eq. 12 of ZOGY (2016))
            - ``S_var`` : the corrected variance image (denominator of eq. 25 of ZOGY (2016))
            - ``Dpsf`` : the PSF of the diffim D, likely never to be used.
        """
        # Some masked regions are NaN or infinite!, and FFTs no likey.
        def fix_nans(im):
            """Replace any NaNs or Infs with the mean of the image."""
            isbad = ~np.isfinite(im)
            if np.any(isbad):
                im[isbad] = np.nan
                im[isbad] = np.nanmean(im)
            return im

        self.im1 = fix_nans(self.im1)
        self.im2 = fix_nans(self.im2)
        self.im1_var = fix_nans(self.im1_var)
        self.im2_var = fix_nans(self.im2_var)

        # Do all in fourier space (needs image-sized PSFs)
        psf1 = ZogyTask._padPsfToSize(self.im1_psf, self.im1.shape)
        psf2 = ZogyTask._padPsfToSize(self.im2_psf, self.im2.shape)

        preqs = self.computePrereqs(psf1, psf2, padSize=0)  # already padded the PSFs

        # Compute D_hat here (don't need D then, for speed)
        R_hat = np.fft.fft2(self.im1)
        N_hat = np.fft.fft2(self.im2)
        D_hat = self.Fr * preqs.Pr_hat * N_hat - self.Fn * preqs.Pn_hat * R_hat
        D_hat /= preqs.denom

        Pd_hat = self.computeDiffimPsf(padSize=0, keepFourier=True, psf1=psf1, psf2=psf2)
        Pd_bar = np.conj(Pd_hat)
        S = np.fft.ifft2(D_hat * Pd_bar)

        # Adjust the variance planes of the two images to contribute to the final detection
        # (eq's 26-29).
        Pn_hat2 = np.conj(preqs.Pn_hat) * preqs.Pn_hat
        Kr_hat = self.Fr * self.Fn**2. * np.conj(preqs.Pr_hat) * Pn_hat2 / preqs.denom**2.
        Pr_hat2 = np.conj(preqs.Pr_hat) * preqs.Pr_hat
        Kn_hat = self.Fn * self.Fr**2. * np.conj(preqs.Pn_hat) * Pr_hat2 / preqs.denom**2.

        Kr_hat2 = np.fft.fft2(np.fft.ifft2(Kr_hat)**2.)
        Kn_hat2 = np.fft.fft2(np.fft.ifft2(Kn_hat)**2.)
        var1c_hat = Kr_hat2 * np.fft.fft2(self.im1_var)
        var2c_hat = Kn_hat2 * np.fft.fft2(self.im2_var)

        # Do the astrometric variance correction
        fGradR, fGradN = self._computeVarAstGradients(xVarAst, yVarAst, inImageSpace=False,
                                                      R_hat=R_hat, Kr_hat=Kr_hat,
                                                      N_hat=N_hat, Kn_hat=Kn_hat)

        S_var = np.sqrt(np.fft.ifftshift(np.fft.ifft2(var1c_hat + var2c_hat)) + fGradR + fGradN)
        S_var *= preqs.Fd

        S = np.fft.ifftshift(np.fft.ifft2(Kn_hat * N_hat - Kr_hat * R_hat))
        S *= preqs.Fd

        Pd = self.computeDiffimPsf(padSize=0)
        return pipeBase.Struct(S=S.real, S_var=S_var.real, Dpsf=Pd)

    def computeScorrImageSpace(self, xVarAst=0., yVarAst=0., padSize=None, **kwargs):
        """Compute corrected likelihood image, optimal for source detection

        Compute ZOGY S_corr image. This image can be thresholded for
        detection without optimal filtering, and the variance image is
        corrected to account for astrometric noise (errors in
        astrometric registration whether systematic or due to effects
        such as DCR). The calculations here are all performed in
        Real (image) space.

        Parameters
        ----------
        xVarAst, yVarAst : `float`
           estimated astrometric noise (variance of astrometric registration errors)

        Returns
        -------
        A `lsst.pipe.base.Struct` containing:
        - S : `lsst.afw.image.Exposure`, the likelihood exposure S (eq. 12 of ZOGY (2016)),
            including corrected variance, masks, and PSF
        - D : `lsst.afw.image.Exposure`, the proper image difference, including correct
            variance, masks, and PSF
        """
        # Do convolutions in image space
        preqs = self.computePrereqs(padSize=0)

        padSize = self.padSize if padSize is None else padSize
        D = self.computeDiffimImageSpace(padSize=padSize).D
        Pd = self.computeDiffimPsf()
        D = self._setNewPsf(D, Pd)
        Pd_bar = np.fliplr(np.flipud(Pd))
        S, _ = self._doConvolve(D, Pd_bar)
        tmp = S.getMaskedImage()
        tmp *= preqs.Fd

        # Adjust the variance planes of the two images to contribute to the final detection
        # (eq's 26-29).
        Pn_hat2 = np.conj(preqs.Pn_hat) * preqs.Pn_hat
        Kr_hat = self.Fr * self.Fn**2. * np.conj(preqs.Pr_hat) * Pn_hat2 / preqs.denom**2.
        Pr_hat2 = np.conj(preqs.Pr_hat) * preqs.Pr_hat
        Kn_hat = self.Fn * self.Fr**2. * np.conj(preqs.Pn_hat) * Pr_hat2 / preqs.denom**2.

        Kr = np.fft.ifft2(Kr_hat).real
        Kr = np.roll(np.roll(Kr, -1, 0), -1, 1)
        Kn = np.fft.ifft2(Kn_hat).real
        Kn = np.roll(np.roll(Kn, -1, 0), -1, 1)
        var1c, _ = self._doConvolve(self.template.getMaskedImage().getVariance(), Kr**2.)
        var2c, _ = self._doConvolve(self.science.getMaskedImage().getVariance(), Kn**2.)

        # Do the astrometric variance correction
        fGradR, fGradN = self._computeVarAstGradients(xVarAst, yVarAst, inImageSpace=True,
                                                      Kr=Kr, Kn=Kn)

        Smi = S.getMaskedImage()
        Smi *= preqs.Fd
        S_var = np.sqrt(var1c.getArray() + var2c.getArray() + fGradR + fGradN)
        S.getMaskedImage().getVariance().getArray()[:, :] = S_var
        S = self._setNewPsf(S, Pd)

        # also return diffim since it was calculated and might be desired
        return pipeBase.Struct(S=S, D=D)

    def computeScorr(self, xVarAst=0., yVarAst=0., inImageSpace=None, padSize=0, **kwargs):
        """Wrapper method to compute ZOGY corrected likelihood image, optimal for
        source detection

        This method should be used as the public interface for
        computing the ZOGY S_corr.

        Parameters
        ----------
        xVarAst, yVarAst : `float`
           estimated astrometric noise (variance of astrometric registration errors)
        inImageSpace : `bool`
           Override config `inImageSpace` parameter
        padSize : `int`
           Override config `padSize` parameter

        Returns
        -------
        S : `lsst.afw.image.Exposure`
            The likelihood exposure S (eq. 12 of ZOGY (2016)),
            including corrected variance, masks, and PSF
        """
        inImageSpace = self.config.inImageSpace if inImageSpace is None else inImageSpace
        if inImageSpace:
            res = self.computeScorrImageSpace(xVarAst=xVarAst, yVarAst=yVarAst, padSize=padSize)
            S = res.S
        else:
            res = self.computeScorrFourierSpace(xVarAst=xVarAst, yVarAst=yVarAst)

            S = self.science.clone()
            S.getMaskedImage().getImage().getArray()[:, :] = res.S
            S.getMaskedImage().getVariance().getArray()[:, :] = res.S_var
            S = self._setNewPsf(S, res.Dpsf)

        return pipeBase.Struct(S=S)

    # DM-23855 modifications
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

    @staticmethod
    def computePsfAtCenter(exposure):
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

        # Secure positive limit to avoid floating point operations resulting in exact zero
        tiny = np.finfo(psf1.dtype).tiny * 100
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
        fltZero = sDenom < tiny
        nZero = np.sum(fltZero)
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
            Pd[fltZero] = 0

        Fd = F1*F2/FdDenom  # Flux scaling of D eq. (15)
        if calculateS:
            c1 = F1*F2*F2*np.conj(psf1)*psfAbsSq2/sDenom
            c2 = F2*F1*F1*np.conj(psf2)*psfAbsSq1/sDenom
            if nZero > 0:
                c1[fltZero] = 0
                c2[fltZero] = 0
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
        psfImg = self.subExpPsf1.Factory(self.subExpPsf1.getDimensions())
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
            psfSImg = self.subExpPsf1.Factory(self.subExpPsf1.getDimensions())
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


class ZogyMapper(ZogyTask, ImageMapper):
    """Task to be used as an ImageMapper for performing
    ZOGY image subtraction on a grid of subimages.
    """
    ConfigClass = ZogyConfig
    _DefaultName = 'ip_diffim_ZogyMapper'

    def __init__(self, *args, **kwargs):
        ImageMapper.__init__(self, *args, **kwargs)

    def run(self, subExposure, expandedSubExposure, fullBBox, template,
            **kwargs):
        """Perform ZOGY proper image subtraction on sub-images

        This method performs ZOGY proper image subtraction on
        `subExposure` using local measures for image variances and
        PSF. `subExposure` is a sub-exposure of the science image. It
        also requires the corresponding sub-exposures of the template
        (`template`). The operations are actually performed on
        `expandedSubExposure` to allow for invalid edge pixels arising
        from convolutions, which are then removed.

        Parameters
        ----------
        subExposure : `lsst.afw.image.Exposure`
            the sub-exposure of the diffim
        expandedSubExposure : `lsst.afw.image.Exposure`
            the expanded sub-exposure upon which to operate
        fullBBox : `lsst.geom.Box2I`
            the bounding box of the original exposure
        template : `lsst.afw.image.Exposure`
            the template exposure, from which a corresponding sub-exposure
            is extracted
        **kwargs
            additional keyword arguments propagated from
            `ImageMapReduceTask.run`. These include:

            ``doScorr`` : `bool`
                Compute and return the corrected likelihood image S_corr
                rather than the proper image difference
            ``inImageSpace`` : `bool`
                Perform all convolutions in real (image) space rather than
                in Fourier space. This option currently leads to artifacts
                when using real (measured and noisy) PSFs, thus it is set
                to `False` by default.
                These kwargs may also include arguments to be propagated to
                `ZogyTask.computeDiffim` and `ZogyTask.computeScorr`.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

                ``subExposure``: Either the subExposure of the proper image difference ``D``,
                    or (if `doScorr==True`) the corrected likelihood exposure ``S``.

        Notes
        -----
        This `run` method accepts parameters identical to those of
        `ImageMapper.run`, since it is called from the
        `ImageMapperTask`. See that class for more information.
        """
        bbox = subExposure.getBBox()
        center = ((bbox.getBeginX() + bbox.getEndX()) // 2., (bbox.getBeginY() + bbox.getEndY()) // 2.)
        center = geom.Point2D(center[0], center[1])

        imageSpace = kwargs.pop('inImageSpace', False)
        doScorr = kwargs.pop('doScorr', False)
        sigmas = kwargs.pop('sigmas', None)
        padSize = kwargs.pop('padSize', 7)

        # Psf and image for science img (index 2)
        subExp2 = expandedSubExposure

        # Psf and image for template img (index 1)
        subExp1 = template.Factory(template, expandedSubExposure.getBBox())

        if sigmas is None:
            sig1 = sig2 = None
        else:
            # for testing, can use the input sigma (e.g., global value for entire exposure)
            sig1, sig2 = sigmas[0], sigmas[1]

        def _makePsfSquare(psf):
            # Sometimes CoaddPsf does this. Make it square.
            if psf.shape[0] < psf.shape[1]:
                psf = np.pad(psf, ((0, psf.shape[1] - psf.shape[0]), (0, 0)), mode='constant',
                             constant_values=0.)
            elif psf.shape[0] > psf.shape[1]:
                psf = np.pad(psf, ((0, 0), (0, psf.shape[0] - psf.shape[1])), mode='constant',
                             constant_values=0.)
            return psf

        psf2 = subExp2.getPsf().computeKernelImage(center).getArray()
        psf2 = _makePsfSquare(psf2)

        psf1 = template.getPsf().computeKernelImage(center).getArray()
        psf1 = _makePsfSquare(psf1)

        # from diffimTests.diffimTests ...
        if subExp1.getDimensions()[0] < psf1.shape[0] or subExp1.getDimensions()[1] < psf1.shape[1]:
            return pipeBase.Struct(subExposure=subExposure)

        def _filterPsf(psf):
            """Filter a noisy Psf to remove artifacts. Subject of future research."""
            # only necessary if it's from a measured psf and PsfEx seems to always make PSFs of size 41x41
            if psf.shape[0] == 41:  # its from a measured psf
                psf = psf.copy()
                psf[psf < 0] = 0
                psf[0:10, :] = psf[:, 0:10] = psf[31:41, :] = psf[:, 31:41] = 0
                psf /= psf.sum()

            return psf

        psf1b = psf2b = None
        if self.config.doFilterPsfs:  # default True
            # Note this *really* helps for measured psfs.
            psf1b = _filterPsf(psf1)
            psf2b = _filterPsf(psf2)

        config = ZogyConfig()
        if imageSpace is True:
            config.inImageSpace = imageSpace
            config.padSize = padSize  # Don't need padding if doing all in fourier space
        task = ZogyTask(templateExposure=subExp1, scienceExposure=subExp2,
                        sig1=sig1, sig2=sig2, psf1=psf1b, psf2=psf2b, config=config)

        if not doScorr:
            res = task.computeDiffim(**kwargs)
            D = res.D
        else:
            res = task.computeScorr(**kwargs)
            D = res.S

        outExp = D.Factory(D, subExposure.getBBox())
        out = pipeBase.Struct(subExposure=outExp)
        return out


class ZogyMapReduceConfig(ImageMapReduceConfig):
    """Config to be passed to ImageMapReduceTask

    This config targets the imageMapper to use the ZogyMapper.
    """
    mapper = pexConfig.ConfigurableField(
        doc='Zogy task to run on each sub-image',
        target=ZogyMapper
    )


class ZogyImagePsfMatchConfig(ImagePsfMatchConfig):
    """Config for the ZogyImagePsfMatchTask"""

    zogyConfig = pexConfig.ConfigField(
        dtype=ZogyConfig,
        doc='ZogyTask config to use when running on complete exposure (non spatially-varying)',
    )

    zogyMapReduceConfig = pexConfig.ConfigField(
        dtype=ZogyMapReduceConfig,
        doc='ZogyMapReduce config to use when running Zogy on each sub-image (spatially-varying)',
    )

    def setDefaults(self):
        self.zogyMapReduceConfig.gridStepX = self.zogyMapReduceConfig.gridStepY = 40
        self.zogyMapReduceConfig.cellSizeX = self.zogyMapReduceConfig.cellSizeY = 41
        self.zogyMapReduceConfig.borderSizeX = self.zogyMapReduceConfig.borderSizeY = 8
        self.zogyMapReduceConfig.reducer.reduceOperation = 'average'
        self.zogyConfig.inImageSpace = False


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
                exposure2 after warping to match exposure1
        """

        if spatiallyVarying:
            raise NotImplementedError(
                "DM-25115 Spatially varying zogy subtraction is not implemented.")

        if not self._validateWcs(exposure1, exposure2):
            if doWarping:
                if warpExposure2:
                    self.log.info("Warping exposure2 to exposure1")
                    xyTransform = afwGeom.makeWcsPairTransform(exposure2.getWcs(),
                                                               exposure1.getWcs())
                    psfWarped = measAlg.WarpedPsf(exposure2.getPsf(), xyTransform)
                    exposure2 = self._warper.warpExposure(
                        exposure1.getWcs(), exposure2, destBBox=exposure1.getBBox())
                    exposure2.setPsf(psfWarped)
                else:
                    self.log.info("Warping exposure1 to exposure2")
                    xyTransform = afwGeom.makeWcsPairTransform(exposure1.getWcs(),
                                                               exposure2.getWcs())
                    psfWarped = measAlg.WarpedPsf(exposure1.getPsf(), xyTransform)
                    exposure1 = self._warper.warpExposure(
                        exposure2.getWcs(), exposure1, destBBox=exposure2.getBBox())
                    exposure1.setPsf(psfWarped)
            else:
                self.log.error("ERROR: Input images not registered")
                raise RuntimeError("Input images not registered")

        config = self.config.zogyConfig
        task = ZogyTask(config=config)
        results = task.run(exposure1, exposure2)
        if warpExposure2:
            results.warpedExposure = exposure2
        else:
            results.warpedExposure = exposure1
        return results

    def subtractExposures(self, templateExposure, scienceExposure,
                          doWarping=True, spatiallyVarying=True, inImageSpace=False,
                          doPreConvolve=False):
        """Register, PSF-match, and subtract two Exposures using the ZOGY algorithm.

        Do the following, in order:
        - Warp templateExposure to match scienceExposure, if their WCSs do not already match
        - Compute subtracted exposure ZOGY image subtraction algorithm on the two exposures

        Parameters
        ----------
        templateExposure : `lsst.afw.image.Exposure`
            exposure to PSF-match to scienceExposure. The exposure's mean value is subtracted
            in-place.
        scienceExposure : `lsst.afw.image.Exposure`
            reference Exposure. The exposure's mean value is subtracted in-place.
        doWarping : `bool`
            what to do if templateExposure's and scienceExposure's WCSs do not match:
            - if True then warp templateExposure to match scienceExposure
            - if False then raise an Exception
        spatiallyVarying : `bool`
            If True, perform the operation over a grid of patches across the two exposures
        inImageSpace : `bool`
            If True, perform the Zogy convolutions in image space rather than in frequency space.
        doPreConvolve : `bool`
            ***Currently not implemented.*** If True assume we are to compute the match filter-convolved
            exposure which can be thresholded for detection. In the case of Zogy this would mean
            we compute the Scorr image.

        Returns
        -------
        A `lsst.pipe.base.Struct` containing these fields:
        - subtractedExposure: subtracted Exposure
        - warpedExposure: templateExposure after warping to match scienceExposure (if doWarping true)
        """

        mn1 = self._computeImageMean(templateExposure)
        mn2 = self._computeImageMean(scienceExposure)
        self.log.info("Exposure means=%f, %f; median=%f, %f:" % (mn1[0], mn2[0], mn1[1], mn2[1]))
        if not np.isnan(mn1[0]) and np.abs(mn1[0]) > 1:
            mi = templateExposure.getMaskedImage()
            mi -= mn1[0]
        if not np.isnan(mn2[0]) and np.abs(mn2[0]) > 1:
            mi = scienceExposure.getMaskedImage()
            mi -= mn2[0]

        self.log.info('Running Zogy algorithm: spatiallyVarying=%r' % spatiallyVarying)

        if not self._validateWcs(templateExposure, scienceExposure):
            if doWarping:
                self.log.info("Astrometrically registering template to science image")
                # Also warp the PSF
                xyTransform = afwGeom.makeWcsPairTransform(templateExposure.getWcs(),
                                                           scienceExposure.getWcs())
                psfWarped = measAlg.WarpedPsf(templateExposure.getPsf(), xyTransform)
                templateExposure = self._warper.warpExposure(scienceExposure.getWcs(),
                                                             templateExposure,
                                                             destBBox=scienceExposure.getBBox())

                templateExposure.setPsf(psfWarped)
            else:
                self.log.error("ERROR: Input images not registered")
                raise RuntimeError("Input images not registered")

        def gm(exp):
            return exp.getMaskedImage().getMask()

        def ga(exp):
            return exp.getMaskedImage().getImage().getArray()

        if self.config.zogyConfig.inImageSpace:
            inImageSpace = True  # Override
        self.log.info('Running Zogy algorithm: inImageSpace=%r' % inImageSpace)
        if spatiallyVarying:
            config = self.config.zogyMapReduceConfig
            task = ImageMapReduceTask(config=config)
            results = task.run(scienceExposure, template=templateExposure, inImageSpace=inImageSpace,
                               doScorr=doPreConvolve, forceEvenSized=False)
            results.D = results.exposure
            # The CoaddPsf, when used for detection does not utilize its spatially-varying
            #   properties; it simply computes the PSF at its getAveragePosition().
            # TODO: we need to get it to return the matchedExposure (convolved template)
            #   too, for dipole fitting; but the imageMapReduce task might need to be engineered
            #   for this purpose.
        else:
            config = self.config.zogyConfig
            task = ZogyTask(scienceExposure=scienceExposure, templateExposure=templateExposure,
                            config=config)
            if not doPreConvolve:
                results = task.computeDiffim(inImageSpace=inImageSpace)
                results.matchedExposure = results.R
            else:
                results = task.computeScorr(inImageSpace=inImageSpace)
                results.D = results.S

        # Make sure masks of input images are propagated to diffim
        mask = results.D.getMaskedImage().getMask()
        mask |= scienceExposure.getMaskedImage().getMask()
        mask |= templateExposure.getMaskedImage().getMask()
        results.D.getMaskedImage().getMask()[:, :] = mask
        badBitsNan = mask.addMaskPlane('UNMASKEDNAN')
        resultsArr = results.D.getMaskedImage().getMask().getArray()
        resultsArr[np.isnan(resultsArr)] |= badBitsNan
        resultsArr[np.isnan(scienceExposure.getMaskedImage().getImage().getArray())] |= badBitsNan
        resultsArr[np.isnan(templateExposure.getMaskedImage().getImage().getArray())] |= badBitsNan

        results.subtractedExposure = results.D
        results.warpedExposure = templateExposure
        return results

    def subtractMaskedImages(self, templateExposure, scienceExposure,
                             doWarping=True, spatiallyVarying=True, inImageSpace=False,
                             doPreConvolve=False):
        raise NotImplementedError


subtractAlgorithmRegistry.register('zogy', ZogyImagePsfMatchTask)
