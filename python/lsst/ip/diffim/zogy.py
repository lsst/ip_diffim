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

    ignoreMaskPlanes = pexConfig.ListField(
        dtype=str,
        default=("INTRP", "EDGE", "DETECTED", "SAT", "CR", "BAD", "NO_DATA", "DETECTED_NEGATIVE"),
        doc="Mask planes to ignore for statistics"
    )


MIN_KERNEL = 1.0e-4


class ZogyTask(pipeBase.Task):
    """Task to perform ZOGY proper image subtraction. See module-level documentation for
    additional details.

    In all methods, im1 is R (reference, or template) and im2 is N (new, or science).
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

        self.Fr = self.config.templateFluxScaling  # default is 1
        self.Fn = self.config.scienceFluxScaling  # default is 1
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
        offset = [size[0]//2 - psf.shape[0]//2 - 1, size[1]//2 - psf.shape[1]//2 - 1]
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
            kern.setCtrX(maxloc[0])
            kern.setCtrY(maxloc[1])
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

    def _computeImageMean(self, exposure):
        """Compute the sigma-clipped mean of the pixels image of `exposure`.
        """
        statsControl = afwMath.StatisticsControl()
        statsControl.setNumSigmaClip(3.)
        statsControl.setNumIter(3)
        ignoreMaskPlanes = ("INTRP", "EDGE", "DETECTED", "SAT", "CR", "BAD", "NO_DATA", "DETECTED_NEGATIVE")
        statsControl.setAndMask(afwImage.Mask.getPlaneBitMask(ignoreMaskPlanes))
        statObj = afwMath.makeStatistics(exposure.getMaskedImage().getImage(),
                                         exposure.getMaskedImage().getMask(),
                                         afwMath.MEANCLIP | afwMath.MEDIAN, statsControl)
        mn = statObj.getValue(afwMath.MEANCLIP)
        med = statObj.getValue(afwMath.MEDIAN)
        return mn, med

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
