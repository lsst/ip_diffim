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

import numpy as np
from scipy.fftpack import fft2, ifft2, ifftshift
from scipy.ndimage.filters import convolve as scipy_convolve

import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

__all__ = ("DecorrelateALKernelTask", "DecorrelateALKernelConfig")


class DecorrelateALKernelConfig(pexConfig.Config):
    """!
    \anchor DecorrelateALKernelConfig_

    \brief Configuration parameters for the DecorrelateALKernelTask

    Currently there are no parameters for DecorrelateALKernelTask.
    """
    pass

## \addtogroup LSST_task_documentation
## \{
## \page DecorrelateALKernelTask
## \ref DecorrelateALKernelTask_ "DecorrelateALKernelTask"
##      Decorrelate the effect of convolution by Alard-Lupton matching kernel in image difference
## \}

class DecorrelateALKernelTask(pipeBase.Task):
    """!
    \anchor DecorrelateALKernelTask_

    \brief Decorrelate the effect of convolution by Alard-Lupton matching kernel in image difference

    \section pipe_tasks_multiBand_Contents Contents

      - \ref ip_diffim_imageDecorrelation_DecorrelateALKernelTask_Purpose
      - \ref ip_diffim_imageDecorrelation_DecorrelateALKernelTask_Config
      - \ref ip_diffim_imageDecorrelation_DecorrelateALKernelTask_Run
      - \ref ip_diffim_imageDecorrelation_DecorrelateALKernelTask_Debug
      - \ref ip_diffim_imageDecorrelation_DecorrelateALKernelTask_Example

    \section ip_diffim_imageDecorrelation_DecorrelateALKernelTask_Purpose	Description

    Pipe-task that removes the neighboring-pixel covariance in an
    image difference that are added when the template image is
    convolved with the Alard-Lupton PSF matching kernel.

    The image differencing pipeline task \link
    ip.diffim.psfMatch.PsfMatchTask PSFMatchTask\endlink and \link
    ip.diffim.psfMatch.PsfMatchConfigAL PSFMatchConfigAL\endlink uses
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

    \section ip_diffim_imageDecorrelation_DecorrelateALKernelTask_Initialize       Task initialization

    \copydoc \_\_init\_\_

    \section ip_diffim_imageDecorrelation_DecorrelateALKernelTask_Run       Invoking the Task

    \copydoc run

    \section ip_diffim_imageDecorrelation_DecorrelateALKernelTask_Config       Configuration parameters

    This task currently has no relevant configuration parameters.
    See \ref DecorrelateALKernelConfig

    \section ip_diffim_imageDecorrelation_DecorrelateALKernelTask_Debug		Debug variables

    This task has no debug variables

    \section ip_diffim_imageDecorrelation_DecorrelateALKernelTask_Example	Example of using DecorrelateALKernelTask

    This task has no standalone example, however it is applied as a
    subtask of \link pipe.tasks.imageDifference.ImageDifferenceTask ImageDifferenceTask\endlink .

    """
    ConfigClass = DecorrelateALKernelConfig
    _DefaultName = "ip_diffim_decorrelateALKernel"

    def __init__(self, *args, **kwargs):
        """! Create the image decorrelation Task
        @param *args arguments to be passed to lsst.pipe.base.task.Task.__init__
        @param **kwargs keyword arguments to be passed to lsst.pipe.base.task.Task.__init__
        """
        pipeBase.Task.__init__(self, *args, **kwargs)

        self.statsControl = afwMath.StatisticsControl()
        self.statsControl.setNumSigmaClip(3.)
        self.statsControl.setNumIter(3)
        self.statsControl.setAndMask(afwImage.MaskU.getPlaneBitMask(["INTRP", "EDGE",
                                                                     "DETECTED", "BAD",
                                                                     "NO_DATA", "DETECTED_NEGATIVE"]))

    @pipeBase.timeMethod
    def run(self, templateExposure, exposure, subtractedExposure, psfMatchingKernel):
        """! Perform decorrelation of an image difference exposure.

        Decorrelates the diffim due to the convolution of the templateExposure with the
        A&L PSF matching kernel. Currently can accept a spatially varying matching kernel but in
        this case it simply uses a static kernel from the center of the exposure.

        @param[in] templateExposure the template afwImage.Exposure used for PSF matching
        @param[in] exposure the science afwImage.Exposure used for PSF matching
        @param[in] subtractedExposure the subtracted exposure produced by
        `ip_diffim.ImagePsfMatchTask.subtractExposures()`
        @param[in] psfMatchingKernel an (optionally spatially-varying) PSF matching kernel produced
        by `ip_diffim.ImagePsfMatchTask.subtractExposures()`

        @return a `pipeBase.Struct` containing:
            * `correctedExposure`: the decorrelated diffim
            * `correctionKernel`: the decorrelation correction kernel (which may be ignored)

        @note The `subtractedExposure` is NOT updated
        @note The returned `correctedExposure` has an updated PSF as well.
        @note Here we currently convert a spatially-varying matching kernel into a constant kernel,
        just by computing it at the center of the image (tickets DM-6243, DM-6244).
        @note We are also using a constant accross-the-image measure of sigma (sqrt(variance)) to compute
        the decorrelation kernel.
        @note Still TBD (ticket DM-6580): understand whether the convolution is correctly modifying
        the variance plane of the new subtractedExposure.
        """
        self.log.info("Starting.")
        kimg = None
        try:
            if psfMatchingKernel.isSpatiallyVarying():
                spatialKernel = psfMatchingKernel
                kimg = afwImage.ImageD(spatialKernel.getDimensions())
                bbox = subtractedExposure.getBBox()
                xcen = (bbox.getBeginX() + bbox.getEndX()) / 2.
                ycen = (bbox.getBeginY() + bbox.getEndY()) / 2.
                spatialKernel.computeImage(kimg, True, xcen, ycen)
        except:  # Not a spatially-varying kernel (usually not true)
            kimg = psfMatchingKernel.computeImage()

        # Compute the images' sigmas (sqrt of variance)
        statObj = afwMath.makeStatistics(templateExposure.getMaskedImage().getVariance(), afwMath.MEANCLIP,
                                         self.statsControl)
        sig1 = np.sqrt(statObj.getValue(afwMath.MEANCLIP))

        statObj = afwMath.makeStatistics(exposure.getMaskedImage().getVariance(), afwMath.MEANCLIP,
                                         self.statsControl)
        sig2 = np.sqrt(statObj.getValue(afwMath.MEANCLIP))

        corrKernel = DecorrelateALKernelTask._computeDecorrelationKernel(kimg.getArray(), sig1=sig1, sig2=sig2)
        fcorrKernel = DecorrelateALKernelTask._fixEvenKernel(corrKernel)
        self.log.info("Convolving.")
        correctedExposure, corrKern = DecorrelateALKernelTask._doConvolve(subtractedExposure, fcorrKernel)
        self.log.info("Updating correctedExposure and its PSF.")

        # Compute the subtracted exposure's updated psf
        psf = subtractedExposure.getPsf().computeImage().getArray()
        psfc = DecorrelateALKernelTask.computeCorrectedDiffimPsf(fcorrKernel, psf, sig1=sig1, sig2=sig2)
        psfcI = afwImage.ImageD(psfc.shape[0], psfc.shape[1])
        psfcI.getArray()[:, :] = psfc
        psfcK = afwMath.FixedKernel(psfcI)
        psfNew = measAlg.KernelPsf(psfcK)
        correctedExposure.setPsf(psfNew)
        self.log.info("Complete.")

        return pipeBase.Struct(correctedExposure=correctedExposure, correctionKernel=corrKern)

    @staticmethod
    def _computeDecorrelationKernel(kappa, sig1=0.2, sig2=0.2):
        """! Compute the Lupton/ZOGY post-conv. kernel for decorrelating an
        image difference, based on the PSF-matching kernel.
        @param kappa  A matching kernel 2-d numpy.array derived from Alard & Lupton PSF matching
        @param sig1   Average sqrt(variance) of template image used for PSF matching
        @param sig2   Average sqrt(variance) of science image used for PSF matching
        @return a 2-d numpy.array containing the correction kernel

        @note As currently implemented, kappa is a static (single) kernel.
        """
        kft = fft2(kappa)
        kft = np.sqrt((sig1**2 + sig2**2) / (sig1**2 + sig2**2 * kft**2))
        pck = ifft2(kft)
        pck = ifftshift(pck.real)

        # I think we may need to "reverse" the PSF, as in the ZOGY (and Kaiser) papers...
        # This is the same as taking the complex conjugate in Fourier space before FFT-ing back to real space.
        if False:  # TBD: figure this out. For now, we are turning it off.
            pck = pck[::-1, :]

        return pck

    @staticmethod
    def computeCorrectedDiffimPsf(kappa, psf, sig1=0.2, sig2=0.2):
        """! Compute the (decorrelated) difference image's new PSF.
        new_psf = psf(k) * sqrt((sig1**2 + sig2**2) / (sig1**2 + sig2**2 * kappa_ft(k)**2))

        @param kappa  A matching kernel array derived from Alard & Lupton PSF matching
        @param psf    The uncorrected psf array of the science image (and also of the diffim)
        @param sig1   Average sqrt(variance) of template image used for PSF matching
        @param sig2   Average sqrt(variance) of science image used for PSF matching
        @return a 2-d numpy.array containing the new PSF
        """
        def post_conv_psf_ft2(psf, kernel, sig1=1., sig2=1.):
            # Pad psf or kernel symmetrically to make them the same size!
            # Note this assumes they are both square (width == height)
            if psf.shape[0] < kernel.shape[0]:
                diff = (kernel.shape[0] - psf.shape[0]) // 2
                psf = np.pad(psf, (diff, diff), mode='constant')
            elif psf.shape[0] > kernel.shape[0]:
                diff = (psf.shape[0] - kernel.shape[0]) // 2
                kernel = np.pad(kernel, (diff, diff), mode='constant')
            psf_ft = fft2(psf)
            kft = fft2(kernel)
            out = psf_ft * np.sqrt((sig1**2 + sig2**2) / (sig1**2 + sig2**2 * kft**2))
            return out

        def post_conv_psf(psf, kernel, sig1=1., sig2=1.):
            kft = post_conv_psf_ft2(psf, kernel, sig1, sig2)
            out = ifft2(kft)
            return out

        pcf = post_conv_psf(psf=psf, kernel=kappa, sig1=sig2, sig2=sig1)
        pcf = pcf.real / pcf.real.sum()
        return pcf

    @staticmethod
    def _fixEvenKernel(kernel):
        """! Take a kernel with even dimensions and make them odd, centered correctly.
        @param kernel a numpy.array
        @return a fixed kernel numpy.array
        """
        # Make sure the peak (close to a delta-function) is in the center!
        maxloc = np.unravel_index(np.argmax(kernel), kernel.shape)
        out = np.roll(kernel, kernel.shape[0]//2 - maxloc[0], axis=0)
        out = np.roll(out, out.shape[1]//2 - maxloc[1], axis=1)
        # Make sure it is odd-dimensioned by trimming it.
        if (out.shape[0] % 2) == 0:
            maxloc = np.unravel_index(np.argmax(out), out.shape)
            if out.shape[0] - maxloc[0] > maxloc[0]:
                out = out[:-1, :]
            else:
                out = out[1:, :]
            if out.shape[1] - maxloc[1] > maxloc[1]:
                out = out[:, :-1]
            else:
                out = out[:, 1:]
        return out

    @staticmethod
    def _doConvolve(exposure, kernel):
        """! Convolve an Exposure with a decorrelation convolution kernel.
        @param exposure Input afw.image.Exposure to be convolved.
        @param kernel Input 2-d numpy.array to convolve the image with
        @return a new Exposure with the convolved pixels and the (possibly
        re-centered) kernel.

        @note We use afwMath.convolve() but keep scipy.convolve for debugging.
        @note We re-center the kernel if necessary and return the possibly re-centered kernel
        """
        fkernel = DecorrelateALKernelTask._fixEvenKernel(kernel)

        kernelImg = afwImage.ImageD(fkernel.shape[0], fkernel.shape[1])
        kernelImg.getArray()[:, :] = fkernel
        kern = afwMath.FixedKernel(kernelImg)
        maxloc = np.unravel_index(np.argmax(fkernel), fkernel.shape)
        kern.setCtrX(maxloc[0])
        kern.setCtrY(maxloc[1])
        outExp = exposure.clone()  # Do this to keep WCS, PSF, masks, etc.
        convCntrl = afwMath.ConvolutionControl(False, True, 0)
        afwMath.convolve(outExp.getMaskedImage(), exposure.getMaskedImage(), kern, convCntrl)

        return outExp, kern
