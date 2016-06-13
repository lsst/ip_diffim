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

__all__ = ("decorrelateExposure", "computeDecorrelationKernel", "doConvolve",
           "computeClippedImageStats")


def computeClippedImageStats(im):
    """! Utility function for sigma-clipped array statistics on an image or exposure.
    @param im An afw.Exposure, masked image, or image.
    @return sigma-clipped mean, std, and variance of input array
    """
    statsControl = afwMath.StatisticsControl()
    statsControl.setNumSigmaClip(3.)
    statsControl.setNumIter(3)
    statsControl.setAndMask(afwImage.MaskU.getPlaneBitMask(["INTRP", "EDGE",
                                                            "DETECTED", "BAD",
                                                            "NO_DATA", "DETECTED_NEGATIVE"]))
    statObj = afwMath.makeStatistics(im, afwMath.MEANCLIP | afwMath.STDEVCLIP | afwMath.VARIANCECLIP,
                                     statsControl)
    mean = statObj.getValue(afwMath.MEANCLIP)
    std = statObj.getValue(afwMath.STDEVCLIP)
    var = statObj.getValue(afwMath.VARIANCECLIP)
    return mean, std, var


def getImageGrid(im):
    """! Convenience function for creating a grid of pixel coordinates covering a 2-d array,
    using `numpy.meshgrid()`
    @param im input 2-d numpy array
    @return two numpy arraysm each with same shape as `im`, containing x- and y- coordinate
    grid.
    """
    xim = np.arange(np.int(-np.floor(im.shape[0]/2.)), np.int(np.floor(im.shape[0]/2)))
    yim = np.arange(np.int(-np.floor(im.shape[1]/2.)), np.int(np.floor(im.shape[1]/2)))
    x0im, y0im = np.meshgrid(xim, yim)
    return x0im, y0im


def computeDecorrelationKernel(kappa, sig1=0.2, sig2=0.2):
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


def computeCorrectedDiffimPsf(kappa, im2_psf, sig1=0.2, sig2=0.2):
    """! Compute the (decorrelated) difference image's new PSF.
    new_psf = phi_1(k) * sqrt((sig1**2 + sig2**2) / (sig1**2 + sig2**2 * kappa_ft(k)**2))

    @param kappa  A matching kernel array derived from Alard & Lupton PSF matching
    @param im2_psf The uncorrected psf array of the science image (and also of the diffim)
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

    pcf = post_conv_psf(psf=im2_psf, kernel=kappa, sig1=sig2, sig2=sig1)
    pcf = pcf.real / pcf.real.sum()
    return pcf


def fixEvenKernel(kernel):
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


def doConvolve(exposure, kernel, use_scipy=False):
    """! Convolve an Exposure with a convolution kernel.
    @param exposure Input afw.image.Exposure to be convolved.
    @param kernel Input 2-d numpy.array to convolve the image with
    @param use_scipy Use scipy to do convolution instead of afwMath
    @return a new Exposure with the convolved pixels.

    @note we use afwMath.convolve() but keep scipy.convolve for debugging.
    """
    outExp = kern = None
    fkernel = fixEvenKernel(kernel)
    if use_scipy:
        pci = scipy_convolve(exposure.getMaskedImage().getImage().getArray(),
                             fkernel, mode='constant', cval=np.nan)
        outExp = exposure.clone()
        outExp.getMaskedImage().getImage().getArray()[:, :] = pci
        kern = fkernel

    else:
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


def decorrelateExposure(templateExposure, exposure, subtractedExposure, psfMatchingKernel, log):
    """! Compute decorrelation correction on an image difference exposure.

    Currently can accept a spatially varying kernel but in this case it simply uses a static
    kernel from the center of the exposure.

    @param templateExposure[in] the template afwImage.Exposure used for PSF matching
    @param exposure[in] the science afwImage.Exposure used for PSF matching
    @param subtractedExposure[in] the subtracted exposure produced by 
    ip_diffim.ImagePsfMatchTask.subtractExposures()
    @param psfMatchingKernel an (optionally spatially-varying) PSF matching kernel produced 
    by ip_diffim.ImagePsfMatchTask.subtractExposures()
    @param log an input pexLog for outputting messages

    @return the the decorrelated diffim and decorrelation correction kernel (which may be ignored)
    The returned decorrelated diffim has an updated PSF as well.

    @note the subtractedExposure is NOT updated
    @note Here we currently convert a spatially-varying matching kernel into a constant kernel, 
    just by computing it at the center of the image (tickets DM-6243, DM-6244).
    @note We are also using a constant accross-the-image measure of sigma (sqrt(variance)) to compute
    the decorrelation kernel.
    @note Still TBD (ticket DM-6580): understand whether the convolution is correctly modifying 
    the variance plane of the new subtractedExposure.
    """
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
    _, _, sig1squared = computeClippedImageStats(templateExposure.getMaskedImage())
    sig1 = np.sqrt(sig1squared)

    _, _, sig2squared = computeClippedImageStats(exposure.getMaskedImage())
    sig2 = np.sqrt(sig2squared)

    corrKernel = computeDecorrelationKernel(kimg.getArray(), sig1=sig1, sig2=sig2)
    fcorrKernel = fixEvenKernel(corrKernel)
    log.info("Decorrelation: Convolving.")
    correctedExposure, corrKern = doConvolve(subtractedExposure, fcorrKernel)
    log.info("Decorrelation: Finished with convolution.")

    # Compute the subtracted exposure's updated psf
    psf = subtractedExposure.getPsf().computeImage().getArray()
    psfc = computeCorrectedDiffimPsf(fcorrKernel, psf, sig1=sig1, sig2=sig2)
    psfcI = afwImage.ImageD(psfc.shape[0], psfc.shape[1])
    psfcI.getArray()[:, :] = psfc
    psfcK = afwMath.FixedKernel(psfcI)
    psfNew = measAlg.KernelPsf(psfcK)
    correctedExposure.setPsf(psfNew)

    return correctedExposure, corrKern
