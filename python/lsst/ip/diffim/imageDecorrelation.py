from __future__ import absolute_import, division, print_function
#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
from scipy.stats import sigmaclip
from scipy.fftpack import fft2, ifft2, ifftshift
from scipy.ndimage.filters import convolve

import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.afw.math as afwMath

__all__ = ("performExposureDecorrelation")


def computeClippedImageStats(im, low=3, high=3):
    """! Utility function for sigma-clipped array statistics.
    @param im   A numpy array. Should work on any dimensions/shape.
    @param low  Lower bound factor of sigma clipping. Default is 3.
    @param high  Upper bound factor of sigma clipping. Default is 3.
    @return sigma-clipped mean and std of input array
    """
    _, low, upp = sigmaclip(im, low=low, high=high)
    tmp = im[(im > low) & (im < upp)]
    mean1 = np.nanmean(tmp)
    sig1 = np.nanstd(tmp)
    return mean1, sig1


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


def computeDecorrelationCorrectionKernel(kappa, sig1=0.2, sig2=0.2):
    """! Compute the Lupton/ZOGY post-conv. kernel for decorrelating an
    image difference, based on the PSF-matching kernel.
    @param kappa  A matching kernel array derived from Alard & Lupton PSF matching
    @param sig1   Average sqrt(variance) of template image used for PSF matching
    @param sig2   Average sqrt(variance) of science image used for PSF matching
    @return a 2-d numpy.array containing the correction kernel

    @note As currently implemented, kappa is a static (single) kernel.
    """
    def post_conv_kernel_ft2(kernel, sig1=1., sig2=1.):
        kft = fft2(kernel)
        return np.sqrt((sig1**2 + sig2**2) / (sig1**2 + sig2**2 * kft**2))

    def post_conv_kernel2(kernel, sig1=1., sig2=1.):
        kft = post_conv_kernel_ft2(kernel, sig1, sig2)
        out = ifft2(kft)
        return out

    pck = post_conv_kernel2(kappa, sig1=sig2, sig2=sig1)
    pck = ifftshift(pck.real)

    # I think we may need to "reverse" the PSF, as in the ZOGY (and Kaiser) papers...
    # This is the same as taking the complex conjugate in Fourier space before FFT-ing back to real space.
    if False:  # TBD: figure this out. For now, we are turning it off.
        pck = pck[::-1, :]

    return pck


def computeCorrectedDiffimPsf(kappa, im2_psf, sig1=0.2, sig2=0.2):
    """! Compute the (decorrelated) difference image's new PSF:
    post_conv_psf = phi_1(k) * sym.sqrt((sig1**2 + sig2**2) / (sig1**2 + sig2**2 * kappa_ft(k)**2))

    @param kappa  A matching kernel array derived from Alard & Lupton PSF matching
    @param im2_psf The uncorrected psf array of the science image (and also of the diffim)
    @param sig1   Average sqrt(variance) of template image used for PSF matching
    @param sig2   Average sqrt(variance) of science image used for PSF matching
    @return a 2-d numpy.array containing the new PSF
    """
    def post_conv_psf_ft2(psf, kernel, sig1=1., sig2=1.):
        # Pad psf or kernel symmetrically to make them the same size!
        if psf.shape[0] < kernel.shape[0]:
            while psf.shape[0] < kernel.shape[0]:
                psf = np.pad(psf, (1, 1), mode='constant')
        elif psf.shape[0] > kernel.shape[0]:
            while psf.shape[0] > kernel.shape[0]:
                kernel = np.pad(kernel, (1, 1), mode='constant')
        psf_ft = fft2(psf)
        kft = fft2(kernel)
        out = psf_ft * np.sqrt((sig1**2 + sig2**2) / (sig1**2 + sig2**2 * kft**2))
        return out

    def post_conv_psf(psf, kernel, sig1=1., sig2=1.):
        kft = post_conv_psf_ft2(psf, kernel, sig1, sig2)
        out = ifft2(kft)
        return out

    im2_psf_small = im2_psf
    # First compute the science image's (im2's) psf, subset on -16:15 coords
    if im2_psf.shape[0] > 50:
        x0im, y0im = getImageGrid(im2_psf)
        x = np.arange(-16, 16, 1)
        y = x.copy()
        x0, y0 = np.meshgrid(x, y)
        im2_psf_small = im2_psf[(x0im.max()+x.min()+1):(x0im.max()-x.min()+1),
                                (y0im.max()+y.min()+1):(y0im.max()-y.min()+1)]
    pcf = post_conv_psf(psf=im2_psf_small, kernel=kappa, sig1=sig2, sig2=sig1)
    pcf = pcf.real / pcf.real.sum()
    return pcf


def performExposureDecorrelation(templateExposure, exposure, subtractedExposure, psfMatchingKernel, log):
    """! Compute decorrelation correction on an image difference exposure, given the
    input template/science exposures and the (potentially spatially-varying) PSF matching kernel.

    @param templateExposure the template afwImage.Exposure used for PSF matching
    @param exposure the science afwImage.Exposure used for PSF matching
    @param subtractedExposure the resulting subtracted exposure produced by 
    ip_diffim.ImagePsfMatchTask.subtractExposures()
    @param psfMatchingKernel an (optionally spatially-varying) PSF matching kernel produced 
    by ip_diffim.ImagePsfMatchTask.subtractExposures()
    @param log an input pexLog for outputting messages

    @return the corrected subtractedExposure (updated in-place) with an updated PSF, and
    the decorrelation correction kernel (which may be ignored)

    @note Here we currently convert a spatially-varying matching kernel into a constant kernel, 
    just compute it for the center of the image.
    @note Still TBD: do we need to modify the new subtractedExposure's variance plane?
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
    sig1 = templateExposure.getMaskedImage().getVariance().getArray()
    sig1squared, _ = computeClippedImageStats(sig1)
    sig1 = np.sqrt(sig1squared)

    sig2 = exposure.getMaskedImage().getVariance().getArray()
    sig2squared, _ = computeClippedImageStats(sig2)
    sig2 = np.sqrt(sig2squared)

    corrKernel = computeDecorrelationCorrectionKernel(kimg.getArray(), sig1=sig1, sig2=sig2)
    # Eventually, use afwMath.convolve(), but for now we just use scipy.
    log.info("Decorrelation: Convolving.")
    pci = convolve(subtractedExposure.getMaskedImage().getImage().getArray(),
                   corrKernel, mode='constant')
    subtractedExposure.getMaskedImage().getImage().getArray()[:, :] = pci
    log.info("Decorrelation: Finished with convolution.")

    # Compute the subtracted exposure's updated psf
    psf = subtractedExposure.getPsf().computeImage().getArray()
    psfc = computeCorrectedDiffimPsf(corrKernel, psf, sig1=sig1, sig2=sig2)
    psfcI = afwImage.ImageD(subtractedExposure.getPsf().computeImage().getBBox())
    psfcI.getArray()[:, :] = psfc
    psfcK = afwMath.FixedKernel(psfcI)
    psfNew = measAlg.KernelPsf(psfcK)
    subtractedExposure.setPsf(psfNew)
    return subtractedExposure, corrKernel
