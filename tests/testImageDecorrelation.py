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

import unittest

import numpy as np

import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.pex.logging as pexLog

import lsst.ip.diffim.imageDecorrelation as ipDiffim_id

def singleGaussian2d(x, y, xc, yc, sigma_x=1., sigma_y=1., theta=0.):
    """! Generate a 2-d Gaussian, possibly elongated and rotated, on a grid of pixel
    coordinates given by x,y.
    @param x,y each a 1-d numpy.array containing x- and y- coordinates for independent variables,
    for example `np.arange(-16, 15)`.
    @param xc,yc each a float giving the centroid of the gaussian
    @param sigma_x,sigma_y each a float giving the sigma of the gaussian
    @param theta a float giving the rotation of the gaussian (degrees)
    @return a 2-d numpy.array containing the normalized 2-d Gaussian
    """
    theta = (theta/180.) * np.pi
    cos_theta2, sin_theta2 = np.cos(theta)**2., np.sin(theta)**2.
    sigma_x2, sigma_y2 = sigma_x**2., sigma_y**2.
    a = cos_theta2/(2.*sigma_x2) + sin_theta2/(2.*sigma_y2)
    b = -(np.sin(2.*theta))/(4.*sigma_x2) + (np.sin(2.*theta))/(4.*sigma_y2)
    c = sin_theta2/(2.*sigma_x2) + cos_theta2/(2.*sigma_y2)
    xxc, yyc = x-xc, y-yc
    out = np.exp(-(a*(xxc**2.) + 2.*b*xxc*yyc + c*(yyc**2.)))
    out /= out.sum()
    return out


def makeFakeImages(xim=None, yim=None, sig1=0.2, sig2=0.2, psf1=2.2, psf2=3.3, offset=None,
                   psf_yvary_factor=0., varSourceChange=1/50., theta1=0., theta2=0., im2background=0.,
                   n_sources=500, seed=66, verbose=False):
    """! Make two exposures: a template and a science exposure.
    Add random sources of identical flux, with randomly-distributed fluxes and a given PSF, then add noise.
    @param xim,yim image pixel coordinates on which to generate the image grid. Default is (-256:256).
    @param sig1,sig2 std. dev. of noise to be generated on input images. Defalt is 0.2 for both.
    @param psf1,psf2 std. dev. of (Gaussian) PSFs for the two images in x,y direction. Default is
    [2.2, 2.2] and [3.3, 3.3] for im1 and im2 respectively.
    @param offset add a constant (pixel) astrometric offset between the two images
    @param psf_yvary_factor vary the PSF of the science image by this much across the image
    @param varSourceChange add this amount of fractional flux to a single source closest to
    the center of the science image
    @param im2background add a constant value to the science image
    @param n_sources the number of sources to add to the images
    @param seed the numpy random seed to set prior to image generation
    @param verbose be verbose

    @return im1, im2: the template and science afwImage.Exposures

    @note having sources near the edges really messes up the
    fitting (probably because of the convolution). So we make sure no
    sources are near the edge.
    @note also it seems that having the variable source with a large
    flux increase also messes up the fitting (seems to lead to
    overfitting -- perhaps to the source itself). This might be fixed by
    adding more constant sources.
    """
    np.random.seed(seed)

    psf1 = [2.2, 2.2] if psf1 is None else psf1
    if not hasattr(psf1, "__len__"):
        psf1 = [psf1, psf1]
    psf2 = [3.3, 3.3] if psf2 is None else psf2
    if not hasattr(psf2, "__len__"):
        psf2 = [psf2, psf2]
    offset = [0., 0.] if offset is None else offset   # astrometric offset (pixels) between the two images
    if verbose:
        print('Template PSF:', psf1, theta1)
        print('Science PSF:', psf2, theta2)
        print(np.sqrt(psf2[0]**2 - psf1[0]**2))
        print('Offset:', offset)

    xim = np.arange(-256, 256, 1) if xim is None else xim
    yim = xim.copy() if yim is None else yim
    x0im, y0im = np.meshgrid(xim, yim)
    fluxes = np.random.uniform(50, 30000, n_sources)
    xposns = np.random.uniform(xim.min()+16, xim.max()-5, n_sources)
    yposns = np.random.uniform(yim.min()+16, yim.max()-5, n_sources)

    # Make the source closest to the center of the image the one that increases in flux
    ind = np.argmin(xposns**2. + yposns**2.)

    im1 = np.random.normal(scale=sig1, size=x0im.shape)  # sigma of template
    im2 = np.random.normal(scale=sig2, size=x0im.shape)  # sigma of science image

    # variation in y-width of psf in science image across (x-dim of) image:
    psf2_yvary = psf_yvary_factor * (yim.mean() - yposns) / yim.max()
    if verbose:
        print('PSF y spatial-variation:', psf2_yvary.min(), psf2_yvary.max())

    for i in range(n_sources):
        flux = fluxes[i]
        tmp1 = flux * singleGaussian2d(x0im, y0im, xposns[i], yposns[i], psf1[0], psf1[1], theta=theta1)
        im1 += tmp1
        if i == ind:
            flux += flux * varSourceChange
        tmp2 = flux * singleGaussian2d(x0im, y0im, xposns[i]+offset[0], yposns[i]+offset[1],
                                       psf2[0], psf2[1]+psf2_yvary[i], theta=theta2)
        im2 += tmp2

    # Add a (constant, for now) background offset to im2
    if im2background != 0.:  # im2background = 10.
        if verbose:
            print('Background:', im2background)
        im2 += im2background

    im1_psf = singleGaussian2d(x0im, y0im, 0, 0, psf1[0], psf1[1], theta=theta1)
    im2_psf = singleGaussian2d(x0im, y0im, offset[0], offset[1], psf2[0], psf2[1], theta=theta2)

    def makeExposure(imgArray, psfArray, imgSigma):
        """! Convert an image numpy.array and corresponding PSF numpy.array into an exposure.

        Add the (constant) variance plane equal to imgSigma**2.

        @param imgArray 2-d numpy.array containing the image
        @param psfArray 2-d numpy.array containing the PSF image
        @param imgSigma std. deviation of input image (equal to sqrt(variance))
        @return a new exposure containing the image, PSF and desired variance plane
        """
        # All this code to convert the template image array/psf array into an exposure.
        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Point2I(imgArray.shape[0]-1, imgArray.shape[1]-1))
        im1ex = afwImage.ExposureD(bbox)
        im1ex.getMaskedImage().getImage().getArray()[:, :] = imgArray
        im1ex.getMaskedImage().getVariance().getArray()[:, :] = imgSigma**2.
        psfBox = afwGeom.Box2I(afwGeom.Point2I(-20, -20), afwGeom.Point2I(20, 20))  # a 41x41 pixel psf
        psf = afwImage.ImageD(psfBox)
        psfBox.shift(afwGeom.Extent2I(256, 256))
        im1_psf_sub = psfArray[psfBox.getMinX():psfBox.getMaxX()+1, psfBox.getMinY():psfBox.getMaxY()+1]
        psf.getArray()[:, :] = im1_psf_sub
        psfK = afwMath.FixedKernel(psf)
        psfNew = measAlg.KernelPsf(psfK)
        im1ex.setPsf(psfNew)
        return im1ex

    im1ex = makeExposure(im1, im1_psf, sig1)
    im2ex = makeExposure(im2, im2_psf, sig2)

    return im1ex, im2ex


class DiffimCorrectionTest(lsst.utils.tests.TestCase):
    """!A test case for the diffim image decorrelation algorithm.
    """

    def setUp(self):
        """!Generate a fake aligned template and science image and analyse the noise in the
        resulting diffim. First, use a non-spatially-varying psf.
        """

        self.psf1_sigma = 2.2  # sigma of psf of template image
        self.psf2_sigma = 3.3  # sigma of psf of science image
        self.sig1 = 0.2  # std.dev of noise in template image
        self.sig2 = 0.2  # std.dev of noise in science image

        self.im1ex, self.im2ex \
            = makeFakeImages(sig1=self.sig1, sig2=self.sig2, psf1=self.psf1_sigma, psf2=self.psf2_sigma,
                             n_sources=50, verbose=True)

    def tearDown(self):
        del self.im1ex
        del self.im2ex

    def testDiffimCorrection(self):
        """! Check that the variance of the corrected diffim matches the theoretical value
        (to within a 1% tolerance).
        """

        # Create the matching kernel. We used Gaussian PSFs for im1 and im2, so we can compute the "expected"
        # matching kernel sigma.
        psf1_sig = self.im1ex.getPsf().computeShape().getDeterminantRadius()
        psf2_sig = self.im2ex.getPsf().computeShape().getDeterminantRadius()
        sig_match = np.sqrt((psf2_sig**2. - psf1_sig**2.))
        self.assertClose(sig_match, np.sqrt((self.psf2_sigma**2. - self.psf1_sigma**2.)), rtol=1e-5)
        #matchingKernel = measAlg.SingleGaussianPsf(31, 31, sig_match)
        x0 = np.arange(-16, 16, 1)
        y0 = x0.copy()
        x0im, y0im = np.meshgrid(x0, y0)
        matchingKernel = singleGaussian2d(x0im, y0im, 0., 0., sigma_x=sig_match, sigma_y=sig_match)
        kernelImg = afwImage.ImageD(matchingKernel.shape[0], matchingKernel.shape[1])
        kernelImg.getArray()[:, :] = matchingKernel
        mKernel = afwMath.FixedKernel(kernelImg)
        mKernel = measAlg.KernelPsf(mKernel)

        # Create the matched template by convolving the template with the matchingKernel
        from scipy.ndimage.filters import convolve
        im1 = self.im1ex.getMaskedImage().getImage().getArray()
        matched_im1 = convolve(im1, matchingKernel, mode='constant')
        matched_im1ex = self.im1ex.clone()
        matched_im1ex.getMaskedImage().getImage().getArray()[:, :] = matched_im1

        # Expected (ideal) variance of difference image
        expected_var = self.sig1**2 + self.sig2**2
        print('Expected variance:', expected_var)

        im2 = self.im2ex.getMaskedImage().getImage().getArray()
        print(np.nan_to_num(matched_im1 - im2).var())
        # Uncorrected diffim - variance is wrong (too low)
        self.assertNotClose(np.nan_to_num(matched_im1 - im2).var(), expected_var, rtol=0.1)
        # In fact, it should have an (incorrect) variance that is close to the variance of the science image.
        self.assertClose((matched_im1 - im2)[~np.isnan(matched_im1 - im2)].var(), self.sig2**2., rtol=0.1)

        diffExp = matched_im1ex.clone()
        tmpArr = diffExp.getMaskedImage().getImage().getArray()
        tmpArr -= self.im2ex.getMaskedImage().getImage().getArray()
        # Uncorrected diffim exposure - variance is wrong (too low)
        self.assertNotClose(tmpArr[~np.isnan(tmpArr)].var(), expected_var, rtol=0.1)

        log = pexLog.Log(pexLog.Log.getDefaultLog(), 'testImageDecorrelation', pexLog.Log.INFO)
        # This corrects the diffExp in-place:
        corrected_diffExp, corrKernel = ipDiffim_id.decorrelateExposure(self.im1ex, self.im2ex,
                                                                        diffExp, mKernel, log)

        # Corrected diffim - variance should be close to expected.
        corrected_diffArr = corrected_diffExp.getMaskedImage().getImage().getArray()
        corrected_diffArr = corrected_diffArr[15:-15, 15:-15]
        var = corrected_diffArr[~np.isnan(corrected_diffArr)].var()
        print(var)
        _, _, var = ipDiffim_id.computeClippedImageStats(corrected_diffExp.getMaskedImage())
        print(var)  # why is this not ignoring the "edge" mask?
        self.assertClose(var, expected_var, rtol=0.02)


def suite():
    """!Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(DiffimCorrectionTest)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """!Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)

