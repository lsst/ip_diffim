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

import unittest
import numpy as np

import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
# import lsst.ip.diffim as ipDiffim
# import lsst.daf.base as dafBase
import lsst.pex.logging as pexLog

from lsst.ip.diffim.imageDecorrelation import performExposureDecorrelation as doDecorr

def singleGaussian2d(x, y, xc, yc, sigma_x=1., sigma_y=1., theta=0.):
    """! Generate a 2-d Gaussian, possibly elongated and rotated, on a grid of pixel
    coordinates given by x,y.
    @param x,y each a 1-d numpy.array containing x- and y- coordinates for independent variables.
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
    """! Make the two "images". im1 is the template, im2 is the science image.
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

    @return im1, im2, im1_psf, im2_psf: the template and science images (as 2-d numpy.arrays) and
    their corresponding PSFs (also as 2-d numpy.arrays)

    @note having sources near the edges really messes up the
    fitting (probably because of the convolution). So make sure no
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
    return im1, im2, im1_psf, im2_psf


class DiffimCorrectionTest(lsst.utils.tests.TestCase):
    """!A test case for the diffim image decorrelation algorithm.
    """

    def setUp(self):
        """!Generate a fake aligned template and science image and analyse the noise in the
        resulting diffim. First, use a non-spatially-varying psf.
        """

        self.psf1 = 2.2  # sigma of psf of template image
        self.psf2 = 3.3  # sigma of psf of science image
        self.sig1 = 0.2  # std.dev of noise in template image
        self.sig2 = 0.2  # std.dev of noise in science image

        self.im1, self.im2, self.im1_psf, self.im2_psf \
            = makeFakeImages(sig1=self.sig1, sig2=self.sig2, psf1=self.psf1, psf2=self.psf2,
                             n_sources=50, verbose=True)
        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Point2I(self.im1.shape[0]-1, self.im1.shape[1]-1))

        # All this code to convert the template image array/psf array into an exposure.
        self.im1ex = afwImage.ExposureF(bbox)
        self.im1ex.getMaskedImage().getImage().getArray()[:, :] = self.im1
        self.im1ex.getMaskedImage().getVariance().getArray()[:, :] = self.sig1**2.
        psfBox = afwGeom.Box2I(afwGeom.Point2I(-20, -20), afwGeom.Point2I(20, 20))
        psf = afwImage.ImageD(psfBox)
        psfBox.shift(afwGeom.Extent2I(256, 256))
        im1_psf_sub = self.im1_psf[psfBox.getMinX():psfBox.getMaxX()+1, psfBox.getMinY():psfBox.getMaxY()+1]
        psf.getArray()[:, :] = im1_psf_sub
        psfK = afwMath.FixedKernel(psf)
        psfNew = measAlg.KernelPsf(psfK)
        self.im1ex.setPsf(psfNew)

        # Now for the 2nd (science) image
        self.im2ex = afwImage.ExposureF(bbox)
        self.im2ex.getMaskedImage().getImage().getArray()[:, :] = self.im2
        self.im2ex.getMaskedImage().getVariance().getArray()[:, :] = self.sig2**2.
        psf = afwImage.ImageD(psfBox)
        im2_psf_sub = self.im2_psf[psfBox.getMinX():psfBox.getMaxX()+1, psfBox.getMinY():psfBox.getMaxY()+1]
        psf.getArray()[:, :] = im2_psf_sub
        psfK = afwMath.FixedKernel(psf)
        psfNew = measAlg.KernelPsf(psfK)
        self.im2ex.setPsf(psfNew)

    def tearDown(self):
        del self.im1
        del self.im1_psf
        del self.im2
        del self.im2_psf
        del self.im1ex
        del self.im2ex

    def testDiffimCorrection(self):
        """! Check that the variance of the corrected diffim matches the theoretical value
        (to within a 1% tolerance).
        """

        # Create the matching kernel. We use gaussians, so it's just:
        sig_match = np.sqrt((self.psf2**2. - self.psf1**2.))
        matchingKernel = measAlg.SingleGaussianPsf(25, 25, sig_match)

        # Create the matched template by convolving the template with the matchingKernel
        from scipy.ndimage.filters import convolve
        matched_im1 = convolve(self.im1, matchingKernel.computeImage().getArray(), mode='constant')
        matched_im1ex = self.im1ex.clone()
        matched_im1ex.getMaskedImage().getImage().getArray()[:, :] = matched_im1

        # Expected (ideal) variance of difference image
        expected_var = self.sig1**2 + self.sig2**2
        print('Expected variance:', expected_var)

        #print(self.im1.var())
        #print(self.im2.var())
        #print((self.im1 - self.im2).var())
        print((matched_im1 - self.im2).var())
        # Uncorrected diffim - variance is wrong (too low)
        self.assertNotClose((matched_im1 - self.im2).var(), expected_var, rtol=0.1)
        # In fact, it should have a variance close to the variance of the science image.
        self.assertClose((matched_im1 - self.im2).var(), self.sig2**2., rtol=0.1)

        diffExp = matched_im1ex.clone()
        tmpArr = diffExp.getMaskedImage().getImage().getArray()
        tmpArr -= self.im2ex.getMaskedImage().getImage().getArray()
        #print(diffExp.getMaskedImage().getImage().getArray().var())
        # Uncorrected diffim exposure - variance is wrong (too low)
        self.assertNotClose(diffExp.getMaskedImage().getImage().getArray().var(), expected_var, rtol=0.1)

        log = pexLog.Log(pexLog.Log.getDefaultLog(), 'testImageDecorrelation', pexLog.Log.INFO)
        corrected_diffExp, _ = doDecorr(self.im1ex, self.im2ex, diffExp, matchingKernel, log)

        corrected_diffArr = corrected_diffExp.getMaskedImage().getImage().getArray()
        print(corrected_diffArr.var())
        # Corrected diffim - variance should be close to expected.
        self.assertClose(corrected_diffArr.var(), expected_var, rtol=0.01)


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

