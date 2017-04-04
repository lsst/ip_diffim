#
# LSST Data Management System
# Copyright 2016-2017 AURA/LSST.
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
from __future__ import absolute_import, division, print_function
import unittest

from builtins import range
from past.builtins import basestring
import numpy as np

import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.daf.base as dafBase

from lsst.ip.diffim.imageDecorrelation import (DecorrelateALKernelTask,
                                               DecorrelateALKernelMapReduceConfig)
from lsst.ip.diffim.imageMapReduce import ImageMapReduceTask

try:
    type(verbose)
except NameError:
    verbose = False


def setup_module(module):
    lsst.utils.tests.init()


def singleGaussian2d(x, y, xc, yc, sigma_x=1., sigma_y=1., theta=0., ampl=1.):
    """! Generate a 2-d Gaussian, possibly elongated and rotated, on a grid of pixel
    coordinates given by x,y.
    @param x,y each a 1-d numpy.array containing x- and y- coordinates for independent variables,
    for example `np.arange(-16, 15)`.
    @param xc,yc each a float giving the centroid of the gaussian
    @param sigma_x,sigma_y each a float giving the sigma of the gaussian
    @param theta a float giving the rotation of the gaussian (degrees)
    @param ampl a float giving the amplitude of the gaussian
    @return a 2-d numpy.array containing the normalized 2-d Gaussian

    @Note this can be done in `astropy.modeling` but for now we have it explicitly here.
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


def makeFakeImages(xim=None, yim=None, svar=0.04, tvar=0.04, psf1=3.3, psf2=2.2, offset=None,
                   psf_yvary_factor=0., varSourceChange=1/50., theta1=0., theta2=0.,
                   n_sources=500, seed=66, verbose=False):
    """! Make two exposures: a template and a science exposure.
    Add random sources of identical flux, with randomly-distributed fluxes and a given PSF, then add noise.
    In all cases below, index (1) is the science image, and (2) is the template.
    @param xim,yim image pixel coordinates on which to generate the image grid. Default is (-256:256).
    @param svar,tar variance of noise to be generated on science/template images. Default is 0.04 for both.
    @param psf1,psf2 std. dev. of (Gaussian) PSFs for the two images in x,y direction. Default is
    [3.3, 3.3] and [2.2, 2.2] for im1 and im2 respectively.
    @param offset add a constant (pixel) astrometric offset between the two images
    @param psf_yvary_factor vary the y-width of the PSF across the x-axis of the science image (zero,
    the default, means no variation)
    @param varSourceChange add this amount of fractional flux to a single source closest to
    the center of the science image
    @param n_sources the number of sources to add to the images
    @param seed the numpy random seed to set prior to image generation
    @param verbose be verbose

    @return im1, im2: the science and template afwImage.Exposures

    @note having sources near the edges really messes up the
    fitting (probably because of the convolution). So we make sure no
    sources are near the edge.
    @note also it seems that having the variable source with a large
    flux increase also messes up the fitting (seems to lead to
    overfitting -- perhaps to the source itself). This might be fixed by
    adding more constant sources.
    """
    np.random.seed(seed)

    psf1 = [3.3, 3.3] if psf1 is None else psf1
    if not hasattr(psf1, "__len__") and not isinstance(psf1, basestring):
        psf1 = [psf1, psf1]
    psf2 = [2.2, 2.2] if psf2 is None else psf2
    if not hasattr(psf2, "__len__") and not isinstance(psf2, basestring):
        psf2 = [psf2, psf2]
    offset = [0., 0.] if offset is None else offset   # astrometric offset (pixels) between the two images
    if verbose:
        print('Science PSF:', psf1, theta1)
        print('Template PSF:', psf2, theta2)
        print(np.sqrt(psf1[0]**2 - psf2[0]**2))
        print('Offset:', offset)

    xim = np.arange(-256, 256, 1) if xim is None else xim
    yim = xim.copy() if yim is None else yim
    x0im, y0im = np.meshgrid(xim, yim)
    fluxes = np.random.uniform(50, 30000, n_sources)
    xposns = np.random.uniform(xim.min()+16, xim.max()-5, n_sources)
    yposns = np.random.uniform(yim.min()+16, yim.max()-5, n_sources)

    # Make the source closest to the center of the image the one that increases in flux
    ind = np.argmin(xposns**2. + yposns**2.)

    im1 = np.random.normal(scale=np.sqrt(svar), size=x0im.shape)  # variance of science image
    im2 = np.random.normal(scale=np.sqrt(tvar), size=x0im.shape)  # variance of template

    # vary the y-width of psf across x-axis of science image (zero means no variation):
    psf1_yvary = psf_yvary_factor * (yim.mean() - yposns) / yim.max()
    if verbose:
        print('PSF y spatial-variation:', psf1_yvary.min(), psf1_yvary.max())

    for i in range(n_sources):
        flux = fluxes[i]
        tmp = flux * singleGaussian2d(x0im, y0im, xposns[i], yposns[i], psf2[0], psf2[1], theta=theta2)
        im2 += tmp
        if i == ind:
            flux += flux * varSourceChange
        tmp = flux * singleGaussian2d(x0im, y0im, xposns[i]+offset[0], yposns[i]+offset[1],
                                      psf1[0], psf1[1]+psf1_yvary[i], theta=theta1)
        im1 += tmp

    im1_psf = singleGaussian2d(x0im, y0im, 0, 0, psf1[0], psf1[1], theta=theta1)
    im2_psf = singleGaussian2d(x0im, y0im, offset[0], offset[1], psf2[0], psf2[1], theta=theta2)

    def makeWcs(offset=0):
        """ Make a fake Wcs

        Parameters
        ----------
        offset : float
          offset the Wcs by this many pixels.
        """
        # taken from $AFW_DIR/tests/testMakeWcs.py
        metadata = dafBase.PropertySet()
        metadata.set("SIMPLE", "T")
        metadata.set("BITPIX", -32)
        metadata.set("NAXIS", 2)
        metadata.set("NAXIS1", 1024)
        metadata.set("NAXIS2", 1153)
        metadata.set("RADECSYS", 'FK5')
        metadata.set("EQUINOX", 2000.)
        metadata.setDouble("CRVAL1", 215.604025685476)
        metadata.setDouble("CRVAL2", 53.1595451514076)
        metadata.setDouble("CRPIX1", 1109.99981456774 + offset)
        metadata.setDouble("CRPIX2", 560.018167811613 + offset)
        metadata.set("CTYPE1", 'RA---SIN')
        metadata.set("CTYPE2", 'DEC--SIN')
        metadata.setDouble("CD1_1", 5.10808596133527E-05)
        metadata.setDouble("CD1_2", 1.85579539217196E-07)
        metadata.setDouble("CD2_2", -5.10281493481982E-05)
        metadata.setDouble("CD2_1", -8.27440751733828E-07)
        return afwImage.makeWcs(metadata)

    def makeExposure(imgArray, psfArray, imgVariance):
        """! Convert an image numpy.array and corresponding PSF numpy.array into an exposure.

        Add the (constant) variance plane equal to `imgVariance`.

        @param imgArray 2-d numpy.array containing the image
        @param psfArray 2-d numpy.array containing the PSF image
        @param imgVariance variance of input image
        @return a new exposure containing the image, PSF and desired variance plane
        """
        # All this code to convert the template image array/psf array into an exposure.
        bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Point2I(imgArray.shape[0]-1, imgArray.shape[1]-1))
        im1ex = afwImage.ExposureD(bbox)
        im1ex.getMaskedImage().getImage().getArray()[:, :] = imgArray
        im1ex.getMaskedImage().getVariance().getArray()[:, :] = imgVariance
        psfBox = afwGeom.Box2I(afwGeom.Point2I(-20, -20), afwGeom.Point2I(20, 20))  # a 41x41 pixel psf
        psf = afwImage.ImageD(psfBox)
        psfBox.shift(afwGeom.Extent2I(256, 256))
        im1_psf_sub = psfArray[psfBox.getMinX():psfBox.getMaxX()+1, psfBox.getMinY():psfBox.getMaxY()+1]
        psf.getArray()[:, :] = im1_psf_sub
        psfK = afwMath.FixedKernel(psf)
        psfNew = measAlg.KernelPsf(psfK)
        im1ex.setPsf(psfNew)
        wcs = makeWcs()
        im1ex.setWcs(wcs)
        return im1ex

    im1ex = makeExposure(im1, im1_psf, svar)  # Science image
    im2ex = makeExposure(im2, im2_psf, tvar)  # Template

    return im1ex, im2ex


class DiffimCorrectionTest(lsst.utils.tests.TestCase):
    """!A test case for the diffim image decorrelation algorithm.
    """

    def setUp(self):
        self.psf1_sigma = 3.3  # sigma of psf of science image
        self.psf2_sigma = 2.2  # sigma of psf of template image

        self.statsControl = afwMath.StatisticsControl()
        self.statsControl.setNumSigmaClip(3.)
        self.statsControl.setNumIter(3)
        self.statsControl.setAndMask(afwImage.MaskU.getPlaneBitMask(["INTRP", "EDGE", "SAT", "CR",
                                                                     "DETECTED", "BAD",
                                                                     "NO_DATA", "DETECTED_NEGATIVE"]))

    def _setUpImages(self, svar=0.04, tvar=0.04, varyPsf=0.):
        """!Generate a fake aligned template and science image.
        """

        self.svar = svar  # variance of noise in science image
        self.tvar = tvar  # variance of noise in template image

        self.im1ex, self.im2ex \
            = makeFakeImages(svar=self.svar, tvar=self.tvar, psf1=self.psf1_sigma, psf2=self.psf2_sigma,
                             n_sources=50, psf_yvary_factor=varyPsf, verbose=False)

    def _computeVarianceMean(self, maskedIm):
        statObj = afwMath.makeStatistics(maskedIm.getVariance(),
                                         maskedIm.getMask(), afwMath.MEANCLIP,
                                         self.statsControl)
        mn = statObj.getValue(afwMath.MEANCLIP)
        return mn

    def _computePixelVariance(self, maskedIm):
        statObj = afwMath.makeStatistics(maskedIm, afwMath.VARIANCECLIP,
                                         self.statsControl)
        var = statObj.getValue(afwMath.VARIANCECLIP)
        return var

    def tearDown(self):
        del self.im1ex
        del self.im2ex

    def _makeAndTestUncorrectedDiffim(self):
        """Create the (un-decorrelated) diffim, and verify that its variance is too low.
        """
        # Create the matching kernel. We used Gaussian PSFs for im1 and im2, so we can compute the "expected"
        # matching kernel sigma.
        psf1_sig = self.im1ex.getPsf().computeShape().getDeterminantRadius()
        psf2_sig = self.im2ex.getPsf().computeShape().getDeterminantRadius()
        sig_match = np.sqrt((psf1_sig**2. - psf2_sig**2.))
        # Sanity check - make sure PSFs are correct.
        self.assertClose(sig_match, np.sqrt((self.psf1_sigma**2. - self.psf2_sigma**2.)), rtol=1e-5)
        # mKernel = measAlg.SingleGaussianPsf(31, 31, sig_match)
        x0 = np.arange(-16, 16, 1)
        y0 = x0.copy()
        x0im, y0im = np.meshgrid(x0, y0)
        matchingKernel = singleGaussian2d(x0im, y0im, -1., -1., sigma_x=sig_match, sigma_y=sig_match)
        kernelImg = afwImage.ImageD(matchingKernel.shape[0], matchingKernel.shape[1])
        kernelImg.getArray()[:, :] = matchingKernel
        mKernel = afwMath.FixedKernel(kernelImg)

        # Create the matched template by convolving the template with the matchingKernel
        matched_im2ex = self.im2ex.clone()
        convCntrl = afwMath.ConvolutionControl(False, True, 0)
        afwMath.convolve(matched_im2ex.getMaskedImage(), self.im2ex.getMaskedImage(), mKernel, convCntrl)

        # Expected (ideal) variance of difference image
        expected_var = self.svar + self.tvar
        if verbose:
            print('EXPECTED VARIANCE:', expected_var)

        # Create the diffim (uncorrected)
        # Uncorrected diffim exposure - variance plane is wrong (too low)
        tmp_diffExp = self.im1ex.getMaskedImage().clone()
        tmp_diffExp -= matched_im2ex.getMaskedImage()
        var = self._computeVarianceMean(tmp_diffExp)
        self.assertLess(var, expected_var)

        # Uncorrected diffim exposure - variance is wrong (too low) - same as above but on pixels
        diffExp = self.im1ex.clone()
        tmp = diffExp.getMaskedImage()
        tmp -= matched_im2ex.getMaskedImage()
        var = self._computePixelVariance(diffExp.getMaskedImage())
        self.assertLess(var, expected_var)

        # Uncorrected diffim exposure - variance plane is wrong (too low)
        mn = self._computeVarianceMean(diffExp.getMaskedImage())
        self.assertLess(mn, expected_var)
        if verbose:
            print('UNCORRECTED VARIANCE:', var, mn)

        return diffExp, mKernel, expected_var

    def _runDecorrelationTask(self, diffExp, mKernel):
        """ Run the decorrelation task on the given diffim with the given matching kernel
        """
        task = DecorrelateALKernelTask()
        decorrResult = task.run(self.im1ex, self.im2ex, diffExp, mKernel)
        corrected_diffExp = decorrResult.correctedExposure
        return corrected_diffExp

    def _testDecorrelation(self, expected_var, corrected_diffExp):
        """ Check that the variance of the corrected diffim matches the theoretical value.
        """
        # Corrected diffim - variance should be close to expected.
        # We set the tolerance a bit higher here since the simulated images have many bright stars
        var = self._computePixelVariance(corrected_diffExp.getMaskedImage())
        self.assertClose(var, expected_var, rtol=0.05)

        # Check statistics of variance plane in corrected diffim
        mn = self._computeVarianceMean(corrected_diffExp.getMaskedImage())
        if verbose:
            print('CORRECTED VARIANCE:', var, mn)
        self.assertClose(mn, expected_var, rtol=0.02)
        self.assertClose(var, mn, rtol=0.05)

    def _testDiffimCorrection(self, svar, tvar):
        """ Run decorrelation and check the variance of the corrected diffim.
        """
        self._setUpImages(svar=svar, tvar=tvar)
        diffExp, mKernel, expected_var = self._makeAndTestUncorrectedDiffim()
        corrected_diffExp = self._runDecorrelationTask(diffExp, mKernel)
        self._testDecorrelation(expected_var, corrected_diffExp)

    def testDiffimCorrection(self):
        """Test decorrelated diffim from images with different combinations of variances.
        """
        # Same variance
        self._testDiffimCorrection(svar=0.04, tvar=0.04)
        # Science image variance is higher than that of the template.
        self._testDiffimCorrection(svar=0.08, tvar=0.04)
        # Template variance is higher than that of the science img.
        self._testDiffimCorrection(svar=0.04, tvar=0.08)

    def _runDecorrelationTaskMapReduced(self, diffExp, mKernel):
        """ Run decorrelation using the imageMapReducer.
        """
        config = DecorrelateALKernelMapReduceConfig()
        config.borderSizeX = config.borderSizeY = 3
        config.reducerSubtask.reduceOperation = 'average'
        task = ImageMapReduceTask(config=config)
        decorrResult = task.run(diffExp, template=self.im2ex, science=self.im1ex,
                                psfMatchingKernel=mKernel, forceEvenSized=True)
        corrected_diffExp = decorrResult.exposure
        return corrected_diffExp

    def _testDiffimCorrection_mapReduced(self, svar, tvar, varyPsf=0.0):
        """ Run decorrelation using the imageMapReduce task, and check the variance of
            the corrected diffim.
        """
        self._setUpImages(svar=svar, tvar=tvar, varyPsf=varyPsf)
        diffExp, mKernel, expected_var = self._makeAndTestUncorrectedDiffim()
        corrected_diffExp = self._runDecorrelationTaskMapReduced(diffExp, mKernel)
        self._testDecorrelation(expected_var, corrected_diffExp)
        # Also compare the diffim generated here vs. the non-ImageMapReduce one
        corrected_diffExp_OLD = self._runDecorrelationTask(diffExp, mKernel)
        self.assertMaskedImagesAlmostEqual(corrected_diffExp.getMaskedImage(),
                                           corrected_diffExp_OLD.getMaskedImage())

    def testDiffimCorrection_mapReduced(self):
        """Test decorrelated diffim when using the imageMapReduce task.
           Compare results with those from the original DecorrelateALKernelTask.
        """
        # Same variance
        self._testDiffimCorrection_mapReduced(svar=0.04, tvar=0.04)
        # Science image variance is higher than that of the template.
        self._testDiffimCorrection_mapReduced(svar=0.04, tvar=0.08)
        # Template variance is higher than that of the science img.
        self._testDiffimCorrection_mapReduced(svar=0.08, tvar=0.04)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
