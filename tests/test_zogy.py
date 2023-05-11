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
import numpy as np
import unittest

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.daf.base as dafBase
import lsst.geom as geom
from lsst.ip.diffim.zogy import ZogyTask, ZogyConfig
import lsst.meas.algorithms as measAlg
import lsst.utils.tests
from test_imageDecorrelation import singleGaussian2d

try:
    type(verbose)
except NameError:
    verbose = False


def setup_module(module):
    lsst.utils.tests.init()


def makeFakeImages(size=(256, 256), svar=0.04, tvar=0.04, psf1=3.3, psf2=2.2, offset=None,
                   psf_yvary_factor=0., varSourceChange=1/50., theta1=0., theta2=0.,
                   n_sources=50, seed=66, verbose=False):
    """Make two exposures: science and template pair with flux sources and random noise.
    In all cases below, index (1) is the science image, and (2) is the template.

    Parameters
    ----------
    size : `tuple` of `int`
        Image pixel size (x,y). Pixel coordinates are set to
        (-size[0]//2:size[0]//2, -size[1]//2:size[1]//2)
    svar, tvar : `float`, optional
        Per pixel variance of the added noise.
    psf1, psf2 : `float`, optional
        std. dev. of (Gaussian) PSFs for the two images in x,y direction. Default is
        [3.3, 3.3] and [2.2, 2.2] for im1 and im2 respectively.
    offset : `float`, optional
        add a constant (pixel) astrometric offset between the two images.
    psf_yvary_factor : `float`, optional
        psf_yvary_factor vary the y-width of the PSF across the x-axis of the science image (zero,
        the default, means no variation)
    varSourceChange : `float`, optional
        varSourceChange add this amount of fractional flux to a single source closest to
        the center of the science image.
    theta1, theta2: `float`, optional
        PSF Gaussian rotation angles in degrees.
    n_sources : `int`, optional
        The number of sources to add to the images. If zero, no sources are
        generated just background noise.
    seed : `int`, optional
        Random number generator seed.
    verbose : `bool`, optional
        Print some actual values.

    Returns
    -------
    im1, im2 : `lsst.afw.image.Exposure`
        The science and template exposures.

    Notes
    -----
    If ``n_sources > 0`` and ``varSourceChange > 0.`` exactly one source,
    that is closest to the center, will have different fluxes in the two
    generated images. The flux on the science image will be higher by
    ``varSourceChange`` fraction.

    Having sources near the edges really messes up the
    fitting (probably because of the convolution). So we make sure no
    sources are near the edge.

    Also it seems that having the variable source with a large
    flux increase also messes up the fitting (seems to lead to
    overfitting -- perhaps to the source itself). This might be fixed by
    adding more constant sources.
    """
    rng = np.random.default_rng(seed)

    psf1 = [3.3, 3.3] if psf1 is None else psf1
    if not hasattr(psf1, "__len__"):
        psf1 = [psf1, psf1]
    psf2 = [2.2, 2.2] if psf2 is None else psf2
    if not hasattr(psf2, "__len__"):
        psf2 = [psf2, psf2]
    offset = [0., 0.] if offset is None else offset   # astrometric offset (pixels) between the two images
    if verbose:
        print('Science PSF:', psf1, theta1)
        print('Template PSF:', psf2, theta2)
        print(np.sqrt(psf1[0]**2 - psf2[0]**2))
        print('Offset:', offset)

    xim = np.arange(-size[0]//2, size[0]//2, 1)  # Beware that -N//2 != -1*(N//2) for odd numbers
    yim = np.arange(-size[1]//2, size[1]//2, 1)
    x0im, y0im = np.meshgrid(xim, yim)

    im1 = rng.normal(scale=np.sqrt(svar), size=x0im.shape)  # variance of science image
    im2 = rng.normal(scale=np.sqrt(tvar), size=x0im.shape)  # variance of template

    if n_sources > 0:
        fluxes = rng.uniform(50, 30000, n_sources)
        xposns = rng.uniform(xim.min() + 16, xim.max() - 5, n_sources)
        yposns = rng.uniform(yim.min() + 16, yim.max() - 5, n_sources)

        # Make the source closest to the center of the image the one that increases in flux
        ind = np.argmin(xposns**2. + yposns**2.)

        # vary the y-width of psf across x-axis of science image (zero means no variation):
        psf1_yvary = psf_yvary_factor*(yim.mean() - yposns)/yim.max()
        if verbose:
            print('PSF y spatial-variation:', psf1_yvary.min(), psf1_yvary.max())

    for i in range(n_sources):
        flux = fluxes[i]
        tmp = flux*singleGaussian2d(x0im, y0im, xposns[i], yposns[i], psf2[0], psf2[1], theta=theta2)
        im2 += tmp
        if i == ind:
            flux += flux*varSourceChange
        tmp = flux*singleGaussian2d(x0im, y0im, xposns[i] + offset[0], yposns[i] + offset[1],
                                    psf1[0], psf1[1] + psf1_yvary[i], theta=theta1)
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
        metadata.set("RADESYS", 'FK5')
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
        return afwGeom.makeSkyWcs(metadata)

    def makeExposure(imgArray, psfArray, imgVariance):
        """Convert an image and corresponding PSF into an exposure.

        Set the (constant) variance plane equal to ``imgVariance``.

        Parameters
        ----------
        imgArray : `numpy.ndarray`
            2D array containing the image.
        psfArray : `numpy.ndarray`
            2D array containing the PSF image.
        imgVariance : `float` or `numpy.ndarray`
            Set the variance plane to this value. If an array, must be broadcastable to ``imgArray.shape``.

        Returns
        -------
        im1ex : `lsst.afw.image.Exposure`
            The new exposure.
        """
        # All this code to convert the template image array/psf array into an exposure.
        bbox = geom.Box2I(geom.Point2I(0, 0), geom.Point2I(imgArray.shape[1] - 1, imgArray.shape[0] - 1))
        im1ex = afwImage.ExposureD(bbox)
        im1ex.image.array[:, :] = imgArray
        im1ex.variance.array[:, :] = imgVariance
        psfBox = geom.Box2I(geom.Point2I(-12, -12), geom.Point2I(12, 12))  # a 25x25 pixel psf
        psf = afwImage.ImageD(psfBox)
        psfBox.shift(geom.Extent2I(-(-size[0]//2), -(-size[1]//2)))  # -N//2 != -(N//2) for odd numbers
        im1_psf_sub = psfArray[psfBox.getMinY():psfBox.getMaxY() + 1, psfBox.getMinX():psfBox.getMaxX() + 1]
        psf.array[:, :] = im1_psf_sub
        psfK = afwMath.FixedKernel(psf)
        psfNew = measAlg.KernelPsf(psfK)
        im1ex.setPsf(psfNew)
        wcs = makeWcs()
        im1ex.setWcs(wcs)
        return im1ex

    im1ex = makeExposure(im1, im1_psf, svar)  # Science image
    im2ex = makeExposure(im2, im2_psf, tvar)  # Template

    return im1ex, im2ex


def isPowerOfTwo(x):
    """Returns True if x is a power of 2"""
    while x > 1:
        if x & 1 != 0:
            return False
        x >>= 1
    return True


class ZogyTest(lsst.utils.tests.TestCase):
    """A test case for the Zogy task.
    """

    def setUp(self):
        self.psf1_sigma = 3.3  # sigma of psf of science image
        self.psf2_sigma = 2.2  # sigma of psf of template image

        self.statsControl = afwMath.StatisticsControl()
        self.statsControl.setNumSigmaClip(3.)
        self.statsControl.setNumIter(3)
        self.statsControl.setAndMask(afwImage.Mask
                                     .getPlaneBitMask(["INTRP", "EDGE", "SAT", "CR",
                                                       "DETECTED", "BAD",
                                                       "NO_DATA", "DETECTED_NEGATIVE"]))

    def _computeVarianceMean(self, maskedIm):
        statObj = afwMath.makeStatistics(maskedIm.variance,
                                         maskedIm.mask, afwMath.MEANCLIP,
                                         self.statsControl)
        mn = statObj.getValue(afwMath.MEANCLIP)
        return mn

    def _computePixelVariance(self, maskedIm):
        statObj = afwMath.makeStatistics(maskedIm, afwMath.VARIANCECLIP,
                                         self.statsControl)
        var = statObj.getValue(afwMath.VARIANCECLIP)
        return var

    def _computePixelMean(self, maskedIm):
        statObj = afwMath.makeStatistics(maskedIm, afwMath.MEANCLIP,
                                         self.statsControl)
        var = statObj.getValue(afwMath.MEANCLIP)
        return var

    def testFourierTransformConvention(self):
        """Test numpy FFT normalization factor convention matches our assumption."""
        D = np.arange(16).reshape(4, 4)
        fD = np.real(np.fft.fft2(D))
        self.assertFloatsAlmostEqual(
            fD[0, 0], 120., rtol=None,
            msg="Numpy FFT does not use expected default normalization"
            " convention (1 in forward, 1/Npix in inverse operation).")

    def testSplitBorder(self):
        """Test outer border box splitting around an inner box"""
        config = ZogyConfig()
        task = ZogyTask(config=config)

        bb = geom.Box2I(geom.Point2I(5, 10), geom.Extent2I(20, 30))
        D = afwImage.ImageI(bb)
        innerbox = bb.erodedBy(geom.Extent2I(3, 4))
        D[innerbox] = 1

        borderboxes = task.splitBorder(innerbox, bb)
        for x in borderboxes:
            D[x] += 1
        # The splitting should cover all border pixels exactly once
        self.assertTrue(np.all(D.array == 1), "Border does not cover all pixels exactly once.")

    def testGenerateGrid(self):
        """Test that the generated grid covers the whole image"""
        config = ZogyConfig()
        task = ZogyTask(config=config)
        bb = geom.Box2I(geom.Point2I(5, 10), geom.Extent2I(200, 300))
        D = afwImage.ImageI(bb)
        grid = task.generateGrid(bb, geom.Extent2I(15, 15), geom.Extent2I(20, 30), powerOfTwo=True)
        for x in grid:
            h = x.outerBox.getHeight()
            w = x.outerBox.getWidth()
            self.assertTrue(isPowerOfTwo(h), "Box height is not power of two")
            self.assertTrue(isPowerOfTwo(w), "Box width is not power of two")
            D[x.innerBox] += 1
        self.assertTrue(np.all(D.array == 1), "Grid inner boxes do not cover all pixels exactly once.")

    def testWholeImageGrid(self):
        """Test that a 1-cell `grid` is actually the whole image"""
        config = ZogyConfig()
        task = ZogyTask(config=config)
        bb = geom.Box2I(geom.Point2I(5, 10), geom.Extent2I(200, 300))
        D = afwImage.ImageI(bb)
        grid = task.generateGrid(bb, geom.Extent2I(15, 15), bb.getDimensions())
        self.assertTrue(len(grid) == 1, "Grid length is not 1")
        x = grid[0]
        D[x.innerBox] += 1
        self.assertTrue(np.all(D.array == 1), "Single cell does not cover the original image.")

    def testZogyNewImplementation(self):
        """DM-25115 implementation test.

        Notes
        -----
        See diffimTests: tickets/DM-25115_zogy_implementation/DM-25115_zogy_unit_test_development.ipynb
        """

        # self.svar = svar  # variance of noise in science image
        # self.tvar = tvar  # variance of noise in template image

        # Sourceless case
        self.im1ex, self.im2ex \
            = makeFakeImages(size=(256, 256), svar=100., tvar=100.,
                             psf1=self.psf1_sigma, psf2=self.psf2_sigma,
                             n_sources=0, psf_yvary_factor=0, varSourceChange=0.1,
                             seed=1, verbose=False)

        config = ZogyConfig()
        config.scaleByCalibration = False
        task = ZogyTask(config=config)
        res = task.run(self.im1ex, self.im2ex)

        bbox = res.diffExp.getBBox()
        subBbox = bbox.erodedBy(lsst.geom.Extent2I(25, 25))
        subExp = res.diffExp[subBbox]
        pixvar = self._computePixelVariance(subExp.maskedImage)
        varmean = self._computeVarianceMean(subExp.maskedImage)
        # Due to 3 sigma clipping, this is not so precise
        self.assertFloatsAlmostEqual(pixvar, 200, rtol=0.1, atol=None)
        self.assertFloatsAlmostEqual(varmean, 200, rtol=0.05, atol=None)
        S = res.scoreExp.image.array / np.sqrt(res.scoreExp.variance.array)
        self.assertLess(np.amax(S), 5.)  # Source not detected

        # ==========
        self.im1ex, self.im2ex \
            = makeFakeImages(size=(256, 256), svar=10., tvar=10.,
                             psf1=self.psf1_sigma, psf2=self.psf2_sigma,
                             n_sources=10, psf_yvary_factor=0, varSourceChange=0.1,
                             seed=1, verbose=False)
        task = ZogyTask(config=config)
        res = task.run(self.im1ex, self.im2ex)
        S = res.scoreExp.image.array / np.sqrt(res.scoreExp.variance.array)
        self.assertGreater(np.amax(S), 5.)  # Source detected


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
