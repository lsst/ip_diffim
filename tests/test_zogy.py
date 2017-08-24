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

import numpy as np

import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from test_imageDecorrelation import makeFakeImages

from lsst.ip.diffim.zogy import ZogyTask, ZogyConfig, ZogyMapReduceConfig
from lsst.ip.diffim.imageMapReduce import ImageMapReduceTask

try:
    type(verbose)
except NameError:
    verbose = False


def setup_module(module):
    lsst.utils.tests.init()


class ZogyTest(lsst.utils.tests.TestCase):
    """A test case for the Zogy task.
    """

    def setUp(self):
        self.psf1_sigma = 3.3  # sigma of psf of science image
        self.psf2_sigma = 2.2  # sigma of psf of template image

        self.statsControl = afwMath.StatisticsControl()
        self.statsControl.setNumSigmaClip(3.)
        self.statsControl.setNumIter(3)
        self.statsControl.setAndMask(afwImage.Mask\
                                     .getPlaneBitMask(["INTRP", "EDGE", "SAT", "CR",
                                                       "DETECTED", "BAD",
                                                       "NO_DATA", "DETECTED_NEGATIVE"]))

    def _setUpImages(self, svar=100., tvar=100., varyPsf=0.):
        """Generate a fake aligned template and science image.
        """
        self.svar = svar  # variance of noise in science image
        self.tvar = tvar  # variance of noise in template image

        seed = 666
        self.im1ex, self.im2ex \
            = makeFakeImages(svar=self.svar, tvar=self.tvar,
                             psf1=self.psf1_sigma, psf2=self.psf2_sigma,
                             n_sources=10, psf_yvary_factor=varyPsf,
                             seed=seed, verbose=False)
        # Create an array corresponding to the "expected" subtraction (noise only)
        np.random.seed(seed)
        self.expectedSubtraction = np.random.normal(scale=np.sqrt(svar), size=self.im1ex.getDimensions())
        self.expectedSubtraction -= np.random.normal(scale=np.sqrt(tvar), size=self.im2ex.getDimensions())
        self.expectedVar = np.var(self.expectedSubtraction)
        self.expectedMean = np.mean(self.expectedSubtraction)

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

    def _computePixelMean(self, maskedIm):
        statObj = afwMath.makeStatistics(maskedIm, afwMath.MEANCLIP,
                                         self.statsControl)
        var = statObj.getValue(afwMath.MEANCLIP)
        return var

    def tearDown(self):
        del self.im1ex
        del self.im2ex

    def _compareExposures(self, D_F, D_R, Scorr=False, tol=0.02):
        """Tests to compare the two images (diffim's or Scorr's).

        See below.  Also compare the diffim pixels with the "expected"
        pixels statistics.  Only do the latter if Scorr==False.
        """
        D_F.getMaskedImage().getMask()[:, :] = D_R.getMaskedImage().getMask()
        varMean_F = self._computeVarianceMean(D_F.getMaskedImage())
        varMean_R = self._computeVarianceMean(D_R.getMaskedImage())
        pixMean_F = self._computePixelMean(D_F.getMaskedImage())
        pixMean_R = self._computePixelMean(D_R.getMaskedImage())
        pixVar_F = self._computePixelVariance(D_F.getMaskedImage())
        pixVar_R = self._computePixelVariance(D_R.getMaskedImage())

        if not Scorr:
            self.assertFloatsAlmostEqual(varMean_F, varMean_R, rtol=tol)
            self.assertFloatsAlmostEqual(pixMean_F, self.expectedMean, atol=tol*2.)
            self.assertFloatsAlmostEqual(pixMean_R, self.expectedMean, atol=tol*2.)
            self.assertFloatsAlmostEqual(pixVar_F, pixVar_R, rtol=tol)
            self.assertFloatsAlmostEqual(pixVar_F, self.expectedVar, rtol=tol*2.)
            self.assertFloatsAlmostEqual(pixVar_R, self.expectedVar, rtol=tol*2.)
        else:
            self.assertFloatsAlmostEqual(varMean_F, varMean_R, atol=tol)  # nearly zero so need to use atol
            self.assertFloatsAlmostEqual(pixVar_F, pixVar_R, atol=tol)

        self.assertFloatsAlmostEqual(pixMean_F, pixMean_R, atol=tol*2.)  # nearly zero so need to use atol

    def testZogyDiffim(self):
        """Compute Zogy diffims using Fourier- and Real-space methods.

        Compare the images.  They are not identical but should be
        similar (within ~2%).
        """
        self._setUpImages()
        config = ZogyConfig()
        task = ZogyTask(templateExposure=self.im2ex, scienceExposure=self.im1ex, config=config)
        D_F = task.computeDiffim(inImageSpace=False)
        D_R = task.computeDiffim(inImageSpace=True)
        self._compareExposures(D_F, D_R)

    def _testZogyScorr(self, varAst=0.):
        """Compute Zogy likelihood images (Scorr) using Fourier- and Real-space methods.

        Compare the images. They are not identical but should be similar (within ~2%).
        """
        config = ZogyConfig()
        task = ZogyTask(templateExposure=self.im2ex, scienceExposure=self.im1ex, config=config)
        D_F = task.computeScorr(inImageSpace=False, xVarAst=varAst, yVarAst=varAst)
        D_R = task.computeScorr(inImageSpace=True, xVarAst=varAst, yVarAst=varAst)
        self._compareExposures(D_F, D_R, Scorr=True)

    def testZogyScorr(self):
        """Compute Zogy likelihood images (Scorr) using Fourier- and Real-space methods.

        Do the computation with "astrometric variance" both zero and non-zero.
        Compare the images. They are not identical but should be similar (within ~2%).
        """
        self._setUpImages()
        self._testZogyScorr()
        self._testZogyScorr(varAst=0.1)

    def _testZogyDiffimMapReduced(self, inImageSpace=False, doScorr=False, **kwargs):
        """Test running Zogy using ImageMapReduceTask framework.

        Compare map-reduced version with non-map-reduced version.
        Do it for pure Fourier-based calc. and also for real-space.
        Also for computing pure diffim D and corrected likelihood image Scorr.
        """
        config = ZogyMapReduceConfig()
        config.gridStepX = config.gridStepY = 9
        config.borderSizeX = config.borderSizeY = 3
        if inImageSpace:
            config.gridStepX = config.gridStepY = 8
            config.borderSizeX = config.borderSizeY = 6  # need larger border size for image-space run
        config.reducerSubtask.reduceOperation = 'average'
        task = ImageMapReduceTask(config=config)
        D_mapReduced = task.run(self.im1ex, template=self.im2ex, inImageSpace=inImageSpace,
                                doScorr=doScorr, forceEvenSized=True, **kwargs).exposure

        config = ZogyConfig()
        task = ZogyTask(templateExposure=self.im2ex, scienceExposure=self.im1ex, config=config)
        if not doScorr:
            D = task.computeDiffim(inImageSpace=inImageSpace, **kwargs)
        else:
            D = task.computeScorr(inImageSpace=inImageSpace, **kwargs)

        self._compareExposures(D_mapReduced, D, tol=0.04, Scorr=doScorr)

    def testZogyDiffimMapReduced(self):
        """Test running Zogy using ImageMapReduceTask framework.

        Compare map-reduced version with non-map-reduced version.
        Do it for pure Fourier-based calc. and also for real-space.
        Do it for ZOGY diffim and corrected likelihood image Scorr.
        For Scorr, do it for zero and non-zero astrometric variance.
        """
        self._setUpImages()
        self._testZogyDiffimMapReduced(inImageSpace=False)
        self._testZogyDiffimMapReduced(inImageSpace=True)
        self._testZogyDiffimMapReduced(inImageSpace=False, doScorr=True)
        self._testZogyDiffimMapReduced(inImageSpace=True, doScorr=True)
        self._testZogyDiffimMapReduced(inImageSpace=False, doScorr=True, xVarAst=0.1, yVarAst=0.1)
        self._testZogyDiffimMapReduced(inImageSpace=True, doScorr=True, xVarAst=0.1, yVarAst=0.1)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
