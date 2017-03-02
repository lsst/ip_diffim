#!/usr/bin/env python
#
# LSST Data Management System
#
# Copyright 2008-2016  AURA/LSST.
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
from __future__ import print_function
from builtins import zip
from builtins import range
import math
import numpy as np
import unittest

import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.afw.coord as afwCoord
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.display.ds9 as ds9
import lsst.daf.base as dafBase
import lsst.meas.algorithms as measAlg
# register the PSF determiner
import lsst.meas.extensions.psfex.psfexPsfDeterminer
assert lsst.meas.extensions.psfex.psfexPsfDeterminer  # make pyflakes happy
from lsst.meas.base import SingleFrameMeasurementTask

try:
    type(verbose)
except NameError:
    verbose = 0
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def psfVal(ix, iy, x, y, sigma1, sigma2, b):
    """Return the value at (ix, iy) of a double Gaussian
       (N(0, sigma1^2) + b*N(0, sigma2^2))/(1 + b)
       centered at (x, y)
    """
    dx, dy = x - ix, y - iy
    theta = np.radians(30)
    ab = 1.0/0.75                       # axis ratio
    c, s = np.cos(theta), np.sin(theta)
    u, v = c*dx - s*dy, s*dx + c*dy

    return (math.exp(-0.5*(u**2 + (v*ab)**2)/sigma1**2) +
            b*math.exp(-0.5*(u**2 + (v*ab)**2)/sigma2**2))/(1 + b)


class SpatialModelPsfTestCase(object):
    """A test case for SpatialModelPsf"""

    def measure(self, footprintSet, exposure):
        """Measure a set of Footprints, returning a SourceCatalog"""
        catalog = afwTable.SourceCatalog(self.schema)
        if display:
            ds9.mtv(exposure, title="Original", frame=0)

        footprintSet.makeSources(catalog)
        print(len(catalog))
        catalog = catalog.copy(deep=True)

        self.measureSources.run(catalog, exposure)
        return catalog

    def makeExposure(self):
        self.width, self.height = 110, 301
        self.mi = afwImage.MaskedImageF(afwGeom.ExtentI(self.width, self.height))
        self.mi.set(0)
        sd = 3                          # standard deviation of image
        self.mi.getVariance().set(sd*sd)
        self.mi.getMask().addMaskPlane("DETECTED")

        self.ksize = 31                      # size of desired kernel

        sigma1 = 1.75
        sigma2 = 2*sigma1

        self.exposure = afwImage.makeExposure(self.mi)
        self.exposure.setPsf(measAlg.DoubleGaussianPsf(self.ksize, self.ksize,
                                                       1.5*sigma1, 1, 0.1))
        crval = afwCoord.makeCoord(afwCoord.ICRS, 0.0*afwGeom.degrees, 0.0*afwGeom.degrees)
        wcs = afwImage.makeWcs(crval, afwGeom.PointD(0, 0), 1.0, 0, 0, 1.0)
        self.exposure.setWcs(wcs)

        #
        # Make a kernel with the exactly correct basis functions.  Useful for debugging
        #
        basisKernelList = afwMath.KernelList()
        for sigma in (sigma1, sigma2):
            basisKernel = afwMath.AnalyticKernel(self.ksize, self.ksize,
                                                 afwMath.GaussianFunction2D(sigma, sigma))
            basisImage = afwImage.ImageD(basisKernel.getDimensions())
            basisKernel.computeImage(basisImage, True)
            basisImage /= np.sum(basisImage.getArray())

            if sigma == sigma1:
                basisImage0 = basisImage
            else:
                basisImage -= basisImage0

            basisKernelList.append(afwMath.FixedKernel(basisImage))

        order = 1                                # 1 => up to linear
        spFunc = afwMath.PolynomialFunction2D(order)

        exactKernel = afwMath.LinearCombinationKernel(basisKernelList, spFunc)
        exactKernel.setSpatialParameters([[1.0, 0, 0],
                                          [0.0, 0.5*1e-2, 0.2e-2]])

        rand = afwMath.Random()               # make these tests repeatable by setting seed

        addNoise = True

        if addNoise:
            im = self.mi.getImage()
            afwMath.randomGaussianImage(im, rand)  # N(0, 1)
            im *= sd                              # N(0, sd^2)
            del im

        xarr, yarr = [], []

        for x, y in [(20, 20), (60, 20),
                     (30, 35),
                     (50, 50),
                     (20, 90), (70, 160), (25, 265), (75, 275), (85, 30),
                     (50, 120), (70, 80),
                     (60, 210), (20, 210),
                     ]:
            xarr.append(x)
            yarr.append(y)

        for x, y in zip(xarr, yarr):
            dx = rand.uniform() - 0.5   # random (centered) offsets
            dy = rand.uniform() - 0.5

            k = exactKernel.getSpatialFunction(1)(x, y)  # functional variation of Kernel ...
            b = (k*sigma1**2/((1 - k)*sigma2**2))       # ... converted double Gaussian's "b"

            #flux = 80000 - 20*x - 10*(y/float(height))**2
            flux = 80000*(1 + 0.1*(rand.uniform() - 0.5))
            I0 = flux*(1 + b)/(2*np.pi*(sigma1**2 + b*sigma2**2))
            for iy in range(y - self.ksize//2, y + self.ksize//2 + 1):
                if iy < 0 or iy >= self.mi.getHeight():
                    continue

                for ix in range(x - self.ksize//2, x + self.ksize//2 + 1):
                    if ix < 0 or ix >= self.mi.getWidth():
                        continue

                    I = I0*psfVal(ix, iy, x + dx, y + dy, sigma1, sigma2, b)
                    Isample = rand.poisson(I) if addNoise else I
                    self.mi.getImage().set(ix, iy, self.mi.getImage().get(ix, iy) + Isample)
                    self.mi.getVariance().set(ix, iy, self.mi.getVariance().get(ix, iy) + I)

    def setExposure(self, exposure):
        self.exposure = exposure
        self.mi = exposure.getMaskedImage()
        self.width, self.height = self.mi.getDimensions()

    def setUp(self):
        config = SingleFrameMeasurementTask.ConfigClass()
        config.slots.apFlux = 'base_CircularApertureFlux_12_0'
        self.schema = afwTable.SourceTable.makeMinimalSchema()

        self.measureSources = SingleFrameMeasurementTask(self.schema, config=config)

        bbox = afwGeom.BoxI(afwGeom.PointI(0, 0), afwGeom.ExtentI(self.width, self.height))
        self.cellSet = afwMath.SpatialCellSet(bbox, 100)

        self.footprintSet = afwDetection.FootprintSet(self.mi, afwDetection.Threshold(100), "DETECTED")

        self.catalog = self.measure(self.footprintSet, self.exposure)

        for source in self.catalog:
            try:
                cand = measAlg.makePsfCandidate(source, self.exposure)
                self.cellSet.insertCandidate(cand)

            except Exception as e:
                print(e)
                continue

    def tearDown(self):
        del self.cellSet
        del self.exposure
        del self.mi
        del self.footprintSet
        del self.catalog
        del self.schema
        del self.measureSources

    def setupDeterminer(self, exposure):
        """Setup the starSelector and psfDeterminer"""
        starSelectorClass = measAlg.starSelectorRegistry["objectSize"]
        starSelectorConfig = starSelectorClass.ConfigClass()
        starSelectorConfig.sourceFluxField = "base_GaussianFlux_flux"
        starSelectorConfig.badFlags = ["base_PixelFlags_flag_edge",
                                       "base_PixelFlags_flag_interpolatedCenter",
                                       "base_PixelFlags_flag_saturatedCenter",
                                       "base_PixelFlags_flag_crCenter",
                                       ]
        starSelectorConfig.widthStdAllowed = 0.5  # Set to match when the tolerance of the test was set

        starSelector = starSelectorClass(schema=self.schema, config=starSelectorConfig)

        psfDeterminerClass = measAlg.psfDeterminerRegistry["psfex"]
        psfDeterminerConfig = psfDeterminerClass.ConfigClass()
        width, height = exposure.getMaskedImage().getDimensions()
        psfDeterminerConfig.sizeCellX = width
        psfDeterminerConfig.sizeCellY = height//3
        psfDeterminerConfig.spatialOrder = 1

        psfDeterminer = psfDeterminerClass(psfDeterminerConfig)

        return starSelector, psfDeterminer

    def subtractStars(self, exposure, catalog, chi_lim=-1):
        """Subtract the exposure's PSF from all the sources in catalog"""
        mi, psf = exposure.getMaskedImage(), exposure.getPsf()

        subtracted = mi.Factory(mi, True)
        for s in catalog:
            xc, yc = s.getX(), s.getY()
            bbox = subtracted.getBBox(afwImage.PARENT)
            if bbox.contains(afwGeom.PointI(int(xc), int(yc))):
                try:
                    measAlg.subtractPsf(psf, subtracted, xc, yc)
                except:
                    pass
        self.subtracted = subtracted
        chi = subtracted.Factory(subtracted, True)
        var = subtracted.getVariance()
        np.sqrt(var.getArray(), var.getArray())  # inplace sqrt
        chi /= var

        if display:
            ds9.mtv(subtracted, title="Subtracted", frame=1)
            ds9.mtv(chi, title="Chi", frame=2)
            xc, yc = exposure.getWidth()//2, exposure.getHeight()//2
            ds9.mtv(psf.computeImage(afwGeom.Point2D(xc, yc)), title="Psf %.1f,%.1f" % (xc, yc), frame=3)

        chi_min, chi_max = np.min(chi.getImage().getArray()), np.max(chi.getImage().getArray())
        if False:
            print(chi_min, chi_max)

        #if chi_lim > 0:
        #    self.assertGreater(chi_min, -chi_lim)
        #    self.assertLess(chi_max, chi_lim)

    def testPsfexDeterminer(self):
        """Test the (Psfex) psfDeterminer on subImages"""

        starSelector, psfDeterminer = self.setupDeterminer(self.exposure)
        metadata = dafBase.PropertyList()

        psfCandidateList = starSelector.run(self.exposure, self.catalog).psfCandidates
        psf, cellSet = psfDeterminer.determinePsf(self.exposure, psfCandidateList, metadata)
        self.exposure.setPsf(psf)

        # Test how well we can subtract the PSF model
        self.subtractStars(self.exposure, self.catalog, chi_lim=4.6)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
