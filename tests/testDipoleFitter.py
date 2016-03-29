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

from numpy import (random as np_random,
                   array as np_array,
                   mgrid as np_mgrid)

from lsst.afw.geom import (Box2I, Point2I, Point2D)

## UTILITY CLASS WITH STATIC METHODS FOR DIPOLE TESTING ###


class DipoleTestUtils(object):

    @staticmethod
    def makeStarImage(w=101, h=101, xc=[15.3], yc=[18.6], flux=[2500], psfSigma=2., noise=10.0,
                      gradientParams=None, schema=None):

        from lsst.meas.base.tests import TestDataset
        bbox = Box2I(Point2I(0, 0), Point2I(w-1, h-1))
        dataset = TestDataset(bbox, psfSigma=psfSigma, threshold=1.)

        for i in xrange(len(xc)):
            dataset.addSource(flux=flux[i], centroid=Point2D(xc[i], yc[i]))

        if schema is None:
            schema = TestDataset.makeMinimalSchema()
        exposure, catalog = dataset.realize(noise=noise, schema=schema)

        if gradientParams is not None:
            y, x = np_mgrid[:w, :h]
            gp = gradientParams
            gradient = gp[0] + gp[1] * x + gp[2] * y
            if len(gradientParams) > 3:  # it includes a set of 2nd-order polynomial params
                gradient += gp[3] * x*y + gp[4] * x*x + gp[5] * y*y
            imgArr = exposure.getMaskedImage().getArrays()[0]
            imgArr += gradient

        return exposure, catalog

    @staticmethod
    def makeDipoleImage(w=101, h=101, xcenPos=[27.], ycenPos=[25.], xcenNeg=[23.], ycenNeg=[25.],
                        psfSigma=2., flux=[30000.], fluxNeg=None, noise=10., gradientParams=None):

        posImage, posCatalog = DipoleTestUtils.makeStarImage(
            w, h, xcenPos, ycenPos, flux=flux, psfSigma=psfSigma,
            gradientParams=gradientParams, noise=noise)

        if fluxNeg is None:
            fluxNeg = flux
        negImage, negCatalog = DipoleTestUtils.makeStarImage(
            w, h, xcenNeg, ycenNeg, flux=fluxNeg, psfSigma=psfSigma,
            gradientParams=gradientParams, noise=noise)

        dipole = posImage.clone()
        di = dipole.getMaskedImage()
        di -= negImage.getMaskedImage()

        # Carry through pos/neg detection masks to new planes in diffim image
        dm = di.getMask()
        posDetectedBits = posImage.getMaskedImage().getMask().getArray() == dm.getPlaneBitMask("DETECTED")
        negDetectedBits = negImage.getMaskedImage().getMask().getArray() == dm.getPlaneBitMask("DETECTED")
        pos_det = dm.addMaskPlane("DETECTED_POS")  # new mask plane -- different from "DETECTED"
        neg_det = dm.addMaskPlane("DETECTED_NEG")  # new mask plane -- different from "DETECTED_NEGATIVE"
        dma = dm.getArray()
        # set the two custom mask planes to these new masks
        dma[:, :] = posDetectedBits * pos_det + negDetectedBits * neg_det
        return dipole, (posImage, posCatalog), (negImage, negCatalog)
