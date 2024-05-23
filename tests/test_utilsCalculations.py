# This file is part of ip_diffim.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np

import lsst.geom
from lsst.ip.diffim.utils import angleMean, getPsfFwhm
import lsst.meas.algorithms as measAlg
import lsst.utils.tests


class UtilsCalculationsTest(lsst.utils.tests.TestCase):

    """Unit tests for calculations in utils.
    """

    def test_angleMean(self):
        """Function for averages of angles.
        """
        seed = 5
        nSrc = 100
        angleOffset = 30
        rng = np.random.RandomState(seed)
        angles = np.radians((rng.rand(nSrc) - 0.5)*20 + angleOffset)
        self.assertFloatsAlmostEqual(angleOffset, angleMean(angles).asDegrees(), rtol=0.01)

    def test_getPsfFwhm(self):
        """Calculation of FWHM from a realization of the PSF
        """
        sigmaToFwhm = 2*np.log(2*np.sqrt(2))

        def make_and_check_psf(xKsize, yKsize, sigma):
            psf = measAlg.SingleGaussianPsf(xKsize, yKsize, sigma)
            psfSize = getPsfFwhm(psf)
            psfSize2d = getPsfFwhm(psf, average=False)
            self.assertFloatsAlmostEqual(sigma*sigmaToFwhm, psfSize, rtol=0.01)
            self.assertFloatsAlmostEqual(psfSize, psfSize2d[0], rtol=0.01)
            self.assertFloatsAlmostEqual(psfSize, psfSize2d[1], rtol=0.01)

        # Test equal and unequal axes with a narrow PSF
        make_and_check_psf(23, 23, 1)
        make_and_check_psf(23, 25, 1)
        make_and_check_psf(23, 21, 1)

        # Test equal and unequal axes with a narrow PSF
        make_and_check_psf(23, 23, 1.23456)
        make_and_check_psf(23, 25, 1.23456)
        make_and_check_psf(23, 21, 1.23456)

        # Test equal and unequal axes with a wide PSF
        make_and_check_psf(23, 23, 2.1)
        make_and_check_psf(23, 25, 2.1)
        make_and_check_psf(23, 21, 2.1)
