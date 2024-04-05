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
from lsst.ip.diffim.utils import angleMean
import lsst.utils.tests


class UtilsCalculationsTest(lsst.utils.tests.TestCase):

    """Unit tests for calculations in utils.
    """

    def test_angleMean(self):
        seed = 5
        nSrc = 100
        angleOffset = 30
        rng = np.random.RandomState(seed)
        angles = (rng.rand(nSrc) - 0.5)*20 + angleOffset
        self.assertFloatsAlmostEqual(angleOffset, angleMean(angles).asDegrees(), rtol=0.01)
