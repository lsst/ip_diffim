# This file is part of meas_transiNet.
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

import unittest

import numpy as np

from lsst.ip.diffim.transiNetInterface import TransiNetInterface
from lsst.afw.image import ExposureF
from lsst.ip.diffim.utils import makeTestImage


class TestTNInterface(unittest.TestCase):
    def setUp(self):
        self.interface = TransiNetInterface('temp_var', 'neighbor')

    def test_single_cutout(self):
        """Test running infer on a single empty cutout.
        """

        # Create a pair of empty lsst expsure images
        template = ExposureF(256,256)
        science = template.clone()
        result = self.interface.infer(template, science)
        
        # Test that the result has the correct shape.
        self.assertEqual(result.getDimensions(), science.getDimensions())
        
    def test_single_visit(self):
        """Test running infer on a pair of images with as large as a single visit.
        """
        # Create a pair of empty lsst expsure images
        template = ExposureF(4000,4000)
        science = template.clone()
        result = self.interface.infer(template, science)
        
        # Test that the result has the correct shape.
        self.assertEqual(result.getDimensions(), science.getDimensions())
        

