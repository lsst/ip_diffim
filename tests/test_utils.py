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


import requests
from unittest import mock

import lsst.afw.cameraGeom.testUtils
import lsst.daf.butler
from lsst.images import obs_info_from_legacy
from lsst.images.tests import get_dp2_exposure_record
from lsst.ip.diffim.utils import populate_sattle_visit_cache, record_from_obs_info
import lsst.utils.tests

from test_detectAndMeasure import makeVisitInfo, MockResponse
from utils import makeTestImage


class ExposureRecordTest(lsst.utils.tests.TestCase):
    """Tests of rebuilding an exposure record from an observation info."""

    def setUp(self):
        self.universe = lsst.daf.butler.DimensionUniverse()
        self.record = get_dp2_exposure_record(self.universe)
        # lsst.images needs per-amplifier raw geometry and a field angle
        # transform, which the trivial test detectors do not have.
        self.detector = list(lsst.afw.cameraGeom.testUtils.CameraWrapper().camera)[0]
        self.obsInfo = self._makeObsInfo(self.record)

    def _makeObsInfo(self, record):
        """Return the observation info that a converted image built from
        ``record`` would carry.
        """
        return obs_info_from_legacy(makeVisitInfo(), record, self.detector, detector_exposure_id=12345)

    def _rebuild(self):
        """Return the record rebuilt from the observation info, with additional
        identifiers added.
        """
        return record_from_obs_info(self.obsInfo, self.record.instrument, self.record.id, self.universe)

    def test_fields(self):
        """Every field the conversion reads comes back with the value the
        original record had.
        """
        record = self._rebuild()
        for field in ("instrument", "id", "obs_id", "group", "day_obs", "physical_filter",
                      "exposure_time", "seq_num", "seq_start", "seq_end", "can_see_sky", "timespan"):
            self.assertEqual(getattr(record, field), getattr(self.record, field), msg=field)
        self.assertFloatsAlmostEqual(record.azimuth, self.record.azimuth, rtol=1e-14)
        self.assertFloatsAlmostEqual(record.zenith_angle, self.record.zenith_angle, rtol=1e-14)

    def test_conversion_unchanged(self):
        """Converting with the rebuilt record gives the same observation
        info as converting with the original one.
        """
        self.assertEqual(self._makeObsInfo(self._rebuild()), self.obsInfo)

    def test_missing_obs_info(self):
        with self.assertRaises(ValueError):
            record_from_obs_info(None, self.record.instrument, self.record.id, self.universe)


class UtilsTest(lsst.utils.tests.TestCase):

    def test_populate_sattle(self):
        response = MockResponse({}, 200, "success")
        visit_info = makeVisitInfo()
        with mock.patch('requests.put', return_value=response):
            populate_sattle_visit_cache(visit_info)

    def test_populate_sattle_raises(self):
        response = MockResponse({}, 500, "failure")
        visit_info = makeVisitInfo()
        with mock.patch('requests.put', return_value=response):
            with self.assertRaises(requests.exceptions.HTTPError):
                populate_sattle_visit_cache(visit_info)

    def test_raise_on_even_kernel(self):
        with self.assertRaises(ValueError):
            makeTestImage(kernelSize=32)
