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

from lsst.ip.diffim.utils import populate_sattle_visit_cache
import lsst.utils.tests

from test_detectAndMeasure import makeVisitInfo, MockResponse


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
