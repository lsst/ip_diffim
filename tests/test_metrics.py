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

import unittest

import astropy.units as u

from lsst.afw.table import SourceCatalog
import lsst.utils.tests
from lsst.verify import Name
from lsst.verify.gen2tasks.testUtils import MetricTaskTestCase
from lsst.verify.tasks import MetricComputationError

from lsst.ip.diffim.metrics import \
    NumberSciSourcesMetricTask, \
    FractionDiaSourcesToSciSourcesMetricTask


def _makeDummyCatalog(size):
    catalog = SourceCatalog(SourceCatalog.Table.makeMinimalSchema())
    for i in range(size):
        catalog.addNew()
    return catalog


class TestNumSciSources(MetricTaskTestCase):

    @classmethod
    def makeTask(cls):
        return NumberSciSourcesMetricTask()

    def testValid(self):
        catalog = _makeDummyCatalog(3)
        result = self.task.run(catalog)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ip_diffim.numSciSources"))
        self.assertEqual(meas.quantity, len(catalog) * u.count)

    def testEmptyCatalog(self):
        catalog = _makeDummyCatalog(0)
        result = self.task.run(catalog)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ip_diffim.numSciSources"))
        self.assertEqual(meas.quantity, 0 * u.count)

    def testMissingData(self):
        result = self.task.run(None)
        meas = result.measurement
        self.assertIsNone(meas)

    def testGetInputDatasetTypes(self):
        config = self.taskClass.ConfigClass()
        types = self.taskClass.getInputDatasetTypes(config)
        # dict.keys() is a collections.abc.Set, which has a narrower interface than __builtins__.set...
        self.assertSetEqual(set(types.keys()), {"sources"})
        self.assertEqual(types["sources"], "src")


class TestFractionDiaSources(MetricTaskTestCase):

    @classmethod
    def makeTask(cls):
        return FractionDiaSourcesToSciSourcesMetricTask()

    def testValid(self):
        sciCatalog = _makeDummyCatalog(5)
        diaCatalog = _makeDummyCatalog(3)
        result = self.task.run(sciCatalog, diaCatalog)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ip_diffim.fracDiaSourcesToSciSources"))
        self.assertEqual(meas.quantity, len(diaCatalog) / len(sciCatalog) * u.dimensionless_unscaled)

    def testEmptyDiaCatalog(self):
        sciCatalog = _makeDummyCatalog(5)
        diaCatalog = _makeDummyCatalog(0)
        result = self.task.run(sciCatalog, diaCatalog)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ip_diffim.fracDiaSourcesToSciSources"))
        self.assertEqual(meas.quantity, 0.0 * u.dimensionless_unscaled)

    def testEmptySciCatalog(self):
        sciCatalog = _makeDummyCatalog(0)
        diaCatalog = _makeDummyCatalog(3)
        with self.assertRaises(MetricComputationError):
            self.task.run(sciCatalog, diaCatalog)

    def testEmptyCatalogs(self):
        sciCatalog = _makeDummyCatalog(0)
        diaCatalog = _makeDummyCatalog(0)
        with self.assertRaises(MetricComputationError):
            self.task.run(sciCatalog, diaCatalog)

    def testMissingData(self):
        result = self.task.run(None, None)
        meas = result.measurement
        self.assertIsNone(meas)

    def testSemiMissingData(self):
        result = self.task.run(sciSources=_makeDummyCatalog(3), diaSources=None)
        meas = result.measurement
        self.assertIsNone(meas)

    def testGetInputDatasetTypes(self):
        config = self.taskClass.ConfigClass()
        types = self.taskClass.getInputDatasetTypes(config)
        # dict.keys() is a collections.abc.Set, which has a narrower interface than __builtins__.set...
        self.assertSetEqual(set(types.keys()), {"sciSources", "diaSources"})
        self.assertEqual(types["sciSources"], "src")
        self.assertEqual(types["diaSources"], "deepDiff_diaSrc")

    # TODO: add a test for the templating in FractionDiaSourcesToSciSourcesMetricConfig
    # once we've migrated to PipelineTaskConfig


# Hack around unittest's hacky test setup system
del MetricTaskTestCase


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
