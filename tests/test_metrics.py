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
import lsst.pipe.base.testUtils
from lsst.verify import Name
from lsst.verify.gen2tasks.testUtils import MetricTaskTestCase
from lsst.verify.tasks import MetricComputationError

from lsst.ip.diffim.metrics import \
    NumberSciSourcesMetricTask, \
    FractionDiaSourcesToSciSourcesMetricTask


def _makeDummyCatalog(size, skyFlag=False):
    """Create a trivial catalog for testing source counts.

    Parameters
    ----------
    size : `int`
        The number of entries in the catalog.
    skyFlag : `bool`
        If set, the schema is guaranteed to have the ``sky_source`` flag, and
        one row has it set to `True`. If not set, the ``sky_source`` flag is
        not present.

    Returns
    -------
    catalog : `lsst.afw.table.SourceCatalog`
        A new catalog with ``size`` rows.
    """
    schema = SourceCatalog.Table.makeMinimalSchema()
    if skyFlag:
        schema.addField("sky_source", type="Flag", doc="Sky objects.")
    catalog = SourceCatalog(schema)
    for i in range(size):
        record = catalog.addNew()
    if skyFlag and size > 0:
        record["sky_source"] = True
    return catalog


class TestNumSciSources(MetricTaskTestCase):

    @classmethod
    def makeTask(cls):
        return NumberSciSourcesMetricTask()

    def testValid(self):
        catalog = _makeDummyCatalog(3)
        result = self.task.run(catalog)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ip_diffim.numSciSources"))
        self.assertEqual(meas.quantity, len(catalog) * u.count)

    def testEmptyCatalog(self):
        catalog = _makeDummyCatalog(0)
        result = self.task.run(catalog)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ip_diffim.numSciSources"))
        self.assertEqual(meas.quantity, 0 * u.count)

    def testSkySources(self):
        catalog = _makeDummyCatalog(3, skyFlag=True)
        result = self.task.run(catalog)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ip_diffim.numSciSources"))
        self.assertEqual(meas.quantity, (len(catalog) - 1) * u.count)

    def testMissingData(self):
        result = self.task.run(None)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)


class TestFractionDiaSources(MetricTaskTestCase):

    @classmethod
    def makeTask(cls):
        return FractionDiaSourcesToSciSourcesMetricTask()

    def testValid(self):
        sciCatalog = _makeDummyCatalog(5)
        diaCatalog = _makeDummyCatalog(3)
        result = self.task.run(sciCatalog, diaCatalog)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ip_diffim.fracDiaSourcesToSciSources"))
        self.assertEqual(meas.quantity, len(diaCatalog) / len(sciCatalog) * u.dimensionless_unscaled)

    def testEmptyDiaCatalog(self):
        sciCatalog = _makeDummyCatalog(5)
        diaCatalog = _makeDummyCatalog(0)
        result = self.task.run(sciCatalog, diaCatalog)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
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
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)

    def testSemiMissingData(self):
        result = self.task.run(sciSources=_makeDummyCatalog(3), diaSources=None)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)

    def testSkySources(self):
        sciCatalog = _makeDummyCatalog(5, skyFlag=True)
        diaCatalog = _makeDummyCatalog(3)
        result = self.task.run(sciCatalog, diaCatalog)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ip_diffim.fracDiaSourcesToSciSources"))
        self.assertEqual(meas.quantity, len(diaCatalog) / (len(sciCatalog) - 1) * u.dimensionless_unscaled)


# Hack around unittest's hacky test setup system
del MetricTaskTestCase


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
