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

import math
import unittest
import uuid

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose
import numpy as np
import pandas

from lsst.geom import SpherePoint
from lsst.afw.table import SourceCatalog
import lsst.utils.tests
import lsst.pipe.base.testUtils
from lsst.verify import Name
from lsst.verify.gen2tasks.testUtils import MetricTaskTestCase
from lsst.verify.tasks import MetricComputationError

from lsst.ip.diffim.metrics import \
    NumberSciSourcesMetricTask, \
    FractionDiaSourcesToSciSourcesMetricTask


def _makeDummyCatalog(size, skyFlag=False, priFlag=False):
    """Create a trivial catalog for testing source counts.

    Parameters
    ----------
    size : `int`
        The number of entries in the catalog.
    skyFlag : `bool`
        If set, the schema is guaranteed to have the ``sky_source`` flag, and
        one row has it set to `True`. If not set, the ``sky_source`` flag is
        not present.
    priFlag : `bool`
        As ``skyFlag``, but for a ``detect_isPrimary`` flag.

    Returns
    -------
    catalog : `lsst.afw.table.SourceCatalog`
        A new catalog with ``size`` rows.
    """
    schema = SourceCatalog.Table.makeMinimalSchema()
    if skyFlag:
        schema.addField("sky_source", type="Flag", doc="Sky source.")
    if priFlag:
        schema.addField("detect_isPrimary", type="Flag", doc="Primary source.")
    catalog = SourceCatalog(schema)
    rng = np.random.Generator(np.random.PCG64(42))
    for i in range(size):
        record = catalog.addNew()
        record[SourceCatalog.Table.getIdKey()] = i + 1  # source ID 0 not allowed
        record[SourceCatalog.Table.getCoordKey()] = SpherePoint(rng.random() * 2 * math.pi,
                                                                (rng.random() - 0.5) * math.pi,
                                                                lsst.geom.radians)
    if priFlag and size > 0:
        record["detect_isPrimary"] = True
    if skyFlag and size > 0:
        record["sky_source"] = True
    return catalog


def _makeDummyFakes(size):
    """Create a trivial fakes catalog for testing fakes exclusion.

    Parameters
    ----------
    size : `int`
        The number of entries in the catalog.

    Returns
    -------
    fakes : `pandas.DataFrame`
        A cross-matched fakes catalog containing at least an ``id`` column,
        with values consistent with the output of ``_makeDummyCatalog``.
    """
    rng = np.random.Generator(np.random.PCG64(43))
    data = {
        "fakeId": [uuid.uuid4().int & (1 << 64) - 1 for n in range(size)],
        "raJ2000": rng.random(size) * 2 * math.pi,
        "decJ2000": (rng.random(size) - 0.5) * math.pi,
        "id": range(1, size + 1),  # source ID 0 not allowed
    }
    return pandas.DataFrame(data)


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
        assert_quantity_allclose(meas.quantity, len(catalog) * u.count)

    def testEmptyCatalog(self):
        catalog = _makeDummyCatalog(0)
        result = self.task.run(catalog)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ip_diffim.numSciSources"))
        assert_quantity_allclose(meas.quantity, 0 * u.count)

    def testSkySources(self):
        catalog = _makeDummyCatalog(3, skyFlag=True)
        result = self.task.run(catalog)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ip_diffim.numSciSources"))
        assert_quantity_allclose(meas.quantity, (len(catalog) - 1) * u.count)

    def testPrimarySources(self):
        catalog = _makeDummyCatalog(3, priFlag=True)
        result = self.task.run(catalog)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ip_diffim.numSciSources"))
        assert_quantity_allclose(meas.quantity, 1 * u.count)

    def testMissingData(self):
        result = self.task.run(None)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement
        self.assertIsNone(meas)

    def testFakesRemoval(self):
        catalog = _makeDummyCatalog(3)
        config = NumberSciSourcesMetricTask.ConfigClass()
        config.removeFakes = True
        config.fakesSourceIdColumn = "id"
        task = NumberSciSourcesMetricTask(config=config)

        result1 = task.run(catalog, _makeDummyFakes(1))
        assert_quantity_allclose(result1.measurement.quantity, (len(catalog) - 1) * u.count)

        resultAll = task.run(catalog, _makeDummyFakes(len(catalog)))
        assert_quantity_allclose(resultAll.measurement.quantity, 0 * u.count)


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
        assert_quantity_allclose(meas.quantity, len(diaCatalog) / len(sciCatalog) * u.dimensionless_unscaled)

    def testEmptyDiaCatalog(self):
        sciCatalog = _makeDummyCatalog(5)
        diaCatalog = _makeDummyCatalog(0)
        result = self.task.run(sciCatalog, diaCatalog)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ip_diffim.fracDiaSourcesToSciSources"))
        assert_quantity_allclose(meas.quantity, 0.0 * u.dimensionless_unscaled)

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
        assert_quantity_allclose(meas.quantity,
                                 len(diaCatalog) / (len(sciCatalog) - 1) * u.dimensionless_unscaled)

    def testPrimarySources(self):
        sciCatalog = _makeDummyCatalog(5, skyFlag=True, priFlag=True)
        diaCatalog = _makeDummyCatalog(3)
        result = self.task.run(sciCatalog, diaCatalog)
        lsst.pipe.base.testUtils.assertValidOutput(self.task, result)
        meas = result.measurement

        self.assertEqual(meas.metric_name, Name(metric="ip_diffim.fracDiaSourcesToSciSources"))
        assert_quantity_allclose(meas.quantity, len(diaCatalog) * u.dimensionless_unscaled)

    def testFakesRemoval(self):
        sciCatalog = _makeDummyCatalog(5)
        diaCatalog = _makeDummyCatalog(3)
        config = FractionDiaSourcesToSciSourcesMetricTask.ConfigClass()
        config.removeFakes = True
        config.fakesSourceIdColumn = "id"
        task = FractionDiaSourcesToSciSourcesMetricTask(config=config)

        result1 = task.run(sciCatalog, diaCatalog, _makeDummyFakes(1))
        assert_quantity_allclose(result1.measurement.quantity,
                                 (len(diaCatalog) - 1) / (len(sciCatalog) - 1) * u.dimensionless_unscaled)

        resultAll = task.run(sciCatalog, diaCatalog, _makeDummyFakes(len(diaCatalog)))
        assert_quantity_allclose(resultAll.measurement.quantity, 0.0 * u.dimensionless_unscaled)


# Hack around unittest's hacky test setup system
del MetricTaskTestCase


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
