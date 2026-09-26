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

import numpy as np

import lsst.afw.image
import lsst.afw.table
import lsst.ip.diffim
import lsst.utils.tests

from utils import makeTestImage


class ComputeSpatiallySampledMetricsTest(lsst.utils.tests.TestCase):
    """Test the metrics that are sampled across the difference image."""

    def setUp(self):
        noiseLevel = 1.
        science, sources = makeTestImage(psfSize=2.4, noiseLevel=noiseLevel, noiseSeed=6)
        template, _ = makeTestImage(psfSize=2.0, noiseLevel=noiseLevel, noiseSeed=7,
                                    templateBorderSize=20, doApplyCalibration=True)
        config = lsst.ip.diffim.AlardLuptonSubtractTask.ConfigClass()
        config.doSubtractBackground = False
        config.sourceSelector.signalToNoise.fluxField = "truth_instFlux"
        config.sourceSelector.signalToNoise.errField = "truth_instFluxErr"
        subtraction = lsst.ip.diffim.AlardLuptonSubtractTask(config=config).run(template, science, sources)
        self.science = subtraction.matchedScience
        self.template = template
        self.difference = subtraction.difference
        self.psfMatchingKernel = subtraction.psfMatchingKernel
        self.diaSources = self._makeDiaSources()

    @staticmethod
    def _makeDiaSources():
        """Return an empty catalog with the fields the metrics read from the
        detected sources.
        """
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        schema.addField("x", "D", "Centroid x.")
        schema.addField("y", "D", "Centroid y.")
        schema.addField("isDipole", "Flag", "Is this source a dipole?")
        schema.addField("dipoleAngle", "D", "Dipole orientation.")
        schema.addField("dipoleLength", "D", "Dipole separation.")
        return lsst.afw.table.SourceCatalog(schema)

    def _run(self, **kwargs):
        """Run the metrics task on the images made in `setUp`."""
        config = lsst.ip.diffim.SpatiallySampledMetricsTask.ConfigClass()
        config.update(**kwargs)
        task = lsst.ip.diffim.SpatiallySampledMetricsTask(config=config)
        return task.run(self.science, self.template, self.difference, self.diaSources,
                        self.psfMatchingKernel).spatiallySampledMetrics

    def testMaskPlaneNotInImage(self):
        """A configured mask plane that is not registered has a fraction of
        zero, rather than raising.

        An image converted from `lsst.images.DifferenceImage` only defines the
        planes that have pixels set, so INJECTED is missing from an ordinary
        AP run.
        """
        present = "DETECTED"
        missing = "NOT_A_MASK_PLANE"
        self.assertIn(present, self.difference.mask.getMaskPlaneDict())
        self.assertNotIn(missing, self.difference.mask.getMaskPlaneDict())

        metrics = self._run(metricsMaskPlanes=[present, missing])

        self.assertGreater(len(metrics), 0)
        self.assertTrue(np.all(metrics[f"{missing.lower()}_mask_fraction"] == 0))
        self.assertTrue(np.all(np.isfinite(metrics[f"{present.lower()}_mask_fraction"])))


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
