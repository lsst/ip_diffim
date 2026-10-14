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

import types
import unittest

import astropy.units as u
import numpy as np

import lsst.afw.cameraGeom.testUtils
import lsst.afw.image
import lsst.afw.table
from lsst.daf.butler import DimensionUniverse
import lsst.ip.diffim
import lsst.utils.tests
from lsst.images import DifferenceImage, VisitImage
from lsst.images.convolution_kernels import ImageBasisConvolutionKernel

from utils import makeTestExposureRecord, makeTestImage, makeTestVisitInfo


class ComputeSpatiallySampledMetricsTest(lsst.utils.tests.TestCase):
    """Test the metrics that are sampled across the difference image."""

    def setUp(self):
        noiseLevel = 1.
        science, sources = makeTestImage(psfSize=2.4, noiseLevel=noiseLevel, noiseSeed=6)
        template, _ = makeTestImage(psfSize=2.0, noiseLevel=noiseLevel, noiseSeed=7,
                                    templateBorderSize=20, doApplyCalibration=True)
        config = lsst.ip.diffim.AlardLuptonSubtractTask.ConfigClass()
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

    def testConnections(self):
        """``image_type`` sets the storage class of the image and kernel
        inputs.
        """
        expected = {
            "legacy": {"science": "ExposureF", "template": "ExposureF", "difference": "ExposureF",
                       "psfMatchingKernel": "MatchingKernel"},
            "future": {"science": "VisitImage", "template": "DifferenceImage",
                       "difference": "DifferenceImage", "psfMatchingKernel": "ConvolutionKernel"},
        }
        for imageType, storageClasses in expected.items():
            with self.subTest(imageType=imageType):
                config = lsst.ip.diffim.SpatiallySampledMetricsTask.ConfigClass()
                config.image_type = imageType
                connections = config.ConnectionsClass(config=config)
                for name, storageClass in storageClasses.items():
                    self.assertEqual(getattr(connections, name).storageClass, storageClass)

    def _makeFutureInputs(self):
        """Convert the images and kernel made in `setUp` to the types the
        ``future`` mode reads.

        Returns
        -------
        inputs : `dict` [`str`, `object`]
            The converted inputs, keyed by connection name.
        """
        detector = list(lsst.afw.cameraGeom.testUtils.CameraWrapper().camera)[0]
        visit = 12345
        exposureRecord = makeTestExposureRecord(DimensionUniverse(), visit=visit)
        inputs = {}
        for name, imageType in (("science", VisitImage), ("template", DifferenceImage),
                                ("difference", DifferenceImage)):
            exposure = getattr(self, name).clone()
            exposure.setDetector(detector)
            exposure.info.setVisitInfo(makeTestVisitInfo(visit))
            inputs[name] = imageType.from_legacy(exposure, unit=u.nJy, exposure_record=exposureRecord)
        inputs["psfMatchingKernel"] = ImageBasisConvolutionKernel.from_legacy(self.psfMatchingKernel)
        return inputs

    def testRunQuantumFuture(self):
        """``runQuantum`` converts the future inputs back to the legacy types,
        giving the same metrics as a legacy run.
        """
        legacyInputs = {"science": self.science, "template": self.template, "difference": self.difference,
                        "psfMatchingKernel": self.psfMatchingKernel}
        metrics = {}
        for imageType, inputs in (("legacy", legacyInputs), ("future", self._makeFutureInputs())):
            config = lsst.ip.diffim.SpatiallySampledMetricsTask.ConfigClass()
            config.image_type = imageType
            task = lsst.ip.diffim.SpatiallySampledMetricsTask(config=config)
            butlerQC = _RecordingQuantumContext()
            task.runQuantum(butlerQC, types.SimpleNamespace(**inputs, diaSources=self.diaSources), None)
            metrics[imageType] = butlerQC.put_values.spatiallySampledMetrics

        self.assertGreater(len(metrics["legacy"]), 0)
        self.assertEqual(metrics["future"].colnames, metrics["legacy"].colnames)
        for name in metrics["legacy"].colnames:
            np.testing.assert_array_equal(metrics["future"][name], metrics["legacy"][name], err_msg=name)


class _RecordingQuantumContext:
    """Minimal `~lsst.pipe.base.QuantumContext` that records what was put."""
    def __init__(self):
        self.put_values = None

    def get(self, refs):
        return dict(refs.__dict__)

    def put(self, values, refs):
        self.put_values = values


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
