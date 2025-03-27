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

import lsst.afw.image
import lsst.afw.math
import lsst.afw.table
import lsst.geom
import lsst.ip.diffim
from lsst.meas.algorithms import SourceDetectionTask
from lsst.meas.base.pluginRegistry import register
from lsst.meas.base.tests import TestDataset
import lsst.utils.tests


@register("base_scoreTest")
class ScoreTestPlugin(lsst.meas.base.SingleFramePlugin):

    def __init__(self, config, name, schema, metadata):
        super().__init__(config, name, schema, metadata)
        self.key = schema.addField(name + "_value", type="F", doc="test value")
        # self.keyY = schema.addField(name + "_y", type="D", doc="peak centroid", units="pixel")

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def measure(self, record, exposure, kernel=None):
        center = record.getCentroid()
        flux = exposure.image[center]
        record.set(self.key, flux)


import lsst.afw.display


class SingleFrameScoreMeasurementTest(lsst.utils.tests.TestCase):
    """Test that SingleFrameScoreMeasurementTask properly deconvolves
    the image.
    """
    def setUp(self):
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(5, 4), lsst.geom.Point2I(155, 125))
        # TODO: psfDim here shouldn't have to depend on the image cutout size.
        dataset = TestDataset(bbox, psfSigma=4.0, psfDim=41)
        # two sources, separated by about twice our PSF image size.
        # TODO: have to separate them more for now, because the convolve below
        # grows them and detection doesn't deblend them.
        dataset.addSource(1e4, lsst.geom.Point2D(50, 50))
        dataset.addSource(1e5, lsst.geom.Point2D(100, 50))
        self.exposure, self.catalog = dataset.realize(2.0, dataset.makeMinimalSchema())

        convolutionControl = lsst.afw.math.ConvolutionControl()
        # convolutionControl.setDoNormalize(False)
        # convolutionControl.setDoCopyEdge(True)
        # import os; print(os.getpid()); import ipdb; ipdb.set_trace();
        # TODO: should this be a double or float, for reduction of numerical issues?
        self.score = lsst.afw.image.ExposureF(self.exposure, deep=True)
        lsst.afw.math.convolve(self.score.maskedImage,
                               self.exposure.maskedImage,
                               self.exposure.psf.getLocalKernel(self.catalog[0].getCentroid()),
                               convolutionControl)
        self.score.mask.clearMaskPlane(self.score.mask.getMaskPlane("DETECTED"))

    def test_callMeasure(self):
        schema = lsst.afw.table.SourceTable.makeMinimalSchema()

        config = lsst.ip.diffim.SingleFrameScoreMeasurementTask.ConfigClass()
        config.plugins = ["base_scoreTest", "base_SdssCentroid"]
        config.slots.apFlux = None
        config.slots.calibFlux = None
        config.slots.gaussianFlux = None
        config.slots.modelFlux = None
        config.slots.psfFlux = None
        config.slots.shape = None
        config.slots.psfShape = None
        task = lsst.ip.diffim.SingleFrameScoreMeasurementTask(schema=schema,
                                                              config=config)

        config = SourceDetectionTask.ConfigClass()
        config.nSigmaToGrow = 0
        detection = SourceDetectionTask(config=config)
        result = detection.run(lsst.afw.table.SourceCatalog(schema), self.score, doSmooth=False)
        catalog = result.sources

        # catalog = lsst.afw.table.SourceCatalog(schema)
        # catalog.addNew()
        # catalog[-1].setFootprint(self.catalog[0].getFootprint())
        # catalog.addNew()
        # catalog[-1].setFootprint(self.catalog[1].getFootprint())

        display = lsst.afw.display.Display()
        display.frame = 1
        display.image(self.exposure, title="exposure")
        display.frame = 2
        display.image(self.score, title="score")
        # display.centroids(catalog, size=10, ctype="red", symbol="x")
        # for x in catalog:
        #     display.dot(x['id'], x.getX(), x.getY(), size=10, ctype="cyan")

        for record in catalog:
            task.callMeasure(record, self.score, kernel=self.exposure.psf)

        for record in catalog:
            with self.subTest(f"record {record['id']} at {record.getCentroid()}"):
                self.assertFloatsAlmostEqual(record["base_scoreTest_value"],
                                             self.exposure.image[record.getCentroid()])

        # import ipdb; ipdb.set_trace();


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
