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

import itertools
import unittest

import numpy as np

import lsst.afw.geom
import lsst.afw.image
import lsst.afw.math
import lsst.geom
import lsst.ip.diffim
import lsst.meas.algorithms
import lsst.meas.base.tests
import lsst.skymap
import lsst.utils.tests

import lsst.afw.display
display = lsst.afw.display.Display()


class GetTemplateTaskTestCase(lsst.utils.tests.TestCase):
    """Checks that GetTemplateTask works on both one tract and multiple tract
    input coadd exposures.

    Makes a synthetic exposure large enough to fit four small tracts with 2x2
    (300x300 pixel) patches each, extracts pixels for thos patches by warping,
    and tests GetTemplateTask against a box that fits in one tract and a box
    that overlaps multiple.
    """
    def setUp(self):
        self.scale = 0.2
        self.skymap = self._makeSkymap()
        self.patches = {}
        self.dataIds = {}
        self.exposure = self._makeExposure()
        self.copy = self.exposure.clone()
        # display.image(self.exposure)
        # print("!!!!!!!")
        # for tract_id in range(4):
        tract = self.skymap.generateTract(0)
        self._makePatches(tract)
        # import os; print(os.getpid()); import ipdb; ipdb.set_trace();

        # display.frame = 0
        # for patch in self.patches:
        #     print(patch)
        #     display.image(self.patches[patch], title=patch)
        #     display.frame += 1
        # display.image(self.exposure)
        # import os; print(os.getpid()); import ipdb; ipdb.set_trace();

    def _makeSkymap(self):
        """Make a Skymap with 4 tracts with 4 patches each.
        """
        tractScale = 0.01
        coords = [(0, 0),
                  (0, tractScale),
                  (tractScale, 0),
                  (tractScale, tractScale),
                  ]
        config = lsst.skymap.DiscreteSkyMap.ConfigClass()
        config.raList = [c[0] for c in coords]
        config.decList = [c[1] for c in coords]
        config.radiusList = [tractScale for c in coords]
        config.projection = "TAN"
        config.pixelScale = self.scale
        config.tractOverlap = 0.001
        config.tractBuilder = "legacy"
        config.tractBuilder["legacy"].patchInnerDimensions = (300, 300)
        config.tractBuilder["legacy"].patchBorder = 10
        # config.tractBuilder = "cells"
        # config.tractBuilder["cells"].numCellsPerPatchInner = 11
        # import os; print(os.getpid()); import ipdb; ipdb.set_trace();
        return lsst.skymap.DiscreteSkyMap(config=config)

    def _makeExposure(self):
        """Create a large image to break up into tracts and patches.

        The image will have a source every 100 pixels in x and y, and a
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(-100, -100), lsst.geom.Point2I(700, 700))
        # Use a WCS close to that of tract 0, so that all the tracts fit.
        cd_matrix = lsst.afw.geom.makeCdMatrix(self.scale*1.05*lsst.geom.arcseconds, 50*lsst.geom.arcseconds)
        wcs = lsst.afw.geom.makeSkyWcs(lsst.geom.Point2D(200, 200),
                                       lsst.geom.SpherePoint(0, 0, lsst.geom.radians),
                                       cd_matrix)
        dataset = lsst.meas.base.tests.TestDataset(box, wcs=wcs)
        for x, y in itertools.product(np.arange(0, 500, 100), np.arange(0, 500, 100)):
            dataset.addSource(1e5, lsst.geom.Point2D(x, y))
        exposure, _ = dataset.realize(1, dataset.makeMinimalSchema())
        exposure.setFilter(lsst.afw.image.FilterLabel("a", "a_test"))
        return exposure

    def _makePatches(self, tract):
        """Fill the dicts of (tract_id, patch_id)->exposure, with the
        exposures being deep copied subsets of the main exposure.
        """
        # print(tract.tract_id, tract.bbox)
        # print(tract.wcs)
        config = lsst.afw.math.Warper.ConfigClass()
        # Use 5th order to minimize artifacts on the templates.
        config.warpingKernelName = "lanczos5"
        warper = lsst.afw.math.Warper.fromConfig(config)
        # import os; print(os.getpid()); import ipdb; ipdb.set_trace();
        # color = ['red', 'green', 'cyan', 'yellow'][tract.tract_id]
        for patch_id in range(tract.num_patches.x*tract.num_patches.y):
            patch = tract.getPatchInfo(patch_id)
            box = patch.getOuterBBox()
            # print(box)
            # points = self.exposure.wcs.skyToPixel(patch.wcs.pixelToSky([lsst.geom.Point2D(x) for x in box.getCorners()]))
            # for p in points:
            #     display.dot(patch_id, p.x, p.y, ctype=color)
            #     print(p)
            # subset = self.exposure[box].clone()

            # This is mostly stolen from pipe_tasks warpAndPsfMatch, but
            # ip_diffim cannot depend on pipe_tasks.
            xyTransform = lsst.afw.geom.makeWcsPairTransform(self.exposure.wcs, patch.wcs)
            warpedPsf = lsst.meas.algorithms.WarpedPsf(self.exposure.psf, xyTransform)
            warped = warper.warpExposure(patch.wcs, self.exposure, destBBox=box)
            warped.setPsf(warpedPsf)
            # psfMatch = lsst.ip.diffim.ModelPsfMatchTask()
            # warpedAndMatched = psfMatch.run(warped, self.exposure.psf).psfMatchedExposure
            self.patches[(tract.tract_id, patch_id)] = warped
            self.dataIds[(tract.tract_id, patch_id)] = {"tract": tract.tract_id,
                                                        "patch": patch_id}

    def _checkMetadata(self, template, config, box, wcs):
        """Check that the various metadata components were set correctly.
        """
        expectedBox = lsst.geom.Box2I(box)
        expectedBox.grow(config.templateBorderSize)
        self.assertEqual(template.getBBox(), expectedBox)
        # WCS should be the distorted one above, not the input exposure.
        # self.assertEqual(template.wcs, wcs)
        # self.assertNotEqual(template.wcs, self.exposure.wcs)
        self.assertEqual(template.photoCalib, self.exposure.photoCalib)
        self.assertEqual(template.getXY0(), expectedBox.getMin())
        self.assertEqual(template.filter.bandLabel, "a")
        self.assertEqual(template.filter.physicalLabel, "a_test")
        # TODO: can I check something better on the psf?
        self.assertIsInstance(template.psf, lsst.meas.algorithms.CoaddPsf)
        # TOOD: need other things to test here!

    def _checkPixels(self, template, config, box):
        # All pixels should have real values!
        expectedBox = lsst.geom.Box2I(box)
        expectedBox.grow(config.templateBorderSize)
        display.frame = 0
        display.image(self.exposure, title="exposure")
        display.frame += 1
        display.image(template, title="template")
        image = template.clone()
        image.maskedImage -= self.exposure.maskedImage[expectedBox]
        display.frame += 1
        display.image(image, title="difference")
        import os; print(os.getpid()); import ipdb; ipdb.set_trace();
        # Check that we fully filled the template from the patches.
        self.assertTrue(np.all(np.isfinite(template.image.array)))
        # Because of the scale changes, there will be small ringing in the
        # difference between the template and the original image.
        self.assertImagesAlmostEqual(template.image, self.exposure[expectedBox].image, atol=5)
        # Variance plane ==1 in the original image, but the warped images will
        # have some structure due to the warping.
        self.assertImagesAlmostEqual(template.variance, self.exposure[expectedBox].variance, atol=0.5)
        # Not checking the mask, as warping changes the sizes of the masks.

    def testRunSameTract(self):
        """Test a bounding box that fully fits inside one tract.
        """
        # A WCS with a slightly larger pixel scale to distort the coadds to.
        # distortion = lsst.afw.geom.makeRadialTransform([0, 1.1])
        # wcs = lsst.afw.geom.makeModifiedWcs(distortion, self.skymap.generateTract(0).wcs, False)
        box = lsst.geom.Box2I(lsst.geom.Point2I(50, 50), lsst.geom.Point2I(200, 200))
        task = lsst.ip.diffim.GetTemplateTask()
        # Task modifies the input bbox, so pass a copy.
        result = task.run(list(self.patches.values()), lsst.geom.Box2I(box), self.exposure.wcs,
                          list(self.dataIds.values()))

        self._checkMetadata(result.template, task.config, box, self.exposure.wcs)
        self._checkPixels(result.template, task.config, box)

    def testRunTwoTracts(self):
        """Test a bounding box that crosses one tract boundary.
        """
        # A WCS with a slightly larger pixel scale to distort the coadds to.
        # distortion = lsst.afw.geom.makeRadialTransform([0, 1.1])
        # wcs = lsst.afw.geom.makeModifiedWcs(distortion, self.skymap.generateTract(0).wcs, False)
        box = lsst.geom.Box2I(lsst.geom.Point2I(100, 200), lsst.geom.Point2I(400, 450))
        task = lsst.ip.diffim.GetTemplateTask()
        # TODO: something should fail here!
        result = task.run(list(self.patches.values()), lsst.geom.Box2I(box), self.exposure.wcs,
                          list(self.dataIds.values()))

        self._checkMetadata(result.template, task.config, box, self.exposure.wcs)
        self._checkPixels(result.template, task.config, box)


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
