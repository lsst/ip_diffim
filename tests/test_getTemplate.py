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

import collections
import itertools
import unittest

import numpy as np

import lsst.afw.geom
import lsst.afw.image
import lsst.afw.math
from lsst.daf.butler import DataCoordinate, DimensionUniverse
import lsst.geom
import lsst.ip.diffim
import lsst.meas.algorithms
import lsst.meas.base.tests
import lsst.pipe.base as pipeBase
import lsst.skymap
import lsst.utils.tests

from utils import generate_data_id

# Change this to True, `setup display_ds9`, and open ds9 (or use another afw
# display backend) to show the tract/patch layouts on the image.
debug = False
if debug:
    import lsst.afw.display
    display = lsst.afw.display.Display()
    display.frame = 1


def _showTemplate(box, template):
    """Show the corners of the template we made in this test."""
    for point in box.getCorners():
        display.dot("+", point.x, point.y, ctype="orange", size=40)
    display.frame = 2
    display.image(template, "warped template")
    display.frame = 3
    display.image(template.variance, "warped variance")


class GetTemplateTaskTestCase(lsst.utils.tests.TestCase):
    """Test that GetTemplateTask works on both one tract and multiple tract
    input coadd exposures.

    Makes a synthetic exposure large enough to fit four small tracts with 2x2
    (300x300 pixel) patches each, extracts pixels for those patches by warping,
    and tests GetTemplateTask's output against boxes that overlap various
    combinations of one or multiple tracts.
    """
    def setUp(self):
        self.scale = 0.2  # arcsec/pixel
        self.skymap = self._makeSkymap()
        self.patches = collections.defaultdict(list)
        self.dataIds = collections.defaultdict(list)
        self.exposure = self._makeExposure()

        if debug:
            display.image(self.exposure, "base exposure")

        for tract_id in range(4):
            tract = self.skymap.generateTract(tract_id)
            self._makePatches(tract)

    def _makeSkymap(self):
        """Make a Skymap with 4 tracts with 4 patches each.
        """
        tractScale = 0.02  # degrees
        # On-sky coordinates of the tract centers.
        coords = [(0, 0),
                  (0, tractScale),
                  (tractScale, 0),
                  (tractScale, tractScale),
                  ]
        config = lsst.skymap.DiscreteSkyMap.ConfigClass()
        config.raList = [c[0] for c in coords]
        config.decList = [c[1] for c in coords]
        # Half the tract center step size, to keep the tract overlap small.
        config.radiusList = [tractScale/2 for c in coords]
        config.projection = "TAN"
        config.pixelScale = self.scale
        config.tractOverlap = 0.0005
        config.tractBuilder = "legacy"
        config.tractBuilder["legacy"].patchInnerDimensions = (300, 300)
        config.tractBuilder["legacy"].patchBorder = 10
        return lsst.skymap.DiscreteSkyMap(config=config)

    def _makeExposure(self):
        """Create a large image to break up into tracts and patches.

        The image will have a source every 100 pixels in x and y, and a WCS
        that results in the tracts all fitting in the image, with tract=0
        in the lower left, tract=1 to the right, tract=2 above, and tract=3
        to the upper right.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(-200, -200), lsst.geom.Point2I(800, 800))
        # This WCS was constructed so that tract 0 mostly fills the lower left
        # quadrant of the image, and the other tracts fill the rest; slight
        # extra rotation as a check on the final warp layout, scaled by 5%
        # from the patch pixel scale.
        cd_matrix = lsst.afw.geom.makeCdMatrix(1.05*self.scale*lsst.geom.arcseconds, 93*lsst.geom.degrees)
        wcs = lsst.afw.geom.makeSkyWcs(lsst.geom.Point2D(120, 150),
                                       lsst.geom.SpherePoint(0, 0, lsst.geom.radians),
                                       cd_matrix)
        dataset = lsst.meas.base.tests.TestDataset(box, wcs=wcs)
        for x, y in itertools.product(np.arange(0, 500, 100), np.arange(0, 500, 100)):
            dataset.addSource(1e5, lsst.geom.Point2D(x, y))
        exposure, _ = dataset.realize(2, dataset.makeMinimalSchema())
        exposure.setFilter(lsst.afw.image.FilterLabel("a", "a_test"))
        return exposure

    def _makePatches(self, tract):
        """Populate the patches and dataId dicts, keyed on tract id, with the
        warps of the main exposure and minimal dataIds, respectively.
        """
        if debug:
            color = ['red', 'green', 'cyan', 'yellow'][tract.tract_id]
            point = self.exposure.wcs.skyToPixel(tract.ctr_coord)
            # Show the tract center, colored by tract id.
            display.dot("x", point.x, point.y, ctype=color, size=30)

        # Use 5th order to minimize artifacts on the templates.
        config = lsst.afw.math.Warper.ConfigClass()
        config.warpingKernelName = "lanczos5"
        warper = lsst.afw.math.Warper.fromConfig(config)
        for patchId in range(tract.num_patches.x*tract.num_patches.y):
            patch = tract.getPatchInfo(patchId)
            box = patch.getOuterBBox()

            if debug:
                # Show the patch corners as patch ids, colored by tract id.
                points = self.exposure.wcs.skyToPixel(patch.wcs.pixelToSky([lsst.geom.Point2D(x)
                                                                           for x in box.getCorners()]))
                for p in points:
                    display.dot(patchId, p.x, p.y, ctype=color)

            # This is mostly taken from drp_tasks makePsfMatchedWarp, but
            # ip_diffim cannot depend on drp_tasks.
            xyTransform = lsst.afw.geom.makeWcsPairTransform(self.exposure.wcs, patch.wcs)
            warpedPsf = lsst.meas.algorithms.WarpedPsf(self.exposure.psf, xyTransform)
            warped = warper.warpExposure(patch.wcs, self.exposure, destBBox=box)
            warped.setPsf(warpedPsf)
            dataRef = pipeBase.InMemoryDatasetHandle(
                warped,
                storageClass="ExposureF",
                copy=True,
                dataId=generate_data_id(
                    tract=tract,
                    patch=patch,
                )
            )
            self.patches[tract.tract_id].append(dataRef)
            dataCoordinate = DataCoordinate.standardize({"tract": tract.tract_id,
                                                         "patch": patchId,
                                                         "band": "a",
                                                         "skymap": "skymap"},
                                                        universe=DimensionUniverse())
            self.dataIds[tract.tract_id].append(dataCoordinate)

    def _checkMetadata(self, template, config, box, wcs, nPsfs):
        """Check that the various metadata components were set correctly.
        """
        expectedBox = lsst.geom.Box2I(box)
        expectedBox.grow(config.templateBorderSize)
        self.assertEqual(template.getBBox(), expectedBox)
        # WCS should match our exposure, not any of the coadd tracts.
        for tract in self.patches:
            self.assertNotEqual(template.wcs, self.patches[tract][0].get().wcs)
        self.assertEqual(template.wcs, self.exposure.wcs)
        self.assertEqual(template.photoCalib, self.exposure.photoCalib)
        self.assertEqual(template.getXY0(), expectedBox.getMin())
        self.assertEqual(template.filter.bandLabel, "a")
        self.assertEqual(template.filter.physicalLabel, "a_test")
        self.assertEqual(template.psf.getComponentCount(), nPsfs)
        self.assertTrue(template.getInfo().hasCoaddInputs())
        self.assertEqual(len(template.getInfo().getCoaddInputs().ccds), nPsfs)

    def _checkPixels(self, template, config, box):
        """Check that the pixel values in the template are close to the
        original image.
        """
        # All pixels should have real values!
        expectedBox = lsst.geom.Box2I(box)
        expectedBox.grow(config.templateBorderSize)

        if debug:
            _showTemplate(expectedBox, template)

        # Check that we fully filled the template from the patches.
        self.assertTrue(np.all(np.isfinite(template.image.array)))
        # Because of the scale changes, there will be some ringing in the
        # difference between the template and the original image; pick
        # tolerances large enough to account for that.
        self.assertImagesAlmostEqual(template.image, self.exposure[expectedBox].image,
                                     rtol=.1, atol=4)
        # Variance plane ==2 in the original image, but the warped images will
        # have some structure due to the warping.
        self.assertImagesAlmostEqual(template.variance, self.exposure[expectedBox].variance,
                                     rtol=0.55, msg="variance planes differ")
        # Not checking the mask, as warping changes the sizes of the masks.

    def testRunOneTractInput(self):
        """Test a bounding box that fully fits inside one tract, with only
        that tract passed as input. This checks that the code handles a single
        tract input correctly.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        task = lsst.ip.diffim.GetTemplateTask()
        # Restrict to tract 0, since the box fits in just that tract.
        # Task modifies the input bbox, so pass a copy.
        result = task.run(coaddExposureHandles={0: self.patches[0]},
                          bbox=lsst.geom.Box2I(box),
                          wcs=self.exposure.wcs,
                          dataIds={0: self.dataIds[0]},
                          physical_filter="a_test")

        # All 4 patches from tract 0 are included in this template.
        self._checkMetadata(result.template, task.config, box, self.exposure.wcs, 4)
        self._checkPixels(result.template, task.config, box)

    def testRunOneTractMultipleInputs(self):
        """Test a bounding box that fully fits inside one tract but where
        multiple tracts were passed in. This checks that patches that are
        mostly NaN after warping are merged correctly in the output.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        task = lsst.ip.diffim.GetTemplateTask()
        # Task modifies the input bbox, so pass a copy.
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")

        # All 4 patches from two tracts are included in this template.
        self._checkMetadata(result.template, task.config, box, self.exposure.wcs, 6)
        self._checkPixels(result.template, task.config, box)

    def testRunTwoTracts(self):
        """Test a bounding box that crosses tract boundaries.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(200, 200), lsst.geom.Point2I(600, 600))
        task = lsst.ip.diffim.GetTemplateTask()
        # Task modifies the input bbox, so pass a copy.
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")

        # All 4 patches from all 4 tracts are included in this template
        self._checkMetadata(result.template, task.config, box, self.exposure.wcs, 9)
        self._checkPixels(result.template, task.config, box)

    def testRunNoTemplate(self):
        """A bounding box that doesn't overlap the patches will raise.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(1200, 1200), lsst.geom.Point2I(1600, 1600))
        task = lsst.ip.diffim.GetTemplateTask()
        with self.assertRaisesRegex(lsst.pipe.base.NoWorkFound, "No patches found"):
            task.run(coaddExposureHandles=self.patches,
                     bbox=lsst.geom.Box2I(box),
                     wcs=self.exposure.wcs,
                     dataIds=self.dataIds,
                     physical_filter="a_test")

    def testMissingPatches(self):
        """Test that a missing patch results in an appropriate mask.

        This fixes the bug reported on DM-44997 (image and variance were NaN
        but the mask was not set to NO_DATA for those pixels).
        """
        # tract=0, patch=1 is the lower-left corner, as displayed in DS9.
        self.patches[0].pop(1)
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))
        task = lsst.ip.diffim.GetTemplateTask()
        # Task modifies the input bbox, so pass a copy.
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")
        no_data = (result.template.mask.array & result.template.mask.getPlaneBitMask("NO_DATA")) != 0
        self.assertTrue(np.isfinite(result.template.image.array).all())
        self.assertTrue(np.isfinite(result.template.variance.array).all())
        self.assertEqual(no_data.sum(), 20990)

    @lsst.utils.tests.methodParameters(
        box=[
            lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180)),
            lsst.geom.Box2I(lsst.geom.Point2I(200, 200), lsst.geom.Point2I(600, 600)),
        ],
        nInput=[8, 16],
    )
    def testNanInputs(self, box=None, nInput=None):
        """Test that the template has finite values when some of the input
        pixels have NaN as variance.
        """
        for tract, patchRefs in self.patches.items():
            for patchRef in patchRefs:
                patchCoadd = patchRef.get()
                bbox = lsst.geom.Box2I()
                bbox.include(lsst.geom.Point2I(patchCoadd.getBBox().getCenter()))
                bbox.grow(3)
                patchCoadd.variance[bbox].array *= np.nan

        box = lsst.geom.Box2I(lsst.geom.Point2I(200, 200), lsst.geom.Point2I(600, 600))
        task = lsst.ip.diffim.GetTemplateTask()
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")
        if debug:
            _showTemplate(box, result.template)
        self._checkMetadata(result.template, task.config, box, self.exposure.wcs, 9)
        # We just check that the pixel values are all finite. We cannot check that pixel values
        # in the template are closer to the original anymore.
        self.assertTrue(np.isfinite(result.template.image.array).all())

    def _scaleInputVariance(self, tract, factor):
        """Return fresh handles for one tract's patches, with their variance
        planes multiplied by ``factor``.

        Parameters
        ----------
        tract : `int`
            Id of the tract whose patches should be copied.
        factor : `float`
            Factor to multiply the input variance planes by.

        Returns
        -------
        handles : `list` [`lsst.pipe.base.InMemoryDatasetHandle`]
            Handles to the modified patches.
        """
        handles = []
        for handle in self.patches[tract]:
            # ``copy=True`` on the original handles means this is a copy, so
            # the patches shared with the other tests are left untouched.
            patch = handle.get()
            patch.variance.array *= factor
            handles.append(pipeBase.InMemoryDatasetHandle(patch,
                                                          storageClass="ExposureF",
                                                          copy=True,
                                                          dataId=handle.dataId))
        return handles

    def testScaleVariance(self):
        """Test that the template variance plane is rescaled to match the
        empirical pixel noise, and that the factor used is recorded in the
        task metadata.
        """
        scaleFactor = 1.345
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Point2I(180, 180))

        def _configureAndRunTask(doScaleVariance, varianceScale=1.):
            """Build a template from tract 0, optionally rescaling the input
            variance planes by ``varianceScale`` first.
            """
            config = lsst.ip.diffim.GetTemplateTask.ConfigClass()
            config.doScaleVariance = doScaleVariance
            task = lsst.ip.diffim.GetTemplateTask(config=config)
            # Task modifies the input bbox, so pass a copy.
            result = task.run(coaddExposureHandles={0: self._scaleInputVariance(0, varianceScale)},
                              bbox=lsst.geom.Box2I(box),
                              wcs=self.exposure.wcs,
                              dataIds={0: self.dataIds[0]},
                              physical_filter="a_test")
            return task, result.template

        # With scaling disabled the subtask is never constructed, and nothing
        # is recorded in the metadata.
        taskOff, templateOff = _configureAndRunTask(False)
        self.assertFalse(hasattr(taskOff, "scaleVariance"))
        self.assertNotIn("scaleTemplateVarianceFactor", taskOff.metadata)

        # Both warps -- lanczos5 in ``_makePatches`` and lanczos3 in the
        # task -- correlate the noise. The variance plane tracks only the
        # per-pixel diagonal, which the second warp leaves too low, so
        # ``scaleVariance`` measures a factor well above 1 even though the
        # input variance planes are correct.
        #
        taskOn, templateOn = _configureAndRunTask(True)
        factor = taskOn.metadata["scaleTemplateVarianceFactor"]
        # TODO DM-55879: this value is pinned on purpose. The lanczos warping
        # kernels introduce small correlations that artificially suppress the
        # image pixel stddev and inflate the variance scaling factor. This
        # should be changed to 1.0 after DM-55879 is merged.
        self.assertFloatsAlmostEqual(factor, 1.1465, atol=0.01,
                                     msg="Measured template variance scaling changed; see the"
                                         " comment above if the correlation correction landed.")
        # The only difference from the unscaled template is the constant
        # factor applied to the variance plane.
        self.assertFloatsAlmostEqual(templateOn.variance.array,
                                     templateOff.variance.array*factor, rtol=1e-5)
        # Tolerance here is float32 round-off: repeated runs of the task are
        # not bitwise identical.
        self.assertImagesAlmostEqual(templateOn.image, templateOff.image, rtol=1e-5, atol=1e-5)

        # If the input variance planes under-estimate the noise by a known
        # factor, the measured factor grows by that amount and the same
        # output variance plane is recovered.
        taskLow, templateLow = _configureAndRunTask(True, varianceScale=1/scaleFactor)
        self.assertFloatsAlmostEqual(taskLow.metadata["scaleTemplateVarianceFactor"],
                                     factor*scaleFactor, rtol=1e-5)
        self.assertImagesAlmostEqual(templateLow.variance, templateOn.variance, rtol=1e-5)


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
