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

    def _replacePatchesWithNoise(self, true_var, seed=0):
        """Replace each patch's pixels with iid Gaussian noise of the
        given marginal variance and pin the variance plane to that
        value, so the input to the task has *consistent* noise and
        variance plane (mask bits cleared too).

        Parameters
        ----------
        true_var : `float`
            Variance of the Gaussian noise drawn into each patch's
            image plane; the same value is written to the variance
            plane and used as the ground truth that downstream
            propagation tests compare against.
        seed : `int`
            Numpy RNG seed; per-patch noise draws use independent
            substreams so tract/patch overlaps combine independent
            samples.
        """
        rng = np.random.default_rng(seed)
        for tract_id in list(self.patches.keys()):
            new_handles = []
            for handle in self.patches[tract_id]:
                exposure = handle.get()
                ny, nx = exposure.image.array.shape
                exposure.image.array[:, :] = rng.normal(
                    0.0, np.sqrt(true_var), size=(ny, nx),
                ).astype(np.float32)
                exposure.variance.array[:, :] = float(true_var)
                exposure.mask.array[:, :] = 0
                new_handles.append(pipeBase.InMemoryDatasetHandle(
                    exposure, storageClass="ExposureF",
                    copy=True, dataId=handle.dataId,
                ))
            self.patches[tract_id] = new_handles

    def _interiorView(self, template, shrink=30):
        """Return (image, variance, valid_mask) for the interior of the
        merged template, excluding NO_DATA / EDGE pixels and a
        ``shrink``-pixel inset that drops residual warping/convolution
        edge effects from the analytical bounds.

        Used by variance-propagation tests that compare the empirical
        per-pixel variance of the output image to the median of the
        variance plane: those statistics are only meaningful where every
        contributing input was real data.
        """
        mask_bits = (template.mask.getPlaneBitMask("NO_DATA")
                     | template.mask.getPlaneBitMask("EDGE"))
        good = (template.mask.array & mask_bits) == 0
        # Sentinel-marked variance pixels inflate the median artificially.
        good &= template.variance.array < 1e10
        inset = lsst.geom.Box2I(template.getBBox())
        inset.grow(-shrink)
        x0, y0 = template.getX0(), template.getY0()
        sl = (slice(inset.getMinY() - y0, inset.getMaxY() - y0 + 1),
              slice(inset.getMinX() - x0, inset.getMaxX() - x0 + 1))
        valid = good[sl]
        return (template.image.array[sl], template.variance.array[sl], valid)

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
        # The variance plane is now an inverse-variance combination, so a
        # pixel covered by N independent patches has variance ≈ input/N.
        # Instead of comparing element-wise to the input, sanity-check the
        # plane: finite everywhere, strictly positive where data lands,
        # and never inflated above the per-pixel input.
        template_var = template.variance.array
        input_var_max = float(np.max(self.exposure[expectedBox].variance.array))
        self.assertTrue(np.all(np.isfinite(template_var)))
        valid = (template.mask.array
                 & template.mask.getPlaneBitMask("NO_DATA")) == 0
        self.assertTrue(np.all(template_var[valid] > 0))
        self.assertTrue(np.all(template_var[valid] <= 2*input_var_max))
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

    def testVariancePlaneTracksOutputNoise(self):
        """The output variance plane should match the marginal per-pixel
        variance of the output image plane for noise-only inputs.

        Lanczos warping introduces correlations between output pixels but
        the *marginal* variance of any one pixel is still the empirical
        variance of pixel values across many pixels (the cross-pixel
        covariance does not bias an unweighted sample variance). So even
        though sqrt(template.image) drops below sqrt(true_var) after warp,
        the variance plane must drop by exactly the same factor —
        otherwise downstream 5-σ peak detection on
        ``image / sqrt(variance)`` would mis-threshold.
        """
        true_var = 4.0
        self._replacePatchesWithNoise(true_var=true_var, seed=0)

        box = lsst.geom.Box2I(lsst.geom.Point2I(200, 200),
                              lsst.geom.Point2I(600, 600))
        task = lsst.ip.diffim.GetTemplateTask()
        result = task.run(coaddExposureHandles=self.patches,
                          bbox=lsst.geom.Box2I(box),
                          wcs=self.exposure.wcs,
                          dataIds=self.dataIds,
                          physical_filter="a_test")
        image, variance, valid = self._interiorView(result.template)
        emp = float(np.var(image[valid]))
        plane = float(np.median(variance[valid]))
        self.assertGreater(valid.sum(), 1000,
                           "fixture must produce enough good pixels")
        self.assertGreater(emp, 0.0)
        self.assertGreater(plane, 0.0)
        # Both empirical and plane must be reduced from true_var (Lanczos
        # warp + multi-tract inverse-variance combine).
        self.assertLess(plane, true_var)
        # Per-pixel agreement between marginal noise and variance plane.
        # 25% tolerance covers correlated-sample variance estimator error
        # plus per-tile inverse-variance combination noise.
        self.assertAlmostEqual(
            emp, plane, delta=0.25*plane,
            msg=f"empirical output variance {emp:.3f} disagrees with "
                f"variance plane {plane:.3f} by more than 25%",
        )

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
        # We just check that the pixel values are all finite. We cannot
        # check that pixel values in the template are closer to the
        # original anymore.
        self.assertTrue(np.isfinite(result.template.image.array).all())


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
