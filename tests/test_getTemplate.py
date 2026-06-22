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
import lsst.afw.table
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

    def testComputeDepthMap(self):
        """Test that _computeDepthMap counts the distinct input images covering
        each pixel and deduplicates inputs on (visit, ccd).
        """
        wcs = lsst.afw.geom.makeSkyWcs(
            lsst.geom.Point2D(50, 50),
            lsst.geom.SpherePoint(0, 0, lsst.geom.degrees),
            lsst.afw.geom.makeCdMatrix(0.2*lsst.geom.arcseconds, 0*lsst.geom.degrees))
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(100, 100))
        fullBox = lsst.geom.Box2D(box)
        leftBox = lsst.geom.Box2D(
            lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(50, 100)))

        # The left half is covered by visits {1, 2, 3} and the right half by {1}.
        # The inputs share the template WCS, so their footprints map directly
        # into the template frame.
        schema = lsst.afw.table.ExposureTable.makeMinimalSchema()
        schema.addField("visit", type=np.int64, doc="visit id")
        schema.addField("ccd", type=np.int32, doc="ccd id")
        ccds = lsst.afw.table.ExposureCatalog(schema)

        def addCcd(visit, ccd, polyBox):
            record = ccds.addNew()
            record.setWcs(wcs)
            record.setBBox(lsst.geom.Box2I(polyBox))
            record.setValidPolygon(lsst.afw.geom.Polygon(polyBox.getCorners()))
            record["visit"] = visit
            record["ccd"] = ccd

        addCcd(1, 1, fullBox)
        addCcd(2, 1, leftBox)
        addCcd(3, 1, leftBox)
        # A duplicate physical image (same visit, ccd) must not inflate the depth.
        addCcd(1, 1, fullBox)

        task = lsst.ip.diffim.GetTemplateTask()
        depth = task._computeDepthMap([ccds], wcs, box)
        self.assertTrue((depth[:, :50] == 3).all())
        self.assertTrue((depth[:, 50:] == 1).all())

    def testMaskShallowCoverage(self):
        """Test the cross-tract AND combination and the grow: a pixel is flagged
        only where every contributing tract is shallow.
        """
        box = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(100, 100))
        height, width = box.getHeight(), box.getWidth()
        # One tract covers the whole image but is shallow everywhere; a second,
        # deep tract covers only the right half (anyDeep). So the left half
        # (shallow, no deep tract) should be flagged, and the right half (filled
        # by the deep tract) should not.
        anyData = np.ones((height, width), dtype=bool)
        anyDeep = np.zeros((height, width), dtype=bool)
        anyDeep[:, 50:] = True

        def flaggedMask(growRadius):
            template = lsst.afw.image.ExposureF(box)
            config = lsst.ip.diffim.GetTemplateConfig()
            config.minNumberOfInputImages = 2
            config.coverageGrowRadius = growRadius
            task = lsst.ip.diffim.GetTemplateTask(config=config)
            task.maskShallowCoverage(template, anyData, anyDeep)
            return (template.mask.array
                    & template.mask.getPlaneBitMask("HIGH_VARIANCE")) != 0

        # With no growth, the shallow-only (left) half is flagged and the half a
        # deep tract fills is left clean.
        flagged = flaggedMask(0)
        self.assertTrue(flagged[:, :50].all())
        self.assertEqual(flagged[:, 50:].sum(), 0)

        # Growing by 5 expands the masked region by 5 columns into the deep half.
        self.assertEqual(flaggedMask(5).any(axis=0).sum(), 55)

    def testLabelConstantDepthRegions(self):
        """Connected regions of constant depth get distinct labels, including
        two disconnected regions that share the same depth value.
        """
        depth = np.zeros((10, 30), dtype=np.int32)
        depth[:, 0:10] = 3
        depth[:, 10:20] = 5
        depth[:, 20:30] = 3
        hasData = np.ones(depth.shape, dtype=bool)
        labels = lsst.ip.diffim.GetTemplateTask._labelConstantDepthRegions(depth, hasData)
        # Two depth-3 blocks (disconnected) and one depth-5 block -> 3 regions.
        self.assertEqual(len(np.unique(labels)), 3)
        self.assertNotEqual(labels[5, 5], labels[5, 25])   # same depth, disconnected
        self.assertNotEqual(labels[5, 5], labels[5, 15])   # different depth

    def testResolveNarrowRegions(self):
        """Narrow shallow regions are masked; narrow deep regions are merged
        into their largest neighbor.
        """
        task = lsst.ip.diffim.GetTemplateTask()  # minRegionWidth=10, minNumberOfInputImages=2
        depth = np.full((100, 100), 5, dtype=np.int32)
        depth[:, 0:4] = 8       # narrow deep strip at the left edge (one neighbor)
        depth[:, 50:54] = 1     # narrow shallow strip, splits the depth-5 background
        hasData = np.ones(depth.shape, dtype=bool)
        labels = task._labelConstantDepthRegions(depth, hasData)
        maskNarrow = task._resolveNarrowRegions(labels, depth)

        # The narrow shallow strip is masked and removed from the labels.
        self.assertEqual(maskNarrow.sum(), 4*100)
        self.assertTrue(maskNarrow[:, 50:54].all())
        self.assertTrue((labels[:, 50:54] == 0).all())
        # The narrow deep strip is merged into its neighbor and not masked.
        self.assertFalse(maskNarrow[:, 0:4].any())
        self.assertTrue((labels[:, 0:4] == labels[5, 10]).all())

    def testFlagDiscontinuousRegions(self):
        """A region whose PSF size deviates from the dominant model is flagged,
        while a smooth gradient consistent with the polynomial order is not.
        """
        task = lsst.ip.diffim.GetTemplateTask()  # order 1, threshold 0.05

        # Two large regions at radius 5 and a small deviant region at radius 8.
        x1 = np.array([-0.8, -0.6, -0.4, -0.8, -0.6, -0.4])
        y1 = np.array([-0.5, -0.5, -0.5, 0.5, 0.5, 0.5])
        x2 = np.array([0.4, 0.6, 0.8, 0.4, 0.6, 0.8])
        y2 = y1
        x3 = np.array([0.9, 0.85])
        y3 = np.array([0.9, 0.95])
        samples = {
            "x": np.concatenate([x1, x2, x3]),
            "y": np.concatenate([y1, y2, y3]),
            "radius": np.concatenate([np.full(6, 5.0), np.full(6, 5.0), np.full(2, 8.0)]),
            "weight": np.concatenate([np.full(6, 4e4), np.full(6, 4e4), np.full(2, 1e3)]),
            "label": np.concatenate([np.full(6, 1), np.full(6, 2), np.full(2, 3)]),
        }
        self.assertEqual(task._flagDiscontinuousRegions(samples), {3})

        # A linear gradient in PSF size is fully absorbed by the order-1 fit.
        xs = np.linspace(-1, 1, 6)
        gradient = {
            "x": xs,
            "y": np.zeros(6),
            "radius": 5 + 2*xs,
            "weight": np.full(6, 4e4),
            "label": np.arange(1, 7),
        }
        self.assertEqual(task._flagDiscontinuousRegions(gradient), set())

    def testMaskPsfDiscontinuity(self):
        """End-to-end: a spatially constant PSF over two depth regions is not
        flagged, while masking that leaves too few usable pixels raises.
        """
        wcs = lsst.afw.geom.makeSkyWcs(
            lsst.geom.Point2D(100, 100),
            lsst.geom.SpherePoint(0, 0, lsst.geom.degrees),
            lsst.afw.geom.makeCdMatrix(0.2*lsst.geom.arcseconds, 0*lsst.geom.degrees))
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(200, 200))

        def makeCcds(boxes):
            schema = lsst.afw.table.ExposureTable.makeMinimalSchema()
            schema.addField("visit", type=np.int64, doc="visit id")
            schema.addField("ccd", type=np.int32, doc="ccd id")
            ccds = lsst.afw.table.ExposureCatalog(schema)
            for visit, ccd, box in boxes:
                record = ccds.addNew()
                record.setWcs(wcs)
                record.setBBox(box)
                record.setValidPolygon(lsst.afw.geom.Polygon(lsst.geom.Box2D(box).getCorners()))
                record["visit"] = visit
                record["ccd"] = ccd
            return ccds

        def makeTemplate():
            template = lsst.afw.image.ExposureF(bbox)
            template.setWcs(wcs)
            template.setPsf(lsst.meas.algorithms.SingleGaussianPsf(25, 25, 2.0))
            return template

        leftHalf = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(100, 200))
        # Depth 2 on the left, 1 on the right, but a spatially constant PSF: no
        # discontinuity, so nothing should be flagged.
        template = makeTemplate()
        lsst.ip.diffim.GetTemplateTask().maskPsfDiscontinuity(
            template, [makeCcds([(1, 1, bbox), (2, 1, leftHalf)])])
        if "PSF_DISCONTINUITY" in template.mask.getMaskPlaneDict():
            bit = template.mask.getPlaneBitMask("PSF_DISCONTINUITY")
            self.assertEqual(np.count_nonzero(template.mask.array & bit), 0)

        # A narrow shallow strip (cols 196-199, depth 1) is masked; requiring a
        # high usable fraction then rejects the template.
        almostFull = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(196, 200))
        config = lsst.ip.diffim.GetTemplateConfig()
        config.psfDiscontinuityMinUsableFraction = 0.99
        with self.assertRaises(lsst.ip.diffim.PsfDiscontinuityError):
            lsst.ip.diffim.GetTemplateTask(config=config).maskPsfDiscontinuity(
                makeTemplate(), [makeCcds([(1, 1, almostFull), (2, 1, bbox)])])


def setup_module(module):
    lsst.utils.tests.init()


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
