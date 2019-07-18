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

import os
import unittest

import lsst.utils.tests
import lsst.utils
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.geom as geom
import lsst.ip.diffim as ipDiffim
import lsst.pex.config as pexConfig
import lsst.log.utils as logUtils
import lsst.afw.table as afwTable

logUtils.traceSetAt("ip.diffim", 4)

# known input images
try:
    defDataDir = lsst.utils.getPackageDir('afwdata')
except Exception:
    defDataDir = None

try:
    display
    defDataDir
except NameError:
    display = False
else:
    import lsst.afw.display as afwDisplay
    afwDisplay.setDefaultMaskTransparency(75)


class DiffimTestCases(lsst.utils.tests.TestCase):

    def setUp(self):
        schema = afwTable.SourceTable.makeMinimalSchema()
        afwTable.Point2DKey.addFields(schema, "Centroid", "input centroid", "pixel")
        schema.addField("PsfFlux_instFlux", type=float)
        schema.addField("PsfFlux_instFluxErr", type=float)
        schema.addField("PsfFlux_flag", type="Flag")
        self.table = afwTable.SourceTable.make(schema)
        self.table.definePsfFlux("PsfFlux")
        self.table.defineCentroid("Centroid")
        self.ss = afwTable.SourceCatalog(self.table)

        self.config = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.config.kernel.name = "DF"
        self.subconfig = self.config.kernel.active

        self.policy = pexConfig.makePolicy(self.subconfig)
        self.policy.set('fitForBackground', True)  # we are testing known background recovery here
        self.policy.set('checkConditionNumber', False)  # just in case
        self.policy.set("useRegularization", False)

        if defDataDir:
            defSciencePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v26-e0",
                                          "v26-e0-c011-a10.sci.fits")
            defTemplatePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                           "v5-e0-c011-a10.sci.fits")

            scienceExposure = afwImage.ExposureF(defSciencePath)
            templateExposure = afwImage.ExposureF(defTemplatePath)
            # set XY0 = 0
            scienceExposure.setXY0(geom.Point2I(0, 0))
            templateExposure.setXY0(geom.Point2I(0, 0))
            # do the warping first so we don't have any masked pixels in the postage stamps
            warper = afwMath.Warper.fromConfig(self.subconfig.warpingConfig)
            templateExposure = warper.warpExposure(scienceExposure.getWcs(), templateExposure,
                                                   destBBox=scienceExposure.getBBox())

            # Change xy0
            # Nice star at position 276, 717
            # And should be at index 40, 40
            # No masked pixels in this one
            self.x02 = 276
            self.y02 = 717
            size = 40
            bbox2 = geom.Box2I(geom.Point2I(self.x02 - size, self.y02 - size),
                               geom.Point2I(self.x02 + size, self.y02 + size))
            self.scienceImage2 = afwImage.ExposureF(scienceExposure, bbox2, origin=afwImage.LOCAL)
            self.templateExposure2 = afwImage.ExposureF(templateExposure, bbox2, origin=afwImage.LOCAL)

    def addNoise(self, mi):
        img = mi.getImage()
        seed = int(afwMath.makeStatistics(mi.getVariance(), afwMath.MEDIAN).getValue())
        rdm = afwMath.Random(afwMath.Random.MT19937, seed)
        rdmImage = img.Factory(img.getDimensions())
        afwMath.randomGaussianImage(rdmImage, rdm)
        img += rdmImage
        return afwMath.makeStatistics(rdmImage, afwMath.MEAN).getValue(afwMath.MEAN)

    def verifyDeltaFunctionSolution(self, solution, kSum=1.0, bg=0.0):
        # when kSum = 1.0, this agrees to the default precision.  when
        # kSum != 1.0 I need to go to only 4 digits.
        #
        # -5.4640810225678728e-06 != 0.0 within 7 places
        #
        bgSolution = solution.getBackground()
        self.assertAlmostEqual(bgSolution, bg, 4)

        # again when kSum = 1.0 this agrees.  otherwise
        #
        # 2.7000000605594079 != 2.7000000000000002 within 7 places
        #
        kSumSolution = solution.getKsum()
        self.assertAlmostEqual(kSumSolution, kSum, 5)

        kImage = solution.makeKernelImage()
        for j in range(kImage.getHeight()):
            for i in range(kImage.getWidth()):

                if (i == kImage.getWidth() // 2) and (j == kImage.getHeight() // 2):
                    self.assertAlmostEqual(kImage[i, j, afwImage.LOCAL], kSum, 5)
                else:
                    self.assertAlmostEqual(kImage[i, j, afwImage.LOCAL], 0., 5)

    @unittest.skipIf(not defDataDir, "Warning: afwdata is not set up")
    def testConstructor(self):
        # Original and uninitialized
        kc = ipDiffim.KernelCandidateF(self.x02, self.y02,
                                       self.templateExposure2.getMaskedImage(),
                                       self.scienceImage2.getMaskedImage(),
                                       self.policy)

        # Kernel not initialized
        self.assertEqual(kc.isInitialized(), False)

        # But this should be set on construction
        try:
            kc.getCandidateRating()
        except Exception as e:
            print(e)
            self.fail()

        # And these should be filled
        try:
            kc.getTemplateMaskedImage()
            kc.getScienceMaskedImage()
        except Exception as e:
            print(e)
            self.fail()

        # And of the right type
        self.assertEqual(type(kc.getTemplateMaskedImage()), afwImage.MaskedImageF)
        self.assertEqual(type(kc.getScienceMaskedImage()), afwImage.MaskedImageF)

        # None of these should work
        for kType in (ipDiffim.KernelCandidateF.ORIG,
                      ipDiffim.KernelCandidateF.PCA,
                      ipDiffim.KernelCandidateF.RECENT):
            for kMethod in (kc.getKernelSolution,
                            kc.getKernel,
                            kc.getBackground,
                            kc.getKsum,
                            kc.getKernelImage,
                            kc.getDifferenceImage):
                try:
                    kMethod(kType)
                except Exception:
                    pass
                else:
                    self.fail()
        try:
            kc.getImage()
        except Exception:
            pass
        else:
            self.fail()

    @unittest.skipIf(not defDataDir, "Warning: afwdata is not set up")
    def testSourceStats(self):
        source = self.ss.addNew()
        source.setId(1)
        source.set(self.table.getCentroidKey().getX(), 276)
        source.set(self.table.getCentroidKey().getY(), 717)
        source.set("slot_PsfFlux_instFlux", 1.)

        kc = ipDiffim.KernelCandidateF(source,
                                       self.templateExposure2.getMaskedImage(),
                                       self.scienceImage2.getMaskedImage(),
                                       self.policy)
        kList = ipDiffim.makeKernelBasisList(self.subconfig)

        kc.build(kList)
        self.assertEqual(kc.isInitialized(), True)

    @unittest.skipIf(not defDataDir, "Warning: afwdata is not set up")
    def testSourceConstructor(self):
        source = self.ss.addNew()
        source.setId(1)
        source.set(self.table.getCentroidKey().getX(), 276)
        source.set(self.table.getCentroidKey().getY(), 717)
        source.set("slot_PsfFlux_instFlux", 1.)

        kc = ipDiffim.KernelCandidateF(source,
                                       self.templateExposure2.getMaskedImage(),
                                       self.scienceImage2.getMaskedImage(),
                                       self.policy)

        # Kernel not initialized
        self.assertEqual(kc.isInitialized(), False)

        # Check that the source is set
        self.assertEqual(kc.getSource(), source)
        self.assertEqual(kc.getCandidateRating(), source.getPsfInstFlux())

        # But this should be set on construction
        try:
            kc.getCandidateRating()
        except Exception as e:
            print(e)
            self.fail()

        # And these should be filled
        try:
            kc.getTemplateMaskedImage()
            kc.getScienceMaskedImage()
        except Exception as e:
            print(e)
            self.fail()

        # And of the right type
        self.assertEqual(type(kc.getTemplateMaskedImage()), afwImage.MaskedImageF)
        self.assertEqual(type(kc.getScienceMaskedImage()), afwImage.MaskedImageF)

        # None of these should work
        for kType in (ipDiffim.KernelCandidateF.ORIG,
                      ipDiffim.KernelCandidateF.PCA,
                      ipDiffim.KernelCandidateF.RECENT):
            for kMethod in (kc.getKernelSolution,
                            kc.getKernel,
                            kc.getBackground,
                            kc.getKsum,
                            kc.getKernelImage,
                            kc.getDifferenceImage):
                try:
                    kMethod(kType)
                except Exception:
                    pass
                else:
                    self.fail()
        try:
            kc.getImage()
        except Exception:
            pass
        else:
            self.fail()

        kList = ipDiffim.makeKernelBasisList(self.subconfig)

        kc.build(kList)
        self.assertEqual(kc.isInitialized(), True)

    @unittest.skipIf(not defDataDir, "Warning: afwdata is not set up")
    def testDeltaFunctionScaled(self, scaling=2.7, bg=11.3):
        sIm = afwImage.MaskedImageF(self.templateExposure2.getMaskedImage(), deep=True)
        sIm *= scaling
        kc = ipDiffim.KernelCandidateF(self.x02, self.y02,
                                       self.templateExposure2.getMaskedImage(),
                                       sIm,
                                       self.policy)

        kList = ipDiffim.makeKernelBasisList(self.subconfig)
        kc.build(kList)
        self.verifyDeltaFunctionSolution(kc.getKernelSolution(ipDiffim.KernelCandidateF.RECENT),
                                         kSum=scaling)

        sIm = afwImage.MaskedImageF(self.templateExposure2.getMaskedImage(), deep=True)
        sIm += bg
        kc = ipDiffim.KernelCandidateF(self.x02, self.y02,
                                       self.templateExposure2.getMaskedImage(),
                                       sIm,
                                       self.policy)

        kList = ipDiffim.makeKernelBasisList(self.subconfig)
        kc.build(kList)
        self.verifyDeltaFunctionSolution(kc.getKernelSolution(ipDiffim.KernelCandidateF.RECENT),
                                         bg=bg)

    @unittest.skipIf(not defDataDir, "Warning: afwdata is not set up")
    def testDeltaFunction(self):
        # Match an image to itself, with delta-function basis set
        # No regularization
        kc = ipDiffim.KernelCandidateF(self.x02, self.y02,
                                       self.templateExposure2.getMaskedImage(),
                                       self.templateExposure2.getMaskedImage(),
                                       self.policy)

        kList = ipDiffim.makeKernelBasisList(self.subconfig)

        kc.build(kList)
        self.assertEqual(kc.isInitialized(), True)

        # These should work
        for kType in (ipDiffim.KernelCandidateF.ORIG,
                      ipDiffim.KernelCandidateF.RECENT):
            for kMethod in (kc.getKernelSolution,
                            kc.getKernel,
                            kc.getBackground,
                            kc.getKsum,
                            kc.getKernelImage,
                            kc.getDifferenceImage):
                try:
                    kMethod(kType)
                except Exception as e:
                    print(kMethod, e)
                    self.fail()
                else:
                    pass
        try:
            kc.getImage()
        except Exception as e:
            print(kMethod, e)
            self.fail()
        else:
            pass

        # None of these should work
        for kType in (ipDiffim.KernelCandidateF.PCA,):
            for kMethod in (kc.getKernelSolution,
                            kc.getKernel,
                            kc.getBackground,
                            kc.getKsum,
                            kc.getKernelImage,
                            kc.getImage,
                            kc.getDifferenceImage):
                try:
                    kMethod(kType)
                except Exception:
                    pass
                else:
                    print(kMethod)
                    self.fail()

        self.verifyDeltaFunctionSolution(kc.getKernelSolution(ipDiffim.KernelCandidateF.RECENT))

    @unittest.skipIf(not defDataDir, "Warning: afwdata is not set up")
    def testGaussianWithNoise(self):
        # Convolve a real image with a gaussian and try and recover
        # it.  Add noise and perform the same test.

        gsize = self.policy.getInt("kernelSize")
        gaussFunction = afwMath.GaussianFunction2D(2, 3)
        gaussKernel = afwMath.AnalyticKernel(gsize, gsize, gaussFunction)
        kImageIn = afwImage.ImageD(geom.Extent2I(gsize, gsize))
        kSumIn = gaussKernel.computeImage(kImageIn, False)

        imX, imY = self.templateExposure2.getMaskedImage().getDimensions()
        smi = afwImage.MaskedImageF(geom.Extent2I(imX, imY))
        afwMath.convolve(smi, self.templateExposure2.getMaskedImage(), gaussKernel, False)

        bbox = gaussKernel.shrinkBBox(smi.getBBox(afwImage.LOCAL))

        tmi2 = afwImage.MaskedImageF(self.templateExposure2.getMaskedImage(), bbox, origin=afwImage.LOCAL)
        smi2 = afwImage.MaskedImageF(smi, bbox, origin=afwImage.LOCAL)

        kc = ipDiffim.KernelCandidateF(self.x02, self.y02, tmi2, smi2, self.policy)
        kList = ipDiffim.makeKernelBasisList(self.subconfig)
        kc.build(kList)
        self.assertEqual(kc.isInitialized(), True)
        kImageOut = kc.getImage()

        soln = kc.getKernelSolution(ipDiffim.KernelCandidateF.RECENT)
        self.assertAlmostEqual(soln.getKsum(), kSumIn)
        # 8.7499380640430563e-06 != 0.0 within 7 places
        self.assertAlmostEqual(soln.getBackground(), 0.0, 4)

        for j in range(kImageOut.getHeight()):
            for i in range(kImageOut.getWidth()):

                # in the outskirts of the kernel, the ratio can get screwed because of low S/N
                # e.g. 7.45817359824e-09 vs. 1.18062529402e-08
                # in the guts of the kernel it should look closer
                if kImageIn[i, j, afwImage.LOCAL] > 1e-4:
                    # sigh, too bad this sort of thing fails..
                    # 0.99941584433815966 != 1.0 within 3 places
                    self.assertAlmostEqual(kImageOut[i, j, afwImage.LOCAL]/kImageIn[i, j, afwImage.LOCAL],
                                           1.0, 2)

        # now repeat with noise added; decrease precision of comparison
        self.addNoise(smi2)
        kc = ipDiffim.KernelCandidateF(self.x02, self.y02, tmi2, smi2, self.policy)
        kList = ipDiffim.makeKernelBasisList(self.subconfig)
        kc.build(kList)
        self.assertEqual(kc.isInitialized(), True)
        kImageOut = kc.getImage()

        soln = kc.getKernelSolution(ipDiffim.KernelCandidateF.RECENT)
        self.assertAlmostEqual(soln.getKsum(), kSumIn, 3)
        if not (self.policy.get("fitForBackground")):
            self.assertEqual(soln.getBackground(), 0.0)

        for j in range(kImageOut.getHeight()):
            for i in range(kImageOut.getWidth()):
                if kImageIn[i, j, afwImage.LOCAL] > 1e-2:
                    self.assertAlmostEqual(kImageOut[i, j, afwImage.LOCAL],
                                           kImageIn[i, j, afwImage.LOCAL], 2)

    def testGaussian(self, imsize=50):
        # Convolve a delta function with a known gaussian; try to
        # recover using delta-function basis

        gsize = self.policy.getInt("kernelSize")
        tsize = imsize + gsize

        gaussFunction = afwMath.GaussianFunction2D(2, 3)
        gaussKernel = afwMath.AnalyticKernel(gsize, gsize, gaussFunction)
        kImageIn = afwImage.ImageD(geom.Extent2I(gsize, gsize))
        gaussKernel.computeImage(kImageIn, False)

        # template image with a single hot pixel in the exact center
        tmi = afwImage.MaskedImageF(geom.Extent2I(tsize, tsize))
        tmi.set(0, 0x0, 1e-4)
        cpix = tsize // 2
        tmi[cpix, cpix, afwImage.LOCAL] = (1, 0x0, 1)

        # science image
        smi = afwImage.MaskedImageF(tmi.getDimensions())
        afwMath.convolve(smi, tmi, gaussKernel, False)

        # get the actual kernel sum (since the image is not infinite)
        gscaling = afwMath.makeStatistics(smi, afwMath.SUM).getValue(afwMath.SUM)

        # grab only the non-masked subregion
        bbox = gaussKernel.shrinkBBox(smi.getBBox(afwImage.LOCAL))

        tmi2 = afwImage.MaskedImageF(tmi, bbox, origin=afwImage.LOCAL)
        smi2 = afwImage.MaskedImageF(smi, bbox, origin=afwImage.LOCAL)

        # make sure its a valid subregion!
        for j in range(tmi2.getHeight()):
            for i in range(tmi2.getWidth()):
                self.assertEqual(tmi2.mask[i, j, afwImage.LOCAL], 0)
                self.assertEqual(smi2.mask[i, j, afwImage.LOCAL], 0)

        kc = ipDiffim.KernelCandidateF(0.0, 0.0, tmi2, smi2, self.policy)
        kList = ipDiffim.makeKernelBasisList(self.subconfig)
        kc.build(kList)
        self.assertEqual(kc.isInitialized(), True)
        kImageOut = kc.getImage()

        soln = kc.getKernelSolution(ipDiffim.KernelCandidateF.RECENT)
        self.assertAlmostEqual(soln.getKsum(), gscaling)
        self.assertAlmostEqual(soln.getBackground(), 0.0)

        for j in range(kImageOut.getHeight()):
            for i in range(kImageOut.getWidth()):
                self.assertAlmostEqual(kImageOut[i, j, afwImage.LOCAL]/kImageIn[i, j, afwImage.LOCAL],
                                       1.0, 5)

    def testZeroVariance(self, imsize=50):
        gsize = self.policy.getInt("kernelSize")
        tsize = imsize + gsize

        tmi = afwImage.MaskedImageF(geom.Extent2I(tsize, tsize))
        tmi.set(0, 0x0, 1.0)
        cpix = tsize // 2
        tmi[cpix, cpix, afwImage.LOCAL] = (1, 0x0, 0.0)
        smi = afwImage.MaskedImageF(geom.Extent2I(tsize, tsize))
        smi.set(0, 0x0, 1.0)
        smi[cpix, cpix, afwImage.LOCAL] = (1, 0x0, 0.0)

        kList = ipDiffim.makeKernelBasisList(self.subconfig)
        self.policy.set("constantVarianceWeighting", False)
        kc = ipDiffim.KernelCandidateF(0.0, 0.0, tmi, smi, self.policy)
        try:
            kc.build(kList)
        except Exception:
            pass
        else:
            self.fail()

    @unittest.skipIf(not defDataDir, "Warning: afwdata is not set up")
    def testConstantWeighting(self):
        self.policy.set("fitForBackground", False)
        self.testGaussian()
        self.testGaussianWithNoise()

    @unittest.skipIf(not defDataDir, "Warning: afwdata is not set up")
    def testNoBackgroundFit(self):
        self.policy.set("constantVarianceWeighting", True)
        self.testGaussian()

    def testInsert(self):
        mi = afwImage.MaskedImageF(geom.Extent2I(10, 10))
        kc = ipDiffim.makeKernelCandidate(0., 0., mi, mi, self.policy)
        kc.setStatus(afwMath.SpatialCellCandidate.GOOD)

        sizeCellX = self.policy.get("sizeCellX")
        sizeCellY = self.policy.get("sizeCellY")
        kernelCellSet = afwMath.SpatialCellSet(geom.Box2I(geom.Point2I(0, 0), geom.Extent2I(1, 1)),
                                               sizeCellX, sizeCellY)
        kernelCellSet.insertCandidate(kc)
        nSeen = 0
        for cell in kernelCellSet.getCellList():
            for cand in cell.begin(True):
                self.assertEqual(cand.getStatus(), afwMath.SpatialCellCandidate.GOOD)
                nSeen += 1
        self.assertEqual(nSeen, 1)

    @unittest.skipIf(not display, "display is None: skipping testDisp")
    def testDisp(self):
        afwDisplay.Display(frame=1).mtv(self.scienceImage2,
                                        title=self._testMethodName + ": scienceImage2")
        afwDisplay.Display(frame=2).mtv(self.templateExposure2,
                                        title=self._testMethodName + ": templateExposure2")

    def tearDown(self):
        del self.policy
        del self.table
        del self.ss
        if defDataDir:
            del self.scienceImage2
            del self.templateExposure2

#####


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
