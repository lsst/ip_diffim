#!/usr/bin/env python
import unittest
import lsst.utils.tests
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDet
import lsst.ip.diffim as ipDiffim
import lsst.pex.logging as pexLog
import lsst.pex.config as pexConfig

pexLog.Trace_setVerbosity('lsst.ip.diffim', 5)

class DiffimTestCases(lsst.utils.tests.TestCase):

    def setUp(self):
        self.config    = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.config.kernel.name = "AL"
        self.subconfig = self.config.kernel.active

        # Test was put together before the min size went to 21
        self.subconfig.kernelSize = 19

        self.subconfig.scaleByFwhm = False
        self.subconfig.fitForBackground = True
        self.subconfig.spatialModelType = "polynomial"
        self.policy = pexConfig.makePolicy(self.subconfig)

        self.smi = afwImage.MaskedImageF('tests/compareToHotpants/scienceMI.fits')
        self.tmi = afwImage.MaskedImageF('tests/compareToHotpants/templateMI.fits')
        self.smi.setXY0(0,0)
        self.tmi.setXY0(0,0)

        # Run detection
        #detConfig = self.subconfig.detectionConfig
        # Note here regarding detConfig:
        #
        # If I set detThresholdType = "pixel_stdev", I get slightly
        # different centroids than if I use "stdev".  These different
        # centroids screw up the testing since hotpants was hardcoded to
        # use the "stdev" centroids.  For completeness these are:
        #
        # 32 32
        # 96 32
        # 160 32
        # 96 95
        # 31 96
        # 160 96
        # 96 160
        # 160 160
        # 32 160

        # As of Winter2013, KernelCandidateDetectionF does not return
        # these exact centroids anymore, so I need to hardcode them
        # in.
        self.footprints = []
        for xc,yc in [(32, 32), (96, 32), (160, 32),
                      (96, 95), (31, 96), (160, 96),
                      (96, 160), (160, 160), (32, 160)]:
            self.footprints.append(afwDet.Footprint(afwGeom.Box2I(
                        afwGeom.Point2I(xc,yc), afwGeom.Extent2I(1,1))))

        #detConfig.detThresholdType = "stdev"
        #kcDetect = ipDiffim.KernelCandidateDetectionF(pexConfig.makePolicy(detConfig))
        #kcDetect.apply(self.smi, self.tmi)
        #self.footprints = kcDetect.getFootprints()

        # Make a basis list that hotpants has been run with
        nGauss = 1
        sGauss = [3.]
        dGauss = [3]
        self.subconfig.alardNGauss = nGauss
        self.subconfig.alardSigGauss = sGauss
        self.subconfig.alardDegGauss = dGauss
        basisList0 = ipDiffim.makeKernelBasisList(self.subconfig)

        # HP does things in a different order, and with different normalization, so reorder list
        order   = [0, 2, 5, 9, 1, 4, 8, 3, 7, 6]
        scaling = [1.000000e+00,
                   8.866037e-02,
                   1.218095e+01,
                   5.099318e-03,
                   8.866037e-02,
                   4.179772e-02,
                   1.138120e-02,
                   1.218095e+01,
                   1.138120e-02,
                   5.099318e-03]

        self.basisList = afwMath.KernelList()
        for i in range(len(order)):
            im  = afwImage.ImageD(basisList0[order[i]].getDimensions())
            basisList0[order[i]].computeImage(im, False)
            im /= scaling[i]
            #im.writeFits('k%d.fits' % (i))
            k   = afwMath.FixedKernel(im)
            self.basisList.append(k)

        # And a place to put candidates
        self.kernelCellSet = afwMath.SpatialCellSet(afwGeom.Box2I(afwGeom.Point2I(0,0),
                                                                  afwGeom.Extent2I(self.smi.getWidth(),
                                                                                   self.smi.getHeight())),
                                                    self.policy.getInt("sizeCellX"),
                                                    self.policy.getInt("sizeCellY"))

        # There are some -1 factors that come from differences in how
        # convolution is done.  Some resulting convovled images end up
        # being a factor of -1 different, therefore the coefficients
        # need to be a factor of -1 different as well.
        self.parity = [1, -1, 1, -1, -1, 1, -1, 1, -1, -1, 1]

    def tearDown(self):
        del self.policy
        del self.tmi
        del self.smi
        del self.basisList
        del self.footprints
        del self.kernelCellSet

    def testSingleNoVariation(self):
        self.policy.set('constantVarianceWeighting', True)
        self.policy.set('spatialKernelOrder', 0)
        self.policy.set('spatialBgOrder', 0)

        # Place candidate footprints within the spatial grid
        for fp in self.footprints:
            bbox = fp.getBBox()

            # Grab the centers in the parent's coordinate system
            xC   = int(0.5 * ( bbox.getMinX() + bbox.getMaxX() ))
            yC   = int(0.5 * ( bbox.getMinY() + bbox.getMaxY() ))

            bbox = afwGeom.Box2I(afwGeom.Point2I(int(xC)-24, int(yC)-24), afwGeom.Extent2I(49, 49))

            tsmi  = afwImage.MaskedImageF(self.tmi, bbox, afwImage.LOCAL)
            ssmi  = afwImage.MaskedImageF(self.smi, bbox, afwImage.LOCAL)

            # Hotpants centroids go from -1 to 1
            # Only one passes
            if xC > 150 and yC > 150:
                cand = ipDiffim.makeKernelCandidate( (xC - 0.5 * self.smi.getWidth()) /
                                                     (0.5 * self.smi.getWidth()),
                                                     (yC - 0.5 * self.smi.getHeight()) /
                                                     (0.5 * self.smi.getHeight()),
                                                     tsmi, ssmi, self.policy)
                self.kernelCellSet.insertCandidate(cand)

        # Visitors
        bbox  = self.kernelCellSet.getBBox()
        bsikv = ipDiffim.BuildSingleKernelVisitorF(self.basisList, self.policy)
        bspkv = ipDiffim.BuildSpatialKernelVisitorF(self.basisList, bbox, self.policy)

        for cell in self.kernelCellSet.getCellList():
            for cand in cell.begin(False): # False = include bad candidates
                cand  = ipDiffim.cast_KernelCandidateF(cand)
                bsikv.processCandidate(cand)
                bspkv.processCandidate(cand)

        HPsingleSolution = [ 0.959086,
                            -0.000344,
                            -0.197758,
                             0.000019,
                            -0.000172,
                             0.000053,
                             0.000018,
                            -0.192776,
                             0.000000,
                             0.000001,
                             0.602642 ]

        HPspatialSolution = HPsingleSolution

        singleSolution = cand.getKernel(ipDiffim.KernelCandidateF.RECENT).getKernelParameters()
        for i in range(len(singleSolution)):
            self.assertAlmostEqual(HPsingleSolution[i] * self.parity[i], singleSolution[i], 5)

        bspkv.solveLinearEquation()
        sk, sb = bspkv.getSolutionPair()
        spatialSolution = sk.getKernelParameters()
        for i in range(len(spatialSolution)):
            self.assertAlmostEqual(HPspatialSolution[i] * self.parity[i], spatialSolution[i], 6)

        self.assertAlmostEqual(sb.getParameters()[0], HPspatialSolution[-1], 5)

    def testFourNoVariation(self):
        self.policy.set('constantVarianceWeighting', True)
        self.policy.set('spatialKernelOrder', 0)
        self.policy.set('spatialBgOrder', 0)

        # Place candidate footprints within the spatial grid
        for fp in self.footprints:
            bbox = fp.getBBox()

            # Grab the centers in the parent's coordinate system
            xC   = int(0.5 * ( bbox.getMinX() + bbox.getMaxX() ))
            yC   = int(0.5 * ( bbox.getMinY() + bbox.getMaxY() ))

            bbox = afwGeom.Box2I(afwGeom.Point2I(int(xC)-24, int(yC)-24), afwGeom.Extent2I(49, 49))

            tsmi  = afwImage.MaskedImageF(self.tmi, bbox, afwImage.LOCAL)
            ssmi  = afwImage.MaskedImageF(self.smi, bbox, afwImage.LOCAL)

            # Hotpants centroids go from -1 to 1
            if xC > 90 and yC > 90:
                cand = ipDiffim.makeKernelCandidate( (xC - 0.5 * self.smi.getWidth()) /
                                                     (0.5 * self.smi.getWidth()),
                                                     (yC - 0.5 * self.smi.getHeight()) /
                                                     (0.5 * self.smi.getHeight()),
                                                     tsmi, ssmi, self.policy)

                self.kernelCellSet.insertCandidate(cand)

        # Visitors
        bbox  = self.kernelCellSet.getBBox()
        bsikv = ipDiffim.BuildSingleKernelVisitorF(self.basisList, self.policy)
        bspkv = ipDiffim.BuildSpatialKernelVisitorF(self.basisList, bbox, self.policy)

        for cell in self.kernelCellSet.getCellList():
            for cand in cell.begin(False): # False = include bad candidates
                cand  = ipDiffim.cast_KernelCandidateF(cand)
                bsikv.processCandidate(cand)
                bspkv.processCandidate(cand)

        HPspatialSolution = [ 0.969559,
                             -0.000223,
                             -0.198374,
                              0.000012,
                             -0.000010,
                              0.000036,
                             -0.000004,
                             -0.206751,
                              0.000012,
                              0.000004,
                              0.452304 ]

        bspkv.solveLinearEquation()
        sk, sb = bspkv.getSolutionPair()
        spatialSolution = sk.getKernelParameters()
        for i in range(len(spatialSolution)):
            self.assertAlmostEqual(HPspatialSolution[i] * self.parity[i], spatialSolution[i], 5)

        self.assertAlmostEqual(sb.getParameters()[0], HPspatialSolution[-1], 5)

    def testFourKernelVariation(self):
        self.policy.set('constantVarianceWeighting', True)
        self.policy.set('spatialKernelOrder', 1)
        self.policy.set('spatialBgOrder', 0)

        # Place candidate footprints within the spatial grid
        for fp in self.footprints:
            bbox = fp.getBBox()

            # Grab the centers in the parent's coordinate system
            xC   = int(0.5 * ( bbox.getMinX() + bbox.getMaxX() ))
            yC   = int(0.5 * ( bbox.getMinY() + bbox.getMaxY() ))

            bbox = afwGeom.Box2I(afwGeom.Point2I(int(xC)-24, int(yC)-24), afwGeom.Extent2I(49, 49))

            tsmi  = afwImage.MaskedImageF(self.tmi, bbox, afwImage.LOCAL)
            ssmi  = afwImage.MaskedImageF(self.smi, bbox, afwImage.LOCAL)

            # Hotpants centroids go from -1 to 1
            if xC > 90 and yC > 90:
                cand = ipDiffim.makeKernelCandidate( (xC - 0.5 * self.smi.getWidth()) /
                                                     (0.5 * self.smi.getWidth()),
                                                     (yC - 0.5 * self.smi.getHeight()) /
                                                     (0.5 * self.smi.getHeight()),
                                                     tsmi, ssmi, self.policy)
                self.kernelCellSet.insertCandidate(cand)

        # Visitors
        bbox  = self.kernelCellSet.getBBox()
        bsikv = ipDiffim.BuildSingleKernelVisitorF(self.basisList, self.policy)
        bspkv = ipDiffim.BuildSpatialKernelVisitorF(self.basisList, bbox, self.policy)

        for cell in self.kernelCellSet.getCellList():
            for cand in cell.begin(False): # False = include bad candidates
                cand  = ipDiffim.cast_KernelCandidateF(cand)
                bsikv.processCandidate(cand)
                bspkv.processCandidate(cand)

        HPspatialSolution = [ [  0.969559,
                                 0.,
                                 0. ],
                              [ -0.000082,
                                -0.000620,
                                 0.000185 ],
                              [ -0.197749,
                                 0.001418,
                                -0.003321 ],
                              [  0.000002,
                                 0.000049,
                                -0.000016 ],
                              [  0.000211,
                                -0.000283,
                                -0.000397 ],
                              [  0.000034,
                                 0.000002,
                                 0.000006 ],
                              [ -0.000013,
                                 0.000041,
                                -0.000010 ],
                              [ -0.220238,
                                 0.028395,
                                 0.013148 ],
                              [  0.000019,
                                -0.000025,
                                 0.000003 ],
                              [  0.000003,
                                 0.000000,
                                 0.000005 ],
                              0.452304 ]

        bspkv.solveLinearEquation()
        sk, sb = bspkv.getSolutionPair()
        spatialSolution = sk.getSpatialParameters()

        for i in range(len(spatialSolution)):
            # HP and LSST switch the order x<->y
            self.assertAlmostEqual(HPspatialSolution[i][0] * self.parity[i], spatialSolution[i][0], 5)
            self.assertAlmostEqual(HPspatialSolution[i][1] * self.parity[i], spatialSolution[i][2], 5)
            self.assertAlmostEqual(HPspatialSolution[i][2] * self.parity[i], spatialSolution[i][1], 5)

        self.assertAlmostEqual(sb.getParameters()[0], HPspatialSolution[-1], 5)

    def testFourBgVariation(self):

        # OK, so these can end up a bit different due to how HP and
        # LSST represent the background in the matrix math.  HP has
        # each pixel have its own coordinate (which goes from -1 to 1
        # across the entire image, by the way), whereas we give all
        # the LSST pixels within a stamp the same coordinate.

        # To make this comparison, I go ahead and edit the Hotpants
        # code to give each pixel the same weight.  For reference this
        # is in fillStamp() and I replace:
        #
        #        //xf = (i - rPixX2) / rPixX2;
        #        xf = (xi - rPixX2) / rPixX2;
        #
        #            //yf = (j - rPixY2) / rPixY2;
        #            yf = (yi - rPixY2) / rPixY2;

        self.policy.set('constantVarianceWeighting', True)
        self.policy.set('spatialKernelOrder', 0)
        self.policy.set('spatialBgOrder', 1)

        # Place candidate footprints within the spatial grid
        for fp in self.footprints:
            bbox = fp.getBBox()

            # Grab the centers in the parent's coordinate system
            xC   = int(0.5 * ( bbox.getMinX() + bbox.getMaxX() ))
            yC   = int(0.5 * ( bbox.getMinY() + bbox.getMaxY() ))

            bbox = afwGeom.Box2I(afwGeom.Point2I(int(xC)-24, int(yC)-24), afwGeom.Extent2I(49, 49))

            tsmi  = afwImage.MaskedImageF(self.tmi, bbox, afwImage.LOCAL)
            ssmi  = afwImage.MaskedImageF(self.smi, bbox, afwImage.LOCAL)

            # Hotpants centroids go from -1 to 1
            if xC > 90 and yC > 90:
                cand = ipDiffim.makeKernelCandidate( (xC - 0.5 * self.smi.getWidth()) /
                                                     (0.5 * self.smi.getWidth()),
                                                     (yC - 0.5 * self.smi.getHeight()) /
                                                     (0.5 * self.smi.getHeight()),
                                                     tsmi, ssmi, self.policy)
                #print 'OBJECT', cand.getId(), 'AT', xC, yC, cand.getXCenter(), cand.getYCenter()
                self.kernelCellSet.insertCandidate(cand)

        # Visitors
        bbox  = self.kernelCellSet.getBBox()
        bsikv = ipDiffim.BuildSingleKernelVisitorF(self.basisList, self.policy)
        bspkv = ipDiffim.BuildSpatialKernelVisitorF(self.basisList, bbox, self.policy)

        for cell in self.kernelCellSet.getCellList():
            for cand in cell.begin(False): # False = include bad candidates
                cand  = ipDiffim.cast_KernelCandidateF(cand)
                bsikv.processCandidate(cand)
                bspkv.processCandidate(cand)

        HPspatialSolution = [ 0.969559,
                             -0.000223,
                             -0.198374,
                              0.000012,
                             -0.000010,
                              0.000036,
                             -0.000004,
                             -0.206751,
                              0.000012,
                              0.000004,
                              [ 0.782113,
                               -0.910963,
                               -0.106636] ]


        bspkv.solveLinearEquation()
        sk, sb = bspkv.getSolutionPair()
        spatialSolution = sk.getKernelParameters()
        for i in range(len(spatialSolution)):
            self.assertAlmostEqual(HPspatialSolution[i] * self.parity[i], spatialSolution[i], 5)

        self.assertAlmostEqual(sb.getParameters()[0], HPspatialSolution[-1][0], 5)
        self.assertAlmostEqual(sb.getParameters()[1], HPspatialSolution[-1][2], 5) # x<->y
        self.assertAlmostEqual(sb.getParameters()[2], HPspatialSolution[-1][1], 5) # x<->y

    def testFourVariation(self):
        self.policy.set('constantVarianceWeighting', True)
        self.policy.set('spatialKernelOrder', 1)
        self.policy.set('spatialBgOrder', 1)

        # Place candidate footprints within the spatial grid
        for fp in self.footprints:
            bbox = fp.getBBox()

            # Grab the centers in the parent's coordinate system
            xC   = int(0.5 * ( bbox.getMinX() + bbox.getMaxX() ))
            yC   = int(0.5 * ( bbox.getMinY() + bbox.getMaxY() ))

            bbox = afwGeom.Box2I(afwGeom.Point2I(int(xC)-24, int(yC)-24), afwGeom.Extent2I(49, 49))

            tsmi  = afwImage.MaskedImageF(self.tmi, bbox, afwImage.LOCAL)
            ssmi  = afwImage.MaskedImageF(self.smi, bbox, afwImage.LOCAL)

            # Hotpants centroids go from -1 to 1
            if xC > 90 and yC > 90:
                cand = ipDiffim.makeKernelCandidate( (xC - 0.5 * self.smi.getWidth()) /
                                                     (0.5 * self.smi.getWidth()),
                                                     (yC - 0.5 * self.smi.getHeight()) /
                                                     (0.5 * self.smi.getHeight()),
                                                     tsmi, ssmi, self.policy)
                self.kernelCellSet.insertCandidate(cand)

        # Visitors
        bbox  = self.kernelCellSet.getBBox()
        bsikv = ipDiffim.BuildSingleKernelVisitorF(self.basisList, self.policy)
        bspkv = ipDiffim.BuildSpatialKernelVisitorF(self.basisList, bbox, self.policy)

        for cell in self.kernelCellSet.getCellList():
            for cand in cell.begin(False): # False = include bad candidates
                cand  = ipDiffim.cast_KernelCandidateF(cand)
                bsikv.processCandidate(cand)
                bspkv.processCandidate(cand)

        HPspatialSolution = [ [  0.969559,
                                 0.,
                                 0. ],
                              [ -0.000082,
                                -0.000620,
                                 0.000185 ],
                              [ -0.197749,
                                 0.001418,
                                -0.003321 ],
                              [  0.000002,
                                 0.000049,
                                -0.000016 ],
                              [  0.000211,
                                -0.000283,
                                -0.000397 ],
                              [  0.000034,
                                 0.000002,
                                 0.000006 ],
                              [ -0.000013,
                                 0.000041,
                                -0.000010 ],
                              [ -0.220238,
                                 0.028395,
                                 0.013148 ],
                              [  0.000019,
                                -0.000025,
                                 0.000003 ],
                              [  0.000003,
                                 0.000000,
                                 0.000005 ],
                              [  0.782113,
                                -0.910963,
                                -0.106636] ]

        bspkv.solveLinearEquation()
        sk, sb = bspkv.getSolutionPair()
        spatialSolution = sk.getSpatialParameters()

        for i in range(len(spatialSolution)):
            # HP and LSST switch the order x<->y
            self.assertAlmostEqual(HPspatialSolution[i][0] * self.parity[i], spatialSolution[i][0], 5)
            self.assertAlmostEqual(HPspatialSolution[i][1] * self.parity[i], spatialSolution[i][2], 5)
            self.assertAlmostEqual(HPspatialSolution[i][2] * self.parity[i], spatialSolution[i][1], 5)

        self.assertAlmostEqual(sb.getParameters()[0], HPspatialSolution[-1][0], 5)
        self.assertAlmostEqual(sb.getParameters()[1], HPspatialSolution[-1][2], 5) # x<->y
        self.assertAlmostEqual(sb.getParameters()[2], HPspatialSolution[-1][1], 5) # x<->y

    def testAllBgVariation2(self):
        # OK, I ran HP on all the things in this image.  Enough for
        # second order spatial variation

        self.policy.set('constantVarianceWeighting', True)
        self.policy.set('spatialKernelOrder', 0)
        self.policy.set('spatialBgOrder', 2)

        # Ignore the whole kernelCellSet thing
        cands = []
        for fp in self.footprints:
            bbox = fp.getBBox()

            # Grab the centers in the parent's coordinate system
            xC   = int(0.5 * ( bbox.getMinX() + bbox.getMaxX() ))
            yC   = int(0.5 * ( bbox.getMinY() + bbox.getMaxY() ))

            bbox = afwGeom.Box2I(afwGeom.Point2I(int(xC)-24, int(yC)-24), afwGeom.Extent2I(49, 49))

            tsmi  = afwImage.MaskedImageF(self.tmi, bbox, afwImage.LOCAL)
            ssmi  = afwImage.MaskedImageF(self.smi, bbox, afwImage.LOCAL)

            # Hotpants centroids go from -1 to 1
            cand = ipDiffim.makeKernelCandidate( (xC - 0.5 * self.smi.getWidth()) /
                                                 (0.5 * self.smi.getWidth()),
                                                 (yC - 0.5 * self.smi.getHeight()) /
                                                 (0.5 * self.smi.getHeight()),
                                                 tsmi, ssmi, self.policy)
            cands.append(cand)

        # Visitors
        bbox  = self.kernelCellSet.getBBox()
        bsikv = ipDiffim.BuildSingleKernelVisitorF(self.basisList, self.policy)
        bspkv = ipDiffim.BuildSpatialKernelVisitorF(self.basisList, bbox, self.policy)

        for cand in cands:
            bsikv.processCandidate(cand)
            bspkv.processCandidate(cand)

        HPspatialSolution = [ 0.968505,
                             -0.000053,
                             -0.206505,
                              0.000005,
                              0.000062,
                              0.000028,
                             -0.000004,
                             -0.206135,
                              0.000005,
                              0.000001,
                             [ 0.812488,
                               0.096456,
                              -1.140900,
                               0.132670,
                              -0.571923,
                              -0.284670] ]

        bspkv.solveLinearEquation()
        sk, sb = bspkv.getSolutionPair()
        spatialSolution = sk.getKernelParameters()

        # Kernel
        for i in range(len(spatialSolution)):
            self.assertAlmostEqual(HPspatialSolution[i] * self.parity[i], spatialSolution[i], 5)

        # Bg
        # Ordering of second order terms is just messy
        spReorder = [0, 3, 1, 5, 4, 2]
        spatialSolution = sb.getParameters()
        for i in range(len(spatialSolution)):
            self.assertAlmostEqual(HPspatialSolution[-1][spReorder[i]], spatialSolution[i], 5)

    def testAllKernelVariation2(self):
        # OK, I ran HP on all the things in this image.  Enough for
        # second order spatial variation

        self.policy.set('constantVarianceWeighting', True)
        self.policy.set('spatialKernelOrder', 2)
        self.policy.set('spatialBgOrder', 0)

        # Ignore the whole kernelCellSet thing
        cands = []
        for fp in self.footprints:
            bbox = fp.getBBox()

            # Grab the centers in the parent's coordinate system
            xC   = int(0.5 * ( bbox.getMinX() + bbox.getMaxX() ))
            yC   = int(0.5 * ( bbox.getMinY() + bbox.getMaxY() ))

            bbox = afwGeom.Box2I(afwGeom.Point2I(int(xC)-24, int(yC)-24), afwGeom.Extent2I(49, 49))

            tsmi  = afwImage.MaskedImageF(self.tmi, bbox, afwImage.LOCAL)
            ssmi  = afwImage.MaskedImageF(self.smi, bbox, afwImage.LOCAL)

            # Hotpants centroids go from -1 to 1
            cand = ipDiffim.makeKernelCandidate( (xC - 0.5 * self.smi.getWidth()) /
                                                 (0.5 * self.smi.getWidth()),
                                                 (yC - 0.5 * self.smi.getHeight()) /
                                                 (0.5 * self.smi.getHeight()),
                                                 tsmi, ssmi, self.policy)
            cands.append(cand)

        # Visitors
        bbox  = self.kernelCellSet.getBBox()
        bsikv = ipDiffim.BuildSingleKernelVisitorF(self.basisList, self.policy)
        bspkv = ipDiffim.BuildSpatialKernelVisitorF(self.basisList, bbox, self.policy)

        for cand in cands:
            bsikv.processCandidate(cand)
            bspkv.processCandidate(cand)

        HPspatialSolution = [ [  0.968505,
                                 0.,
                                 0.,
                                 0.,
                                 0.,
                                 0.,
                                 0.],
                              [ -0.000094,
                                -0.000206,
                                -0.000027,
                                 0.000025,
                                -0.000705,
                                 0.000162],
                              [ -0.188375,
                                 0.001801,
                                -0.027534,
                                 0.008718,
                                 0.016346,
                                -0.033879],
                              [  0.000004,
                                 0.000012,
                                 0.000023,
                                 0.000004,
                                 0.000017,
                                -0.000020],
                              [  0.000128,
                                -0.000218,
                                 0.000304,
                                -0.000038,
                                -0.000151,
                                -0.000531],
                              [ -0.000011,
                                -0.000013,
                                 0.000038,
                                -0.000017,
                                 0.000133,
                                 0.000093],
                              [ -0.000003,
                                 0.000003,
                                 0.000006,
                                 0.000008,
                                 0.000002,
                                -0.000010],
                              [ -0.212235,
                                -0.000856,
                                 0.012246,
                                -0.010893,
                                 0.049302,
                                 0.008249],
                              [  0.000014,
                                -0.000002,
                                -0.000050,
                                -0.000001,
                                 0.000030,
                                 0.000020],
                              [ -0.000001,
                                 0.000010,
                                -0.000012,
                                -0.000007,
                                 0.000015,
                                 0.000019],
                              0.392482]

        bspkv.solveLinearEquation()
        sk, sb = bspkv.getSolutionPair()
        spatialSolution = sk.getSpatialParameters()

        # Kernel
        spReorder = [0, 3, 1, 5, 4, 2]
        for i in range(len(spatialSolution)):
            for j in range(len(spReorder)):
                self.assertAlmostEqual(HPspatialSolution[i][spReorder[j]] * self.parity[i],
                                       spatialSolution[i][j], 5)

        self.assertAlmostEqual(sb.getParameters()[0], HPspatialSolution[-1], 5)

    def testAllVariation2(self):
        # OK, I ran HP on all the things in this image.  Enough for
        # second order spatial variation

        self.policy.set('constantVarianceWeighting', True)
        self.policy.set('spatialKernelOrder', 2)
        self.policy.set('spatialBgOrder', 2)

        # Ignore the whole kernelCellSet thing
        cands = []
        for fp in self.footprints:
            bbox = fp.getBBox()

            # Grab the centers in the parent's coordinate system
            xC   = int(0.5 * ( bbox.getMinX() + bbox.getMaxX() ))
            yC   = int(0.5 * ( bbox.getMinY() + bbox.getMaxY() ))

            bbox = afwGeom.Box2I(afwGeom.Point2I(int(xC)-24, int(yC)-24), afwGeom.Extent2I(49, 49))

            tsmi  = afwImage.MaskedImageF(self.tmi, bbox, afwImage.LOCAL)
            ssmi  = afwImage.MaskedImageF(self.smi, bbox, afwImage.LOCAL)

            # Hotpants centroids go from -1 to 1
            cand = ipDiffim.makeKernelCandidate( (xC - 0.5 * self.smi.getWidth()) /
                                                 (0.5 * self.smi.getWidth()),
                                                 (yC - 0.5 * self.smi.getHeight()) /
                                                 (0.5 * self.smi.getHeight()),
                                                 tsmi, ssmi, self.policy)
            cands.append(cand)

        # Visitors
        bbox  = self.kernelCellSet.getBBox()
        bsikv = ipDiffim.BuildSingleKernelVisitorF(self.basisList, self.policy)
        bspkv = ipDiffim.BuildSpatialKernelVisitorF(self.basisList, bbox, self.policy)

        for cand in cands:
            bsikv.processCandidate(cand)
            bspkv.processCandidate(cand)

        HPspatialSolution =[ [ 0.968505,
                               0.,
                               0.,
                               0.,
                               0.,
                               0.],
                             [-0.000094,
                              -0.000206,
                              -0.000027,
                               0.000025,
                              -0.000705,
                               0.000162],
                             [-0.188375,
                               0.001801,
                              -0.027534,
                               0.008718,
                               0.016346,
                              -0.033879],
                             [ 0.000004,
                               0.000012,
                               0.000023,
                               0.000004,
                               0.000017,
                              -0.000020],
                             [ 0.000128,
                              -0.000218,
                               0.000304,
                              -0.000038,
                              -0.000151,
                              -0.000531],
                             [-0.000011,
                              -0.000013,
                               0.000038,
                              -0.000017,
                               0.000133,
                               0.000093],
                             [-0.000003,
                               0.000003,
                               0.000006,
                               0.000008,
                               0.000002,
                              -0.000010],
                             [-0.212235,
                              -0.000856,
                               0.012246,
                              -0.010893,
                               0.049302,
                               0.008249],
                             [ 0.000014,
                              -0.000002,
                              -0.000050,
                              -0.000001,
                               0.000030,
                               0.000020],
                             [-0.000001,
                               0.000010,
                              -0.000012,
                              -0.000007,
                               0.000015,
                               0.000019],
                             [ 0.812488,
                               0.096456,
                              -1.140900,
                               0.132670,
                              -0.571923,
                              -0.284670] ]

        bspkv.solveLinearEquation()
        sk, sb = bspkv.getSolutionPair()
        spatialSolution = sk.getSpatialParameters()

        # Kernel
        spReorder = [0, 3, 1, 5, 4, 2]
        for i in range(len(spatialSolution)):
            for j in range(len(spReorder)):
                self.assertAlmostEqual(HPspatialSolution[i][spReorder[j]] * self.parity[i],
                                       spatialSolution[i][j], 5)
        # Bg
        spatialSolution = sb.getParameters()
        for i in range(len(spatialSolution)):
            self.assertAlmostEqual(HPspatialSolution[-1][spReorder[i]], spatialSolution[i], 5)

#####

class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass

def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()