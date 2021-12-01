import os
import unittest


import lsst.utils.tests
import lsst.utils
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.geom as geom
import lsst.ip.diffim as ipDiffim
import lsst.ip.diffim.diffimTools as diffimTools
import lsst.utils.logging as logUtils
import lsst.pex.config as pexConfig

logUtils.trace_set_at("lsst.ip.diffim", 2)

# known input images
try:
    defDataDir = lsst.utils.getPackageDir('afwdata')
except Exception:
    defDataDir = None


class DiffimTestCases(lsst.utils.tests.TestCase):

    def setUp(self):
        self.config = ipDiffim.ImagePsfMatchTask.ConfigClass()
        self.subconfig = self.config.kernel.active
        self.ps = pexConfig.makePropertySet(self.subconfig)
        self.kSize = self.ps['kernelSize']

        # gaussian reference kernel
        self.gSize = self.kSize
        self.gaussFunction = afwMath.GaussianFunction2D(2, 3)
        self.gaussKernel = afwMath.AnalyticKernel(self.gSize, self.gSize, self.gaussFunction)

        if defDataDir:
            defImagePath = os.path.join(defDataDir, "DC3a-Sim", "sci", "v5-e0",
                                        "v5-e0-c011-a00.sci.fits")
            self.templateImage = afwImage.MaskedImageF(defImagePath)
            self.scienceImage = self.templateImage.Factory(self.templateImage.getDimensions())

            convolutionControl = afwMath.ConvolutionControl()
            convolutionControl.setDoNormalize(False)
            afwMath.convolve(self.scienceImage, self.templateImage, self.gaussKernel, convolutionControl)

    def tearDown(self):
        del self.config
        del self.ps
        del self.gaussFunction
        del self.gaussKernel
        if defDataDir:
            del self.templateImage
            del self.scienceImage

    @unittest.skipIf(not defDataDir,
                     "Warning: afwdata is not set up; not running KernelCandidateDetection.py")
    def testGetCollection(self):
        # NOTE - you need to subtract off background from the image
        # you run detection on.  Here it is the template.
        bgConfig = self.subconfig.afwBackgroundConfig
        diffimTools.backgroundSubtract(bgConfig, [self.templateImage, ])

        detConfig = self.subconfig.detectionConfig
        maskPlane = detConfig.badMaskPlanes[0]
        maskVal = afwImage.Mask.getPlaneBitMask(maskPlane)

        kcDetect = ipDiffim.KernelCandidateDetectionF(pexConfig.makePropertySet(detConfig))
        kcDetect.apply(self.templateImage, self.scienceImage)
        fpList1 = kcDetect.getFootprints()

        self.assertNotEqual(len(fpList1), 0)

        for fp in fpList1:
            bbox = fp.getBBox()
            tmi = afwImage.MaskedImageF(self.templateImage, bbox, origin=afwImage.LOCAL)
            smi = afwImage.MaskedImageF(self.scienceImage, bbox, origin=afwImage.LOCAL)
            tmask = tmi.getMask()
            smask = smi.getMask()

            for j in range(tmask.getHeight()):
                for i in range(tmask.getWidth()):
                    # No masked pixels in either image
                    self.assertEqual(tmask[i, j, afwImage.LOCAL], 0)
                    self.assertEqual(smask[i, j, afwImage.LOCAL], 0)

        # add a masked pixel to the template image and make sure you don't get it
        tp = geom.Point2I(tmask.getWidth()//2, tmask.getHeight()//2)
        self.templateImage.mask[fpList1[0].getBBox(), afwImage.LOCAL][tp, afwImage.LOCAL] = maskVal
        kcDetect.apply(self.templateImage, self.scienceImage)
        fpList2 = kcDetect.getFootprints()
        self.assertEqual(len(fpList2), (len(fpList1)-1))

        # add a masked pixel to the science image and make sure you don't get it
        sp = geom.Point2I(smask.getWidth()//2, smask.getHeight()//2)
        self.scienceImage.mask[fpList1[1].getBBox(), afwImage.LOCAL][sp, afwImage.LOCAL] = maskVal
        self.scienceImage.mask[fpList1[2].getBBox(), afwImage.LOCAL][sp, afwImage.LOCAL] = maskVal
        kcDetect.apply(self.templateImage, self.scienceImage)
        fpList3 = kcDetect.getFootprints()
        self.assertEqual(len(fpList3), (len(fpList1)-3))

#####


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
