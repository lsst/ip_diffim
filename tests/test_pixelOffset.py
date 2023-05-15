import unittest
import numpy as np

import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
from lsst.ip.diffim.zogy import ZogyTask, ZogyConfig


class PixelOffsetTest(unittest.TestCase):
    """A test case for the pixel offset.
    """

    @staticmethod
    def _makeImage(w, h, psfsize, psfsig, set_peak=False):
        """Make fake images for testing. If set_peak is True, the
        image would have a peak at [300][300].
        """
        # Make image
        image = afwImage.MaskedImageF(w, h)
        image.set(0)
        array = image.image.array
        if set_peak:
            # set the peak value
            array[300][300] = 1000
        var = image.getVariance()
        var.set(1.0)
        # Make PSF
        psf = measAlg.DoubleGaussianPsf(psfsize, psfsize, psfsig, psfsig, 0.)
        # Make exposure
        exp = afwImage.makeExposure(image)
        exp.setPsf(psf)
        return exp

    @staticmethod
    def _find_max(data):
        """ Find the maximum position of the data.
        """
        return np.unravel_index(np.argmax(data), data.shape)

    def _setUpImages(self):
        """ Setup test images.
        """
        self.imrex = PixelOffsetTest._makeImage(1000, 2000, 17, 1, set_peak=False)
        self.imnex = PixelOffsetTest._makeImage(1000, 2000, 17, 2, set_peak=True)

    def testPixelOffset(self):
        """ Test whether the peak position is at [300][300].
        """
        self._setUpImages()
        config = ZogyConfig(scaleByCalibration=False)
        task = ZogyTask(config=config)
        D_F = task.run(self.imnex, self.imrex, calculateScore=False)
        max_loc = PixelOffsetTest._find_max(D_F.diffExp.image.array)

        self.assertEqual(max_loc, (300, 300))
