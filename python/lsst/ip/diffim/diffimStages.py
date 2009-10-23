import re, time

from warpTemplateExposure import *
from subtractExposure import *
from createSdqaRatingVector import *

from lsst.pex.harness.Stage import Stage
import lsst.afw.image as afwImage
import lsst.sdqa as sdqa
import lsst.pex.logging as pexLog
import lsst.afw.display.ds9 as ds9

display = True

class DiffimStage(Stage):

    def process(self):
        self.activeClipboard = self.inputQueue.getNextDataset()
        
        self.log = pexLog.Log(pexLog.Log.getDefaultLog(),
                              "ip.diffim.DiffimStage")

        scienceExposureKey = self._policy.get("scienceExposureKey")
        templateExposureKey = self._policy.get("templateExposureKey")
        
        scienceExposure = self.activeClipboard.get(scienceExposureKey)
        templateExposure = self.activeClipboard.get(templateExposureKey)
        #
        # We may have been passed an Image, but we need an Exposure
        #
        if re.search(r"ImageBase<", templateExposure.repr()):
            # Yes, an Image of some sort
            # N.b. we don't use type() as we don't know what sort of Image it'll be,
            # but repr can be fooled by a Mask
            im = templateExposure
            msk = afwImage.MaskU(im.getDimensions()); msk.set(0x0)
            var = afwImage.ImageF(im.getDimensions()); var.set(0.0)
            maskedImage = afwImage.makeMaskedImage(im, msk, var)
            del im; del msk; del var

            templateExposure = afwImage.makeExposure(maskedImage)

            wcsKey = self._policy.get("templateWcsKey")
            wcs = self.activeClipboard.get(wcsKey)
            bBoxKey = self._policy.get("templateBBoxKey")
            bBox = self.activeClipboard.get(bBoxKey)

            nwcs = wcs.clone()
            nwcs.shiftReferencePixel(bBox.get("llcx"), bBox.get("llcy"))
            templateExposure.setWcs(nwcs)
       
        diffimPolicy = self._policy.get("diffimPolicy")

        # step 1
        self.log.log(pexLog.Log.INFO, "Starting warp : %s" % (time.ctime()))
        remapedTemplateExposure = warpTemplateExposure(templateExposure,
                                                       scienceExposure, 
                                                       diffimPolicy)
        self.log.log(pexLog.Log.INFO, "Ending warp : %s" % (time.ctime()))
        
        if display:
            frame = 1
            ds9.mtv(remapedTemplateExposure, frame=frame)
            ds9.dot("Warped Template", 0, 0, frame=frame)

            frame = 2
            ds9.mtv(scienceExposure, frame=frame)
            ds9.dot("Science Exposure", 0, 0, frame=frame)

        # step 2
        self.log.log(pexLog.Log.INFO, "Starting subtract : %s" % (time.ctime()))
        try:
            result = subtractExposure(remapedTemplateExposure, 
                                      scienceExposure, 
                                      diffimPolicy)
        except:
            pexLog.Trace("lsst.ip.diffim.DiffimStage", 1,
                         "ERROR: Unable to calculate psf matching kernel")
            raise RuntimeException("DiffimStage.subtractExposure failed")
        else:
            differenceExposure, spatialKernel, spatialBg, kernelCellSet = result
            
        self.log.log(pexLog.Log.INFO, "Ending subtract : %s" % (time.ctime()))

        if display:
            frame = 3
            ds9.mtv(differenceExposure, frame=frame)
            ds9.dot("Difference Exposure", 0, 0, frame=frame)

        sdqaVector = createSdqaRatingVector(kernelCellSet, spatialKernel, spatialBg)
        persistableSdqaVector = sdqa.PersistableSdqaRatingVector(sdqaVector)

        exposureKey = self._policy.getString("differenceExposureKey")
        self.activeClipboard.put(exposureKey, differenceExposure)

        sdqaKey = self._policy.getString("sdqaRatingSetKey")
        self.activeClipboard.put(sdqaKey, persistableSdqaVector)

        self.outputQueue.addDataset(self.activeClipboard)
