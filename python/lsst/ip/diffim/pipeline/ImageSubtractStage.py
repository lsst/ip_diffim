import re
import lsst.pex.harness.Stage
import lsst.ip.diffim
import lsst.afw.image as afwImage

__all__ = ["ImageSubtractStage"]

class ImageSubtractStage(lsst.pex.harness.Stage.Stage):
    def process(self):
        print 'Python ImageSubtractStage process : _rank %i stageId %d' % (self._rank, self.stageId)
        activeClipboard = self.inputQueue.getNextDataset()

#        lsst.pex.logging.Trace_setVerbosity("lsst.ip.diffim", 5)

        ###########
        #
        # Get objects from clipboard
        #
        templateExposure = activeClipboard.get('TemplateExposure')
        #
        # We may have been passed an Image, but we need an Exposure
        #
        if re.search(r"ImageBase<", templateExposure.repr()): # Yes, an Image of some sort
            # N.b. we don't use type() as we don't know what sort of Image it'll be, but repr can be fooled by a Mask
            im = templateExposure
            msk = afwImage.MaskU(im.getDimensions()); msk.set(0x0)
            var = afwImage.ImageF(im.getDimensions()); var.set(0.0)
            maskedImage = afwImage.makeMaskedImage(im, msk, var);  del im; del msk; del var

            templateExposure = afwImage.makeExposure(maskedImage)

        scienceExposure = activeClipboard.get('ScienceExposure')
        templateMaskedImage = templateExposure.getMaskedImage()
        scienceMaskedImage = scienceExposure.getMaskedImage()
        try:
            psfMatchBasisKernelSet = activeClipboard.get('PsfMatchBasisKernelSet')
        except KeyError:
            psfMatchBasisKernelSet = None
#        print "ImageSubtractState.process clipboard contains:"
#        for key in activeClipboard.getKeys():
#            print "* %s: %r" % (key, activeClipboard.get(key))

        differenceImage, psfMatchingKernelPtr, backgroundFunction = lsst.ip.diffim.imageSubtract(
            imageToConvolve = templateMaskedImage,
            imageToNotConvolve = scienceMaskedImage,
            policy = self._policy,
            psfMatchBasisKernelSet = psfMatchBasisKernelSet,
            footprintList = None,
        )
        
        # Turn the difference image into an exposure. Since we don't have a factory function for this,
        # we have to know the Exposure type, which is based on the image type.
        # The type of the difference image matches the type of the image that was convolved,
        # which is the template image. It would be safer to use a factory function.,
        differenceExposureClass = type(templateExposure)
        if scienceExposure.hasWcs():
            differenceExposure = differenceExposureClass(differenceImage, scienceExposure.getWcs())
        else:
            differenceExposure = differenceExposureClass(differenceImage)

        ###########
        #
        # Post results to clipboard
        #
        activeClipboard.put('DifferenceExposure', differenceExposure)
        activeClipboard.put('PsfMatchKernel', psfMatchingKernelPtr)
        activeClipboard.put('BackgroundModel', backgroundFunction)
        
        self.outputQueue.addDataset(activeClipboard)
