import lsst.fw.Core.fwLib as fw
import lsst.dps.Stage
import lsst.imageproc
#from lsst.imageproc.imageSubtract import imageSubtract

__all__ = ["ImageSubtractStage"]

class ImageSubtractStage(lsst.dps.Stage.Stage):
    def process(self):
        print 'Python ImageSubtractStage process : _rank %i stageId %d' % (self._rank, self.stageId)
        activeClipboard = self.inputQueue.getNextDataset()

#        lsst.mwi.utils.Trace_setVerbosity("lsst.imageproc", 5)

        ###########
        #
        # Get objects from clipboard
        #
        templateExposure = activeClipboard.get('TemplateExposure')
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

        differenceImage, psfMatchingKernelPtr, backgroundFunctionPtr = lsst.imageproc.imageSubtract(
            imageToConvolve = templateMaskedImage,
            imageToNotConvolve = scienceMaskedImage,
            policy = self._policy,
            psfMatchBasisKernelSet = psfMatchBasisKernelSet,
            footprintList = None,
        )
        
        if scienceExposure.hasWcs():
            differenceExposure = fw.ExposureD(differenceImage, scienceExposure.getWcs())
        else:
            differenceExposure = fw.ExposureD(differenceImage)

        ###########
        #
        # Post results to clipboard
        #
        activeClipboard.put('DifferenceExposure', differenceExposure)
        activeClipboard.put('PsfMatchKernel', psfMatchingKernelPtr)
        activeClipboard.put('BackgroundModel', backgroundFunctionPtr)
        
        self.outputQueue.addDataset(activeClipboard)
