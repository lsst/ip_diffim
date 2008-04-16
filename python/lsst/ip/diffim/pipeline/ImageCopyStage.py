import lsst.pex.harness.Stage
import lsst.pex.logging as pexLog

__all__ = ["ImageCopyStage"]

class ImageCopyStage(lsst.pex.harness.Stage.Stage):
    """A simple stage to test running an image processing pipeline.
    
    This stage simply copies "inputImage" on input to "outputImage" on output
    and prints some diagnostic information on the way.
    The image can be of any type that has getCols and getRows methods.
    """
    def process(self):
        pexLog.Trace("ip.diffim.imageCopy", 3,
            "ImageCopyStage process: _rank %i stageId %d" % (self._rank, self.stageId))
        
        activeClipboard = self.inputQueue.getNextDataset()

        inputImage = activeClipboard.get('InputImage')
        
        pexLog.Trace("ip.diffim.imageCopy", 3,
            "ImageCopyStage read image from clipboard as InputImage; cols=%r, rows=%r" % \
                (inputImage.getCols(), inputImage.getRows()))

        activeClipboard.put('OutputImage', inputImage)

        pexLog.Trace("ip.diffim.imageCopy", 3,
            "ImageCopyStage wrote image to clipboard as OutputImage")

        self.outputQueue.addDataset(activeClipboard)
