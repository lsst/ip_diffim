import os
import lsst.dps.Pipeline
#import lsst.imageproc.ImageSubtractStage
import lsst.fw.Core.fwLib as fw
import lsst.mwi.data as datap

if __name__ == '__main__':
    pyPipeline = lsst.dps.Pipeline.Pipeline()
    pyPipeline.configurePipeline()
    pyPipeline.initializeQueues()

    # Add image subtraction stage
    #isStage = ImageSubtractStage()
    #pyPipeline.stageList.append(isStage)

    # Populate the clipboard
    scienceMaskedImage  = fw.MaskedImageD()
    templateMaskedImage = fw.MaskedImageD()
    scienceMaskedImage.readFits(os.path.join(os.environ['FWDATA_DIR'], '871034p_1_MI'))
    templateMaskedImage.readFits(os.path.join(os.environ['FWDATA_DIR'], '871034p_1_MI'))
    pyPipeline.populateClipboard('TemplateImage', templateMaskedImage)
    pyPipeline.populateClipboard('InputImage', scienceMaskedImage)
    
    pyPipeline.initializeStages()
    pyPipeline.startSlices()
    pyPipeline.startInitQueue()
    pyPipeline.startStagesLoop()
    pyPipeline.shutdown()
