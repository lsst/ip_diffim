import lsst.daf.base as dafBase
import lsst.pex.harness.Stage
import lsst.pex.logging as pexLog
import lsst.ip.diffim as ipDiffim

__all__ = ["TemplateBBoxStage"]

class TemplateBBoxStage(lsst.pex.harness.Stage.Stage):
    """A pipeline stage that computes the bounding box of a template Exposure
    for difference imaging, given a science Exposure.
    
    Inputs: see TemplateBBoxDictionary.paf
    
    Outputs:
    - templateBBoxProperties: bounding box of the template exposure as a PropertySet containing
        four integers named: llcx, llcy, width, height
    """
    def process(self):
        log = pexLog.Log(pexLog.Log.getDefaultLog(), "ip.diffim.TemplateBBoxStage")

        log.log(pexLog.Log.INFO,
            "TemplateBBoxStage process: _rank %i stageId %d" % (self._rank, self.stageId))
        
        activeClipboard = self.inputQueue.getNextDataset()

        scienceExposureName = self._policy.getString("scienceExposureName")
        scienceExposure = activeClipboard.get(scienceExposureName)
        templateDimensionsName = self._policy.getString("templateDimensionsName")
        templateDimensions = activeClipboard.get(templateDimensionsName)
        templateWcsName = self._policy.getString("templateWcsName")
        templateWcs = activeClipboard.get(templateWcsName)
        borderWidth = self._policy.getInt("borderWidth")
        
        scienceDimensions = scienceExposure.getMaskedImage().getDimensions()
        if not scienceExposure.hasWcs():
            raise RuntimeError("science exposure has no Wcs")
        scienceWcs = scienceExposure.getWcs()

        log.log(pexLog.Log.INFO,
            "scienceDimensions=%s, %s; templateDimensions=%s, %s" % \
                (scienceDimensions[0], scienceDimensions[1], templateDimensions[0], templateDimensions[1]))

        templateBBox = lsst.ip.diffim(scienceWcs, scienceDimensions, templateWcs, templateDimensions, borderWidth)

        nameValPairs = (
            ("llcx", templateBBox.getX0()),
            ("llcy", templateBBox.getY0()),
            ("width", templateBBox.getWidth()),
            ("height", templateBBox.getHeight()),
        )
        templateBBoxProperties = dafBase.PropertySet()
        outStrList = []
        for name, val in nameValPairs:
            templateBBoxProperties.add(name, val)
            outStrList.append("%s=%s" % (name, val))

        log.log(pexLog.Log.INFO,
            "templateBBoxProperties: %s" % (";".join(outStrList)))

        activeClipboard.put('templateBBoxProperties', templateBBoxProperties)

        self.outputQueue.addDataset(activeClipboard)
