#config.doPreConvolve=False
#config.doMatchSources=False
#config.doAddMetrics=False
#config.doUseRegister=False
#config.doSelectSources=False
#config.kernelSourcesFromRef=False
config.makeDiffim.doWriteSubtractedExp=True
config.makeDiffim.doWriteMatchedExp=True
config.makeDiffim.doDecorrelation=True
config.makeDiffim.subtract='al'

config.makeDiffim.subtract['zogy'].zogyConfig.inImageSpace=False

from lsst.ip.diffim.getTemplate import GetCalexpAsTemplateTask
config.getTemplate.retarget(GetCalexpAsTemplateTask)
