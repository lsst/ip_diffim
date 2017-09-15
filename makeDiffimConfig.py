#config.doPreConvolve=False
#config.doMatchSources=False
#config.doAddMetrics=False
#config.doUseRegister=False
#config.doSelectSources=False
#config.kernelSourcesFromRef=False
config.doWriteSubtractedExp=True
config.doWriteMatchedExp=True
config.doDecorrelation=True
config.subtract='al'

config.subtract['zogy'].zogyConfig.inImageSpace=False

from lsst.ip.diffim.getTemplate import GetCalexpAsTemplateTask
config.getTemplate.retarget(GetCalexpAsTemplateTask)
