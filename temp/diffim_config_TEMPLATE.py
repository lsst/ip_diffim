from lsst.ip.diffim.getTemplate import GetCalexpAsTemplateTask
config.getTemplate.retarget(GetCalexpAsTemplateTask)
# config.doMeasurement=True
# config.doPreConvolve=False
#config.doWriteMatchedExp=True
#config.doUseRegister=True
config.detection.thresholdValue=5.0
