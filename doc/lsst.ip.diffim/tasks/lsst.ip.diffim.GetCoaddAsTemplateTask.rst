.. lsst-task-topic:: lsst.ip.diffim.GetCoaddAsTemplateTask

##########################
GetCoaddAsTemplateTask
##########################

From the given skymap, the closest tract is selected; multiple tracts are
not supported. The assembled template inherits the WCS of the selected
skymap tract and the resolution of the template exposures. Overlapping box
regions of the input template patches are pixel by pixel copied into the
assembled template image. There is no warping or pixel resampling.

Pixels with no overlap of any available input patches are set to ``nan``
value and ``NO_DATA`` flagged.

.. _lsst.ip.diffim.GetCoaddAsTemplateTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.GetCoaddAsTemplateTask

.. _lsst.ip.diffim.GetCoaddAsTemplateTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.GetCoaddAsTemplateTask

.. _lsst.ip.diffim.GetCoaddAsTemplateTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.GetCoaddAsTemplateTask
