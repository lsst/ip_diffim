.. lsst-task-topic:: lsst.ip.diffim.GetTemplateTask

###############
GetTemplateTask
###############

Build a template from existing coadds, which may span multiple tracts.
The assembled template inherits the WCS of the selected
skymap tract and the resolution of the template exposures. Overlapping box
regions of the input template patches are pixel by pixel copied into the
assembled template image. There is no warping or pixel resampling.

Pixels with no overlap of any available input patches are set to ``nan``
value and ``NO_DATA`` flagged.

.. _lsst.ip.diffim.GetTemplateTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.GetTemplateTask

.. _lsst.ip.diffim.GetTemplateTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.GetTemplateTask

.. _lsst.ip.diffim.GetTemplateTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.GetTemplateTask
