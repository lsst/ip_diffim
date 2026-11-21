.. lsst-task-topic:: lsst.ip.diffim.HomogenizePsfTask

#################
HomogenizePsfTask
#################

HomogenizePsfTask broadens the PSF of a template to a common width, and drops the parts of the template whose PSF is wider.

.. _lsst.ip.diffim.HomogenizePsfTask-description:

Description
===========

A template assembled from coadds has a PSF that steps wherever the set of input images changes, which the smooth spatial model of PSF matching cannot follow.
HomogenizePsfTask measures the width of the template CoaddPsf on a grid of blocks, and takes the width at ``config.quantile`` of the blocks as the common width.
Each narrower block is convolved with the circular Gaussian that brings its PSF up to the common width.
Blocks wider than the common width by more than ``config.tolerance`` are set to ``NO_DATA``, because reaching the common width would need sharpening.
The task is used by `~lsst.ip.diffim.GetTemplateTask` when ``doHomogenizePsf`` is set.

.. _lsst.ip.diffim.HomogenizePsfTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.HomogenizePsfTask

.. _lsst.ip.diffim.HomogenizePsfTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.HomogenizePsfTask

.. _lsst.ip.diffim.HomogenizePsfTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.HomogenizePsfTask
