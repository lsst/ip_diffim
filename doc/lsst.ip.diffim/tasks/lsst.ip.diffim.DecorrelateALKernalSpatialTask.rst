.. lsst-task-topic:: lsst.ip.diffim.DecorrelateALKernelSpatialTask

##########################
DecorrelateALKernelSpatialTask
##########################

Pipe-task that removes the neighboring-pixel covariance in an
image difference that are added when the template image is
convolved with the Alard-Lupton PSF matching kernel.

.. _lsst.ip.diffim.DecorrelateALKernelSpatialTask-description:

Description
==================

This task is a simple wrapper around @ref DecorrelateALKernelTask,
which takes a `spatiallyVarying` parameter in its `run` method. If
it is `False`, then it simply calls the `run` method of @ref
DecorrelateALKernelTask. If it is True, then it uses the @ref
ImageMapReduceTask framework to break the exposures into
subExposures on a grid, and performs the `run` method of @ref
DecorrelateALKernelTask on each subExposure. This enables it to
account for spatially-varying PSFs and noise in the exposures when
performing the decorrelation.

This task has no standalone example, however it is applied as a
subtask of pipe.tasks.imageDifference.ImageDifferenceTask.
There is also an example of its use in `tests/testImageDecorrelation.py`.

.. _lsst.ip.diffim.DecorrelateALKernelSpatialTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.DecorrelateALKernelSpatialTask

.. _lsst.ip.diffim.DecorrelateALKernelSpatialTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.DecorrelateALKernelSpatialTask

.. _lsst.ip.diffim.DecorrelateALKernelSpatialTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.DecorrelateALKernelSpatialTask
