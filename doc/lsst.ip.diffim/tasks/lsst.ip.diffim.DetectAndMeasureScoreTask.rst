.. lsst-task-topic:: lsst.ip.diffim.DetectAndMeasureScoreTask

#########################
DetectAndMeasureScoreTask
#########################

DetectAndMeasureScoreTask performs source detection and execution of measurement plugins in difference images.

.. _lsst.ip.diffim.DetectAndMeasureScoreTask-description:

Description
===========

This task detects sources in score-like difference images, such as the ones obtained with pre-convolution image differencing, or ZOGY "S" images.
This means that cross-convolution is not needed, and the task obtains maximum-likelihood source candidates in a single step. 
Once candidate sources are obtained this task performs the same operations as `~lsst.ip.diffim.DetectAndMeasureTask`, including measurements in the conventional difference image, not the score-like image; for further reference consult the documentation of `~lsst.ip.diffim.DetectAndMeasureTask`.

.. _lsst.ip.diffim.DetectAndMeasureScoreTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.DetectAndMeasureScoreTask

.. _lsst.ip.diffim.DetectAndMeasureScoreTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.DetectAndMeasureScoreTask

.. _lsst.ip.diffim.DetectAndMeasureScoreTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.DetectAndMeasureScoreTask

.. _lsst.ip.diffim.DetectAndMeasureScoreTask-debug:

.. Debugging
.. =========


.. The ``pipetask`` command line interface supports a ``--debug`` flag to import
.. ``debug.py`` from your PYTHONPATH; see :ref:`lsstDebug` for more about ``debug.py``
.. files.
.. The available variables in DetectAndMeasureScoreTask include:


.. display : `bool`
..     Enable debug display output.
.. maskTransparency : `float`
..     Transparency of mask planes in the output display.
.. displayDiaSources : `bool`
..     Show exposure with dipole fitting results.
